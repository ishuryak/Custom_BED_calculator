"""Plan step 5 gate: simulation-based calibration (SBC) of the primary candidate model.

Usage:
    python3 tools/early_tox/step5_sbc.py [--reps 500] [--out DIR]

Primary model: U02 with arm-level overdispersion (Tk 7, chemotherapy term, scale-group offset).
Each replicate draws every parameter from the priors in numpy (independently of the Stan code, so
a prior coded differently in Stan shows up as non-uniform ranks), simulates counts on the real
design, fits with 2 chains x 500 draws, and records the rank of each true value among 100 thinned
posterior draws (0-100). Uniform ranks mean calibrated posteriors. One row is appended to
ranks.csv per replicate, and a rerun skips replicates already present (plan rule 12).
Parameters registered in sim_provenance.yaml (et_prior_*, et_sbc_*).
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import importlib.util
import json
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("step5_fit", os.path.join(ROOT, "tools", "early_tox", "step5_fit.py"))
s5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s5)

SEED = 20261112
PARAMS = ["Lstar", "ab", "r", "gamma", "sigma_arm", "m0", "tau0"]
OUT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5", "sbc")


def prior_draw(rng, d):
    P = s5.PRIOR
    th = dict(m0=rng.normal(0, 2.5), tau0=abs(rng.normal(0, 1)), Lstar=abs(rng.normal(0, P["Lstar_sd"])),
              ab=float(np.exp(rng.normal(P["ab_mu"], P["ab_sd"]))), r=float(np.exp(rng.normal(P["r_mu"], P["r_sd"]))),
              gamma=rng.normal(0, P["gamma_sd"]), sigma_arm=abs(rng.normal(0, P["od_sd"])))
    z, delta, zarm = rng.normal(size=d["J"]), rng.normal(size=d["G"]), rng.normal(size=d["I"])
    eta = th["m0"] + np.array(d["grp"]).reshape(d["J"], d["G"]) @ delta + th["tau0"] * z
    dE = np.array(d["dD"]) + np.array(d["dQ"]) / th["ab"] - th["r"] * np.array(d["dTe"])
    lin = (eta[np.array(d["fam"]) - 1] + th["Lstar"] / 10 * dE + th["gamma"] * np.array(d["chemo"])
           + th["sigma_arm"] * zarm)
    y = rng.binomial(np.array(d["N"]), 1 / (1 + np.exp(-lin)))
    return th, y


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--summarize-only", action="store_true", help="recompute summary.json from ranks.csv; no sampling")
    a = ap.parse_args(argv[1:])
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "ranks.csv")
    if a.summarize_only:
        return summarize(path, a.out)
    done = set()
    if os.path.exists(path):
        with open(path, newline="") as f:
            done = {int(r["rep"]) for r in csv.DictReader(f)}
    rows, fams = s5.load_rows()
    base = s5.stan_data(rows, fams, 7.0, dict(use_od=1))
    sampler = dict(chains=2, parallel_chains=1, iter_warmup=500, iter_sampling=500, adapt_delta=0.95)
    identity = s5.run_identity(base, dict(use_od=1, sbc=True), sampler, extra=dict(seed=SEED, reps=a.reps))
    man = os.path.join(a.out, "manifest.json")
    if done:
        if not os.path.exists(man) or json.load(open(man)).get("identity") != identity:
            raise RuntimeError(f"{path}: existing ranks were made with another (or unrecorded) identity; refusing to append")
    else:
        with open(man, "w") as f:
            json.dump(dict(identity=identity, seed=SEED, reps=a.reps, sampler=sampler), f, indent=1)
    model = s5.compile_model()
    s5.wait_for_memory()
    t0 = time.time()
    for rep in range(a.reps):
        if rep in done:
            continue
        rng = np.random.default_rng([SEED, rep])
        th, y = prior_draw(rng, base)
        d = dict(base, y=y.tolist())
        fit = model.sample(data=d, chains=2, parallel_chains=1, iter_warmup=500, iter_sampling=500,
                           adapt_delta=0.95, seed=SEED + rep, show_progress=False, show_console=False)
        row = dict(rep=rep, divergences=int(np.sum(fit.divergences)))
        for p in PARAMS:
            draws = np.asarray(fit.stan_variable(p)).reshape(-1)[::10]          # 100 thinned draws
            row[p] = int(np.sum(draws < th[p]))
            row[f"true_{p}"] = round(float(th[p]), 5)
        new = not os.path.exists(path)
        with open(path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(row), lineterminator="\n")
            if new:
                w.writeheader()
            w.writerow(row)
        print(f"rep {rep + 1}/{a.reps} ({time.time() - t0:.0f} s) divergences {row['divergences']}", flush=True)
        if rep % 25 == 24:
            s5.wait_for_memory()
    summarize(path, a.out)


RANK_BINS = [(0, 10)] + [(10 * b + 1, 10 * b + 10) for b in range(1, 10)]     # 101 ranks: 11, then 10 each


def bin_probabilities():
    """Exact discrete-rank bin probabilities for ranks 0..100 (review F5): 11/101, then 10/101."""
    return np.array([(hi - lo + 1) / 101 for lo, hi in RANK_BINS])


def summarize(path, out):
    """Chi-square test of ranks against the exact discrete bin probabilities, per parameter, and the
    list of replicates with divergences. 'No detected rank departure' is the strongest statement."""
    from scipy.stats import chi2
    with open(path, newline="") as f:
        R = list(csv.DictReader(f))
    probs = bin_probabilities()
    div = [dict(rep=int(r["rep"]), divergences=int(r["divergences"])) for r in R if int(r["divergences"]) > 0]
    res = dict(reps=len(R), seed=SEED, reps_with_divergences=len(div),
               divergent_transitions=sum(d["divergences"] for d in div), divergent_reps=div,
               bin_probabilities=[round(float(p), 6) for p in probs], params={},
               statement="no detected rank departure at the 5% level" )
    for p in PARAMS:
        x = np.array([int(r[p]) for r in R])
        counts = np.array([np.sum((x >= lo) & (x <= hi)) for lo, hi in RANK_BINS])
        e = len(x) * probs
        stat = float(np.sum((counts - e) ** 2 / e))
        res["params"][p] = dict(bin_counts=counts.tolist(), chisq=round(stat, 2), p=round(float(chi2.sf(stat, 9)), 4),
                                mean_rank=round(float(x.mean()), 1))
        if res["params"][p]["p"] < 0.05:
            res["statement"] = "rank departure detected for at least one parameter"
    tmp = os.path.join(out, f"summary.json.partial-{os.getpid()}")
    with open(tmp, "w") as f:
        json.dump(res, f, indent=1)
    os.replace(tmp, os.path.join(out, "summary.json"))
    for p, v in res["params"].items():
        print(f"{p:10s} chisq {v['chisq']:6.2f} p {v['p']:.3f} mean rank {v['mean_rank']} bins {v['bin_counts']}")
    print(f"replicates with divergences: {res['reps_with_divergences']} / {res['reps']} "
          f"({res['divergent_transitions']} transitions); {res['statement']}")
    return res


if __name__ == "__main__":
    main(sys.argv)
