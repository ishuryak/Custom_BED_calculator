"""Plan step 6: validation without leakage, Fowler benchmark (U01) against the fitted candidate (U02).

Usage:
    python3 tools/early_tox/step6_validate.py [--out DIR] [--inner 100]

Both models use the step 5 primary settings: overdispersed (C3), Tk 7, chemotherapy term (C1),
scale-group offset (M1). For each of the 22 primary families, both models are refitted without that
family (its arms, counts and reference arm all removed; no prior in either model is estimated from
any trial in this data set) and then predict it:

  new-trial task          log p(y_j) with a new trial baseline eta ~ N(m0 + delta_g, tau0^2)
  conditional contrast    log p(y_j,non-ref | y_j,ref): the reference arm informs the baseline,
                          the other arms are predicted (log p(all) - log p(ref))

Each predictive density averages over posterior draws and, within each draw, over `inner` Monte
Carlo draws of the trial baseline and arm effects (seeded). The score per family is the log
predictive density; the paired difference d_j = candidate - benchmark is summarized with a Bayesian
bootstrap of family weights (Dirichlet(1,...,1), seeds 20261113 and 20261114, 20,000 replicates).
Promotion convention (plan step 6): the lower 95% bound of the weighted mean d must be above 0 in
both seeds, and the leave-one-family influence must not change the sign.

External checks, from fits on all 22 primary families: HYPNO (benchmark, R6) and the sensitivity
trials (IAEA-ACC, Pinto, Poddar, Sanchiz, Dwivedi, ELAN-RT), conditional contrast task.

Fits are checkpointed per (model, held-out family) in OUT/cache/ (git-ignored) and skipped on rerun.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import importlib.util
import json
import math
import sys
import time

import numpy as np
from scipy.special import gammaln, logsumexp

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("step5_fit", os.path.join(ROOT, "tools", "early_tox", "step5_fit.py"))
s5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s5)

SEED = 20261113
BOOT_SEEDS = (20261113, 20261114)
N_BOOT = 20000
TK = 7.0
MODELS = {"benchmark": dict(fixed_E=1, use_od=1), "candidate": dict(use_od=1)}
OUT_DEFAULT = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step6")
PARAMS = ["m0", "tau0", "kappa", "ab", "r", "gamma", "sigma_arm"]


def training_set(rows, fams, held_out):
    """Rows and families with every arm of the held-out family removed."""
    return [x for x in rows if x["family"] != held_out], [g for g in fams if g != held_out]


def binom_logpmf(y, n, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1) + y * np.log(p) + (n - y) * np.log1p(-p)


def family_data(fam, rows_all):
    """Contrast terms of one family's arms against its own reference arm."""
    rows = [x for x in rows_all if x["family"] == fam]
    d = s5.stan_data(rows, [fam], TK, {})
    ref = [i for i, x in enumerate(rows) if x["arm_id"] == s5.REF_ARM[fam]]
    assert len(ref) == 1
    return dict(rows=rows, ref=ref[0], y=np.array(d["y"]), N=np.array(d["N"]), dD=np.array(d["dD"]),
                dQ=np.array(d["dQ"]), dTe=np.array(d["dTe"]), chemo=np.array(d["chemo"]),
                group=rows[0]["scale_group"])


def predictive(draws, fd, inner, rng):
    """(log p(all arms), log p(reference arm)) by Monte Carlo over draws x inner baselines."""
    S = len(draws["kappa"])
    I = len(fd["y"])
    dE = fd["dD"][None, :] + fd["dQ"][None, :] / draws["ab"][:, None] - draws["r"][:, None] * fd["dTe"][None, :]
    base = draws["m0"] + draws["delta"]                                    # (S,)
    eta = base[:, None] + draws["tau0"][:, None] * rng.standard_normal((S, inner))       # (S, inner)
    arm = draws["sigma_arm"][:, None, None] * rng.standard_normal((S, inner, I))
    lin = eta[:, :, None] + (draws["kappa"][:, None] * dE + draws["gamma"][:, None] * fd["chemo"][None, :])[:, None, :] + arm
    ll = binom_logpmf(fd["y"][None, None, :], fd["N"][None, None, :], 1 / (1 + np.exp(-lin)))   # (S, inner, I)
    tot = S * inner
    return float(logsumexp(ll.sum(axis=2)) - math.log(tot)), float(logsumexp(ll[:, :, fd["ref"]]) - math.log(tot))


def group_delta(fit_draws, group, train_groups, rng):
    """delta for the held-out family's scale group: fitted if the group is in training, else prior N(0, 1)."""
    S = len(fit_draws["kappa"])
    if group == "G1":
        return np.zeros(S)
    if group in train_groups:
        return fit_draws["delta_all"][:, train_groups.index(group)]
    return rng.standard_normal(S)


def fit_cached(model, name, rows, fams, settings, cache):
    path = os.path.join(cache, f"{name}.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=False)
        return {k: z[k] for k in z.files}
    data = s5.stan_data(rows, fams, TK, settings)
    s5.wait_for_memory()
    fit = model.sample(data=data, chains=4, parallel_chains=1, iter_warmup=1000, iter_sampling=1000,
                       adapt_delta=0.95, seed=SEED, show_progress=False, show_console=False)
    s = fit.summary()
    keep = [n for n in s.index if n.split("[")[0] in ("m0", "tau0", "Lstar", "ab", "r", "gamma", "sigma_arm", "delta", "eta")]
    sub = s.loc[keep]
    gate = (int(np.sum(fit.divergences)) == 0 and float(sub["R_hat"].max()) < 1.01
            and float(sub["ESS_bulk"].min()) > 400 and float(sub["ESS_tail"].min()) > 400)
    out = {p: np.asarray(fit.stan_variable(p)).reshape(len(fit.stan_variable("kappa")), -1)[:, 0]
           if p in ("gamma", "sigma_arm") else np.asarray(fit.stan_variable(p)) for p in PARAMS}
    out["delta_all"] = np.asarray(fit.stan_variable("delta")).reshape(len(out["kappa"]), -1)
    out["families"] = np.array(fams)
    out["groups"] = np.array([g for g in s5.GROUP_COLS if any(x["scale_group"] == g for x in rows)])
    out["gate"] = np.array([gate, int(np.sum(fit.divergences)), float(sub["R_hat"].max()),
                            float(sub["ESS_bulk"].min()), float(sub["ESS_tail"].min())])
    np.savez_compressed(path + ".tmp.npz", **out)
    os.replace(path + ".tmp.npz", path)
    return out


def score_family(fd, draws_fit, inner, seed):
    rng = np.random.default_rng(seed)
    d = dict(draws_fit)
    d["delta"] = group_delta(draws_fit, fd["group"], list(draws_fit["groups"]), rng)
    lp_all, lp_ref = predictive(d, fd, inner, rng)
    return lp_all, lp_all - lp_ref


def bayes_boot(dvals, seed, n=N_BOOT):
    rng = np.random.default_rng(seed)
    w = rng.dirichlet(np.ones(len(dvals)), size=n)
    m = w @ np.asarray(dvals)
    return dict(mean=float(np.mean(dvals)), lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)),
                p_positive=float(np.mean(m > 0)))


def contrast_check(fd, draws_fit, seed, inner):
    """Observed against predicted log odds ratio of each non-reference arm, given the reference arm."""
    rng = np.random.default_rng(seed)
    S = len(draws_fit["kappa"])
    dE = fd["dD"][None, :] + fd["dQ"][None, :] / draws_fit["ab"][:, None] - draws_fit["r"][:, None] * fd["dTe"][None, :]
    lor = draws_fit["kappa"][:, None] * dE + draws_fit["gamma"][:, None] * fd["chemo"][None, :]
    # predictive spread for a new arm pair adds two arm effects
    lor_pred = lor + draws_fit["sigma_arm"][:, None] * math.sqrt(2) * rng.standard_normal(lor.shape)
    y, N, r = fd["y"], fd["N"], fd["ref"]
    out = []
    for i, x in enumerate(fd["rows"]):
        if i == r:
            continue
        obs = math.log((y[i] + 0.5) / (N[i] - y[i] + 0.5)) - math.log((y[r] + 0.5) / (N[r] - y[r] + 0.5))
        q = np.percentile(lor_pred[:, i], [5, 50, 95])
        out.append(dict(arm_id=x["arm_id"], observed_lor=round(obs, 4), predicted_lor_median=round(float(q[1]), 4),
                        predicted_lor_5=round(float(q[0]), 4), predicted_lor_95=round(float(q[2]), 4),
                        inside_90=bool(q[0] <= obs <= q[2]), pit=round(float(np.mean(lor_pred[:, i] <= obs)), 4)))
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--inner", type=int, default=100)
    ap.add_argument("--families", default="", help="comma list to restrict held-out families (smoke)")
    a = ap.parse_args(argv[1:])
    cache = os.path.join(a.out, "cache")
    os.makedirs(cache, exist_ok=True)
    model = s5.compile_model()
    prim, fams = s5.load_rows()
    allrows, allfams = s5.load_rows(m4=True)
    held = [f for f in fams if not a.families or f in a.families.split(",")]
    t0 = time.time()
    per_family, contrasts, gates = [], [], []
    for k, f in enumerate(held):
        rows, tfams = training_set(prim, fams, f)
        fd = family_data(f, prim)
        rec = dict(family=f, n_arms=len(fd["y"]), group=fd["group"])
        for m, st in MODELS.items():
            dr = fit_cached(model, f"{m}__without_{f}", rows, tfams, st, cache)
            assert f not in set(dr["families"].tolist()), f"leak: {f} in its own training fit"
            gates.append(dict(model=m, held_out=f, gate_pass=bool(dr["gate"][0]), divergences=int(dr["gate"][1]),
                              rhat=round(float(dr["gate"][2]), 4), ess_bulk=int(dr["gate"][3]), ess_tail=int(dr["gate"][4])))
            new, cond = score_family(fd, dr, a.inner, SEED + 1000 * k)      # common random numbers
            rec[f"{m}_new_trial"] = round(new, 4)
            rec[f"{m}_conditional"] = round(cond, 4)
            contrasts += [dict(model=m, family=f, **c) for c in contrast_check(fd, dr, SEED + k, a.inner)]
        rec["d_new_trial"] = round(rec["candidate_new_trial"] - rec["benchmark_new_trial"], 4)
        rec["d_conditional"] = round(rec["candidate_conditional"] - rec["benchmark_conditional"], 4)
        per_family.append(rec)
        print(f"held out {f} ({k + 1}/{len(held)}, {time.time() - t0:.0f} s): d_new {rec['d_new_trial']:+.3f}, "
              f"d_cond {rec['d_conditional']:+.3f}", flush=True)

    summary = dict(seed=SEED, inner=a.inner, n_boot=N_BOOT, tk=TK, families=len(per_family),
                   revision=os.popen(f"git -C {ROOT} rev-parse --short HEAD").read().strip(), tasks={})
    for task in ("new_trial", "conditional"):
        d = [r[f"d_{task}"] for r in per_family]
        boots = {str(s): bayes_boot(d, s) for s in BOOT_SEEDS}
        infl = {r["family"]: round(bayes_boot([x for x, rr in zip(d, per_family) if rr is not r], BOOT_SEEDS[0])["mean"], 4)
                for r in per_family}
        summary["tasks"][task] = dict(
            total_benchmark=round(sum(r[f"benchmark_{task}"] for r in per_family), 3),
            total_candidate=round(sum(r[f"candidate_{task}"] for r in per_family), 3),
            sum_d=round(float(np.sum(d)), 3), families_candidate_better=int(np.sum(np.array(d) > 0)),
            bootstrap=boots, leave_one_family_mean=infl,
            influence_sign_changes=int(sum((v > 0) != (np.mean(d) > 0) for v in infl.values())),
            promote=bool(all(b["lo"] > 0 for b in boots.values())
                         and all((v > 0) for v in infl.values())))

    # external checks from fits on all primary families
    ext = []
    for m, st in MODELS.items():
        dr = fit_cached(model, f"{m}__all_primary", prim, fams, st, cache)
        gates.append(dict(model=m, held_out="none (all primary)", gate_pass=bool(dr["gate"][0]),
                          divergences=int(dr["gate"][1]), rhat=round(float(dr["gate"][2]), 4),
                          ess_bulk=int(dr["gate"][3]), ess_tail=int(dr["gate"][4])))
        hy = external_family("HYPNO")
        targets = [("HYPNO", hy)] + [(f, family_data(f, allrows)) for f in sorted(s5.SENSITIVITY_FAMILIES)]
        for j, (f, fd) in enumerate(targets):
            new, cond = score_family(fd, dr, a.inner, SEED + 5000 + j)
            cc = contrast_check(fd, dr, SEED + 7000 + j, a.inner)
            ext.append(dict(model=m, family=f, conditional=round(cond, 4), new_trial=round(new, 4),
                            arms=";".join(f"{c['arm_id']} obs {c['observed_lor']:+.2f} pred {c['predicted_lor_median']:+.2f} "
                                          f"[{c['predicted_lor_5']:+.2f},{c['predicted_lor_95']:+.2f}]" for c in cc)))
    summary["external"] = ext
    summary["gates_all_pass"] = all(g["gate_pass"] for g in gates)
    write(os.path.join(a.out, "per_family.csv"), per_family)
    write(os.path.join(a.out, "contrasts.csv"), contrasts)
    write(os.path.join(a.out, "fit_gates.csv"), gates)
    write(os.path.join(a.out, "external.csv"), ext)
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: summary[k] for k in ("tasks", "gates_all_pass")}, indent=1))


def external_family(fam):
    """HYPNO is the R6 benchmark: never fitted; its reference arm is the 66 Gy / 33 control."""
    rows = [x for x in s5_benchmark_rows() if x["family"] == fam]
    old = s5.REF_ARM.get(fam)
    s5.REF_ARM[fam] = "HYPNO_NFX"
    try:
        return family_data(fam, rows)
    finally:
        if old is None:
            del s5.REF_ARM[fam]


def s5_benchmark_rows():
    """Benchmark arms in the same row format as step5_fit.load_rows."""
    fr = {}
    for r in s5.read("arm_fractions.csv"):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    fr = {k: ([d for _, d, _ in sorted(v)], [t for _, _, t in sorted(v)]) for k, v in fr.items()}
    arms = {x["arm_id"]: x for x in s5.read("arms.csv")}
    outs = {o["outcome_id"]: o for o in s5.read("outcomes.csv")}
    rows = []
    for m in s5.read("endpoint_mapping.csv"):
        if m["analysis_set"] != "benchmark":
            continue
        o, x = outs[m["outcome_id"]], arms[m["arm_id"]]
        variants = [(fr[v.split(":")[0]], float(v.split(":")[1])) for v in x["schedule_variants"].split()]
        rows.append(dict(outcome_id=m["outcome_id"], arm_id=m["arm_id"], family=m["family_id"],
                         scale_group=m["scale_group"], N=int(o["denominator"]), y=int(o["events"]),
                         variants=variants, chemo_text=x["planned_systemic_therapy"]))
    return rows


def write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv)
