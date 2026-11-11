"""Plan step 4: identification and prior predictive checks, before any outcome-based fit.

Usage:
    python3 tools/early_tox/step4_identification.py [--reps N] [--out DIR]

Uses only the trial DESIGN (arms, schedules, patient numbers) from data/early_tox/. Outcome counts
enter only one place: the observed reference-arm rates set plausible baselines for the simulations.
No posterior for the real outcomes is computed here (that is step 5).

Model (plan step 5, arm-only candidate), primary analysis set:
    y_ij ~ Binomial(N_ij, p_ij)
    logit p_ij = eta_j + kappa * (E_ij - E_j,ref)
    E = early BED from tools/early_tox/evaluator.py (repair half-time 3.2 h, elapsed clock)
    theta = (alpha/beta, r, Tk); Tk is fixed within a scenario
Priors (decision D3, analyst-chosen): eta_j ~ N(m0, tau0^2), m0 ~ N(0, 2.5^2), tau0 ~ HalfNormal(1);
    L* = 10 kappa ~ HalfNormal(1); alpha/beta ~ lognormal(log 10, 0.354); r ~ lognormal(log 0.8, 0.5).
For recovery fits the eta_j are fixed effects with a N(0, 2.5^2) prior (a simplification that can
only widen, not narrow, the recovered uncertainty of the shared parameters).

Outputs in OUT (default docs/development/2026_11_10_early_toxicity/step4/):
    contrasts.csv, jacobian_summary.csv, tk_profile.csv, fisher_se.csv,
    recovery_replicates.csv (one row per replicate, written as it completes),
    prior_predictive.csv, summary.json
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
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_spec = importlib.util.spec_from_file_location("evaluator", os.path.join(ROOT, "tools", "early_tox", "evaluator.py"))
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)

SEED = 20261110
HALFTIME_H = 3.2
TRUTH = dict(ab=10.0, r=ev.fowler_rate(), Tk=7.0, kappa=0.1)   # L* = 1 per 10 BED units
PRIOR = dict(eta_sd=2.5, ab_mu=math.log(10), ab_sd=0.354, r_mu=math.log(0.8), r_sd=0.5, Lstar_sd=1.0)

# Reference (comparator) arm of each primary family: the conventional or control arm.
REF_ARM = {
    "GORTEC9902": "GORTEC9902_CRT", "RTOG9003": "RTOG9003_SF", "RTOG0129": "RTOG0129_SFX",
    "EORTC22791": "EORTC22791_CF", "EORTC22851": "EORTC22851_CF", "PMH_HF": "PMH_CF",
    "Ghoshal": "Ghoshal_CF", "CAIR1": "CAIR1_CF", "Vienna": "Vienna_CF", "ORO9301": "Olmi_CF",
    "Bartelink2002": "Bartelink_CF", "Tandon": "Tandon_SIB", "Sikdar": "Sikdar_SEQ",
    "Shivananjappa": "Shivananjappa_CONV", "Mani": "Mani_SEQ",
}


def read(name):
    with open(os.path.join(ROOT, "data", "early_tox", name), newline="") as f:
        return list(csv.DictReader(f))


def load_design():
    fr = {}
    for r in read("arm_fractions.csv"):
        fr.setdefault(r["schedule_id"], []).append((int(r["fraction"]), float(r["dose_gy"]), float(r["time_h"])))
    fr = {k: ([d for _, d, _ in sorted(v)], [t for _, _, t in sorted(v)]) for k, v in fr.items()}
    arms = {a["arm_id"]: a for a in read("arms.csv")}
    out = {o["outcome_id"]: o for o in read("outcomes.csv")}
    rows = []
    for m in read("endpoint_mapping.csv"):
        if m["in_primary"] != "1":
            continue
        o, a = out[m["outcome_id"]], arms[m["arm_id"]]
        variants = []
        for v in a["schedule_variants"].split():
            sid, w = v.split(":")
            variants.append((fr[sid], float(w)))
        rows.append(dict(outcome_id=m["outcome_id"], arm_id=m["arm_id"], family=m["family_id"],
                         window=m["window"], N=int(o["denominator"]), y=int(o["events"]), variants=variants,
                         chemo=a["planned_systemic_therapy"]))
    fams = sorted({r["family"] for r in rows})
    missing = [f for f in fams if f not in REF_ARM or REF_ARM[f] not in {r["arm_id"] for r in rows}]
    if missing:
        raise ValueError(f"no reference arm for {missing}")
    return rows, fams


def arm_terms(row, Tk):
    """Weighted D, Q (3.2 h repair), T_eff = max(0, T - Tk) and elapsed T over schedule variants."""
    tot = sum(w for _, w in row["variants"])
    D = Q = Te = T = 0.0
    for (doses, times), w in row["variants"]:
        D += w * sum(doses)
        Q += w * ev.repair_q(doses, times, HALFTIME_H)
        T += w * ev.elapsed_days(times)
        Te += w * max(0.0, ev.elapsed_days(times) - Tk)
    return D / tot, Q / tot, Te / tot, T / tot


def exposures(rows, ab, r, Tk):
    terms = np.array([arm_terms(x, Tk) for x in rows])
    E = terms[:, 0] + terms[:, 1] / ab - r * terms[:, 2]
    return E, terms


def contrasts(rows, fams, ab, r, Tk):
    E, terms = exposures(rows, ab, r, Tk)
    ref = {f: [i for i, x in enumerate(rows) if x["arm_id"] == REF_ARM[f]][0] for f in fams}
    dE = np.array([E[i] - E[ref[x["family"]]] for i, x in enumerate(rows)])
    dterms = np.array([terms[i] - terms[ref[x["family"]]] for i, x in enumerate(rows)])
    return dE, dterms, ref


def jacobian(rows, fams, ab, r, Tk, h=1e-4):
    """d(Delta E)/d(log ab, log r, Tk) by central differences, contrasts only (reference rows dropped)."""
    def f(x):
        return contrasts(rows, fams, math.exp(x[0]), math.exp(x[1]), x[2])[0]
    x0 = np.array([math.log(ab), math.log(r), Tk])
    J = np.column_stack([(f(x0 + e) - f(x0 - e)) / (2 * h) for e in np.eye(3) * h])
    _, _, ref = contrasts(rows, fams, ab, r, Tk)
    keep = [i for i in range(len(rows)) if i not in ref.values()]
    return J[keep], keep


def neg_log_post(phi, rows, fams, Tk, y, fam_idx):
    k = len(fams)
    eta, logk, logab, logr = phi[:k], phi[k], phi[k + 1], phi[k + 2]
    kappa, ab, r = math.exp(logk), math.exp(logab), math.exp(logr)
    dE, _, _ = contrasts(rows, fams, ab, r, Tk)
    lin = eta[fam_idx] + kappa * dE
    N = np.array([x["N"] for x in rows])
    ll = np.sum(y * lin - N * np.logaddexp(0, lin))
    lp = -0.5 * np.sum((eta / PRIOR["eta_sd"]) ** 2)
    Lstar = 10 * kappa                                          # HalfNormal(1) on L*, sampled on log scale
    lp += -0.5 * (Lstar / PRIOR["Lstar_sd"]) ** 2 + logk
    lp += -0.5 * ((logab - PRIOR["ab_mu"]) / PRIOR["ab_sd"]) ** 2
    lp += -0.5 * ((logr - PRIOR["r_mu"]) / PRIOR["r_sd"]) ** 2
    return -(ll + lp)


def numeric_hessian(fun, x, h=1e-3):
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei, ej = np.eye(n)[i] * h, np.eye(n)[j] * h
            H[i, j] = H[j, i] = (fun(x + ei + ej) - fun(x + ei - ej) - fun(x - ei + ej) + fun(x - ei - ej)) / (4 * h * h)
    return H


def fit_map(rows, fams, Tk, y, fam_idx, start):
    fun = lambda p: neg_log_post(p, rows, fams, Tk, y, fam_idx)
    res = minimize(fun, start, method="L-BFGS-B")
    H = numeric_hessian(fun, res.x)
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
        ok = bool(np.all(np.linalg.eigvalsh(H) > 0))
    except np.linalg.LinAlgError:
        se, ok = np.full(len(start), np.nan), False
    return res, se, ok


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step4"))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    rows, fams = load_design()
    fam_idx = np.array([fams.index(x["family"]) for x in rows])
    N = np.array([x["N"] for x in rows])
    t0 = time.time()

    # 1. contrasts and structural identification at the Fowler values
    dE, dterms, ref = contrasts(rows, fams, TRUTH["ab"], TRUTH["r"], TRUTH["Tk"])
    E, terms = exposures(rows, TRUTH["ab"], TRUTH["r"], TRUTH["Tk"])
    with open(os.path.join(args.out, "contrasts.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["outcome_id", "arm_id", "family", "is_reference", "N", "D_gy", "Q_gy2", "elapsed_days",
                    "E_fowler_rep3p2", "dE_vs_reference", "chemo_planned"])
        for i, x in enumerate(rows):
            w.writerow([x["outcome_id"], x["arm_id"], x["family"], int(i in ref.values()), x["N"],
                        round(terms[i, 0], 4), round(terms[i, 1], 4), round(terms[i, 3], 4),
                        round(E[i], 4), round(dE[i], 4), x["chemo"]])
    keep = [i for i in range(len(rows)) if i not in ref.values()]
    X = np.column_stack([dterms[keep, 0], dterms[keep, 1], dterms[keep, 2]])   # dD, dQ, dT_eff
    sv = np.linalg.svd(X / np.linalg.norm(X, axis=0), compute_uv=False)
    J, _ = jacobian(rows, fams, TRUTH["ab"], TRUTH["r"], TRUTH["Tk"])
    jac = dict(n_contrasts=len(keep), design_rank=int(np.linalg.matrix_rank(X)),
               design_condition_number=float(sv[0] / sv[-1]),
               dE_dlogab_range=[float(J[:, 0].min()), float(J[:, 0].max())],
               dE_dlogr_range=[float(J[:, 1].min()), float(J[:, 1].max())],
               dE_dTk_max_abs=float(np.abs(J[:, 2]).max()),
               corr_dD_dQ=float(np.corrcoef(X[:, 0], X[:, 1])[0, 1]),
               corr_dD_dT=float(np.corrcoef(X[:, 0], X[:, 2])[0, 1]),
               corr_dQ_dT=float(np.corrcoef(X[:, 1], X[:, 2])[0, 1]))
    with open(os.path.join(args.out, "jacobian_summary.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        for k, v in jac.items():
            w.writerow([k, json.dumps(v)])

    # 2. Tk profile: how much the contrasts change as Tk moves (zero if every arm ends after Tk)
    with open(os.path.join(args.out, "tk_profile.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["Tk_days", "max_abs_change_in_dE_vs_Tk7", "n_arms_ending_before_Tk"])
        base = dE[keep]
        T_el = terms[:, 3]
        for Tk in range(0, 41, 2):
            d2, _, _ = contrasts(rows, fams, TRUTH["ab"], TRUTH["r"], float(Tk))
            w.writerow([Tk, round(float(np.abs(d2[keep] - base).max()), 4), int(np.sum(T_el < Tk))])

    # 3. expected posterior SDs (Laplace at the noise-free expected counts) for three slope sizes
    obs_ref = {f: rows[ref[f]]["y"] / rows[ref[f]]["N"] for f in fams}
    eta_true = np.array([math.log(min(max(obs_ref[f], 0.03), 0.97) / (1 - min(max(obs_ref[f], 0.03), 0.97))) for f in fams])
    names = ["log_kappa", "log_ab", "log_r"]
    prior_sd = {"log_kappa": float("nan"), "log_ab": PRIOR["ab_sd"], "log_r": PRIOR["r_sd"]}
    fisher = []
    for Lstar in (0.5, 1.0, 2.0):
        kap = Lstar / 10
        phi_k = np.concatenate([eta_true, [math.log(kap), math.log(TRUTH["ab"]), math.log(TRUTH["r"])]])
        p_k = 1 / (1 + np.exp(-(eta_true[fam_idx] + kap * dE)))
        Hk = numeric_hessian(lambda p: neg_log_post(p, rows, fams, TRUTH["Tk"], N * p_k, fam_idx), phi_k)
        cov = np.linalg.inv(Hk)
        for j, nm in enumerate(names):
            se = float(math.sqrt(cov[len(fams) + j, len(fams) + j]))
            fisher.append(dict(Lstar=Lstar, parameter=nm, expected_posterior_sd=round(se, 4), prior_sd=prior_sd[nm],
                               ratio_to_prior=round(se / prior_sd[nm], 3) if prior_sd[nm] == prior_sd[nm] else ""))
    p_true = 1 / (1 + np.exp(-(eta_true[fam_idx] + TRUTH["kappa"] * dE)))
    with open(os.path.join(args.out, "fisher_se.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fisher[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(fisher)

    # 4. recovery simulations on the real design (one row per replicate, resumable)
    rec_path = os.path.join(args.out, "recovery_replicates.csv")
    done = set()
    if os.path.exists(rec_path):
        with open(rec_path, newline="") as f:
            done = {int(r["rep"]) for r in csv.DictReader(f)}
    fieldnames = ["rep", "converged", "hessian_pd", "kappa_hat", "kappa_se_log", "ab_hat", "ab_se_log",
                  "r_hat", "r_se_log", "covers_kappa", "covers_ab", "covers_r"]
    new_file = not os.path.exists(rec_path)
    with open(rec_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        if new_file:
            w.writeheader()
        for rep in range(args.reps):
            if rep in done:
                continue
            rng = np.random.default_rng([SEED, rep])
            y = rng.binomial(N, p_true)
            start = np.concatenate([np.zeros(len(fams)), [math.log(0.05), math.log(10), math.log(0.8)]])
            res, se, ok = fit_map(rows, fams, TRUTH["Tk"], y, fam_idx, start)
            k = len(fams)
            est = res.x[k:k + 3]
            truth = [math.log(TRUTH["kappa"]), math.log(TRUTH["ab"]), math.log(TRUTH["r"])]
            cov_ = [int(abs(est[j] - truth[j]) <= 1.96 * se[k + j]) if np.isfinite(se[k + j]) else "" for j in range(3)]
            w.writerow(dict(rep=rep, converged=int(res.success), hessian_pd=int(ok),
                            kappa_hat=round(math.exp(est[0]), 5), kappa_se_log=round(float(se[k]), 4),
                            ab_hat=round(math.exp(est[1]), 4), ab_se_log=round(float(se[k + 1]), 4),
                            r_hat=round(math.exp(est[2]), 4), r_se_log=round(float(se[k + 2]), 4),
                            covers_kappa=cov_[0], covers_ab=cov_[1], covers_r=cov_[2]))
            f.flush()
            print(f"recovery {rep + 1}/{args.reps}", flush=True)

    # 5. prior predictive (D3 priors), within-trial risk differences
    rng = np.random.default_rng([SEED, 999999])
    M = 2000
    pp = []
    for b in range(M):
        m0, tau0 = rng.normal(0, 2.5), abs(rng.normal(0, 1))
        eta = rng.normal(m0, tau0, len(fams))
        kappa = abs(rng.normal(0, 1)) / 10
        ab, r = math.exp(rng.normal(PRIOR["ab_mu"], PRIOR["ab_sd"])), math.exp(rng.normal(PRIOR["r_mu"], PRIOR["r_sd"]))
        d, _, _ = contrasts(rows, fams, ab, r, TRUTH["Tk"])
        p = 1 / (1 + np.exp(-(eta[fam_idx] + kappa * d)))
        rd = p[keep] - np.array([p[ref[rows[i]["family"]]] for i in keep])
        pp.append((float(np.max(np.abs(rd))), float(np.median(np.abs(rd))), float(np.mean(p)), math.exp(10 * kappa)))
    pp = np.array(pp)
    with open(os.path.join(args.out, "prior_predictive.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["quantity", "q05", "q50", "q95"])
        for j, nm in enumerate(["max_abs_within_trial_risk_difference", "median_abs_within_trial_risk_difference",
                                "mean_arm_risk", "odds_ratio_per_10_BED"]):
            q = np.quantile(pp[:, j], [0.05, 0.5, 0.95])
            w.writerow([nm] + [round(float(v), 4) for v in q])

    # summary
    with open(rec_path, newline="") as f:
        rec = list(csv.DictReader(f))
    def frac(col):
        v = [int(r[col]) for r in rec if r[col] != ""]
        return round(sum(v) / len(v), 3) if v else None
    summary = dict(seed=SEED, reps=len(rec), truth=TRUTH, halftime_h=HALFTIME_H,
                   n_primary_arms=len(rows), n_families=len(fams), jacobian=jac,
                   recovery=dict(converged=frac("converged"), hessian_pd=frac("hessian_pd"),
                                 coverage_kappa=frac("covers_kappa"), coverage_ab=frac("covers_ab"),
                                 coverage_r=frac("covers_r"),
                                 ab_hat_q05_q50_q95=[round(float(q), 3) for q in np.quantile([float(r["ab_hat"]) for r in rec], [.05, .5, .95])],
                                 r_hat_q05_q50_q95=[round(float(q), 4) for q in np.quantile([float(r["r_hat"]) for r in rec], [.05, .5, .95])],
                                 kappa_hat_q05_q50_q95=[round(float(q), 4) for q in np.quantile([float(r["kappa_hat"]) for r in rec], [.05, .5, .95])]),
                   fisher=fisher, seconds=round(time.time() - t0, 1),
                   versions=dict(python=sys.version.split()[0], numpy=np.__version__,
                                 scipy=__import__("scipy").__version__))
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
