"""Calibrated DD model fitted conditional on the DI refit, and DI-vs-DD regimen agreement.

For each of N_DI evenly spaced DI posterior draws, alpha, alpha/beta and g are taken
from the draw and the DD onset constant is solved so that DD repopulation equals DI
repopulation at the anchor (35 x 2 Gy in 47 days), as the app does (R/bed_models.R,
dd_calibration). Only lambda_dd, the trial intercepts and the overdispersion are
fitted (fit/stan/dd_cal.stan). A smooth exponential penalty keeps every fitted arm's
repopulation onset at or after ONSET_MIN days.

Each DI draw's fit is saved as soon as it finishes and skipped on a rerun.
Then every regimen's tumor-dose change versus the anchor is computed under DI and
calibrated DD on the same pooled draws, and their agreement is summarized.

Run:  BED_FIT_DATA=data/refit_trials_2026_10_22.csv \
      BED_FIT_OUT=fit/output_refit_2026_10_22 python3 -I fit/fit_dd_calibrated.py
Env:  DDCAL_NDI (default 40), DDCAL_ONSET_MIN (default 7).
"""
import json, os, site, sys, time
if site.getusersitepackages() not in sys.path:
    sys.path.append(site.getusersitepackages())
import numpy as np

SMOOTH_K = 3.0
ANCHOR = (35, 2.0, 47.0)            # fractions, Gy per fraction, overall days
ONSET_SCALE = 1.0                   # days; penalty exp((ONSET_MIN - onset) / ONSET_SCALE)
PRIOR_LOG_LAMBDA = (np.log(0.5), 1.5)  # same lambda prior as the free DD fit

# label, fractions, Gy per fraction, overall time (days, first to last fraction inclusive)
REGIMENS = [
    ("HYPNO 55 Gy/20 (26 d)", 20, 2.75, 26),
    ("HYHOPE-type 46.5 Gy/15 (19 d)", 15, 3.1, 19),
    ("17 x 3.05 Gy (23 d)", 17, 3.05, 23),
    ("60 Gy/25 (33 d)", 25, 2.4, 33),
    ("65 Gy/30 (40 d)", 30, 65 / 30, 40),
    ("66 Gy/33 (45 d)", 33, 2.0, 45),
    ("30 x 2 Gy (40 d), last week dropped", 30, 2.0, 40),
    ("33 x 2 Gy, 6 per week (38 d)", 33, 2.0, 38),
    ("70 Gy/35, 6 per week (40 d)", 35, 2.0, 40),
    ("68 x 1.2 Gy twice daily (47 d)", 68, 1.2, 47),
]


def smax(x, k=SMOOTH_K):
    x = np.asarray(x, float)
    return np.where(k * x > 30, x, np.log1p(np.exp(np.minimum(k * x, 30))) / k)


def inv_smax(y, k=SMOOTH_K):
    y = np.asarray(y, float)
    return np.where(k * y > 30, y, np.log(np.expm1(np.minimum(k * y, 30))) / k)


def bed(n, d, ab):
    return n * d * (1 + d / ab)


def calibrate(alpha, ab, lam_di, Tk, lam_dd, anchor=ANCHOR):
    """Onset at the anchor and c_over_alpha (Gy) such that DD repopulation = DI repopulation there."""
    n0, d0, T0 = anchor
    K0 = bed(n0, d0, ab)
    R0 = lam_di * smax(T0 - Tk)
    onset_ref = T0 - inv_smax(R0 / (lam_dd * -np.expm1(-alpha * K0 / T0)))
    return onset_ref, K0 * onset_ref / T0


def rep_di(lam_di, Tk, T):
    return lam_di * smax(T - Tk)


def rep_dd(alpha, ab, lam_dd, c_over_alpha, n, d, T):
    K = bed(n, d, ab)
    return lam_dd * -np.expm1(-alpha * K / T) * smax(T - c_over_alpha * T / K)


def eff_bed(alpha, ab, g, rep, n, d, T):
    return bed(n, d, ab) - (rep + g * T) / alpha


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def main():
    import pandas as pd
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import fit_refit as F
    ndi = int(os.environ.get("DDCAL_NDI", "40"))
    onset_min = float(os.environ.get("DDCAL_ONSET_MIN", "7"))
    out = os.path.join(F.OUT, f"dd_calibrated_onset{onset_min:g}")
    os.makedirs(out, exist_ok=True)
    rows, trials, base = F.load()
    di = pd.read_csv(os.path.join(F.OUT, "di_primary", "draws_params.csv"))
    pick = np.linspace(0, len(di) - 1, ndi).round().astype(int)
    log("data", F.DATA, "arms", base["N"], "units", base["J"], "DI draws", ndi,
        "onset_min", onset_min, "out", out)
    m = F.model("dd_cal")
    n0, d0, T0 = ANCHOR
    for k, idx in enumerate(pick):
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        if os.path.exists(path):
            continue
        r = di.iloc[idx]
        K0 = bed(n0, d0, r.ab)
        data = dict(base, alpha=r.alpha, ab=r.ab, g=r.g, R0=float(rep_di(r["lambda"], r.Tk, T0)),
                    K0=K0, T0=T0, onset_min=onset_min, onset_scale=ONSET_SCALE,
                    prior_log_lambda_mean=PRIOR_LOG_LAMBDA[0], prior_log_lambda_sd=PRIOR_LOG_LAMBDA[1])
        rng = np.random.default_rng(F.SEED + int(idx))
        inits = [dict(lambda_dd=0.5 * np.exp(rng.normal(0, .3)), c=list(rng.normal(0, .2, base["J"])),
                      sigma_raw=0.2, z=list(rng.normal(0, .1, base["N"]))) for _ in range(4)]
        f = m.sample(data=data, chains=4, parallel_chains=4, iter_warmup=600, iter_sampling=250,
                     seed=F.SEED + int(idx), adapt_delta=0.99, max_treedepth=12, show_progress=False,
                     output_dir=os.path.join(F.RAW, "dd_cal", f"onset{onset_min:g}", str(idx)), inits=inits)
        dr = f.draws_pd()
        s = f.summary()
        rh = float(s.loc["lambda_dd", [c for c in s.columns if c.lower().startswith("r_hat")][0]])
        # The calibration must hold exactly in every draw.
        err = float(np.max(np.abs(dr["rep_anchor"] - data["R0"])))
        assert err < 1e-6 * max(1.0, data["R0"]), f"calibration broken: {err}"
        res = pd.DataFrame(dict(di_draw=idx, lambda_dd=dr["lambda_dd"], onset_ref=dr["onset_ref"],
                                c_over_alpha=dr["c_over_alpha"], onset_min_fitted=dr["onset_min_fitted"],
                                loglik=dr[[f"log_lik[{i}]" for i in range(1, base["N"] + 1)]].sum(axis=1),
                                divergent=dr["divergent__"].to_numpy(), rhat_lambda=rh))
        res.to_csv(path, index=False)
        log(f"DI draw {idx} ({k + 1}/{ndi}): lambda_dd median {res.lambda_dd.median():.3g}, "
            f"onset_ref {res.onset_ref.median():.1f} d, divergences {int(res.divergent.sum())}, rhat {rh:.3f}")
    summarize(di, pick, out, onset_min)


def summarize(di, pick, out, onset_min):
    import pandas as pd
    dd = pd.concat([pd.read_csv(os.path.join(out, f"di_draw_{i:05d}.csv")) for i in pick], ignore_index=True)
    p = dd.merge(di.rename(columns={"draw": "di_draw", "lambda": "lambda_di"}), on="di_draw")
    a, ab, g = p.alpha.to_numpy(), p.ab.to_numpy(), p.g.to_numpy()
    n0, d0, T0 = ANCHOR

    def delta(n, d, T):
        e_di = eff_bed(a, ab, g, rep_di(p.lambda_di.to_numpy(), p.Tk.to_numpy(), T), n, d, T) \
            - eff_bed(a, ab, g, rep_di(p.lambda_di.to_numpy(), p.Tk.to_numpy(), T0), n0, d0, T0)
        e_dd = eff_bed(a, ab, g, rep_dd(a, ab, p.lambda_dd.to_numpy(), p.c_over_alpha.to_numpy(), n, d, T), n, d, T) \
            - eff_bed(a, ab, g, rep_dd(a, ab, p.lambda_dd.to_numpy(), p.c_over_alpha.to_numpy(), n0, d0, T0), n0, d0, T0)
        onset = p.c_over_alpha.to_numpy() * T / bed(n, d, ab)
        return e_di, e_dd, onset

    q = lambda x: dict(median=float(np.median(x)), q10=float(np.quantile(x, .1)), q90=float(np.quantile(x, .9)))
    regs, D_di, D_dd = [], [], []
    for lab, n, d, T in REGIMENS:
        e_di, e_dd, on = delta(n, d, T)
        D_di.append(e_di); D_dd.append(e_dd)
        regs.append(dict(regimen=lab, di=q(e_di), dd=q(e_dd), dd_minus_di=q(e_dd - e_di),
                         p_same_sign=float(np.mean(np.sign(e_di) == np.sign(e_dd))),
                         p_dd_onset_before_min=float(np.mean(on < onset_min)), dd_onset=q(on)))
    D_di, D_dd = np.array(D_di), np.array(D_dd)
    # Pairwise ordering agreement between regimens (anchor included as 0), per draw.
    Z = np.zeros((1, D_di.shape[1]))
    A, B = np.vstack([Z, D_di]), np.vstack([Z, D_dd])
    iu = np.triu_indices(A.shape[0], 1)
    agree = (np.sign(A[:, None, :] - A[None, :, :]) == np.sign(B[:, None, :] - B[None, :, :]))[iu[0], iu[1], :]
    lam = p.lambda_dd.to_numpy()
    res = dict(onset_min=onset_min, onset_scale=ONSET_SCALE, n_di_draws=len(pick), n_pooled=len(p),
               divergences=int(dd.divergent.sum()), max_rhat_lambda=float(dd.rhat_lambda.max()),
               lambda_dd=dict(q(lam), q2_5=float(np.quantile(lam, .025)), q97_5=float(np.quantile(lam, .975))),
               lambda_dd_prior_95=[float(np.exp(PRIOR_LOG_LAMBDA[0] - 1.96 * PRIOR_LOG_LAMBDA[1])),
                                   float(np.exp(PRIOR_LOG_LAMBDA[0] + 1.96 * PRIOR_LOG_LAMBDA[1]))],
               onset_ref=q(p.onset_ref), min_fitted_onset=q(p.onset_min_fitted),
               pairwise_order_agreement=float(agree.mean()), regimens=regs)
    json.dump(res, open(os.path.join(out, "summary.json"), "w"), indent=1)
    log("lambda_dd", res["lambda_dd"], "onset_ref", res["onset_ref"], "divergences", res["divergences"],
        "pairwise order agreement", round(res["pairwise_order_agreement"], 3))
    for r in regs:
        log(f"{r['regimen']:38s} DI {r['di']['median']:+5.1f} ({r['di']['q10']:+.1f},{r['di']['q90']:+.1f})  "
            f"DD {r['dd']['median']:+5.1f} ({r['dd']['q10']:+.1f},{r['dd']['q90']:+.1f})  "
            f"same sign {r['p_same_sign']:.2f}  DD onset<{onset_min:g}d {r['p_dd_onset_before_min']:.2f}")


if __name__ == "__main__":
    main()
