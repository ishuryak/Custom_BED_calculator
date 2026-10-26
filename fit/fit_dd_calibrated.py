"""Calibrated DD model fitted conditional on the DI refit, and DI-vs-DD regimen agreement.

For each of N_DI evenly spaced DI posterior draws, alpha, alpha/beta and g are taken from the draw
and the DD onset constant is solved so that DD repopulation equals DI repopulation at the anchor
(35 x 2 Gy in 47 days), as the app does (R/bed_models.R, dd_calibration). Only lambda_dd, the
trial intercepts and the overdispersion are fitted (fit/stan/dd_cal.stan). This is a modular
(cut) analysis, p_DI(theta | y) p_DD(phi | theta, y): DD does not reweight DI, and each DD draw is
paired with the DI draw it was fitted on.

A smooth exponential penalty -sum exp((ONSET_MIN - onset_i) / 1 day) over the fitted arms
discourages onsets before ONSET_MIN days. It is a soft prior preference (Igor Shuryak,
2026-10-24), not a bound: onsets below ONSET_MIN are reported separately for fitted arms and
for the evaluation regimens.

Run identity: the output folder holds identity.json binding the parent DI fit (its identity and
draws file hash), the data, dd_cal.stan, this file and all settings. Each DI draw's fit is saved
as soon as it finishes, carries the DI values it used, and is reused only if those values match
the parent file. A fit that fails the diagnostic gate (fit/gate.py: any divergence or tree-depth
hit, R-hat > 1.01, bulk or tail ESS < 400, E-BFMI < 0.3) is rerun with longer chains, up to two
times. Before 2026-10-26 this loop omitted ESS (review R01); fits made then are reclassified and
repaired by fit/repair_dd_ess.py.

Summaries weight every outer DI draw equally (1/N_DI), whatever the number of rows its fit
produced, so a fit rerun with longer chains does not gain weight (review R06).

Run:  BED_FIT_DATA=data/X.csv BED_FIT_OUT=fit/output_X python3 -I fit/fit_dd_calibrated.py
Env:  DDCAL_NDI (default 200), DDCAL_ONSET_MIN (default 7), DDCAL_DI_VARIANT (default primary).
"""
import json, os, site, sys, time, hashlib
if site.getusersitepackages() not in sys.path:
    sys.path.append(site.getusersitepackages())
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_refit as F

SMOOTH_K = 3.0
ANCHOR = F.ANCHOR
REGIMENS = F.REGIMENS
ONSET_SCALE = 1.0                   # days; penalty exp((ONSET_MIN - onset) / ONSET_SCALE)
PRIOR_LOG_LAMBDA = (np.log(0.5), 1.5)  # same lambda prior as the free DD fit
INNER = dict(chains=4, warmup=600, samples=250)
if os.environ.get("DDCAL_INNER"):   # "warmup,samples", smoke tests only (recorded in the identity)
    INNER["warmup"], INNER["samples"] = (int(x) for x in os.environ["DDCAL_INNER"].split(","))
MAX_ATTEMPTS = 3


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


def penalty_adjusted_prior(r, base, onset_min, n=200_000, seed=1):
    """lambda_dd prior alone versus prior x onset penalty, at one DI draw r (no likelihood)."""
    rng = np.random.default_rng(seed)
    lam = np.exp(rng.normal(PRIOR_LOG_LAMBDA[0], PRIOR_LOG_LAMBDA[1], n))
    _, co = calibrate(r.alpha, r.ab, r["lambda"], r.Tk, lam)
    K = np.array(base["dsum"]) + np.array(base["d2sum"]) / r.ab
    on = co[:, None] * np.array(base["T"])[None, :] / K[None, :]
    with np.errstate(over="ignore"):   # exp overflow -> weight 0, as intended
        w = np.exp(-np.exp((onset_min - on) / ONSET_SCALE).sum(axis=1))
    o = np.argsort(lam); cw = np.cumsum(w[o]) / w.sum()
    wq = lambda q: float(lam[o][np.searchsorted(cw, q)])
    return dict(nominal_95=[float(np.quantile(lam, .025)), float(np.quantile(lam, .975))],
                nominal_median=float(np.median(lam)),
                penalized_95=[wq(.025), wq(.975)], penalized_median=wq(.5))


def main():
    import pandas as pd
    ndi = int(os.environ.get("DDCAL_NDI", "200"))
    onset_min = float(os.environ.get("DDCAL_ONSET_MIN", "7"))
    di_variant = os.environ.get("DDCAL_DI_VARIANT", "primary")
    out = os.path.join(F.OUT, f"dd_calibrated_onset{onset_min:g}_n{ndi}" + ("" if di_variant == "primary" else f"_{di_variant}"))
    os.makedirs(out, exist_ok=True)
    rows, trials, base = F.load()
    di_dir = os.path.join(F.OUT, f"di_{di_variant}")
    di_ident = json.load(open(os.path.join(di_dir, "identity.json")))
    di_file = os.path.join(di_dir, "draws_params.csv")
    di = pd.read_csv(di_file)
    pick = np.linspace(0, len(di) - 1, ndi).round().astype(int)
    ident = dict(parent_di_identity=di_ident["identity_sha256"], parent_di_draws_sha256=F.sha_file(di_file),
                 data_file_sha256=F.sha_file(F.DATA), stan_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "dd_cal.stan")),
                 common_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "common_functions.stanfunctions")),
                 runner_sha256=F.sha_file(os.path.abspath(__file__)), fit_refit_sha256=F.sha_file(F.__file__),
                 di_variant=di_variant, n_di=ndi, picks_sha256=hashlib.sha256(pick.tobytes()).hexdigest(), onset_min=onset_min,
                 onset_scale=ONSET_SCALE, prior_log_lambda=list(PRIOR_LOG_LAMBDA), prior_c_sd=F.PRIOR_C_SD,
                 inner=INNER, max_attempts=MAX_ATTEMPTS, seed=F.SEED, software=F.versions())
    ident["identity_sha256"] = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()
    idp = os.path.join(out, "identity.json")
    if os.path.exists(idp) and json.load(open(idp)).get("identity_sha256") != ident["identity_sha256"]:
        sys.exit(f"REFUSED: {out} was started with a different identity. Use a new output folder.")
    json.dump(ident, open(idp, "w"), indent=1)
    log("data", F.DATA, "arms", base["N"], "units", base["J"], "DI draws", ndi, "onset_min", onset_min, "out", out)
    m = F.model("dd_cal")
    n0, d0, T0 = ANCHOR
    R = len(REGIMENS)
    for k, idx in enumerate(pick):
        r = di.iloc[idx]
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        if os.path.exists(path):
            old = pd.read_csv(path, nrows=1)
            same = all(np.isclose(old[c].iloc[0], r[c], rtol=0, atol=1e-12) for c in ("alpha", "ab", "g", "lambda", "Tk"))
            if not same:
                sys.exit(f"REFUSED: {path} was fitted on other DI values than draw {idx} of the parent file.")
            continue
        K0 = bed(n0, d0, r.ab)
        data = dict(base, alpha=r.alpha, ab=r.ab, g=r.g, R0=float(rep_di(r["lambda"], r.Tk, T0)),
                    K0=K0, T0=T0, onset_min=onset_min, onset_scale=ONSET_SCALE,
                    prior_log_lambda_mean=PRIOR_LOG_LAMBDA[0], prior_log_lambda_sd=PRIOR_LOG_LAMBDA[1])
        for attempt in range(1, MAX_ATTEMPTS + 1):
            mult = 2 ** (attempt - 1)
            seed = F.SEED + int(idx) + 1_000_000 * (attempt - 1)
            rng = np.random.default_rng(seed)
            inits = [dict(lambda_dd=0.5 * np.exp(rng.normal(0, .3)), c=list(rng.normal(0, .2, base["J"])),
                          sigma_raw=0.2, z=list(rng.normal(0, .1, base["N"]))) for _ in range(INNER["chains"])]
            f = m.sample(data=data, chains=INNER["chains"], parallel_chains=INNER["chains"],
                         iter_warmup=INNER["warmup"] * mult, iter_sampling=INNER["samples"] * mult,
                         seed=seed, adapt_delta=0.99, max_treedepth=12, show_progress=False,
                         output_dir=os.path.join(F.RAW, "dd_cal", os.path.basename(F.OUT), os.path.basename(out), f"{idx}_{attempt}"),
                         inits=inits)
            dg = F.diagnostics(f, ["lambda_dd", "sigma_raw"])
            ok = F.G.passes(dg)
            if ok:
                break
            log(f"  DI draw {idx} attempt {attempt} failed the gate: {'; '.join(F.G.failures(dg))}")
        dr = f.draws_pd()
        # The calibration must hold exactly in every draw.
        err = float(np.max(np.abs(dr["rep_anchor"] - data["R0"])))
        assert err < 1e-6 * max(1.0, data["R0"]), f"calibration broken: {err}"
        res = pd.DataFrame(dict(di_draw=idx, alpha=r.alpha, ab=r.ab, g=r.g, **{"lambda": r["lambda"]}, Tk=r.Tk,
                                chain=dr["chain__"], iter=dr["iter__"], lambda_dd=dr["lambda_dd"],
                                onset_ref=dr["onset_ref"], c_over_alpha=dr["c_over_alpha"],
                                onset_min_fitted=dr["onset_min_fitted"],
                                loglik=dr[[f"log_lik[{i}]" for i in range(1, base["N"] + 1)]].sum(axis=1),
                                divergent=dr["divergent__"].to_numpy(), attempt=attempt, gate_pass=ok,
                                max_rhat=dg["max_rhat"], min_ess_bulk=dg["min_ess_bulk"],
                                min_ess_tail=dg["min_ess_tail"], min_ebfmi=dg["min_ebfmi"],
                                treedepth_hits=dg["treedepth_hits"]))
        for j in range(1, R + 1):
            res[f"dd_delta_{j}"] = dr[f"delta_eval[{j}]"].to_numpy()
            res[f"dd_onset_{j}"] = dr[f"onset_eval[{j}]"].to_numpy()
        res.to_csv(path, index=False, float_format="%.10g")
        log(f"DI draw {idx} ({k + 1}/{ndi}): lambda_dd median {res.lambda_dd.median():.3g}, onset_ref "
            f"{res.onset_ref.median():.1f} d, attempt {attempt}, gate {'PASS' if ok else 'FAIL'}, rhat {dg['max_rhat']:.4f}")
    summarize(di, pick, out, onset_min, base)


def outer_se(x, groups):
    """Standard error of a pooled mean when rows cluster within outer DI draws."""
    g = np.unique(groups, return_inverse=True)[1]
    m = np.bincount(g, weights=x) / np.bincount(g)
    return float(np.std(m, ddof=1) / np.sqrt(len(m)))


def outer_weights(groups):
    """Weight 1 / (N_outer x rows of that outer draw): each outer DI draw sums to 1/N_outer."""
    g = np.unique(groups, return_inverse=True)[1]
    cnt = np.bincount(g)
    return 1.0 / (len(cnt) * cnt[g])


def wquantile(x, w, p):
    """Weighted quantile: smallest x whose weighted cumulative share reaches p (inverse of the
    weighted empirical CDF). Exactly unchanged when one outer draw's rows are duplicated."""
    x = np.asarray(x, float); o = np.argsort(x, kind="stable"); x, w = x[o], np.asarray(w, float)[o]
    c = np.cumsum(w) / w.sum()
    return float(x[min(np.searchsorted(c, p - 1e-12), len(x) - 1)])


def fit_diagnostics(rows):
    """One fit's diagnostics from its saved rows: the per-fit columns plus the divergence count."""
    r = rows.iloc[0]
    return dict(divergences=int(rows.divergent.astype(bool).sum()), treedepth_hits=r.get("treedepth_hits"),
                max_rhat=r.get("max_rhat"), min_ess_bulk=r.get("min_ess_bulk"), min_ess_tail=r.get("min_ess_tail"),
                min_ebfmi=r.get("min_ebfmi"))


def per_fit_gate(dd):
    """Gate verdict per outer fit, recomputed from saved diagnostics (never from a stored boolean)."""
    cols = ["divergent", "treedepth_hits", "max_rhat", "min_ess_bulk", "min_ess_tail", "min_ebfmi"]
    return dd.groupby("di_draw")[cols].apply(lambda x: F.G.passes(fit_diagnostics(x)))


def summarize(di, pick, out, onset_min, base):
    import pandas as pd
    dd = pd.concat([pd.read_csv(os.path.join(out, f"di_draw_{i:05d}.csv")) for i in pick], ignore_index=True)
    R = len(REGIMENS)
    dic = di.set_index("draw")
    D_di = np.vstack([dic.loc[dd.di_draw, f"delta_eval[{j}]"].to_numpy() for j in range(1, R + 1)])
    D_dd = np.vstack([dd[f"dd_delta_{j}"].to_numpy() for j in range(1, R + 1)])
    # Cross-check the Stan contrasts against the Python formulas on the same draws.
    a, ab, g = dd.alpha.to_numpy(), dd.ab.to_numpy(), dd.g.to_numpy()
    n0, d0, T0 = ANCHOR
    for j, (lab, n, d, T) in enumerate(REGIMENS):
        e_di = eff_bed(a, ab, g, rep_di(dd["lambda"].to_numpy(), dd.Tk.to_numpy(), T), n, d, T) \
            - eff_bed(a, ab, g, rep_di(dd["lambda"].to_numpy(), dd.Tk.to_numpy(), T0), n0, d0, T0)
        e_dd = eff_bed(a, ab, g, rep_dd(a, ab, dd.lambda_dd.to_numpy(), dd.c_over_alpha.to_numpy(), n, d, T), n, d, T) \
            - eff_bed(a, ab, g, rep_dd(a, ab, dd.lambda_dd.to_numpy(), dd.c_over_alpha.to_numpy(), n0, d0, T0), n0, d0, T0)
        assert np.max(np.abs(e_di - D_di[j])) < 1e-4 and np.max(np.abs(e_dd - D_dd[j])) < 1e-4, lab
    grp = dd.di_draw.to_numpy()
    w = outer_weights(grp)
    q = lambda x: dict(median=wquantile(x, w, .5), q10=wquantile(x, w, .1), q90=wquantile(x, w, .9))
    wmean = lambda x: float(np.sum(w * x) / np.sum(w))
    regs = []
    for j, (lab, n, d, T) in enumerate(REGIMENS):
        same = (np.sign(D_di[j]) == np.sign(D_dd[j])).astype(float)
        on = dd[f"dd_onset_{j + 1}"].to_numpy()
        regs.append(dict(regimen=lab, di=q(D_di[j]), dd=q(D_dd[j]), dd_minus_di=q(D_dd[j] - D_di[j]),
                         p_same_sign=wmean(same), p_same_sign_outer_se=outer_se(same, grp),
                         n_dd_onset_before_min=int((on < onset_min).sum()), share_dd_onset_before_min=wmean((on < onset_min).astype(float)), dd_onset=q(on)))
    Z = np.zeros((1, D_di.shape[1]))
    A, B = np.vstack([Z, D_di]), np.vstack([Z, D_dd])
    iu = np.triu_indices(A.shape[0], 1)
    agree = (np.sign(A[:, None, :] - A[None, :, :]) == np.sign(B[:, None, :] - B[None, :, :]))[iu[0], iu[1], :]
    agree_draw = agree.mean(axis=0)   # per pooled row, averaged over regimen pairs
    lam = dd.lambda_dd.to_numpy()
    per = dd.groupby("di_draw").first()
    gate_ok = per_fit_gate(dd).reindex(per.index)
    mid = di.iloc[[int(np.argsort(di.alpha.to_numpy())[len(di) // 2])]].iloc[0]
    res = dict(analysis="calibrated DD conditional on DI (modular); DD draws paired with their DI draw",
               onset_min=onset_min, onset_scale=ONSET_SCALE, onset_penalty="soft prior preference, not a bound",
               n_di_draws=len(pick), n_pooled=len(dd), divergences=int(dd.divergent.sum()),
               weighting="each outer DI draw weight 1/N_DI, split equally over its rows (review R06)",
               gate_policy=F.G.POLICY,
               outer_fits_failing_gate=int((~gate_ok).sum()),
               outer_fits_failing_gate_ids=[int(i) for i in per.index[~gate_ok.to_numpy()]],
               outer_fits_retried=int((per.attempt > 1).sum()),
               max_rhat=float(per.max_rhat.max()), min_ess_bulk_per_fit=float(per.min_ess_bulk.min()),
               min_ess_tail_per_fit=float(per.min_ess_tail.min()), min_ebfmi=float(per.min_ebfmi.min()),
               treedepth_hits=int(per.treedepth_hits.sum()),
               lambda_dd=dict(q(lam), q2_5=wquantile(lam, w, .025), q97_5=wquantile(lam, w, .975)),
               lambda_dd_prior=penalty_adjusted_prior(mid, base, onset_min),
               onset_ref=q(dd.onset_ref),
               fitted_onset_below_min=dict(n_draws=int((dd.onset_min_fitted < onset_min).sum()),
                                           share=wmean((dd.onset_min_fitted < onset_min).to_numpy().astype(float)),
                                           n_outer_fits=int(dd[dd.onset_min_fitted < onset_min].di_draw.nunique()),
                                           minimum=float(dd.onset_min_fitted.min())),
               min_fitted_onset=q(dd.onset_min_fitted),
               pairwise_order_agreement=wmean(agree_draw),
               pairwise_order_agreement_outer_se=outer_se(agree_draw, grp),
               n_regimen_pairs=int(len(iu[0])),
               regimens=regs)
    json.dump(res, open(os.path.join(out, "summary.json"), "w"), indent=1)
    log("lambda_dd", res["lambda_dd"], "divergences", res["divergences"], "failing fits", res["outer_fits_failing_gate"],
        "order agreement", round(res["pairwise_order_agreement"], 4), "+-", round(res["pairwise_order_agreement_outer_se"], 4))
    for r in regs:
        log(f"{r['regimen']:38s} DI {r['di']['median']:+5.1f} ({r['di']['q10']:+.1f},{r['di']['q90']:+.1f})  "
            f"DD {r['dd']['median']:+5.1f} ({r['dd']['q10']:+.1f},{r['dd']['q90']:+.1f})  "
            f"same sign {r['p_same_sign']:.3f} +- {r['p_same_sign_outer_se']:.3f}  onset<{onset_min:g}d {r['n_dd_onset_before_min']}")


if __name__ == "__main__":
    main()
