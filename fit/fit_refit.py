"""Bayesian refit of the DI and DD repopulation models to arm-level trial data, offline,
with CmdStan.

Default data: Shuryak, Hall and Brenner (2018) Table 1 (data/shuryak2018_table1.csv), the fit the
app shipped in version 3.0. BED_FIT_DATA / BED_FIT_OUT point a refit at another data file and
output folder.

Likelihood: binomial on y = round(N x reported control rate). Most reported rates are Kaplan-Meier
or cumulative-incidence estimates, not counts, so this is an approximation (plan decision 5);
the effective-N sensitivity set halves N.

Run identity: each fit folder holds identity.json (data bytes, Stan sources, this file, priors,
sampler settings, seed, software versions). A finished fit is reused only when its identity
matches; otherwise the run stops rather than mixing or relabeling results. Each fit is written as
soon as it finishes, so a killed run resumes where it stopped. Progress lines go to stdout.

Diagnostics cover every sampled parameter and the regimen contrasts: divergences, tree-depth hits,
R-hat, bulk and tail ESS and E-BFMI per chain. A fit passes the release gate (plan decision 9)
only with 0 divergences, 0 tree-depth hits, R-hat <= 1.01, bulk and tail ESS >= 400 and
E-BFMI >= 0.3.

Run from the repository root:  python3 fit/fit_refit.py [--loto]
Env: BED_FIT_PRIOR_C_SD (trial intercept prior sd, default 3), BED_FIT_VARIANTS (comma list),
     BED_FIT_ITER ("warmup,samples", smoke tests only).
"""
import csv, json, math, os, sys, time, hashlib
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.environ.get("BED_FIT_OUT", os.path.join(ROOT, "fit", "output"))
# Raw chain CSVs go to a local work folder: writing them inside a OneDrive
# folder failed with "Operation not permitted" (file locking) on 2026-10-21.
RAW = os.environ.get("BED_FIT_WORKDIR", os.path.expanduser("~/.cache/bed_fit_raw"))
DATA = os.environ.get("BED_FIT_DATA", os.path.join(ROOT, "data", "shuryak2018_table1.csv"))
PRIOR_C_SD = float(os.environ.get("BED_FIT_PRIOR_C_SD", "3"))
SEED = 20261021
CHAINS, WARMUP, SAMPLES = 4, 1500, 2000
# BED_FIT_ITER="warmup,samples" shortens chains for a smoke test (recorded in each fit identity).
if os.environ.get("BED_FIT_ITER"):
    WARMUP, SAMPLES = (int(x) for x in os.environ["BED_FIT_ITER"].split(","))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate as G   # the one diagnostic gate (review R01, 2026-10-26)
GATE = {k: G.POLICY[k] for k in ("max_rhat", "min_ess_bulk", "min_ess_tail", "min_ebfmi")}

# Prior on tumor alpha/beta: lognormal, median 10 Gy; sd 0.354 puts 95% within
# 5-20 Gy (the plausible range stated in the 2018 paper).
VARIANTS = {
    "primary":   dict(prior_log_ab_sd=0.354, fix_ab=0, fix_g=0, overdispersion=1),
    "binomial":  dict(prior_log_ab_sd=0.354, fix_ab=0, fix_g=0, overdispersion=0),
    "paper_restrictions": dict(prior_log_ab_sd=0.354, fix_ab=1, fix_g=1, overdispersion=1),
    "wide_ab":   dict(prior_log_ab_sd=0.707, fix_ab=0, fix_g=0, overdispersion=1),
}
PARAMS = {"di": ["alpha", "ab", "lambda", "Tk", "g", "sigma"],
          "dd": ["alpha", "ab", "lambda", "C", "c_over_alpha", "g", "sigma"]}
# Parameters declared in each model's parameters block (diagnosed in full, with c[] and z[]).
SAMPLED = {"di": ["log_ab_raw", "alpha", "lambda", "Tk", "g_free", "sigma_raw"],
           "dd": ["log_ab_raw", "alpha", "lambda", "C", "g_free", "sigma_raw"]}

ANCHOR = (35, 2.0, 47.0)            # fractions, Gy per fraction, overall days (inclusive)
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


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def picks_sha256(pick):
    """Hash of the selected outer-draw indices as explicit little-endian int64 bytes.

    The same bytes NumPy 2 on 64-bit Linux wrote natively when the 2026-10-25c identities were
    recorded, so those identities verify unchanged; the explicit dtype makes the hash portable to
    platforms whose default integer is 32-bit (review F07 of 2026-10-26: Windows NumPy 1.26).
    """
    return hashlib.sha256(np.asarray(pick).astype("<i8").tobytes()).hexdigest()


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load(path=None):
    """Arm rows, fitting units, and the Stan data shared by all variants.
    BED_LQ uses fraction components when the file has them (total_dose_gy, sum_d2_gy2),
    else n x d and n x d^2."""
    rows = list(csv.DictReader(open(path or DATA, encoding="utf-8")))
    refs = [r["ref_no"] for r in rows]
    trials = sorted(set(refs), key=refs.index)
    npat = [int(r["n_patients"]) for r in rows]
    tcp = [float(r["tcp_actual_pct"]) / 100 for r in rows]
    d = [float(r["dose_per_fraction_gy"]) for r in rows]
    n = [float(r["n_fractions"]) for r in rows]
    has = all(r.get("total_dose_gy") not in (None, "") and r.get("sum_d2_gy2") not in (None, "") for r in rows)
    dsum = [float(r["total_dose_gy"]) for r in rows] if has else [a * b for a, b in zip(n, d)]
    d2sum = [float(r["sum_d2_gy2"]) for r in rows] if has else [a * b * b for a, b in zip(n, d)]
    ev = [(k * x, k * x * x, T) for _, k, x, T in REGIMENS]
    n0, d0, T0 = ANCHOR
    return rows, trials, dict(
        N=len(rows), J=len(trials), trial=[trials.index(t) + 1 for t in refs],
        dsum=dsum, d2sum=d2sum, T=[float(r["total_time_days"]) for r in rows],
        npat=npat, y=[int(round(k * p)) for k, p in zip(npat, tcp)],
        R=len(REGIMENS), eval_dsum=[e[0] for e in ev], eval_d2sum=[e[1] for e in ev], eval_T=[e[2] for e in ev],
        anchor_dsum=n0 * d0, anchor_d2sum=n0 * d0 * d0, anchor_T=T0, prior_c_sd=PRIOR_C_SD)


def stan_data(base, variant):
    v = VARIANTS[variant]
    return dict(base, prior_log_ab_mean=math.log(10), prior_log_ab_sd=v["prior_log_ab_sd"],
                fix_ab=v["fix_ab"], ab_fixed=10.0, fix_g=v["fix_g"], g_fixed=0.0154,
                overdispersion=v["overdispersion"])


def versions():
    import cmdstanpy
    return dict(cmdstan=".".join(map(str, cmdstanpy.cmdstan_version() or ())), cmdstanpy=cmdstanpy.__version__,
                numpy=np.__version__, python=sys.version.split()[0])


def identity(name, data, seed, chains, warmup, samples, extra=None):
    """Everything a fit's result depends on. Two fits with equal identity hashes are interchangeable."""
    stan = {f: sha_file(os.path.join(ROOT, "fit", "stan", f))
            for f in (name + ".stan", "common_functions.stanfunctions")}
    data_hash = hashlib.sha256(json.dumps(data, sort_keys=True, default=float).encode()).hexdigest()
    ident = dict(model=name, data_file=os.path.basename(DATA), data_file_sha256=sha_file(DATA),
                 stan_data_sha256=data_hash, stan_sha256=stan, runner_sha256=sha_file(os.path.abspath(__file__)),
                 seed=seed, chains=chains, warmup=warmup, samples=samples, adapt_delta=0.99, max_treedepth=12,
                 software=versions(), **(extra or {}))
    ident["identity_sha256"] = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()
    return ident


def check_identity(outdir, ident, done_file="summary.json"):
    """True if a finished fit with this identity exists; stops on a finished fit with another identity."""
    idp = os.path.join(outdir, "identity.json")
    if not os.path.exists(os.path.join(outdir, done_file)):
        return False
    old = json.load(open(idp)) if os.path.exists(idp) else {}
    if old.get("identity_sha256") != ident["identity_sha256"]:
        diff = sorted(k for k in set(old) | set(ident) if old.get(k) != ident.get(k))
        sys.exit(f"REFUSED: {outdir} holds a finished fit with a different identity ({', '.join(diff)}). "
                 "Use a new output folder.")
    return True


def model(name):
    # Compile and run from a local copy: launching a freshly built executable
    # inside the OneDrive folder failed with "Operation not permitted".
    import shutil, cmdstanpy, logging
    logging.getLogger("cmdstanpy").setLevel(logging.WARNING)
    build = os.path.join(RAW, "build")
    os.makedirs(build, exist_ok=True)
    for f in (name + ".stan", "common_functions.stanfunctions"):
        src, dst = os.path.join(ROOT, "fit", "stan", f), os.path.join(build, f)
        if not os.path.exists(dst) or open(src, "rb").read() != open(dst, "rb").read():
            shutil.copyfile(src, dst)
    return cmdstanpy.CmdStanModel(stan_file=os.path.join(build, name + ".stan"),
                                  stanc_options={"include-paths": [build]})


def fit_one(m, name, data, seed, outdir, chains=CHAINS, warmup=WARMUP, samples=SAMPLES):
    os.makedirs(outdir, exist_ok=True)
    raw = os.path.join(RAW, os.path.relpath(outdir, os.path.dirname(OUT)))
    os.makedirs(raw, exist_ok=True)
    # Plausible, jittered starting values per chain: Stan's default random
    # inits (alpha up to e^2 per Gy) give TCP of exactly 0 or 1 and fail.
    rng = np.random.default_rng(seed)
    def init():
        x = dict(log_ab_raw=0.0, alpha=0.05 * math.exp(rng.normal(0, .3)),
                 g_free=0.012, c=list(rng.normal(0, .2, data["J"])),
                 sigma_raw=0.2, z=list(rng.normal(0, .1, data["N"])))
        x["lambda"] = (0.02 if name == "di" else 0.5) * math.exp(rng.normal(0, .3))
        if name == "di": x["Tk"] = 28 + rng.normal(0, 3)
        else: x["C"] = 10 * math.exp(rng.normal(0, .3))
        return x
    return m.sample(data=data, chains=chains, parallel_chains=chains, iter_warmup=warmup,
                    iter_sampling=samples, seed=seed, adapt_delta=0.99, max_treedepth=12,
                    show_progress=False, output_dir=raw, inits=[init() for _ in range(chains)])


def diagnostics(f, sampled, max_treedepth=12):
    """Full diagnostics: every sampled parameter (vectors included) and every delta_eval contrast."""
    s = f.summary()
    col = lambda key: [c for c in s.columns if c.lower().replace("-", "_") == key][0]
    rh, eb, et = col("r_hat"), col("ess_bulk"), col("ess_tail")
    names = [i for i in s.index if i.split("[")[0] in set(sampled) | {"c", "z", "delta_eval"}]
    sub = s.loc[names]
    # Every quantity is checked on its own before min/max, which skip NaN (review F02). A quantity
    # constant across draws is a failure, not an exclusion: none of these models samples a fixed one.
    import pandas as pd
    tab = pd.DataFrame(dict(rhat=sub[rh], ess_bulk=sub[eb], ess_tail=sub[et],
                            sd=sub["StdDev"] if "StdDev" in sub.columns else np.nan), index=sub.index)
    bad = G.invalid_quantities(tab, list(sampled) + ["delta_eval"])
    mv = f.method_variables()
    div = int(mv["divergent__"].sum())
    tdh = int((mv["treedepth__"] >= max_treedepth).sum())
    e = f.draws(concat_chains=False)[:, :, list(f.column_names).index("energy__")]
    ebfmi = [float(np.sum(np.diff(e[:, k]) ** 2) / np.sum((e[:, k] - e[:, k].mean()) ** 2)) for k in range(e.shape[1])]
    res = dict(divergences=div, treedepth_hits=tdh, max_rhat=float(sub[rh].max()), worst_rhat=str(sub[rh].idxmax()),
               min_ess_bulk=float(sub[eb].min()), worst_ess_bulk=str(sub[eb].idxmin()),
               min_ess_tail=float(sub[et].min()), worst_ess_tail=str(sub[et].idxmin()),
               min_ebfmi=float(min(ebfmi)), n_quantities=int(len(sub)), invalid_quantities=bad)
    res["gate_pass"] = G.passes(res)
    res["gate_policy"] = G.POLICY["version"]
    return res


def summarize(f, name, outdir, rows):
    keep = PARAMS[name]
    dg = diagnostics(f, SAMPLED[name])
    dr = f.draws_pd()
    post = {p: dict(median=float(dr[p].median()), q2_5=float(dr[p].quantile(.025)),
                    q97_5=float(dr[p].quantile(.975)), q10=float(dr[p].quantile(.10)),
                    q90=float(dr[p].quantile(.90))) for p in keep}
    contrasts = []
    for r, (lab, *_rest) in enumerate(REGIMENS, start=1):
        x = dr[f"delta_eval[{r}]"]
        contrasts.append(dict(regimen=lab, median=float(x.median()), q10=float(x.quantile(.1)), q90=float(x.quantile(.9))))
    # Posterior predictive check per arm: pseudo-count within the 90% interval of y_rep.
    # obs_source_pct is the reported rate; obs_pct is the rounded pseudo-count the likelihood uses.
    ppc = []
    for i, r in enumerate(rows, start=1):
        yr = dr[f"y_rep[{i}]"].to_numpy()
        n = int(r["n_patients"]); yobs = int(round(n * float(r["tcp_actual_pct"]) / 100))
        lo, hi = np.quantile(yr, [.05, .95])
        ppc.append(dict(row=i, author=r["first_author"], arm=r["arm"], n=n,
                        obs_source_pct=float(r["tcp_actual_pct"]), obs_pct=100 * yobs / n,
                        pred_median_pct=100 * float(np.median(yr)) / n, pred_lo_pct=100 * lo / n,
                        pred_hi_pct=100 * hi / n, inside_90=bool(lo <= yobs <= hi)))
    ll = dr[[f"log_lik[{i}]" for i in range(1, len(rows) + 1)]].to_numpy()
    res = dict(model=name, **dg, n_draws=len(dr), posterior=post, contrasts_vs_anchor=contrasts, ppc=ppc,
               ppc_inside_90=sum(p["inside_90"] for p in ppc), likelihood="binomial on rounded N x rate (approximation)",
               mean_log_lik=float(ll.sum(axis=1).mean()),
               correlations={f"{a}~{b}": float(np.corrcoef(np.log(dr[a]), np.log(dr[b]))[0, 1])
                             for k, a in enumerate(keep) for b in keep[k + 1:]
                             if dr[a].std() > 0 and dr[b].std() > 0 and (dr[a] > 0).all() and (dr[b] > 0).all()})
    json.dump(res, open(os.path.join(outdir, "summary.json"), "w"), indent=1)
    cols = ["chain__", "iter__"] + keep + ["lp__"] + [f"delta_eval[{r}]" for r in range(1, len(REGIMENS) + 1)]
    dr[cols].to_csv(os.path.join(outdir, "draws_params.csv"), index_label="draw")
    return res


def core_eta(name, dr, dsum, d2sum, T):
    """eta without the trial intercept and the arm effect, for every draw (rows) and arm (cols)."""
    a = dr["alpha"].to_numpy()[:, None]; ab = dr["ab"].to_numpy()[:, None]
    lam = dr["lambda"].to_numpy()[:, None]; g = dr["g"].to_numpy()[:, None]
    bed = dsum[None, :] + d2sum[None, :] / ab
    sm = lambda x: np.where(3 * x > 30, x, np.log1p(np.exp(np.minimum(3 * x, 30))) / 3)
    if name == "di":
        rep = lam * sm(T[None, :] - dr["Tk"].to_numpy()[:, None])
    else:
        pk = -np.expm1(-a * bed / T[None, :])
        onset = dr["c_over_alpha"].to_numpy()[:, None] * T[None, :] / bed
        rep = lam * pk * sm(T[None, :] - onset)
    return -a * bed + rep + g * T[None, :]


def loto(name, m, base, rows, trials):
    """Leave one STUDY out (all of its fitting units, e.g. both CAIR periods). For each held-out
    unit, predict its other arms given its first-listed arm, which fixes that unit's intercept.
    This is a conditional plug-in score (reference-arm sampling uncertainty is not integrated),
    not a full held-out predictive density. Returns elpd per study with fold diagnostics."""
    from scipy.stats import binom
    rng = np.random.default_rng(SEED)
    res_path = os.path.join(OUT, f"loso_{name}.json")
    dsum = np.array(base["dsum"]); d2sum = np.array(base["d2sum"]); T = np.array(base["T"])
    tr = np.array(base["trial"]); y = np.array(base["y"]); npat = np.array(base["npat"])
    study = np.array([r["first_author"] for r in rows])
    studies = sorted(set(study), key=list(study).index)
    # Cache identity (review R02): saved study scores are reused only for the same data, model,
    # runner, priors, sampler, study list and units; anything else is refused.
    ident = identity(name, stan_data(base, "primary"), SEED, 4, min(800, WARMUP), min(500, SAMPLES),
                     extra=dict(analysis="leave one study out", studies=studies, units=list(trials),
                                prior_c_sd=PRIOR_C_SD, score="conditional plug-in, 20 normal draws per posterior draw"))
    res = json.load(open(res_path)) if os.path.exists(res_path) else {"_identity": ident}
    if res.get("_identity", {}).get("identity_sha256") != ident["identity_sha256"]:
        sys.exit(f"REFUSED: {res_path} holds scores from another data set, model or runner. Use a new output folder.")
    for s_i, st in enumerate(studies, start=1):
        if st in res:
            continue
        keep = study != st
        old = sorted(set(tr[keep])); mp = {o: k + 1 for k, o in enumerate(old)}
        sub = dict(base, N=int(keep.sum()), J=len(old), trial=[mp[t] for t in tr[keep]],
                   dsum=list(dsum[keep]), d2sum=list(d2sum[keep]), T=list(T[keep]),
                   npat=[int(x) for x in npat[keep]], y=[int(x) for x in y[keep]])
        f = fit_one(m, name, stan_data(sub, "primary"), SEED + s_i, os.path.join(OUT, f"loso_{name}", st),
                    chains=4, warmup=min(800, WARMUP), samples=min(500, SAMPLES))
        dr = f.draws_pd()
        sig = dr["sigma"].to_numpy()[:, None]
        elpd, npred, units = 0.0, 0, {}
        for j in sorted(set(tr[~keep])):
            idx = np.where(tr == j)[0]
            h = core_eta(name, dr, dsum[idx], d2sum[idx], T[idx])
            r0 = idx[0]
            p0 = np.clip(y[r0] / npat[r0], 0.5 / npat[r0], 1 - 0.5 / npat[r0])
            target = math.log(-math.log(p0))
            e_u = 0.0
            for k in range(1, len(idx)):
                eps = rng.normal(size=(len(dr), 20)) * sig * math.sqrt(2)
                p = np.exp(-np.exp(target + (h[:, k] - h[:, 0])[:, None] + eps))
                e_u += math.log(max(binom.pmf(y[idx[k]], npat[idx[k]], p).mean(), 1e-300))
            units[trials[j - 1]] = e_u; elpd += e_u; npred += len(idx) - 1
        mv = f.method_variables()
        res[st] = dict(elpd=elpd, n_pred_arms=npred, units=units,
                       divergences=int(mv["divergent__"].sum()),
                       score="conditional plug-in (first arm of each held-out unit fixes its intercept)")
        json.dump(res, open(res_path, "w"), indent=1)
        log(f"loso {name} study {st} ({s_i}/{len(studies)}) elpd {elpd:.2f} div {res[st]['divergences']}")
    return res


def main():
    import cmdstanpy
    rows, trials, base = load()
    log("data", DATA, "sha256", sha_file(DATA), "arms", base["N"], "units", base["J"],
        "prior_c_sd", PRIOR_C_SD, "cmdstan", cmdstanpy.cmdstan_version())
    want = os.environ.get("BED_FIT_VARIANTS")
    variants = [v for v in VARIANTS if not want or v in want.split(",")]
    models = {n: model(n) for n in ("di", "dd")}
    idents = {}
    for name, m in models.items():
        for k, variant in enumerate(VARIANTS):
            if variant not in variants:
                continue
            outdir = os.path.join(OUT, f"{name}_{variant}")
            seed = SEED + 100 * k + (0 if name == "di" else 50)
            data = stan_data(base, variant)
            ident = identity(name, data, seed, CHAINS, WARMUP, SAMPLES, dict(variant=variant))
            idents[f"{name}_{variant}"] = ident["identity_sha256"]
            if check_identity(outdir, ident):
                log("skip (done, identity matches)", name, variant); continue
            log("fitting", name, variant)
            f = fit_one(m, name, data, seed, outdir)
            r = summarize(f, name, outdir, rows)
            json.dump(ident, open(os.path.join(outdir, "identity.json"), "w"), indent=1)
            log(f"done {name} {variant}: div {r['divergences']} rhat {r['max_rhat']:.4f} ({r['worst_rhat']}) "
                f"ess bulk {r['min_ess_bulk']:.0f} tail {r['min_ess_tail']:.0f} ebfmi {r['min_ebfmi']:.2f} "
                f"gate {'PASS' if r['gate_pass'] else 'FAIL'} ppc {r['ppc_inside_90']}/{len(rows)}")
    if "--loto" in sys.argv:
        for name, m in models.items():
            loto(name, m, base, rows, trials)
    json.dump(dict(data_file=DATA, data_sha256=sha_file(DATA), prior_c_sd=PRIOR_C_SD, seed=SEED, chains=CHAINS,
                   warmup=WARMUP, samples=SAMPLES, variants=VARIANTS, fit_identities=idents, gate=GATE, **versions()),
              open(os.path.join(OUT, "provenance.json"), "w"), indent=1)
    log("all done")


if __name__ == "__main__":
    main()
