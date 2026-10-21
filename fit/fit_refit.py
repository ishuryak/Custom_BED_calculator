"""Bayesian refit of the DI and DD repopulation models to Shuryak, Hall and
Brenner (2018) Table 1 (data/shuryak2018_table1.csv), offline, with CmdStan.

Each fit is written to output/<model>_<variant>/ as soon as it finishes and is
skipped when its summary already exists, so a killed run resumes where it
stopped. Progress lines go to stdout. Seeds are fixed per fit.

Run from the repository root:  python3 fit/fit_refit.py [--loto]
"""
import csv, json, math, os, sys, time, hashlib
import numpy as np
import cmdstanpy, logging
logging.getLogger("cmdstanpy").setLevel(logging.WARNING)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "fit", "output")
# Raw chain CSVs go to a local work folder: writing them inside a OneDrive
# folder failed with "Operation not permitted" (file locking) on 2026-10-21.
RAW = os.environ.get("BED_FIT_WORKDIR", os.path.expanduser("~/.cache/bed_fit_raw"))
DATA = os.path.join(ROOT, "data", "shuryak2018_table1.csv")
SEED = 20261021
CHAINS, WARMUP, SAMPLES = 4, 1500, 2000

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


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def load():
    rows = list(csv.DictReader(open(DATA, encoding="utf-8")))
    refs = [r["ref_no"] for r in rows]
    trials = sorted(set(refs), key=refs.index)
    npat = [int(r["n_patients"]) for r in rows]
    tcp = [float(r["tcp_actual_pct"]) / 100 for r in rows]
    return rows, trials, dict(
        N=len(rows), J=len(trials), trial=[trials.index(t) + 1 for t in refs],
        d=[float(r["dose_per_fraction_gy"]) for r in rows],
        nfx=[float(r["n_fractions"]) for r in rows],
        T=[float(r["total_time_days"]) for r in rows],
        npat=npat, y=[int(round(n * p)) for n, p in zip(npat, tcp)])


def stan_data(base, variant):
    v = VARIANTS[variant]
    return dict(base, prior_log_ab_mean=math.log(10), prior_log_ab_sd=v["prior_log_ab_sd"],
                fix_ab=v["fix_ab"], ab_fixed=10.0, fix_g=v["fix_g"], g_fixed=0.0154,
                overdispersion=v["overdispersion"])


def model(name):
    # Compile and run from a local copy: launching a freshly built executable
    # inside the OneDrive folder failed with "Operation not permitted".
    import shutil
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
    raw = os.path.join(RAW, os.path.relpath(outdir, OUT))
    os.makedirs(raw, exist_ok=True)
    # Plausible, jittered starting values per chain: Stan's default random
    # inits (alpha up to e^2 per Gy) give TCP of exactly 0 or 1 and fail.
    rng = np.random.default_rng(seed)
    def init():
        x = dict(log_ab_raw=0.0, alpha=0.05 * math.exp(rng.normal(0, .3)),
                 lambda_=0.0, g_free=0.012, c=list(rng.normal(0, .2, data["J"])),
                 sigma_raw=0.2, z=list(rng.normal(0, .1, data["N"])))
        x["lambda"] = (0.02 if name == "di" else 0.5) * math.exp(rng.normal(0, .3)); del x["lambda_"]
        if name == "di": x["Tk"] = 28 + rng.normal(0, 3)
        else: x["C"] = 10 * math.exp(rng.normal(0, .3))
        return x
    f = m.sample(data=data, chains=chains, parallel_chains=chains, iter_warmup=warmup,
                 iter_sampling=samples, seed=seed, adapt_delta=0.99, max_treedepth=12,
                 show_progress=False, output_dir=raw, inits=[init() for _ in range(chains)])
    return f


def summarize(f, name, outdir, rows):
    s = f.summary()
    keep = PARAMS[name]
    diag = f.method_variables()
    div = int(diag["divergent__"].sum())
    tdmax = int((diag["treedepth__"] >= 12).sum())
    rhat = float(s.loc[keep, [c for c in s.columns if c.lower().startswith("r_hat")][0]].max())
    esscol = [c for c in s.columns if "ess_bulk" in c.lower()][0]
    ess = float(s.loc[keep, esscol].min())
    dr = f.draws_pd()
    post = {p: dict(median=float(dr[p].median()), q2_5=float(dr[p].quantile(.025)),
                    q97_5=float(dr[p].quantile(.975)), q10=float(dr[p].quantile(.10)),
                    q90=float(dr[p].quantile(.90))) for p in keep}
    # posterior predictive check per arm: observed count within the 90% interval of y_rep
    ppc = []
    for i, r in enumerate(rows, start=1):
        yr = dr[f"y_rep[{i}]"].to_numpy()
        n = int(r["n_patients"]); yobs = int(round(n * float(r["tcp_actual_pct"]) / 100))
        lo, hi = np.quantile(yr, [.05, .95])
        ppc.append(dict(row=i, author=r["first_author"], arm=r["arm"], n=n, obs_pct=100 * yobs / n,
                        pred_median_pct=100 * float(np.median(yr)) / n, pred_lo_pct=100 * lo / n,
                        pred_hi_pct=100 * hi / n, inside_90=bool(lo <= yobs <= hi)))
    ll = dr[[f"log_lik[{i}]" for i in range(1, len(rows) + 1)]].to_numpy()
    res = dict(model=name, divergences=div, treedepth_hits=tdmax, max_rhat=rhat, min_ess_bulk=ess,
               n_draws=len(dr), posterior=post, ppc=ppc,
               ppc_inside_90=sum(p["inside_90"] for p in ppc),
               mean_log_lik=float(ll.sum(axis=1).mean()),
               correlations={f"{a}~{b}": float(np.corrcoef(np.log(dr[a]), np.log(dr[b]))[0, 1])
                             for k, a in enumerate(keep) for b in keep[k + 1:]
                             if dr[a].std() > 0 and dr[b].std() > 0 and (dr[a] > 0).all() and (dr[b] > 0).all()})
    json.dump(res, open(os.path.join(outdir, "summary.json"), "w"), indent=1)
    cols = keep + ["lp__"]
    dr[cols].to_csv(os.path.join(outdir, "draws_params.csv"), index_label="draw")
    return res


def core_eta(name, dr, d, nfx, T):
    """eta without the trial intercept and the arm effect, for every draw (rows) and arm (cols)."""
    a = dr["alpha"].to_numpy()[:, None]; ab = dr["ab"].to_numpy()[:, None]
    lam = dr["lambda"].to_numpy()[:, None]; g = dr["g"].to_numpy()[:, None]
    bed = nfx[None, :] * d[None, :] * (1 + d[None, :] / ab)
    sm = lambda x: np.where(3 * x > 30, x, np.log1p(np.exp(np.minimum(3 * x, 30))) / 3)
    if name == "di":
        rep = lam * sm(T[None, :] - dr["Tk"].to_numpy()[:, None])
    else:
        pk = -np.expm1(-a * bed / T[None, :])
        onset = dr["c_over_alpha"].to_numpy()[:, None] * T[None, :] / bed
        rep = lam * pk * sm(T[None, :] - onset)
    return -a * bed + rep + g * T[None, :]


def loto(name, m, base, rows, trials):
    """Leave one trial out. Fit without trial j; predict trial j's other arms given
    its first-listed arm (which fixes the trial intercept). Returns elpd per trial."""
    rng = np.random.default_rng(SEED)
    res_path = os.path.join(OUT, f"loto_{name}.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}
    d = np.array(base["d"]); nfx = np.array(base["nfx"]); T = np.array(base["T"])
    tr = np.array(base["trial"]); y = np.array(base["y"]); npat = np.array(base["npat"])
    from scipy.stats import binom
    for j, ref in enumerate(trials, start=1):
        if ref in res:
            continue
        keep = tr != j
        idx = np.where(tr == j)[0]
        # renumber trials
        old = sorted(set(tr[keep])); mp = {o: k + 1 for k, o in enumerate(old)}
        sub = dict(base, N=int(keep.sum()), J=len(old), trial=[mp[t] for t in tr[keep]],
                   d=list(d[keep]), nfx=list(nfx[keep]), T=list(T[keep]),
                   npat=[int(x) for x in npat[keep]], y=[int(x) for x in y[keep]])
        f = fit_one(m, name, stan_data(sub, "primary"), SEED + j,
                    os.path.join(OUT, f"loto_{name}", ref), chains=4, warmup=800, samples=500)
        dr = f.draws_pd()
        h = core_eta(name, dr, d[idx], nfx[idx], T[idx])
        r0 = idx[0]
        p0 = np.clip(y[r0] / npat[r0], 0.5 / npat[r0], 1 - 0.5 / npat[r0])
        target = math.log(-math.log(p0))
        sig = dr["sigma"].to_numpy()[:, None]
        elpd = 0.0
        for k in range(1, len(idx)):
            eps = rng.normal(size=(len(dr), 20)) * sig * math.sqrt(2)
            eta = target + (h[:, k] - h[:, 0])[:, None] + eps
            p = np.exp(-np.exp(eta))
            dens = binom.pmf(y[idx[k]], npat[idx[k]], p).mean()
            elpd += math.log(max(dens, 1e-300))
        res[ref] = dict(elpd=elpd, n_pred_arms=len(idx) - 1,
                        divergences=int(f.method_variables()["divergent__"].sum()))
        json.dump(res, open(res_path, "w"), indent=1)
        log(f"loto {name} trial {ref} ({j}/{len(trials)}) elpd {elpd:.2f}")
    return res


def main():
    rows, trials, base = load()
    sha = hashlib.sha256(open(DATA, "rb").read()).hexdigest()
    log("data", DATA, "sha256", sha, "arms", base["N"], "trials", base["J"],
        "cmdstan", cmdstanpy.cmdstan_version())
    models = {n: model(n) for n in ("di", "dd")}
    for name, m in models.items():
        for k, variant in enumerate(VARIANTS):
            outdir = os.path.join(OUT, f"{name}_{variant}")
            if os.path.exists(os.path.join(outdir, "summary.json")):
                log("skip (done)", name, variant); continue
            log("fitting", name, variant)
            f = fit_one(m, name, stan_data(base, variant), SEED + 100 * k + (0 if name == "di" else 50), outdir)
            r = summarize(f, name, outdir, rows)
            log(f"done {name} {variant}: div {r['divergences']} rhat {r['max_rhat']:.3f} "
                f"ess {r['min_ess_bulk']:.0f} ppc {r['ppc_inside_90']}/35")
    if "--loto" in sys.argv:
        for name, m in models.items():
            loto(name, m, base, rows, trials)
    json.dump(dict(data_sha256=sha, seed=SEED, chains=CHAINS, warmup=WARMUP, samples=SAMPLES,
                   cmdstan=".".join(map(str, cmdstanpy.cmdstan_version())),
                   cmdstanpy=cmdstanpy.__version__, numpy=np.__version__, variants=VARIANTS),
              open(os.path.join(OUT, "provenance.json"), "w"), indent=1)
    log("all done")


if __name__ == "__main__":
    main()
