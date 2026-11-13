"""Assemble the comparison of plan revision 3, section 4, from VERIFIED results only.

    python3 fit/frequentist/compare.py      -> fit/output_frequentist/comparison.json + figures/*.png

Refuses unless: each maximum-likelihood record and every profile verifies its own digest and carries
the same run identity; each bootstrap generation loads through fit_ml.load_generation (manifest
present, every record and the summary verified). Posterior summaries are computed from the deployed
draws with the same target functions; prior intervals by simulation from the full declared joint
prior through the same functions (seed 20261112, 200,000 draws), never from marginal endpoints.
"""
import csv
import glob
import gzip
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from frequentist import fit_ml as FM  # noqa: E402
from frequentist import likelihood as L  # noqa: E402

ROOT = FM.ROOT
N_PRIOR = 200_000
LOG_SCALE = {"alpha", "ab", "lam", "r"}                 # lognormal priors: widths compared on the log scale
U02 = os.path.join(ROOT, "docs", "development", "2026_11_10_early_toxicity", "step5_review", "overdispersed", "U02")


def verified(path, identity=None):
    r = json.load(open(path))
    if r.get("record_sha256") != FM.record_digest(r) or (identity and r.get("identity") != identity):
        raise RuntimeError(f"{path}: does not verify; refusing to assemble")
    return r


def generation(kind, support):
    m = FM.Model(kind, support)
    ident = FM.run_identity(kind, m.d, support)
    gdir = FM.gen_dir(kind, support, ident)
    mle = verified(os.path.join(gdir, "mle.json"), ident)
    prof = {}
    for p in glob.glob(os.path.join(gdir, "profiles", "*.json")):
        r = verified(p, ident)
        prof[os.path.basename(p)[:-5]] = r
    boot = None
    bdirs = sorted(glob.glob(os.path.join(gdir, "bootstrap_B*")))
    if bdirs:
        man, recs, summ = FM.load_generation(bdirs[-1])
        if man["identity"] != ident:
            raise RuntimeError(f"{bdirs[-1]}: bootstrap of another identity")
        boot = summ
    return m, ident, gdir, mle, prof, boot


# --- posterior and prior samples of every target ------------------------------------------------
def tumor_posterior(m):
    rows = list(csv.DictReader(open(os.path.join(ROOT, "draws", "tumor_di_primary.csv"))))
    first = {}
    for r in rows:                                   # 1000 distinct DI draws, each repeated 4 times
        first.setdefault(r["di_draw"], r)
    th = {k: np.array([float(r[k]) for r in first.values()]) for k in ("alpha", "ab", "Tk", "g")}
    th["lam"] = np.array([float(r["lambda"]) for r in first.values()])
    out = {k: f(th, m.d) for k, f in m.targets.items() if k != "sigma" and not k.startswith("contrast")}
    con = np.array([L.tumor_contrasts({k: v[i] for k, v in th.items()}, m.d) for i in range(len(th["ab"]))])
    out.update({f"contrast[{k}]": con[:, k] for k in range(10)})
    s = json.load(open(os.path.join(ROOT, "fit", "output_refit_2026_10_25c_main", "di_primary", "summary.json")))["posterior"]["sigma"]
    return out, dict(sigma=(s["median"], s["q2_5"], s["q97_5"])), len(first)


def early_posterior(m):
    import pandas as pd
    dr = pd.read_csv(os.path.join(U02, "draws.csv.gz"))
    out = dict(kappa=dr["kappa"].to_numpy(), ab=dr["ab"].to_numpy(), r=dr["r"].to_numpy(), gamma=dr["gamma[1]"].to_numpy(),
               tau=dr["tau0"].to_numpy(), sigma_arm=dr["sigma_arm[1]"].to_numpy())
    th = dict(kappa=out["kappa"][:, None], ab=out["ab"][:, None], r=out["r"][:, None])
    lor = L.early_contrasts(th, m.d)
    out.update({f"lor[{k}]": lor[:, k] for k in range(4)})
    return out, {}, len(dr)


def tumor_prior(m, rng):
    n = N_PRIOR
    Tk = rng.normal(28, 10, 4 * n); Tk = Tk[(Tk > 0) & (Tk < 60)][:n]           # truncated to the Stan bounds
    th = dict(alpha=rng.lognormal(math.log(0.2), 1, n), ab=rng.lognormal(math.log(10), 0.354, n),
              lam=rng.lognormal(math.log(0.1), 1.5, n), Tk=Tk, g=rng.uniform(0.0046, 0.021, n),
              sigma=np.abs(rng.normal(0, 0.5, n)))
    out = {k: f(th, m.d) for k, f in m.targets.items() if not k.startswith("contrast")}
    con = L.tumor_contrasts({k: v[:, None] for k, v in th.items()}, m.d)
    out.update({f"contrast[{k}]": con[:, k] for k in range(10)})
    return out


def early_prior(m, rng):
    n = N_PRIOR
    th = dict(kappa=np.abs(rng.normal(0, 1, n)) / 10, ab=rng.lognormal(math.log(10), 0.354, n),
              r=rng.lognormal(math.log(0.8), 0.5, n), gamma=rng.normal(0, 1, n), tau=np.abs(rng.normal(0, 1, n)),
              sigma_arm=np.abs(rng.normal(0, 0.5, n)))
    out = {k: th[k] for k in ("kappa", "ab", "r", "gamma", "tau", "sigma_arm")}
    lor = L.early_contrasts({k: th[k][:, None] for k in ("kappa", "ab", "r")}, m.d)
    out.update({f"lor[{k}]": lor[:, k] for k in range(4)})
    return out


def q95(v):
    return [float(np.quantile(v, 0.025)), float(np.median(v)), float(np.quantile(v, 0.975))]


def width(lo, hi, name):
    if lo is None or hi is None:
        return None
    if name in LOG_SCALE:
        return math.log(hi) - math.log(lo) if lo > 0 and hi > 0 else None
    return hi - lo


def interval_text(piece):
    lo, hi = piece["lower"], piece["upper"]
    return (lo["value"] if lo["closed"] else None, hi["value"] if hi["closed"] else None,
            "" if lo["closed"] else lo.get("why", ""), "" if hi["closed"] else hi.get("why", ""))


def assemble():
    res = dict(models={})
    rng = np.random.default_rng(FM.SEED)
    for kind, supports in FM.MODEL_SUPPORTS.items():
        for support in supports:
            m, ident, gdir, mle, prof, boot = generation(kind, support)
            post, post_fixed, n_post = (tumor_posterior if kind == "tumor" else early_posterior)(m)
            prior = (tumor_prior if kind == "tumor" else early_prior)(m, rng)
            rows = []
            for t in m.targets:
                p = prof.get(t)
                row = dict(target=t, units=FM.UNITS[t], mle=mle["targets"][t])
                if p is not None:
                    pieces = [interval_text(c) for c in p["confidence_set"]]
                    row.update(profile_set=p["confidence_set"], profile_flat=p["flat"], profile_disconnected=p["disconnected"],
                               profile_better_point=p["better_point"] is not None)
                    if len(pieces) == 1:
                        row.update(profile_lo=pieces[0][0], profile_hi=pieces[0][1], profile_lo_open=pieces[0][2],
                                   profile_hi_open=pieces[0][3])
                if boot and t in boot["targets"]:
                    row["bootstrap"] = boot["targets"][t]
                if t in post:
                    row["posterior"] = q95(post[t])
                elif t in post_fixed:
                    row["posterior"] = [post_fixed[t][1], post_fixed[t][0], post_fixed[t][2]]
                row["prior"] = q95(prior[t]) if t in prior else None
                if row["prior"] and row.get("profile_lo") is not None and row.get("profile_hi") is not None:
                    wl, wp = width(row["profile_lo"], row["profile_hi"], t), width(row["prior"][0], row["prior"][2], t)
                    row["width_ratio"] = None if (wl is None or not wp) else wl / wp
                rows.append(row)
            sens = {k: v for k, v in prof.items() if "__sigma_" in k}
            res["models"][f"{kind}_{support}"] = dict(
                identity=ident, generation=os.path.relpath(gdir, ROOT), mle=dict(ll=mle["ll"], boundary=mle["boundary"],
                                                                                    n_ties=mle["n_ties"], nodes=mle["nodes"]),
                n_posterior_draws=n_post, rows=rows,
                bootstrap=None if not boot else {k: boot[k] for k in ("B", "ok", "unresolved_failures", "failure_rate",
                                                                      "provisional", "boundary_solutions", "seconds_total")},
                sigma_sensitivity={k: dict(target=v["target"], sigma=v["sigma_fixed"], estimate=v["estimate"],
                                           set=[interval_text(c) for c in v["confidence_set"]]) for k, v in sens.items()})
            res["models"][f"{kind}_{support}"]["_prior"] = {t: prior[t] for t in ("ab", "Tk", "g", "lam_plus_g", "r") if t in prior}
            res["models"][f"{kind}_{support}"]["_post"] = {t: post[t] for t in ("ab", "Tk", "g", "lam_plus_g", "r") if t in post}
            res["models"][f"{kind}_{support}"]["_profiles"] = {t: prof[t]["curve"] for t in prof if "__" not in t}
    return res


FIG_TARGETS = [("tumor_matched", "ab", "Tumor α/β (Gy)", True), ("tumor_matched", "Tk", "Tumor onset Tk (days)", False),
               ("tumor_matched", "g", "Background growth g (per day)", False),
               ("tumor_matched", "lam_plus_g", "λ + g (per day)", False),
               ("early_early_matched", "ab", "Mucosal α/β (Gy)", True), ("early_early_matched", "r", "Time-loss rate r (Gy BED per day)", False)]


def figures(res, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    plt.style.use(os.path.expanduser("~/.claude/figure-style/publication.mplstyle"))
    os.makedirs(out, exist_ok=True)
    fig, axes = plt.subplots(3, 2, figsize=(12, 13))
    for ax, (mk, t, label, logx) in zip(axes.flat, FIG_TARGETS):
        mm = res["models"][mk]
        curve = np.array([c[:2] for c in mm["_profiles"][t]], float)
        ll_hat = mm["mle"]["ll"]
        if logx:
            curve = curve[curve[:, 0] > 0]
        ax.plot(curve[:, 0], curve[:, 1] - ll_hat, "-o", color="#0072B2", ms=3, lw=2, label="Profile log-likelihood")
        ax.axhline(-FM.CUTOFF, color="black", ls="--", lw=1.2, label="95% cutoff (−1.92)")
        row = next(r for r in mm["rows"] if r["target"] == t)
        lo, med, hi = row["posterior"]
        ax.axvspan(lo, hi, color="#E69F00", alpha=0.25, label="Posterior 95% interval")
        ax.axvline(med, color="#E69F00", lw=2)
        # Prior density, drawn on the same axes rescaled to the plotted range (0 at -8, peak at 0):
        # exact where the declared prior has a closed form, a kernel estimate of the simulated joint
        # prior otherwise (lambda + g).
        from scipy import stats
        xs = np.geomspace(max(curve[:, 0].min(), 1e-6), curve[:, 0].max(), 600) if logx else np.linspace(curve[:, 0].min(), curve[:, 0].max(), 600)
        exact = {("tumor_matched", "ab"): stats.lognorm(0.354, scale=10).pdf,
                 ("early_early_matched", "ab"): stats.lognorm(0.354, scale=10).pdf,
                 ("early_early_matched", "r"): stats.lognorm(0.5, scale=0.8).pdf,
                 ("tumor_matched", "g"): stats.uniform(0.0046, 0.021 - 0.0046).pdf,
                 ("tumor_matched", "Tk"): stats.truncnorm((0 - 28) / 10, (60 - 28) / 10, loc=28, scale=10).pdf}
        if (mk, t) in exact:
            dens = exact[(mk, t)](xs)
        else:
            dens = stats.gaussian_kde(np.asarray(mm["_prior"][t])[::20])(xs)
        ax.plot(xs, -8 + 8 * dens / dens.max(), color="#009E73", lw=2)
        ax.set_ylim(-8, 0.5)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5)) if not logx else None
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("Log-likelihood drop")
    h1, l1 = axes.flat[0].get_legend_handles_labels()
    fig.legend(h1 + [plt.Line2D([], [], color="#009E73", lw=2)], l1 + ["Prior density (rescaled, peak at 0)"], loc="lower center", ncol=2,
               frameon=False, fontsize=13)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    p = os.path.join(out, "frequentist_profiles.png")
    fig.savefig(p, dpi=600, bbox_inches="tight")
    return p


def main():
    res = assemble()
    out = os.path.join(FM.OUT, "comparison")
    os.makedirs(out, exist_ok=True)
    fig = figures(res, os.path.join(out, "figures"))
    for mm in res["models"].values():
        for k in ("_prior", "_post", "_profiles"):
            mm.pop(k)
    res["figure"] = os.path.relpath(fig, ROOT)
    json.dump(res, open(os.path.join(out, "comparison.json"), "w"), indent=1, default=float)
    print("wrote", os.path.join(out, "comparison.json"), fig)


if __name__ == "__main__":
    main()
