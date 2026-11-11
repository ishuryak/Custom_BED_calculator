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
Score criterion (plan step 6): the lower 95% bound of the weighted mean d is above 0 in both seeds and
no single family reverses the sign. numerical_rule_pass also requires the full family scope and every
fit passing the complete-schema gate (review F4/F5). It is not a promotion decision: `decision`
stays with Igor Shuryak.

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


def family_data(fam, rows_all, settings=None):
    """Contrast terms of one family's arms against its own reference arm (repair as in settings)."""
    rows = [x for x in rows_all if x["family"] == fam]
    d = s5.stan_data(rows, [fam], TK, {k: v for k, v in (settings or {}).items() if k == "halftime_h"})
    ref = [i for i, x in enumerate(rows) if x["arm_id"] == s5.REF_ARM[fam]]
    assert len(ref) == 1
    return dict(rows=rows, ref=ref[0], y=np.array(d["y"]), N=np.array(d["N"]), dD=np.array(d["dD"]),
                dQ=np.array(d["dQ"]), dTe=np.array(d["dTe"]), chemo=np.array(d["chemo"]),
                group=rows[0]["scale_group"])


def arm_probabilities(draws, fd, eta, arm):
    """p (S, inner, I). An arm's chemotherapy may be a probability mixture fd['chemo_mix'][i] =
    [(weight, code), ...], mixed on the probability scale before the binomial (review F1)."""
    dE = fd["dD"][None, :] + fd["dQ"][None, :] / draws["ab"][:, None] - draws["r"][:, None] * fd["dTe"][None, :]
    base = eta[:, :, None] + (draws["kappa"][:, None] * dE)[:, None, :] + arm
    mix = fd.get("chemo_mix") or [[(1.0, c)] for c in fd["chemo"]]
    p = np.zeros_like(base)
    for i, comps in enumerate(mix):
        if abs(sum(w for w, _ in comps) - 1) > 1e-9:
            raise ValueError("mixture weights must sum to 1")
        for w, c in comps:
            p[:, :, i] += w / (1 + np.exp(-(base[:, :, i] + draws["gamma"][:, None] * c)))
    return p


def latent_draws(draws, fd, inner, rng):
    S, I = len(draws["kappa"]), len(fd["y"])
    eta = (draws["m0"] + draws["delta"])[:, None] + draws["tau0"][:, None] * rng.standard_normal((S, inner))
    arm = draws["sigma_arm"][:, None, None] * rng.standard_normal((S, inner, I))
    return eta, arm


def predictive(draws, fd, inner, rng):
    """(log p(all arms), log p(reference arm)) by Monte Carlo over draws x inner baselines."""
    eta, arm = latent_draws(draws, fd, inner, rng)
    p = arm_probabilities(draws, fd, eta, arm)
    ll = binom_logpmf(fd["y"][None, None, :], fd["N"][None, None, :], p)   # (S, inner, I)
    tot = ll.shape[0] * ll.shape[1]
    return float(logsumexp(ll.sum(axis=2)) - math.log(tot)), float(logsumexp(ll[:, :, fd["ref"]]) - math.log(tot))


def group_delta(fit_draws, group, train_groups, rng):
    """delta for the held-out family's scale group: fitted if the group is in training, else prior N(0, 1)."""
    S = len(fit_draws["kappa"])
    if group == "G1":
        return np.zeros(S)
    if group in train_groups:
        return fit_draws["delta_all"][:, train_groups.index(group)]
    return rng.standard_normal(S)


STEP6_SAMPLER = dict(s5.SAMPLER, seed=SEED)


def fit_cached(model, name, rows, fams, settings, cache, rescore_only=False, accept_legacy=False):
    """Posterior draws for one fit, reused only on an identical run identity (review F4). In
    rescore-only mode a missing cache is an error and nothing is sampled."""
    data = s5.stan_data(rows, fams, TK, settings)
    identity = s5.run_identity(data, settings, STEP6_SAMPLER, extra=dict(step6=name, families=fams))
    path = os.path.join(cache, f"{name}.npz")
    if os.path.exists(path):
        try:
            z = np.load(path, allow_pickle=False)
            out = {k: z[k] for k in z.files}
        except Exception as err:                                        # truncated or corrupt file
            raise RuntimeError(f"{path}: unreadable cache ({err}); refusing to reuse")
        if "identity" not in out:
            if not accept_legacy:
                raise RuntimeError(f"{path}: cache has no run identity (legacy); refusing to reuse")
            if sorted(out["families"].tolist()) != sorted(fams):
                raise RuntimeError(f"{path}: legacy cache family set differs")
            out["gate_json"] = np.array(json.dumps(dict(gate_pass=bool(out["gate"][0]), scope="legacy narrow gate")))
            return out
        if str(out["identity"]) != identity:
            raise RuntimeError(f"{path}: run identity differs from the current inputs; refusing to reuse")
        return out
    if rescore_only:
        raise RuntimeError(f"{path}: no cache; rescore-only mode does not sample")
    s5.wait_for_memory()
    fit = model.sample(data=data, show_progress=False, show_console=False, **STEP6_SAMPLER)
    summ = fit.summary(percentiles=(5, 50, 95))
    ebfmi, depth = s5.chain_energy_checks(fit)
    gate = s5.gate_from_summary(summ, data, np.sum(fit.divergences), ebfmi, depth)
    gate["scope"] = "complete schema gate"
    S = len(fit.stan_variable("kappa"))
    out = {p: np.asarray(fit.stan_variable(p)).reshape(S, -1)[:, 0] if p in ("gamma", "sigma_arm")
           else np.asarray(fit.stan_variable(p)) for p in PARAMS}
    out["delta_all"] = np.asarray(fit.stan_variable("delta")).reshape(S, -1)
    out["families"] = np.array(fams)
    out["groups"] = np.array([g for g in s5.GROUP_COLS if any(x["scale_group"] == g for x in rows)])
    out["gate"] = np.array([gate["gate_pass"], gate["divergences"], gate["rhat"], gate["ess_bulk"], gate["ess_tail"]])
    out["gate_json"] = np.array(json.dumps(gate, default=str))
    out["chain_archive"] = np.array(json.dumps(s5.archive_chains(fit, identity)))
    out["identity"] = np.array(identity)                     # written last into the same atomic file
    tmp = path + f".partial-{os.getpid()}.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, path)
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


def cc_lor(y1, n1, y0, n0):
    """Continuity-corrected sample log odds ratio, the statistic compared on both scales."""
    return np.log((y1 + 0.5) / (n1 - y1 + 0.5)) - np.log((y0 + 0.5) / (n0 - y0 + 0.5))


def weighted_quantiles(x, w, qs):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return [float(x[o][min(np.searchsorted(cw, q), len(x) - 1)]) for q in qs]


def contrast_check(fd, draws_fit, seed, inner):
    """Conditional predictive check of each non-reference arm GIVEN the observed reference count
    (review F2). Baselines and arm effects are weighted by the reference-arm likelihood; the other
    arms' counts are simulated with their own denominators; the same continuity-corrected sample
    log OR is computed on replicates and data (observed scale). The latent log OR, also weighted, is
    reported separately. ESS of the weights is recorded."""
    rng = np.random.default_rng(seed)
    if "delta" not in draws_fit:                                           # fitted draws: attach the group offset
        draws_fit = dict(draws_fit, delta=group_delta(draws_fit, fd["group"], list(draws_fit["groups"]), rng))
    eta, arm = latent_draws(draws_fit, fd, inner, rng)
    p = arm_probabilities(draws_fit, fd, eta, arm)                         # (S, inner, I)
    y, N, r = fd["y"], fd["N"], fd["ref"]
    lw = binom_logpmf(y[r], N[r], p[:, :, r]).ravel()
    w = np.exp(lw - lw.max())
    ess = float(w.sum() ** 2 / np.sum(w ** 2))
    out = []
    for i, x in enumerate(fd["rows"]):
        if i == r:
            continue
        yrep = rng.binomial(N[i], p[:, :, i]).ravel()
        stat_rep = cc_lor(yrep, N[i], y[r], N[r])
        obs = float(cc_lor(y[i], N[i], y[r], N[r]))
        latent = (np.log(p[:, :, i] / (1 - p[:, :, i])) - np.log(p[:, :, r] / (1 - p[:, :, r]))).ravel()
        q = weighted_quantiles(stat_rep, w, [0.05, 0.5, 0.95])
        qc = weighted_quantiles(yrep.astype(float), w, [0.05, 0.5, 0.95])
        ql = weighted_quantiles(latent, w, [0.05, 0.5, 0.95])
        below = float(np.sum(w[stat_rep < obs]) / w.sum())
        equal = float(np.sum(w[stat_rep == obs]) / w.sum())
        out.append(dict(arm_id=x["arm_id"], observed_lor=round(obs, 4),
                        predicted_obs_lor_5=round(q[0], 4), predicted_obs_lor_median=round(q[1], 4),
                        predicted_obs_lor_95=round(q[2], 4), inside_90=bool(q[0] <= obs <= q[2]),
                        width_90=round(q[2] - q[0], 4), mid_pit=round(below + 0.5 * equal, 4),
                        observed_count=int(y[i]), predicted_count_5=qc[0], predicted_count_median=qc[1],
                        predicted_count_95=qc[2],
                        latent_lor_5=round(ql[0], 4), latent_lor_median=round(ql[1], 4), latent_lor_95=round(ql[2], 4),
                        target="observed-scale sample log OR given the observed reference count",
                        reference_arm=fd["rows"][r]["arm_id"], seed=seed, weight_ess=round(ess, 1)))
    return out


def gate_row(model, held_out, dr):
    g = json.loads(str(dr["gate_json"]))
    return dict(model=model, held_out=held_out, gate_pass=bool(g["gate_pass"]), scope=g.get("scope", ""),
                divergences=g.get("divergences", ""), rhat=g.get("rhat", ""), ess_bulk=g.get("ess_bulk", ""),
                ess_tail=g.get("ess_tail", ""), ebfmi_min=g.get("ebfmi_min", ""),
                missing=len(g.get("missing", [])), nonfinite=len(g.get("nonfinite", [])),
                identity=str(dr["identity"]) if "identity" in dr else "legacy (none)")


def with_chemo(fd, arm_id, code=None, mix=None):
    """A copy of a family's data with one arm's chemotherapy code or probability mixture replaced."""
    f = dict(fd)
    i = [x["arm_id"] for x in fd["rows"]].index(arm_id)
    if mix is not None:
        f["chemo_mix"] = [[(1.0, c)] for c in fd["chemo"]]
        f["chemo_mix"][i] = list(mix)
    else:
        ch = fd["chemo"].copy()
        ch[i] = code
        f["chemo"] = ch
    return f


def hypno_scenarios(fd):
    out = {}
    for name, spec in s5.EXTERNAL_CHEMO_SCENARIOS["HYPNO_HFX"].items():
        if isinstance(spec, tuple):
            _, share, code = spec
            out[name] = with_chemo(fd, "HYPNO_HFX", mix=[(share, code), (1 - share, 0.0)])
        else:
            out[name] = with_chemo(fd, "HYPNO_HFX", code=spec)
    return out


STRUCTURAL_HALFTIMES = (2.0, 4.0, s5.COMPLETE_REPAIR)


def structural_repair(fc, prim, fams, allrows, gates, inner):
    """Review F8 / plan D2: both models refitted on all primary families at other repair half-times;
    parameters, decision odds ratios and external checks per variant (no held-out refits)."""
    out = []
    for h in (s5.s4.HALFTIME_H,) + STRUCTURAL_HALFTIMES:
        for m, st in MODELS.items():
            stx = dict(st, halftime_h=h) if h != s5.s4.HALFTIME_H else st
            name = f"{m}__all_primary" + ("" if h == s5.s4.HALFTIME_H else f"__h{h}")
            dr = fc(name, prim, fams, stx)
            if h != s5.s4.HALFTIME_H:
                gates.append(gate_row(m, f"none (all primary, h={h})", dr))
            hh = s5.halftime_of(stx)
            ref70 = s5.schedule_terms(*s5.ev.regular_schedule(35, 2.0), TK, hh)
            row = dict(model=m, halftime_h=h, r_median=round(float(np.median(dr["r"])), 3),
                       ab_median=round(float(np.median(dr["ab"])), 2), kappa_median=round(float(np.median(dr["kappa"])), 4))
            for n, dd in s5.DECISION:
                t = np.subtract(s5.schedule_terms(*s5.ev.regular_schedule(n, dd), TK, hh), ref70)
                lor = dr["kappa"] * (t[0] + t[1] / dr["ab"] - dr["r"] * t[2])
                row[f"OR_{n}x{dd}"] = round(float(np.exp(np.median(lor))), 3)
            for j, (f, fd) in enumerate([(f"HYPNO [{k}]", v) for k, v in hypno_scenarios(external_family("HYPNO", stx)).items()]
                                        + [("ELAN", family_data("ELAN", allrows, stx))]):
                row[f"cond_{f}"] = round(score_family(fd, dr, 300, SEED + 9000 + j)[1], 4)
            out.append(row)
    return out


def numerical_rule(score_criterion_met, scope_complete, gates):
    """True only for the full family scope with every fit passing the complete-schema gate."""
    return bool(score_criterion_met and scope_complete and gates
                and all(g["gate_pass"] and g["scope"] == "complete schema gate" for g in gates))


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--inner", type=int, default=100)
    ap.add_argument("--families", default="", help="comma list to restrict held-out families (smoke; never eligible)")
    ap.add_argument("--rescore-only", action="store_true", help="never sample; a missing cache is an error")
    ap.add_argument("--accept-legacy", action="store_true", help="read caches written before run identities (labelled)")
    ap.add_argument("--structural", action="store_true",
                    help="also fit both models on all primary families at repair half-times 2 h, 4 h and complete")
    a = ap.parse_args(argv[1:])
    cache = os.path.join(a.out, "cache")
    os.makedirs(cache, exist_ok=True)
    model = None if a.rescore_only else s5.compile_model()
    fc = lambda name, rows, fams_, st: fit_cached(model, name, rows, fams_, st, cache, a.rescore_only, a.accept_legacy)
    prim, fams = s5.load_rows()
    allrows, allfams = s5.load_rows(m4=True)
    held = [f for f in fams if not a.families or f in a.families.split(",")]
    t0 = time.time()
    per_family, contrasts, gates, chemo_sens = [], [], [], []
    for k, f in enumerate(held):
        rows, tfams = training_set(prim, fams, f)
        fd = family_data(f, prim)
        rec = dict(family=f, n_arms=len(fd["y"]), group=fd["group"])
        for m, st in MODELS.items():
            dr = fc(f"{m}__without_{f}", rows, tfams, st)
            assert f not in set(dr["families"].tolist()), f"leak: {f} in its own training fit"
            gates.append(gate_row(m, f, dr))
            new, cond = score_family(fd, dr, a.inner, SEED + 1000 * k)      # common random numbers
            rec[f"{m}_new_trial"] = round(new, 4)
            rec[f"{m}_conditional"] = round(cond, 4)
            contrasts += [dict(model=m, family=f, **c) for c in contrast_check(fd, dr, SEED + k, a.inner)]
            for arm_id, code in s5.CHEMO_SENSITIVITY.items():                # decision D5 (held-out fold only)
                if any(x["arm_id"] == arm_id for x in fd["rows"]):
                    _, c2 = score_family(with_chemo(fd, arm_id, code=code), dr, a.inner, SEED + 1000 * k)
                    chemo_sens.append(dict(model=m, family=f, arm_id=arm_id, code=code, conditional=round(c2, 4),
                                           conditional_at_0=round(cond, 4)))
        rec["d_new_trial"] = round(rec["candidate_new_trial"] - rec["benchmark_new_trial"], 4)
        rec["d_conditional"] = round(rec["candidate_conditional"] - rec["benchmark_conditional"], 4)
        per_family.append(rec)
        print(f"held out {f} ({k + 1}/{len(held)}, {time.time() - t0:.0f} s): d_new {rec['d_new_trial']:+.3f}, "
              f"d_cond {rec['d_conditional']:+.3f}", flush=True)

    scope_complete = len(per_family) == len(fams) and not a.families
    summary = dict(seed=SEED, inner=a.inner, n_boot=N_BOOT, tk=TK, families=len(per_family),
                   scope="all primary families" if scope_complete else f"SUBSET {held} (smoke; never eligible)",
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
            score_criterion_met=bool(all(b["lo"] > 0 for b in boots.values()) and all(v > 0 for v in infl.values())))

    # external checks from fits on all primary families; HYPNO under declared chemotherapy scenarios
    ext = []
    for m, st in MODELS.items():
        dr = fc(f"{m}__all_primary", prim, fams, st)
        gates.append(gate_row(m, "none (all primary)", dr))
        hy = external_family("HYPNO")
        targets = [(f"HYPNO [{n}]", fdx) for n, fdx in hypno_scenarios(hy).items()]
        targets += [(f, family_data(f, allrows)) for f in sorted(s5.SENSITIVITY_FAMILIES)]
        for j, (f, fd) in enumerate(targets):
            new, cond = score_family(fd, dr, 300 if f.startswith("HYPNO") else a.inner, SEED + 5000 + j)
            cc = contrast_check(fd, dr, SEED + 7000 + j, a.inner)
            ext.append(dict(model=m, family=f, conditional=round(cond, 4), new_trial=round(new, 4),
                            arms=";".join(f"{c['arm_id']} obs {c['observed_lor']:+.2f} pred(obs) {c['predicted_obs_lor_median']:+.2f} "
                                          f"[{c['predicted_obs_lor_5']:+.2f},{c['predicted_obs_lor_95']:+.2f}] "
                                          f"latent {c['latent_lor_median']:+.2f}" for c in cc)))
    structural = []
    if a.structural:
        structural = structural_repair(fc, prim, fams, allrows, gates, a.inner)
        write(os.path.join(a.out, "structural_repair.csv"), structural)
    gates_ok = all(g["gate_pass"] for g in gates)
    gates_complete = all(g["scope"] == "complete schema gate" for g in gates)
    for task in summary["tasks"].values():
        task["numerical_rule_pass"] = numerical_rule(task["score_criterion_met"], scope_complete, gates)
    summary.update(external=ext, chemo_sensitivity=chemo_sens, gates_all_pass=gates_ok,
                   gates_complete_schema=gates_complete,
                   decision="undecided: Igor Shuryak decides after Phases A-D (Fowler stays the app default)")
    write(os.path.join(a.out, "per_family.csv"), per_family)
    write(os.path.join(a.out, "contrasts.csv"), contrasts)
    write(os.path.join(a.out, "fit_gates.csv"), gates)
    write(os.path.join(a.out, "external.csv"), ext)
    write(os.path.join(a.out, "chemo_sensitivity.csv"), chemo_sens)
    tmp = os.path.join(a.out, f"summary.json.partial-{os.getpid()}")
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=1)
    os.replace(tmp, os.path.join(a.out, "summary.json"))
    print(json.dumps({k: summary[k] for k in ("tasks", "gates_all_pass", "gates_complete_schema", "scope")}, indent=1))


def external_family(fam, settings=None):
    """HYPNO is the R6 benchmark: never fitted; its reference arm is the 66 Gy / 33 control."""
    rows = [x for x in s5_benchmark_rows() if x["family"] == fam]
    old = s5.REF_ARM.get(fam)
    s5.REF_ARM[fam] = "HYPNO_NFX"
    try:
        return family_data(fam, rows, settings)
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
