"""RETIRED 2026-10-26 (review R01): this first repair judged fits by a gate without ESS. It is
kept as the record of the repair of 2026-10-26 02:12 (repair.json in each folder) and refuses to
run. Use fit/repair_dd_ess.py, which applies the full gate of fit/gate.py.

Re-fit the calibrated-DD fits that failed the diagnostic gate, with longer chains.

fit_dd_calibrated.py tries each DI draw up to 3 times (chain length x1, x2, x4). A few fits still
end with R-hat just above 1.01 (slow mixing of the overdispersion sigma_raw; no divergences). This
script continues the same attempt sequence for those fits only: attempts 4 and 5 (x8, x16), with
the same seed rule, data, Stan model and gate. The gate is not changed (plan decision 9).

The runner itself is left unchanged so that identity.json of the finished folders stays true. The
repair refuses unless the folder's identity matches the current data, Stan model, runner and parent
DI fit. A replaced file is kept in repaired_originals/, every repair is listed in repair.json, and
summary.json is recomputed.

Run with the same environment as the original step:
  BED_FIT_DATA=data/X.csv BED_FIT_OUT=fit/output_X DDCAL_NDI=1000 DDCAL_DI_VARIANT=primary \
    python3 -I fit/repair_dd_calibrated.py
"""
import hashlib, json, os, shutil, site, sys
if site.getusersitepackages() not in sys.path:
    sys.path.append(site.getusersitepackages())
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_dd_calibrated as C
F = C.F

FIRST_ATTEMPT, LAST_ATTEMPT = C.MAX_ATTEMPTS + 1, C.MAX_ATTEMPTS + 2


def folder():
    ndi = int(os.environ.get("DDCAL_NDI", "200"))
    onset_min = float(os.environ.get("DDCAL_ONSET_MIN", "7"))
    v = os.environ.get("DDCAL_DI_VARIANT", "primary")
    out = os.path.join(F.OUT, f"dd_calibrated_onset{onset_min:g}_n{ndi}" + ("" if v == "primary" else f"_{v}"))
    return out, ndi, onset_min, v


def check_identity(out, di_dir):
    ident = json.load(open(os.path.join(out, "identity.json")))
    di_ident = json.load(open(os.path.join(di_dir, "identity.json")))
    now = dict(parent_di_identity=di_ident["identity_sha256"],
               parent_di_draws_sha256=F.sha_file(os.path.join(di_dir, "draws_params.csv")),
               data_file_sha256=F.sha_file(F.DATA),
               stan_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "dd_cal.stan")),
               common_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "common_functions.stanfunctions")),
               runner_sha256=F.sha_file(C.__file__), fit_refit_sha256=F.sha_file(F.__file__))
    bad = [k for k, v in now.items() if ident.get(k) != v]
    if bad:
        sys.exit(f"REFUSED: {out} does not match the current {', '.join(bad)}.")
    return ident


def fit_draw(m, r, base, onset_min, idx, out):
    """Attempts 4.. of the runner's loop, same data, seeds and gate. Returns (DataFrame, diagnostics, attempt, ok)."""
    n0, d0, T0 = C.ANCHOR
    data = dict(base, alpha=r.alpha, ab=r.ab, g=r.g, R0=float(C.rep_di(r["lambda"], r.Tk, T0)),
                K0=C.bed(n0, d0, r.ab), T0=T0, onset_min=onset_min, onset_scale=C.ONSET_SCALE,
                prior_log_lambda_mean=C.PRIOR_LOG_LAMBDA[0], prior_log_lambda_sd=C.PRIOR_LOG_LAMBDA[1])
    for attempt in range(FIRST_ATTEMPT, LAST_ATTEMPT + 1):
        mult = 2 ** (attempt - 1)
        seed = F.SEED + int(idx) + 1_000_000 * (attempt - 1)
        rng = np.random.default_rng(seed)
        inits = [dict(lambda_dd=0.5 * np.exp(rng.normal(0, .3)), c=list(rng.normal(0, .2, base["J"])),
                      sigma_raw=0.2, z=list(rng.normal(0, .1, base["N"]))) for _ in range(C.INNER["chains"])]
        f = m.sample(data=data, chains=C.INNER["chains"], parallel_chains=C.INNER["chains"],
                     iter_warmup=C.INNER["warmup"] * mult, iter_sampling=C.INNER["samples"] * mult,
                     seed=seed, adapt_delta=0.99, max_treedepth=12, show_progress=False,
                     output_dir=os.path.join(F.RAW, "dd_cal", os.path.basename(F.OUT), os.path.basename(out), f"{idx}_{attempt}"),
                     inits=inits)
        dg = F.diagnostics(f, ["lambda_dd", "sigma_raw"])
        ok = dg["divergences"] == 0 and dg["treedepth_hits"] == 0 and dg["max_rhat"] <= F.GATE["max_rhat"] \
            and dg["min_ebfmi"] >= F.GATE["min_ebfmi"]
        C.log(f"  repair DI draw {idx} attempt {attempt}: div {dg['divergences']} rhat {dg['max_rhat']:.4f} "
              f"({dg['worst_rhat']}) ebfmi {dg['min_ebfmi']:.2f} {'PASS' if ok else 'FAIL'}")
        if ok:
            break
    dr = f.draws_pd()
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
    for j in range(1, len(C.REGIMENS) + 1):
        res[f"dd_delta_{j}"] = dr[f"delta_eval[{j}]"].to_numpy()
        res[f"dd_onset_{j}"] = dr[f"onset_eval[{j}]"].to_numpy()
    return res, dg, attempt, ok


def main():
    sys.exit("REFUSED: fit/repair_dd_calibrated.py is retired (its gate omits ESS). Use fit/repair_dd_ess.py.")
    out, ndi, onset_min, v = folder()
    di_dir = os.path.join(F.OUT, f"di_{v}")
    ident = check_identity(out, di_dir)
    rows, trials, base = F.load()
    di = pd.read_csv(os.path.join(di_dir, "draws_params.csv"))
    pick = np.linspace(0, len(di) - 1, ndi).round().astype(int)
    if hashlib.sha256(pick.tobytes()).hexdigest() != ident["picks_sha256"]:
        sys.exit("REFUSED: DI draw picks differ from the folder's identity.")
    failing = []
    for idx in pick:
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        if not os.path.exists(path):
            sys.exit(f"REFUSED: {path} is missing; finish fit_dd_calibrated.py first.")
        head = pd.read_csv(path, nrows=1)
        if not bool(head.gate_pass.iloc[0]):
            failing.append(int(idx))
    rep_path = os.path.join(out, "repair.json")
    log_ = json.load(open(rep_path)) if os.path.exists(rep_path) else dict(
        script="fit/repair_dd_calibrated.py", gate="unchanged: no divergence or tree-depth hit, R-hat <= 1.01, E-BFMI >= 0.3",
        attempts=[FIRST_ATTEMPT, LAST_ATTEMPT], chain_multipliers=[2 ** (FIRST_ATTEMPT - 1), 2 ** (LAST_ATTEMPT - 1)],
        repairs=[])
    log_["repair_script_sha256"] = F.sha_file(os.path.abspath(__file__))
    C.log(f"repair {out}: {len(failing)} fits fail the gate")
    if failing:
        m = F.model("dd_cal")
        os.makedirs(os.path.join(out, "repaired_originals"), exist_ok=True)
    for idx in failing:
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        before = pd.read_csv(path, nrows=1)
        res, dg, attempt, ok = fit_draw(m, di.iloc[idx], base, onset_min, idx, out)
        keep = os.path.join(out, "repaired_originals", os.path.basename(path))
        if not os.path.exists(keep):
            shutil.copy2(path, keep)
        res.to_csv(path, index=False, float_format="%.10g")
        log_["repairs"].append(dict(di_draw=idx, rhat_before=float(before.max_rhat.iloc[0]),
                                    attempt_before=int(before.attempt.iloc[0]), rhat_after=dg["max_rhat"],
                                    worst_after=dg["worst_rhat"], attempt_after=attempt, gate_pass=bool(ok)))
        json.dump(log_, open(rep_path, "w"), indent=1)
    C.summarize(di, pick, out, onset_min, base)


if __name__ == "__main__":
    main()
