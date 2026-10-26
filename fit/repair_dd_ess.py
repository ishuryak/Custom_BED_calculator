"""Repair generation 2 (review R01, 2026-10-26): refit the calibrated-DD fits that fail the full
diagnostic gate of fit/gate.py, which includes bulk and tail ESS >= 400.

Until 2026-10-26 the calibrated-DD loop (fit_dd_calibrated.py) and the first repair
(repair_dd_calibrated.py) accepted a fit without checking ESS. This script:

1. Verifies the folder against the bytes it was fitted on: data file, parent DI draws, the
   parent DI identity, dd_cal.stan, common functions and the selected outer draws. Runner hashes
   are not compared with the current runner (the runner now carries the full gate); the folder's
   original identity is kept unchanged and referenced from repair_ess.json.
2. Reclassifies every outer fit from its saved diagnostics under fit/gate.py, never from the
   stored gate_pass column.
3. Refits only the failing fits with longer chains (4x, 8x, 16x the base length), seeds offset by
   SEED_OFFSET from every earlier attempt so no earlier chain is repeated. The replaced file is kept
   in repaired_originals_ess/. Each refit row records the gate policy version.
4. Recomputes summary.json with equal outer weights.

Diagnostics cover lambda_dd, sigma_raw, the trial intercepts, the arm effects and every regimen
contrast. Onsets are not listed separately: within one fit each onset is a monotone function of
lambda_dd alone, and rank-based bulk and tail ESS are unchanged by monotone transforms.

Run (one folder per process):
  BED_FIT_DATA=data/X.csv BED_FIT_OUT=fit/output_X DDCAL_NDI=1000 DDCAL_DI_VARIANT=primary \
    python3 -I fit/repair_dd_ess.py
"""
import hashlib, json, os, shutil, site, sys, time
if site.getusersitepackages() not in sys.path:
    sys.path.append(site.getusersitepackages())
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_dd_calibrated as C
import gate as G
F = C.F

MULTIPLIERS = (4, 8, 16)
SEED_OFFSET = 50_000_000
GENERATION = "repair_ess_2026_10_26"
DIAG_COLS = {"divergent", "treedepth_hits", "max_rhat", "min_ess_bulk", "min_ess_tail", "min_ebfmi"}


def folder():
    ndi = int(os.environ.get("DDCAL_NDI", "200"))
    onset_min = float(os.environ.get("DDCAL_ONSET_MIN", "7"))
    v = os.environ.get("DDCAL_DI_VARIANT", "primary")
    out = os.path.join(F.OUT, f"dd_calibrated_onset{onset_min:g}_n{ndi}" + ("" if v == "primary" else f"_{v}"))
    return out, ndi, onset_min, v


def check_inputs(out, di_dir, pick):
    """Refuse unless the folder was fitted on exactly these bytes. Returns the folder identity."""
    ident = json.load(open(os.path.join(out, "identity.json")))
    di_ident = json.load(open(os.path.join(di_dir, "identity.json")))
    now = dict(parent_di_identity=di_ident["identity_sha256"],
               parent_di_draws_sha256=F.sha_file(os.path.join(di_dir, "draws_params.csv")),
               data_file_sha256=F.sha_file(F.DATA),
               stan_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "dd_cal.stan")),
               common_sha256=F.sha_file(os.path.join(F.ROOT, "fit", "stan", "common_functions.stanfunctions")),
               picks_sha256=hashlib.sha256(pick.tobytes()).hexdigest())
    bad = [k for k, v in now.items() if ident.get(k) != v]
    if bad:
        sys.exit(f"REFUSED: {out} does not match the current {', '.join(bad)}.")
    if di_ident.get("data_file_sha256") != now["data_file_sha256"]:
        sys.exit(f"REFUSED: the parent DI fit {di_dir} was fitted on another data file.")
    return ident


def fit_draw(m, r, base, onset_min, idx, out):
    """Refit one outer draw at 4x, 8x, 16x length until it passes. Returns (rows, diagnostics, mult, ok)."""
    n0, d0, T0 = C.ANCHOR
    data = dict(base, alpha=r.alpha, ab=r.ab, g=r.g, R0=float(C.rep_di(r["lambda"], r.Tk, T0)),
                K0=C.bed(n0, d0, r.ab), T0=T0, onset_min=onset_min, onset_scale=C.ONSET_SCALE,
                prior_log_lambda_mean=C.PRIOR_LOG_LAMBDA[0], prior_log_lambda_sd=C.PRIOR_LOG_LAMBDA[1])
    for k, mult in enumerate(MULTIPLIERS):
        seed = F.SEED + int(idx) + SEED_OFFSET + 1_000_000 * k
        rng = np.random.default_rng(seed)
        inits = [dict(lambda_dd=0.5 * np.exp(rng.normal(0, .3)), c=list(rng.normal(0, .2, base["J"])),
                      sigma_raw=0.2, z=list(rng.normal(0, .1, base["N"]))) for _ in range(C.INNER["chains"])]
        f = m.sample(data=data, chains=C.INNER["chains"], parallel_chains=C.INNER["chains"],
                     iter_warmup=C.INNER["warmup"] * mult, iter_sampling=C.INNER["samples"] * mult,
                     seed=seed, adapt_delta=0.99, max_treedepth=12, show_progress=False,
                     output_dir=os.path.join(F.RAW, "dd_cal_ess", os.path.basename(F.OUT), os.path.basename(out), f"{idx}_x{mult}"),
                     inits=inits)
        dg = F.diagnostics(f, ["lambda_dd", "sigma_raw"])
        fails = G.failures(dg)
        C.log(f"  ess-repair DI draw {idx} x{mult}: {'PASS' if not fails else 'FAIL ' + '; '.join(fails)} "
              f"(ess bulk {dg['min_ess_bulk']:.0f}, tail {dg['min_ess_tail']:.0f}, rhat {dg['max_rhat']:.4f})")
        if not fails:
            break
    dr = f.draws_pd()
    err = float(np.max(np.abs(dr["rep_anchor"] - data["R0"])))
    assert err < 1e-6 * max(1.0, data["R0"]), f"calibration broken: {err}"
    res = pd.DataFrame(dict(di_draw=idx, alpha=r.alpha, ab=r.ab, g=r.g, **{"lambda": r["lambda"]}, Tk=r.Tk,
                            chain=dr["chain__"], iter=dr["iter__"], lambda_dd=dr["lambda_dd"],
                            onset_ref=dr["onset_ref"], c_over_alpha=dr["c_over_alpha"],
                            onset_min_fitted=dr["onset_min_fitted"],
                            loglik=dr[[f"log_lik[{i}]" for i in range(1, base["N"] + 1)]].sum(axis=1),
                            divergent=dr["divergent__"].to_numpy(), attempt=10 + k, chain_multiplier=mult, gate_pass=not fails,
                            gate_policy=G.POLICY["version"], seed=seed,
                            max_rhat=dg["max_rhat"], min_ess_bulk=dg["min_ess_bulk"],
                            min_ess_tail=dg["min_ess_tail"], min_ebfmi=dg["min_ebfmi"],
                            treedepth_hits=dg["treedepth_hits"]))
    for j in range(1, len(C.REGIMENS) + 1):
        res[f"dd_delta_{j}"] = dr[f"delta_eval[{j}]"].to_numpy()
        res[f"dd_onset_{j}"] = dr[f"onset_eval[{j}]"].to_numpy()
    return res, dg, mult, not fails


def main():
    out, ndi, onset_min, v = folder()
    di_dir = os.path.join(F.OUT, f"di_{v}")
    di = pd.read_csv(os.path.join(di_dir, "draws_params.csv"))
    pick = np.linspace(0, len(di) - 1, ndi).round().astype(int)
    ident = check_inputs(out, di_dir, pick)
    rows, trials, base = F.load()
    failing = []
    for idx in pick:
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        if not os.path.exists(path):
            sys.exit(f"REFUSED: {path} is missing.")
        reasons = G.failures(C.fit_diagnostics(pd.read_csv(path, usecols=lambda c: c in DIAG_COLS)))
        if reasons:
            failing.append((int(idx), reasons))
    rep_path = os.path.join(out, "repair_ess.json")
    log_ = json.load(open(rep_path)) if os.path.exists(rep_path) else dict(
        generation=GENERATION, gate_policy=G.POLICY, multipliers=list(MULTIPLIERS), seed_offset=SEED_OFFSET,
        folder_identity_sha256=ident["identity_sha256"], original_runner_sha256=ident.get("runner_sha256"),
        n_failing_at_start=len(failing), repairs=[])
    log_["script_sha256"] = F.sha_file(os.path.abspath(__file__))
    log_["runner_sha256_now"] = F.sha_file(C.__file__)
    C.log(f"ess-repair {out}: {len(failing)} of {len(pick)} fits fail the gate")
    if failing:
        m = F.model("dd_cal")
        os.makedirs(os.path.join(out, "repaired_originals_ess"), exist_ok=True)
    limit = int(os.environ.get("ESS_REPAIR_LIMIT", "0"))   # smoke tests only: repair at most this many
    if limit:
        C.log(f"SMOKE: ESS_REPAIR_LIMIT={limit}")
        failing = failing[:limit]
    for k, (idx, reasons) in enumerate(failing):
        path = os.path.join(out, f"di_draw_{idx:05d}.csv")
        before = pd.read_csv(path, nrows=1).iloc[0]
        res, dg, mult, ok = fit_draw(m, di.iloc[idx], base, onset_min, idx, out)
        keep = os.path.join(out, "repaired_originals_ess", os.path.basename(path))
        if not os.path.exists(keep):
            shutil.copy2(path, keep)
        res.to_csv(path, index=False, float_format="%.10g")
        log_["repairs"].append(dict(di_draw=idx, reasons_before=reasons, original_sha256=F.sha_file(keep),
                                    ess_bulk_before=float(before.min_ess_bulk), ess_tail_before=float(before.min_ess_tail),
                                    ess_bulk_after=dg["min_ess_bulk"], ess_tail_after=dg["min_ess_tail"],
                                    rhat_after=dg["max_rhat"], multiplier=mult, gate_pass=bool(ok)))
        json.dump(log_, open(rep_path, "w"), indent=1)
        C.log(f"ess-repair progress {k + 1}/{len(failing)}")
    C.summarize(di, pick, out, onset_min, base)


if __name__ == "__main__":
    main()
