"""Write the tumor draws the app reads (draws/) and their manifest, from the verified-trial refit.

Version 3.1 (plan PLAN_codex_review_2026_10_25.md, stage E). Source: the corrected main fit
(fit/output_refit_2026_10_25c_main). For each variant:
- tumor_di_<variant>.csv: DI posterior draws of the refit.
- tumor_dd_<variant>.csv: calibrated DD, fitted conditional on the DI draw in the SAME ROW
  (modular coupling: DD never reweights DI). Row i of the two files is a matched pair by
  construction, not by file order.
Each of N_OUTER DI draws contributes K inner DD draws (evenly spaced), so the files have
N_OUTER x K rows, and the DI file repeats each DI draw K times.

Variants: "primary" (alpha/beta and g estimated) and "paper_restrictions" (alpha/beta fixed at
10 Gy and g at 0.0154/day, as in the 2018 paper). Units: alpha (1/Gy), ab (Gy), lambda (1/day),
Tk (inclusive days, DI) or c_over_alpha (Gy, DD), g (1/day).

The export refuses fits that fail the release gate (fit/gate.py, plan decision 9): the DI fit
must pass, and every conditional DD fit must pass, judged from each fit's own saved diagnostics.

Identity chain (review R02, 2026-10-26). Before anything is written, the export checks the data
file bytes against the data hash in both fit identities, the parent DI draw bytes against the
hash the DD fits recorded, the parent identity link, the variant, the outer-draw selection and
the Stan sources; every conditional file against all five parent values (alpha, alpha/beta, g,
lambda, Tk); finiteness and domain of every exported value; and DD = DI repopulation at the
anchor for every exported pair (tolerance ANCHOR_TOL after CSV round trip). All files are
written to a staging folder and replace draws/ only when the complete set has passed.

Run from the repository root after the fits:  python3 -I fit/export_draws.py
Env: BED_EXPORT_FROM (fit folder), BED_EXPORT_DD_N (outer draws, default 1000).
"""
import hashlib, json, os, site, sys
if site.getusersitepackages() not in sys.path:   # python3 -I drops the user site, where numpy lives
    sys.path.append(site.getusersitepackages())
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "fit"))
import fit_dd_calibrated as C
import gate as G

SRC = os.environ.get("BED_EXPORT_FROM", os.path.join(ROOT, "fit", "output_refit_2026_10_25c_main"))
N_OUTER = int(os.environ.get("BED_EXPORT_DD_N", "1000"))
K = 4
DRAWS = os.environ.get("BED_EXPORT_TO", os.path.join(ROOT, "draws"))
DRAWS_VERSION = "3.2"
ANCHOR_TOL = 1e-4   # Gy, after CSV round trip (actual residuals are recorded in the manifest)
PARENT_COLS = ("alpha", "ab", "g", "lambda", "Tk")
VARIANTS = {"primary": "alpha/beta and g estimated",
            "paper_restrictions": "alpha/beta fixed at 10 Gy and g at 0.0154 per day, as in the 2018 paper"}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def dd_dir(variant):
    return os.path.join(SRC, f"dd_calibrated_onset7_n{N_OUTER}" + ("" if variant == "primary" else f"_{variant}"))


def refuse(problems):
    sys.exit("REFUSED:\n" + "\n".join(problems))


def check_identity_chain(variant, di_id, dd_id, picks):
    """Every recorded input hash must match the bytes on disk now (review R02)."""
    di_file = os.path.join(SRC, f"di_{variant}", "draws_params.csv")
    data = os.path.join(ROOT, "data", di_id["data_file"])
    stan_dir = os.path.join(ROOT, "fit", "stan")
    p = []
    if not os.path.exists(data):
        p.append(f"data file {di_id['data_file']} is missing")
    else:
        h = sha(data)
        if h != di_id["data_file_sha256"]:
            p.append(f"DI {variant}: data file {di_id['data_file']} differs from the one fitted")
        if h != dd_id.get("data_file_sha256"):
            p.append(f"calibrated DD {variant}: data file differs from the one fitted")
    if dd_id.get("parent_di_identity") != di_id["identity_sha256"]:
        p.append(f"calibrated DD {variant} was fitted on another DI fit")
    if dd_id.get("parent_di_draws_sha256") != sha(di_file):
        p.append(f"calibrated DD {variant}: parent DI draws differ from those it was fitted on")
    if dd_id.get("di_variant") != variant or di_id.get("variant") != variant:
        p.append(f"{variant}: variant label mismatch")
    if dd_id.get("picks_sha256") != hashlib.sha256(picks.tobytes()).hexdigest():
        p.append(f"calibrated DD {variant}: outer DI draw selection differs")
    for f, h in di_id["stan_sha256"].items():
        if sha(os.path.join(stan_dir, f)) != h:
            p.append(f"DI {variant}: Stan source {f} differs from the one fitted")
    if sha(os.path.join(stan_dir, "dd_cal.stan")) != dd_id.get("stan_sha256") or \
            sha(os.path.join(stan_dir, "common_functions.stanfunctions")) != dd_id.get("common_sha256"):
        p.append(f"calibrated DD {variant}: Stan source differs from the one fitted")
    return p


def anchor_residual(di_rows, dd_rows):
    """|DD - DI| repopulation at the anchor for every exported pair (Gy)."""
    n0, d0, T0 = C.ANCHOR
    r_di = C.rep_di(di_rows["lambda"].to_numpy(), di_rows.Tk.to_numpy(), T0)
    r_dd = C.rep_dd(dd_rows.alpha.to_numpy(), dd_rows.ab.to_numpy(), dd_rows["lambda"].to_numpy(),
                    dd_rows.c_over_alpha.to_numpy(), n0, d0, T0)
    return np.abs(r_dd - r_di)


def export_variant(variant, stage):
    di_sum = json.load(open(os.path.join(SRC, f"di_{variant}", "summary.json")))
    di_id = json.load(open(os.path.join(SRC, f"di_{variant}", "identity.json")))
    dd_sum = json.load(open(os.path.join(dd_dir(variant), "summary.json")))
    dd_id = json.load(open(os.path.join(dd_dir(variant), "identity.json")))
    di = pd.read_csv(os.path.join(SRC, f"di_{variant}", "draws_params.csv")).set_index("draw")
    picks = np.linspace(0, len(di) - 1, N_OUTER).round().astype(int)
    problems = check_identity_chain(variant, di_id, dd_id, picks)
    di_fail = G.failures(di_sum)
    if di_fail:
        problems.append(f"DI {variant} fails the gate: {'; '.join(di_fail)}")
    if problems:
        refuse(problems)
    rows_di, rows_dd, dd_fail, fit_diag = [], [], [], []
    for idx in picks:
        f = pd.read_csv(os.path.join(dd_dir(variant), f"di_draw_{idx:05d}.csv"))
        dg = C.fit_diagnostics(f); fit_diag.append(dg)
        reasons = G.failures(dg)
        if reasons:
            dd_fail.append(f"outer draw {idx}: {'; '.join(reasons)}")
        p = di.loc[idx]
        for c in PARENT_COLS:
            if not np.allclose(f[c].to_numpy(), p[c], rtol=0, atol=1e-12):
                problems.append(f"calibrated DD {variant} outer draw {idx}: {c} differs from the parent DI draw")
        inner = f.iloc[np.linspace(0, len(f) - 1, K).round().astype(int)]
        for _, q in inner.iterrows():
            rows_di.append(dict(alpha=p.alpha, ab=p.ab, **{"lambda": p["lambda"]}, Tk=p.Tk, g=p.g, di_draw=idx))
            rows_dd.append(dict(alpha=p.alpha, ab=p.ab, **{"lambda": q.lambda_dd}, c_over_alpha=q.c_over_alpha, g=p.g,
                                di_draw=idx, dd_chain=int(q["chain"]), dd_iter=int(q["iter"])))
    if dd_fail:
        problems.append(f"calibrated DD {variant}: {len(dd_fail)} of {len(picks)} fits fail the gate, e.g. " + " | ".join(dd_fail[:5]))
    if problems:
        refuse(problems)
    out, resid = {}, {}
    frames = dict(di=pd.DataFrame(rows_di), dd=pd.DataFrame(rows_dd))
    for model, d in frames.items():
        num = d.drop(columns=[c for c in ("di_draw", "dd_chain", "dd_iter") if c in d])
        pos = ["alpha", "ab", "lambda"] + (["Tk"] if model == "di" else ["c_over_alpha"])
        if not np.isfinite(num.to_numpy()).all() or (num[pos] <= 0).any().any() or (num["g"] < 0).any():
            refuse([f"{model} {variant}: non-finite, non-positive ({', '.join(pos)}) or negative g values in the export"])
        name = f"tumor_{model}_{variant}.csv"
        d.to_csv(os.path.join(stage, name), index=False, float_format="%.10g", lineterminator="\n")
    # Round trip: re-read what was written and check pairing and the anchor on those bytes.
    rd = {m: pd.read_csv(os.path.join(stage, f"tumor_{m}_{variant}.csv")) for m in ("di", "dd")}
    if len(rd["di"]) != N_OUTER * K or not (rd["di"].di_draw.to_numpy() == rd["dd"].di_draw.to_numpy()).all():
        refuse([f"{variant}: exported DI and DD rows are not paired"])
    res = anchor_residual(rd["di"], rd["dd"])
    if not res.max() <= ANCHOR_TOL:
        refuse([f"{variant}: DD differs from DI at the anchor by up to {res.max():.3g} Gy (tolerance {ANCHOR_TOL:g})"])
    for model in ("di", "dd"):
        name = f"tumor_{model}_{variant}.csv"
        out[name] = dict(sha256=sha(os.path.join(stage, name)), n_draws=N_OUTER * K, n_distinct_di_draws=int(len(picks)),
                         variant=variant, model=model, variant_label=VARIANTS[variant])
    di_diag = {k: di_sum[k] for k in ("divergences", "treedepth_hits", "max_rhat", "worst_rhat", "min_ess_bulk",
                                      "min_ess_tail", "min_ebfmi", "n_quantities", "ppc_inside_90")}
    di_diag.update(gate_pass=True, gate_policy=G.POLICY["version"])
    fd = pd.DataFrame(fit_diag)
    rep = os.path.join(dd_dir(variant), "repair.json"); rep2 = os.path.join(dd_dir(variant), "repair_ess.json")
    dd_diag = dict(divergences=int(fd.divergences.sum()), treedepth_hits=int(fd.treedepth_hits.sum()),
                   max_rhat=float(fd.max_rhat.max()), min_ess_bulk_per_fit=float(fd.min_ess_bulk.min()),
                   min_ess_tail_per_fit=float(fd.min_ess_tail.min()), min_ebfmi=float(fd.min_ebfmi.min()),
                   outer_fits_failing_gate=0, gate_pass=True, gate_policy=G.POLICY["version"],
                   outer_fits_repaired=len(json.load(open(rep))["repairs"]) if os.path.exists(rep) else 0,
                   outer_fits_repaired_ess=len(json.load(open(rep2))["repairs"]) if os.path.exists(rep2) else 0,
                   anchor_max_abs_residual_gy=float(res.max()))
    out[f"tumor_di_{variant}.csv"].update(diagnostics=di_diag, fit_identity=di_id["identity_sha256"],
                                          min_ess_bulk=di_sum["min_ess_bulk"])
    out[f"tumor_dd_{variant}.csv"].update(diagnostics=dd_diag, fit_identity=dd_id["identity_sha256"],
                                          parent_di_identity=dd_id["parent_di_identity"],
                                          min_ess_bulk=dd_diag["min_ess_bulk_per_fit"],
                                          fitted_onset_below_7=dd_sum["fitted_onset_below_min"])
    return out, di_id


def point_values(stage):
    """Deterministic "fitted values" preset: medians of the identified DI quantities (alpha,
    lambda/alpha, g/alpha, Tk) with lambda and g rebuilt from them; alpha/beta stays 10 Gy because
    the data do not identify it (posterior close to prior). DD: the pooled median of lambda_DD,
    calibrated at the anchor with these DI values, which is the calibrated-DD model at point values."""
    dr = pd.read_csv(os.path.join(stage, "tumor_di_primary.csv")).drop_duplicates("di_draw")
    a = float(dr["alpha"].median())
    di = dict(alpha=a, ab=10.0, **{"lambda": float((dr["lambda"] / dr["alpha"]).median()) * a},
              Tk=float(dr["Tk"].median()), g=float((dr["g"] / dr["alpha"]).median()) * a,
              loss_per_day_after_onset=float((dr["lambda"] / dr["alpha"]).median()),
              growth_loss_per_day=float((dr["g"] / dr["alpha"]).median()))
    dd_draws = pd.read_csv(os.path.join(stage, "tumor_dd_primary.csv"))
    lam_dd = float(dd_draws["lambda"].median())
    onset, co = C.calibrate(di["alpha"], di["ab"], di["lambda"], di["Tk"], lam_dd)
    # lambda_over_alpha is lambda_DD / alpha, not a loss per day: the DD daily loss also carries the
    # regimen's kill probability p (review 2026-10-26), so it is not named as one.
    dd = dict(alpha=a, ab=10.0, **{"lambda": lam_dd}, c_over_alpha=float(co), onset_at_anchor_days=float(onset),
              g=di["g"], lambda_over_alpha=lam_dd / a, growth_loss_per_day=di["growth_loss_per_day"])
    return dict(di=di, dd=dd)


def main():
    import shutil
    stage = DRAWS.rstrip("/\\") + ".staging"
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(stage)
    files, ident = {}, None
    for v in VARIANTS:
        f, i = export_variant(v, stage)
        files.update(f)
        ident = ident or i
    man = dict(created_by="fit/export_draws.py", draws_version=DRAWS_VERSION,
               analysis="verified-trial refit 2026-10-25 (corrected main set): DI posterior; DD calibrated to DI "
                        "at 35 x 2 Gy in 47 days, fitted conditional on each DI draw",
               dataset=ident["data_file"], data_sha256=ident["data_file_sha256"],  # verified against the file above
               gate_policy=G.POLICY,
               pairing="row i of tumor_dd_<variant>.csv was fitted conditional on row i of tumor_di_<variant>.csv",
               toxicity_pairing="toxicity prior draws are independent of the tumor draws (assumption)",
               likelihood="binomial on rounded N x reported rate (approximation: most rates are KM or CIF estimates)",
               onset_penalty="soft preference for DD onset at or after day 7, not a bound",
               mcse_note="draws are not independent: each DI draw is repeated in K rows; Monte Carlo error "
                         "is approximate (IID benchmark on distinct DI draws)",
               source_fit=os.path.relpath(SRC, ROOT), files=files)
    json.dump(man, open(os.path.join(stage, "draws_manifest.json"), "w"), indent=1)
    pv = point_values(stage)
    json.dump(dict(source="draws/tumor_*_primary.csv", dataset=ident["data_file"], rule="medians of identified quantities; "
                   "DD lambda median calibrated at the anchor with the DI point values", values=pv),
              open(os.path.join(stage, "refit_point_values.json"), "w"), indent=1)
    # The complete set has passed: replace draws/ in one step (the old folder is removed after).
    old = DRAWS.rstrip("/\\") + ".previous"
    shutil.rmtree(old, ignore_errors=True)
    if os.path.exists(DRAWS):
        os.rename(DRAWS, old)
    os.rename(stage, DRAWS)
    shutil.rmtree(old, ignore_errors=True)
    for k, v in files.items():
        print(k, v["n_draws"], v["sha256"][:12], "anchor residual", v["diagnostics"].get("anchor_max_abs_residual_gy", ""))
    print(json.dumps(pv, indent=1))


if __name__ == "__main__":
    main()
