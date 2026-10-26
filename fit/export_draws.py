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

The export refuses fits that fail the release gate (plan decision 9): the DI fit must pass, and
every conditional DD fit must pass with no divergence.

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

SRC = os.environ.get("BED_EXPORT_FROM", os.path.join(ROOT, "fit", "output_refit_2026_10_25c_main"))
N_OUTER = int(os.environ.get("BED_EXPORT_DD_N", "1000"))
K = 4
DRAWS = os.environ.get("BED_EXPORT_TO", os.path.join(ROOT, "draws"))
VARIANTS = {"primary": "alpha/beta and g estimated",
            "paper_restrictions": "alpha/beta fixed at 10 Gy and g at 0.0154 per day, as in the 2018 paper"}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def dd_dir(variant):
    return os.path.join(SRC, f"dd_calibrated_onset7_n{N_OUTER}" + ("" if variant == "primary" else f"_{variant}"))


def export_variant(variant):
    di_sum = json.load(open(os.path.join(SRC, f"di_{variant}", "summary.json")))
    di_id = json.load(open(os.path.join(SRC, f"di_{variant}", "identity.json")))
    dd_sum = json.load(open(os.path.join(dd_dir(variant), "summary.json")))
    dd_id = json.load(open(os.path.join(dd_dir(variant), "identity.json")))
    problems = []
    if not di_sum["gate_pass"]:
        problems.append(f"DI {variant} fails the gate: {json.dumps({k: di_sum[k] for k in ('divergences', 'max_rhat', 'min_ess_bulk', 'min_ess_tail', 'min_ebfmi')})}")
    if dd_sum["outer_fits_failing_gate"] or dd_sum["divergences"]:
        problems.append(f"calibrated DD {variant}: {dd_sum['outer_fits_failing_gate']} fits fail the gate, {dd_sum['divergences']} divergences")
    if dd_id["parent_di_identity"] != di_id["identity_sha256"]:
        problems.append(f"calibrated DD {variant} was fitted on another DI fit")
    if problems:
        sys.exit("REFUSED:\n" + "\n".join(problems))
    di = pd.read_csv(os.path.join(SRC, f"di_{variant}", "draws_params.csv")).set_index("draw")
    picks = np.linspace(0, len(di) - 1, N_OUTER).round().astype(int)
    rows_di, rows_dd = [], []
    for idx in picks:
        f = pd.read_csv(os.path.join(dd_dir(variant), f"di_draw_{idx:05d}.csv"))
        inner = f.iloc[np.linspace(0, len(f) - 1, K).round().astype(int)]
        p = di.loc[idx]
        for _, q in inner.iterrows():
            assert np.isclose(q.alpha, p.alpha, rtol=0, atol=1e-12) and np.isclose(q.ab, p.ab, rtol=0, atol=1e-12)
            rows_di.append(dict(alpha=p.alpha, ab=p.ab, **{"lambda": p["lambda"]}, Tk=p.Tk, g=p.g, di_draw=idx))
            rows_dd.append(dict(alpha=p.alpha, ab=p.ab, **{"lambda": q.lambda_dd}, c_over_alpha=q.c_over_alpha, g=p.g,
                                di_draw=idx, dd_chain=int(q["chain"]), dd_iter=int(q["iter"])))
    out = {}
    for model, rows in (("di", rows_di), ("dd", rows_dd)):
        name = f"tumor_{model}_{variant}.csv"
        path = os.path.join(DRAWS, name)
        pd.DataFrame(rows).to_csv(path, index=False, float_format="%.10g", lineterminator="\n")
        out[name] = dict(sha256=sha(path), n_draws=len(rows), n_distinct_di_draws=int(len(picks)), variant=variant,
                         model=model, variant_label=VARIANTS[variant])
    di_diag = {k: di_sum[k] for k in ("divergences", "treedepth_hits", "max_rhat", "worst_rhat", "min_ess_bulk",
                                      "min_ess_tail", "min_ebfmi", "n_quantities", "gate_pass", "ppc_inside_90")}
    dd_diag = {k: dd_sum[k] for k in ("divergences", "treedepth_hits", "max_rhat", "min_ess_bulk_per_fit",
                                      "min_ess_tail_per_fit", "min_ebfmi", "outer_fits_failing_gate", "outer_fits_retried")}
    # Reached only when every conditional fit passed (refused above otherwise).
    dd_diag["gate_pass"] = bool(dd_sum["outer_fits_failing_gate"] == 0 and dd_sum["divergences"] == 0
                                and dd_sum["treedepth_hits"] == 0)
    rep = os.path.join(dd_dir(variant), "repair.json")
    dd_diag["outer_fits_repaired"] = len(json.load(open(rep))["repairs"]) if os.path.exists(rep) else 0
    out[f"tumor_di_{variant}.csv"].update(diagnostics=di_diag, fit_identity=di_id["identity_sha256"],
                                          min_ess_bulk=di_sum["min_ess_bulk"])
    out[f"tumor_dd_{variant}.csv"].update(diagnostics=dd_diag, fit_identity=dd_id["identity_sha256"],
                                          parent_di_identity=dd_id["parent_di_identity"],
                                          min_ess_bulk=dd_sum["min_ess_bulk_per_fit"],
                                          fitted_onset_below_7=dd_sum["fitted_onset_below_min"])
    return out, di_id


def point_values():
    """Deterministic "fitted values" preset: medians of the identified DI quantities (alpha,
    lambda/alpha, g/alpha, Tk) with lambda and g rebuilt from them; alpha/beta stays 10 Gy because
    the data do not identify it (posterior close to prior). DD: the pooled median of lambda_DD,
    calibrated at the anchor with these DI values, which is the calibrated-DD model at point values."""
    dr = pd.read_csv(os.path.join(DRAWS, "tumor_di_primary.csv")).drop_duplicates("di_draw")
    a = float(dr["alpha"].median())
    di = dict(alpha=a, ab=10.0, **{"lambda": float((dr["lambda"] / dr["alpha"]).median()) * a},
              Tk=float(dr["Tk"].median()), g=float((dr["g"] / dr["alpha"]).median()) * a,
              loss_per_day_after_onset=float((dr["lambda"] / dr["alpha"]).median()),
              growth_loss_per_day=float((dr["g"] / dr["alpha"]).median()))
    dd_draws = pd.read_csv(os.path.join(DRAWS, "tumor_dd_primary.csv"))
    lam_dd = float(dd_draws["lambda"].median())
    onset, co = C.calibrate(di["alpha"], di["ab"], di["lambda"], di["Tk"], lam_dd)
    dd = dict(alpha=a, ab=10.0, **{"lambda": lam_dd}, c_over_alpha=float(co), onset_at_anchor_days=float(onset),
              g=di["g"], loss_per_day_after_onset=lam_dd / a, growth_loss_per_day=di["growth_loss_per_day"])
    return dict(di=di, dd=dd)


def main():
    os.makedirs(DRAWS, exist_ok=True)
    files, ident = {}, None
    for v in VARIANTS:
        f, i = export_variant(v)
        files.update(f)
        ident = ident or i
    data = os.path.join(ROOT, "data", ident["data_file"])
    man = dict(created_by="fit/export_draws.py", draws_version="3.1",
               analysis="verified-trial refit 2026-10-25 (corrected main set): DI posterior; DD calibrated to DI "
                        "at 35 x 2 Gy in 47 days, fitted conditional on each DI draw",
               dataset=ident["data_file"], data_sha256=sha(data),
               pairing="row i of tumor_dd_<variant>.csv was fitted conditional on row i of tumor_di_<variant>.csv",
               toxicity_pairing="toxicity prior draws are independent of the tumor draws (assumption)",
               likelihood="binomial on rounded N x reported rate (approximation: most rates are KM or CIF estimates)",
               onset_penalty="soft preference for DD onset at or after day 7, not a bound",
               mcse_note="draws are not independent: each DI draw is repeated in K rows; Monte Carlo error "
                         "is approximate (IID benchmark on distinct DI draws)",
               source_fit=os.path.relpath(SRC, ROOT), files=files)
    json.dump(man, open(os.path.join(DRAWS, "draws_manifest.json"), "w"), indent=1)
    pv = point_values()
    json.dump(dict(source="draws/tumor_*_primary.csv", dataset=ident["data_file"], rule="medians of identified quantities; "
                   "DD lambda median calibrated at the anchor with the DI point values", values=pv),
              open(os.path.join(DRAWS, "refit_point_values.json"), "w"), indent=1)
    for k, v in files.items():
        print(k, v["n_draws"], v["sha256"][:12])
    print(json.dumps(pv, indent=1))


if __name__ == "__main__":
    main()
