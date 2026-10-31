# Fitting pipeline: where the retained fit artifacts live

The app never runs Stan. It reads the shipped draws in `draws/`, so nothing on this page affects
how the running app works. This page is for anyone who wants to validate or re-export those draws.

## This checkout is not a complete fit archive

The v4 checkout has the fitting code, the fit identities (`identity.json`), the summaries
(`summary.json`, `repair_ess.json`) and the shipped draws. It does not have the large fit CSVs,
which `.gitignore` excludes: `draws_params.csv` of each DI and DD fit and the per-fit
`di_draw_*.csv` files of the calibrated DD fits. For example,
`fit/output_refit_2026_10_25c_main/di_primary/draws_params.csv` is absent here. A relative
`source_fit` label in `draws/draws_manifest.json` names a folder; it does not show that the
folder's files are present.

The authoritative retained generation is the 2026-10-25c refit, kept in the sibling v3 checkout:

    ../Custom_BED_calculator_v3/fit/output_refit_2026_10_25c_main      (main fit, exported to draws/)
    ../Custom_BED_calculator_v3/fit/output_refit_2026_10_25c_<name>    (the sensitivity fits)

The ten calibrated DD identities there are byte-for-byte the same as the copies in this checkout.
Keep those folders when archived app versions are cleaned up.

## Pointing the export at the retained generation

`fit/export_draws.py` reads its source from `BED_EXPORT_FROM` and writes to `BED_EXPORT_TO`
(default `draws/`, the app's active draws). For validation, always write to a scratch folder:

    BED_EXPORT_FROM=../Custom_BED_calculator_v3/fit/output_refit_2026_10_25c_main \
    BED_EXPORT_TO=/tmp/bed_export_check \
      python3 -I fit/export_draws.py

Then compare the hashes of the four `tumor_*.csv` files in the scratch folder with those in
`draws/draws_manifest.json`. Replace `draws/` only as a deliberate release step.

## Missing or mismatched artifacts must fail, never trigger a refit

The export checks the data file, the parent DI draws, the identities and the Stan sources before
writing anything. A missing or mismatched file stops it with an error, and it never starts a fit.
The fit and repair scripts (`fit_refit.py`, `fit_dd_calibrated.py`, `repair_dd_ess.py`) write into
`BED_FIT_OUT` and do sample Stan. Do not run them to repopulate a checkout, and do not point them
at the retained generation except for a deliberate, planned repair. A changed anchor, likelihood,
prior or evaluation regimen is a new scientific generation with its own output folder, and
`repair_dd_ess.py` refuses it (see the target check in that file).
