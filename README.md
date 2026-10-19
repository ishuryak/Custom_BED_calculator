# Custom BED calculator, version 2.4

R Shiny app for head and neck squamous cell carcinoma (HNSCC). It computes three biologically effective doses (BEDs):

- tumor BED under the dose-independent (DI) and dose-dependent (DD) repopulation models
- early (acute mucosal) toxicity BED
- late toxicity BED.

It also screens shorter, hypofractionated regimens against a comparison regimen (default 35 × 2 Gy in 47 days, or either HYPNO arm): the tumor BED change must meet a chosen minimum, which may be negative as a non-inferiority margin, under the DI or DD model, while early and late toxicity BED stay within chosen limits. A sensitivity tab shows the tumor BED difference across tumor parameter sets. Results download as a ZIP with a scenario file that reruns the search.

It is a research tool for hypothesis generation, not for treatment decisions.

**Not deployed.** Version 2.4 is on branch `v2-toxicity-hypofractionation` and has not been published to shinyapps.io.

**Version 1 warning.** Version 1 (https://ishuryak.shinyapps.io/custom_bed_calculator/, branch `main`) is unchanged and still carries the DD calibration error described below. Do not use version 1 DD results computed at a tumor α/β other than 10 Gy.

## Files

- `R/bed_models.R`: all formulas, as plain functions.
- `R/export.R`: the scenario file and ZIP export, and `rerun_scenario()`.
- `app.R`: the user interface.
- `about.md`: the "Model and sources" page.
- `tests/testthat/`: tests of the formulas against an independent implementation and published values, of the server logic, and of the export round trip.
- `run_app.bat`: double-click launcher for Windows.
- `.rscignore`: keeps `References/`, `tests/` and `docs/` out of a shinyapps.io deployment.
- `docs/FIXES.md`: one row per landed change, with its evidence.
- `docs/development/`: review records and plans. Start with `2026_10_18_synthesis_plan.md`.
- `References/`: sources, with `References/README.md` giving each one's license and the quoted parameter values. Paywalled or non-redistributable PDFs are git-ignored.

## Run locally

```r
shiny::runApp("path/to/Custom_BED_calculator_v2")
testthat::test_dir("tests/testthat")
```

## Changes in version 2.4

- **Model agreement tab:** regimens whose tumor BED change meets a margin under both the DI and DD models, by default against HYPNO 55 Gy/20 fractions with a −3 Gy₁₀ margin, with each model's range over the 15 sensitivity scenarios, the number of scenarios in which both meet the margin, and the DD and DI onset days that explain the result. Defaults give 17 × 3.05, 18 × 2.90 and 19 × 2.80 Gy (23–25 days), within −3 Gy₁₀ in all 15 scenarios.

## Changes in version 2.3

- **New defaults (Igor Shuryak, 2026-10-19):** tumor scenario D18 (α 0.35 Gy⁻¹, α/β 10 Gy, DI rate 0.23/day, DD maximum rate 0.40/day, onset day 28) and the acute BED of HYPNO 55 Gy/20 fractions as the early-toxicity limit. The version 1 tumor values stay available as a preset. Default results change: see `docs/FIXES.md`.

## Changes in version 2.2

- Optional Fenwick et al. 2008 acute boundary (cell-kill BED, up to 50 days), as a fourth early-toxicity limit and as a line on the Results tab.
- Fixes from an independent review of 2.1: a model whose comparison BED is negative now shows no differences (it showed large spurious gains, although it no longer ranked or credited rows); a separate message when both models are negative; the scenario ID depends only on settings that change the search; an empty export records why it is empty; the Withers rule survives a briefly empty α/β field.
- A test-only GitHub Actions workflow (`.github/workflows/tests.yml`). It never deploys.

## Changes in version 2.1

- **Acute BED on Fowler's clock.** Fowler's acute parameters and tolerance zone count elapsed days. Acute BEDs are now 0.79 Gy₁₀ higher than in version 2.0 (35 × 2 Gy: 53.1 instead of 52.3), which makes the fixed Fowler limit stricter. Comparisons with a reference beyond the onset day do not change.
- **Fixed DD calibration anchor.** The DD model is calibrated at an anchor regimen (default 35 × 2 Gy, 47 days) that no longer moves with the comparison regimen. The old coupling is available as a labeled legacy option.
- **Presets** set every value their source gives: version 1 defaults, Withers (follows the tumor α/β), Hartley 2011, the HYHOPE protocol, and scenario D18. Defaults are unchanged.
- **Search:** blocks a model whose comparison BED is negative, shows which model credits each row, flags fraction sizes above 3.1 Gy, allows a negative minimum change, and offers HYPNO 55 Gy/20 fractions as an acute and as a late limit. Default maximum dose per fraction 5 Gy.
- Optional background growth g, late BED at α/β 2 and 4 Gy, HYHOPE benchmark row, ZIP export with a rerunnable scenario, Sensitivity tab.

Details and evidence: `docs/FIXES.md`.

## Change from version 1

Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days, but hard-codes that regimen's value at a tumor α/β of 10 Gy. At any other α/β the calibration is wrong. For example, at α/β = 5 Gy the DD tumor BED is 73.8 Gy instead of 79.0 Gy. Version 2 recalculates the calibration for any α/β. At α/β = 10 Gy the two versions agree to a relative difference below 1e-8 for tumor BED. Version 2 also marks the DD model as unavailable when no calibration with a repopulation onset after the start of treatment exists, instead of showing a meaningless number.

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
