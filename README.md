# Custom BED calculator, version 2.6

R Shiny app for head and neck squamous cell carcinoma (HNSCC). It computes three biologically effective doses (BEDs):

- tumor BED under the dose-independent (DI) and dose-dependent (DD) repopulation models
- early (acute mucosal) toxicity BED
- late toxicity BED.

It also screens shorter, hypofractionated regimens against a comparison regimen (default 35 × 2 Gy in 47 days, or either HYPNO arm): the tumor BED change must meet a chosen minimum, which may be negative (an allowed loss of modeled tumor BED, not a clinical non-inferiority margin), under the DI or DD model, while early and late toxicity BED stay within chosen limits. A sensitivity tab shows the tumor BED difference across tumor parameter sets. Results download as a ZIP with a scenario file that reruns the search.

It is a research tool for hypothesis generation, not for treatment decisions.

**Not deployed.** Version 2.6 is on branch `v2-toxicity-hypofractionation` and has not been published to shinyapps.io.

**Version 1 warning.** Version 1 (https://ishuryak.shinyapps.io/custom_bed_calculator/, branch `main`) is unchanged and still carries the DD calibration error described below. Do not use version 1 DD results computed at a tumor α/β other than 10 Gy.

## Files

- `R/bed_models.R`: all formulas, as plain functions.
- `R/export.R`: the scenario file and ZIP export, and `rerun_scenario()`.
- `app.R`: the user interface.
- `about.md`: the "Model and sources" page.
- `tests/testthat/`: tests of the formulas against an independent implementation and published values, of the server logic, and of the export round trip.
- `run_app.bat`: double-click launcher for Windows. It uses `RSCRIPT` if set, else Rscript on the PATH, else the newest R under Program Files, and starts `tools/launch.R`, which uses the release library when it has been restored.
- `tools/`: `restore_library.R` (installs the versions in `renv.lock`), `run_tests.R` (full suite in that library only), `launch.R`, `env.R`.
- `.rscignore`: keeps `References/`, `tests/` and `docs/` out of a shinyapps.io deployment.
- `docs/FIXES.md`: one row per landed change, with its evidence.
- `docs/development/`: review records and plans. Start with `2026_10_19_v26_review_response.md`.
- `renv.lock`: the release package versions (R 4.5.1), all from the Posit Package Manager snapshot of 2026-10-15, which serves binaries of exactly these versions for Windows, macOS and Linux.
- `References/`: sources, with `References/README.md` giving each one's license and the quoted parameter values. Paywalled or non-redistributable PDFs are git-ignored.

## Run locally

From the app folder, once per computer (about one minute, binaries only, no compiler needed):

```
Rscript tools/restore_library.R
```

This installs the locked versions into a per-user folder outside the project (`tools::R_user_dir("HNSCC_BED_calculator", "cache")`; set `BED_LIBRARY` to choose another). Then:

```
Rscript tools/run_tests.R      # full suite in that library only; exit status 1 on any failure
run_app.bat                     # or: Rscript tools/launch.R
```

Each prints the R and package versions it used. Without the restore step the app still starts with whatever packages are installed in R, and says so; results from such a session are not from the release environment. On 2026-10-19 a restore into an empty library installed 68 packages as binaries, and `tools/run_tests.R` then gave 2541 expectations, 0 failed, 0 errors, 0 warnings.

## Changes in version 2.6

Fixes from the external review of 2.5.1 (`docs/development/2026_10_19_v26_review_response.md`). No equation and no default changed; the default candidate and agreement rows are unchanged.

- **Ranking (result change in edge cases).** A regimen is ranked only when its ranking value can be evaluated. Before, the "smaller change (conservative)" ranking could list a regimen whose DD BED was not positive, scored by DI alone: in the review's example it listed 32 × 2.50 Gy (DD BED −0.60, score +22) instead of 32 × 2.10 Gy (score +3.31). A one-model ranking could list a row with no score. Selection version 2.6.
- **Agreement ranges** are now over the same evaluable scenarios as the counts (before, each model's range used its own evaluable scenarios). Each model's own count is in the download.
- **Replay** is strict: only files of the current versions are rerun; malformed fields (for example a Fenwick flag "true" instead of true), contradictory agreement settings, a wrong saved onset and oversized grids are refused. Older files are reanalyzed only on request with `migrate_scenario()`, which records their original versions and ID.
- **Wording.** Each analysis tab states its comparison, anchor and numerical limits; the sidebar says which tabs use its comparison; the version 1 and 2.0 anchor history is corrected; the HYPNO dose-response illustration states its extra assumptions and uses the exact logistic curve; the default narrative and preset notes match the outputs; source records are reconciled.
- **Candidate table:** fractions saved, and early and late limit minus the regimen's BED.
- **Environment:** a lockfile that restores without a compiler, used by the launcher, the tests and CI.

## Changes in version 2.5.1

Interface only; no result changes (`docs/development/2026_10_19_interface_clarity_and_input_audit.md`).

- The legacy checkbox "calibrate at the comparison regimen" is removed. To reproduce versions 1 and 2.0, enter the comparison regimen as the DD anchor. The anchor note now explains in plain words what the anchor does.
- The preset menu labels the default (scenario D18) and explains the selected preset: its source, its values, and that D18 is not fitted.
- A browser audit changed every menu box one at a time and checked which outputs change and which must not.

## Changes in version 2.5

Fixes from the 2026-10-19 review (`docs/development/2026_10_19_review_and_plan.md`, which verifies the external review in `docs/development/2026_10_19_external_review_codex.md`). No equation and no default changed; default results are unchanged.

- **One validity rule** for the search, the Sensitivity tab and the Model agreement tab: a tumor model's difference is used only when the model is calibrated and both BEDs are positive. The agreement tab counts evaluable scenarios separately and never counts a non-evaluable one as passing.
- **Wording:** no claims of parameter uncertainty, clinical equivalence or bounds that the calculator does not compute. The numerical control is an allowed BED change, not a non-inferiority margin. The D18/HYPNO statement is conditional on a cited dose-response slope.
- **Independent tabs:** the agreement tab no longer depends on candidate-only inputs, the candidate export no longer depends on the Results regimen. The legacy coupling switch was removed on 2026-10-19 (see `docs/FIXES.md`).
- **Export and replay:** the agreement analysis has its own ZIP with every regimen under every scenario. Scenario files (schema 2) record the analysis, selection and scenario-set versions and the code identity. Replay refuses unsupported versions, hinge or clocks, and checks the saved calibration. The scenario ID no longer depends on key order.
- **Labels and records:** BED units follow the α/β in use, the fraction-size flag is named for the benchmark set, and the Model page has a parameter-origin table. Source records are corrected. New `renv.lock` and a launcher that finds R.

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
- **Fixed DD calibration anchor.** The DD model is calibrated at an anchor regimen (default 35 × 2 Gy, 47 days) that no longer moves with the comparison regimen. To reproduce versions 1 and 2.0, which calibrated at the comparison regimen, enter that regimen as the anchor.
- **Presets** set every value their source gives: version 1 defaults, Withers (follows the tumor α/β), Hartley 2011, the HYHOPE protocol, and scenario D18. Defaults are unchanged.
- **Search:** blocks a model whose comparison BED is negative, shows which model credits each row, flags fraction sizes above 3.1 Gy, allows a negative minimum change, and offers HYPNO 55 Gy/20 fractions as an acute and as a late limit. Default maximum dose per fraction 5 Gy.
- Optional background growth g, late BED at α/β 2 and 4 Gy, HYHOPE benchmark row, ZIP export with a rerunnable scenario, Sensitivity tab.

Details and evidence: `docs/FIXES.md`.

## Change from version 1

Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days, but hard-codes that regimen's value at a tumor α/β of 10 Gy. At any other α/β the calibration is wrong. For example, at α/β = 5 Gy the DD tumor BED is 73.8 Gy instead of 79.0 Gy. Version 2 recalculates the calibration for any α/β. At α/β = 10 Gy the two versions agree to a relative difference below 1e-8 for tumor BED. Version 2 also marks the DD model as unavailable when no calibration with a repopulation onset after the start of treatment exists, instead of showing a meaningless number.

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
