# Custom BED calculator, version 4 (branch `v4`, not deployed)

Version 4 keeps the version 3.2.1 engine and changes the layout for clinicians. The calculator opens on one page: a list of shorter, hypofractionated schedules (5 fractions a week) that, compared with the standard 35 × 2 Gy in 47 days, keep the tumor BED (DI or DD model) at least the standard's minus 5 Gy and the early and late toxicity BEDs at most the standard's plus 5 Gy. Each row gives a number of fractions with the range of doses per fraction that qualifies, and says whether both tumor models agree. Every version 3 tab is still available as a page of a side menu. Version 4 is meant for its own site, `custom_bed_calculator_v4`. The version 3 site (https://ishuryak.shinyapps.io/custom_bed_calculator_v3/) runs 3.2.0 and is not touched. Reasoning, decisions and checks: `docs/development/2026_10_28_v4_simplified_layout.md`.

R Shiny app for head and neck squamous cell carcinoma (HNSCC), for designing hypofractionated trial regimens. It computes:

- tumor BED under the dose-independent (DI) and dose-dependent (DD) repopulation models
- early (acute mucosal) and late toxicity BED
- organ-at-risk doses against QUANTEC or HYHOPE constraints: maximum doses converted to any number of fractions, mean doses compared on physical dose (new in version 3)
- the uncertainty of the tumor BED differences, from a Bayesian refit to 45 arms of 19 randomized trials re-extracted from their publications (version 3.1).

Every tab compares the evaluated regimen with one comparison regimen and one tumor threshold, shown in a banner. A search lists shorter regimens that meet the threshold within the toxicity and organ limits; a both-model screen requires the threshold under DI and DD. It is a research tool for hypothesis generation, not for treatment decisions.

**Live:** https://ishuryak.shinyapps.io/custom_bed_calculator_v3/ (version 3.2.0, revision `4876ae4`, deployed 2026-10-26 with `Rscript tools/deploy.R` as its own app; previous release 3.1.0 at `e4ff22f`, the rollback point). Live browser walkthrough `docs/development/2026_10_21_v3/walkthrough_v3.R --live`: 35 of 35 checks on 2026-10-26. `tools/deploy.R` refuses the names of the live version 1 and version 2 apps, an unconfigured account, uncommitted changes, a revision on which the full suite has not passed, a bundle missing any file the app reads, and a bundle holding any file not tracked by git. Version 2.7 (https://ishuryak.shinyapps.io/custom_bed_calculator_v2/) and version 1 (https://ishuryak.shinyapps.io/custom_bed_calculator/) are unchanged.

**The refit is for review, not a default.** Version 3.2 uses the refit to 45 trial arms. The DI model is informed by the data (repopulation loss λ/α 0.63 Gy per day at each draw's α/β, 95% interval 0.21–1.19, after about day 37; background growth is essentially its prior). A DD model fitted on its own does not sample reliably, so DD is calibrated to DI at 35 × 2 Gy in 47 days and its rate fitted given each DI draw (0.41 per day, 0.12–1.04). Every DI and conditional DD fit meets the full diagnostic gate, effective sample size included. For HYPNO's 55 Gy in 20 fractions against 35 × 2 Gy the median tumor BED change is −2.0 Gy under DI and −3.3 Gy under DD. The defaults remain Scenario D18 with the DD model calibrated to DI. The records are `docs/development/2026_10_25_stages_BCDE_refit.md` and `docs/development/2026_10_26_review_v31_fixes.md`.

**Version 1 warning.** Version 1 (branch `main`) still carries the DD calibration error described below. Do not use version 1 DD results computed at a tumor α/β other than 10 Gy.

## Files

- `R/bed_models.R`: all formulas, as plain functions; `dd_setup()` chooses the calibrated or the fitted DD model.
- `R/schedule.R`: treatment calendars (5 or 6 fractions per week, first day, days without treatment).
- `R/uncertainty.R`: tumor and toxicity BEDs over regimens × parameter draws, shares and quantiles, the search under uncertainty.
- `R/oar.R`: organ constraint library, LQ conversion, DVH reading (refuses identifier fields) and the feasibility table.
- `R/mod_uncertainty.R`, `R/mod_oar.R`: the Uncertainty and Organs at risk tabs (Shiny modules).
- `R/export.R`: the scenario file (schema 3) and ZIP exports, `rerun_scenario()` and `rerun_uncertainty()`.
- `data/`: the 2018 Table 1 trial arms (two independent transcriptions, no disagreement) and the organ constraints, each with source, page and quotation.
- `draws/`: 4000 joint posterior draws per model and fit variant, with SHA-256 manifest. Produced offline by `fit/` (CmdStan); Stan never runs in the app.
- `app.R`: the user interface.
- `about.md`: the "Model and sources" page.
- `tests/testthat/`: tests of the formulas against an independent implementation and published values, of the server logic, and of the export round trip.
- `run_app.bat`: double-click launcher for Windows. It uses `RSCRIPT` if set, else Rscript on the PATH, else the newest R under Program Files, and starts `tools/launch.R`, which uses the release library when it has been restored.
- `tools/`: `restore_library.R` (installs the versions in `renv.lock`), `run_tests.R` (full suite in that library only), `launch.R`, `env.R`, `deploy.R` (publishes to shinyapps.io as `custom_bed_calculator_v3` only; refuses the version 1 and 2 names; pinned account; clean tree; suite passed on the same revision).
- `.rscignore`: keeps `References/`, `tests/`, `docs/`, `fit/`, `.Rhistory` and `.Renviron` out of a shinyapps.io deployment.
- `docs/FIXES.md`: one row per landed change, with its evidence.
- `docs/user_guide/`: user guide to version 3.2.1 (every tab and menu item, the model fitting and its literature, and the uncertainty estimation), as Word and Markdown.
- `docs/development/`: review records and plans. For version 3 start with `2026_10_21_v3_build.md` and `2026_10_21_v3/fit_report.md`. For version 3.1, `PLAN_codex_review_2026_10_25.md` and `2026_10_25_stages_BCDE_refit.md`.
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

## Changes in version 4.0.0

- New main page (`R/v4_home.R`, `app.R`): the list described above, built by the unchanged `suggest_regimens()` with comparison 35 × 2 Gy, criterion "DI or DD", tumor margin −5 Gy, early and late margins +5 Gy, 5–34 fractions, 0.05 Gy steps up to 8 Gy. A row click shows each dose with its BED differences.
- Layout: the version 3 tabs are menu pages. The evaluated regimen and the comparison are set on the Report card page, the model parameters on a Model settings page. Visible text names pages instead of the sidebar and tabs.
- Deploy target `custom_bed_calculator_v4`. The preflight now also refuses the version 3 site.
- Stale-constant registry: three new constants (margins, standard overall time, search ceiling) and three `app.R` site pointers repaired (they had moved with PR #1).
- No equation, parameter, draw or selection rule changes. Every version 3 page gives the same numbers for the same inputs.

## Changes in version 3.2.1

Response to the independent review of version 3.2 (2026-10-26). No BED equation, selection rule,
draw file or displayed tumor, acute or late number changes. Record: `docs/FIXES.md` and
`docs/development/2026_10_26_v321_review_fixes.md`.

- **DVH import (F01).** A DVH row with a missing or extra field, or two structures with the same
  name, is now refused with its line number. Before, a three-column row missing one field was
  recycled to fill the row and a repeated name replaced the first block, and both could turn an
  exceeded organ limit into "within".
- **Diagnostic gate (F02).** Every quantity's R-hat and ESS must be finite before minima are taken,
  a sampled quantity may not be constant, and counts must be exactly 0, in Python and in R (shared
  fixture table). Fits from before 3.2.1 kept only per-fit aggregates; only those were re-checked.
- **Anchor residual (F03)** is now in Gy of BED (shipped maxima 9.8e-7 and 7.8e-7 Gy, tolerance
  1e-4 Gy). **Repair (F04)** refuses a changed prior, onset penalty, anchor or variant.
  **Draw publication (F05)** restores the previous release if it fails, and uses a lock.
- **Text (F06, F07).** The Both-model screen lists regimens at the current settings and does not
  require 15 of 15 (17 x 3.05 Gy is listed at 12 of 15). Repair count 752. Monte Carlo errors
  re-measured for named scenarios. Selection hash portable across NumPy integer widths.
- **Maintenance (M01-M04).** Deploy identity from checked git calls; unused DVH fraction input
  removed; toxicity shares near the required share shown to one decimal with "meets" or "below";
  the version 2.6 parts of `References/README.md` marked historical.

## Changes in version 3.2

From the review of 2026-10-26 (plan `docs/development/PLAN_codex_review_2026_10_26.md`, record `docs/development/2026_10_26_review_v31_fixes.md`).

- **What moves.**
  - The default evaluated regimen is 18 × 2.9 Gy (52.2 Gy in 24 days) instead of HYPNO's 55 Gy in 20 fractions, so every tab opens on that regimen. The comparison (35 × 2 Gy) and DD anchor are unchanged.
  - The DD draws are regenerated (the DI draws are unchanged). HYPNO against 35 × 2 Gy: DI −2.0 Gy (80% range −4.3 to +1.0, full DI posterior, unchanged), DD −3.27 Gy (−5.92 to −0.53, conditional fits, was −3.3). The change is within Monte Carlo error.
  - On the Results tab, the HYPNO control and HYHOPE benchmark rows now show their tumor BED differences (+1.11/−2.87 and −10.60/−15.14 Gy₁₀) instead of "not available".
  - An organ maximum that is indeterminate at the DVH resolution is no longer reported as "within" over the α/β range.
  - Uncertainty scenario files of version 3.1 are refused with a reason, because their draws were replaced. Deterministic scenario files replay as before.
- **One diagnostic gate.** `fit/gate.py` is used by every fit, the repair, the export and (for draws version 3.2) the app. All 3600 conditional DD fits now pass it, 752 after refitting with longer chains (`fit/repair_dd_ess.py`): 748 had been below the bulk-ESS criterion and 4 below the tail criterion only (the 748 first reported counted bulk ESS only; corrected 2026-10-26 from `repair_ess.json`).
- **Identity chain.** The export checks every input byte and every pair's anchor equality, and replaces `draws/` only when the whole set passes.
- **Inputs.** One organ-mapping validator for both replay paths. Malformed Eclipse DVH rows are refused with their line number. Replayed searches are bounded at the largest grid the app can request.
- **Provenance.** `data/*_provenance.csv` keep source intervals, estimators, population exceptions and prescription categories beside the unchanged fitting files.
- **Wording.** Cohort, chemotherapy sensitivity, α, tumor control, QUANTEC and mean-dose scope, acute onset prior and predictive-check statements corrected.

## Changes in version 3.1

- **What moves.** The Uncertainty tab, the report card's refit column and the refit preset now come from the verified-trial refit (45 arms, 19 trials) instead of the 35 arms of the 2018 Table 1, and every DD draw is paired with the DI draw it was fitted on. Example: HYPNO 55 Gy/20 against 35 × 2 Gy, first 50 shipped DI draws, mean ΔDI changes from +1.75 to −2.28 Gy. Deterministic results with the D18 defaults are unchanged. Mean-dose organ limits are compared on physical dose, which changes those organ statuses (stage A).
- **Fits gated (corrected in 3.2).** The DI fits passed 0 divergences, R-hat ≤ 1.01, bulk and tail ESS ≥ 400 and E-BFMI ≥ 0.3. The conditional DD fits were not checked for ESS: 177 of 1000 primary and 272 of 1000 restricted fits had bulk ESS below 400 (minimum 151). Version 3.2 refits them.
- **Units.** BED differences over draws are labeled Gy₁₀ only when α/β is 10 Gy in every draw.
- **Records.** Plan `docs/development/PLAN_codex_review_2026_10_25.md`, stage A `docs/development/2026_10_25_stage_A_app_fixes.md`, stages B–E `docs/development/2026_10_25_stages_BCDE_refit.md`, `docs/FIXES.md`.

## Changes in version 3.0

- **What moves.** For the same inputs, every BED, search, sensitivity and both-model result is identical to version 2.7 (the full version 2 test suite passes unchanged except where a default or label changed on purpose, and version 2.6/2.7 scenario files replay exactly). Two defaults changed: the evaluated regimen is now HYPNO's 55 Gy in 20 fractions (was 20 × 2.65 Gy, not a published schedule), and the both-model screen now uses the global comparison and threshold (35 × 2 Gy, +1 Gy), so it opens empty instead of listing 17, 18 and 19 fractions against HYPNO at −3 Gy. Choosing HYPNO and −3 Gy (globally or as a highlighted override) restores that list.
- **One comparison and one threshold** for every tab, in a banner (plan decision 7a).
- **Report card** and **Compare regimens** tabs; each toxicity limit's value beside its control; differences labeled with their comparison and unit; "about the same" below 0.5 Gy; days instead of weeks; colorblind-safe colors.
- **Uncertainty** tab from a Bayesian refit of the 2018 trial data (`fit/`, `draws/`), and an optional **fitted DD** model (decisions 1a and 3a; not the default).
- **Organs at risk** tab with QUANTEC and HYHOPE constraint libraries (decision 4a/5a), and organ limits optionally applied in the searches.
- **Calendars**: 6 fractions per week, a chosen first day and days without treatment (decision 8a).
- **Records**: `docs/development/2026_10_21_v3_build.md`, `docs/development/2026_10_21_v3/`.

## Changes in version 2.7

Interface changes from a colleague's comments of 2026-10-20 and an independent review of the plan (`docs/development/2026_10_20_v27_colleague_review.md`). **A fresh page now shows a different candidate table:** it opens with both models required and the smaller change used for ranking, and at the defaults no regimen qualifies (before: three regimens that met +1 Gy under DI only, with DD losses of 2.05 to 2.86 Gy₁₀). No equation, parameter default, toxicity limit or selection rule changed, and the search function's own defaults are unchanged. Compatible 2.6 scenario files replay with the same numerical results and selection rules.

- **Two-model default.** "Tumor BED criterion" defaults to both models and "Rank by" to the smaller change. The empty result explains itself and how the Model agreement tab differs. "DI or DD model" and the larger-change ranking remain available.
- **Parameter settings.** Agreement and Sensitivity counts are labeled as selected parameter settings ("15 of 15 settings") and say they are not independent confirmations, confidence intervals or probabilities. Export field names are unchanged.
- **Model limits beside the results.** A short line on the Results, Candidate and Model agreement tabs lists what is not modeled.
- **Source link.** The Model and sources page links this branch.
- **Corrected text.** Chemotherapy and treatment-gap limitations, and the history of the DD calibration anchor (below), are stated more precisely.

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

- The legacy checkbox "calibrate at the comparison regimen" is removed. To reproduce version 2.0, enter the comparison regimen as the DD anchor. Version 1 calibrated at a fixed 35 × 2 Gy, 47-day course (the default anchor), with an α/β error described below. The anchor note now explains in plain words what the anchor does.
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
- **Fixed DD calibration anchor.** The DD model is calibrated at an anchor regimen (default 35 × 2 Gy, 47 days) that no longer moves with the comparison regimen. To reproduce version 2.0, which calibrated at the comparison regimen, enter that regimen as the anchor. Version 1 calibrated at a fixed 35 × 2 Gy, 47-day course.
- **Presets** set every value their source gives: version 1 defaults, Withers (follows the tumor α/β), Hartley 2011, the HYHOPE protocol, and scenario D18. Defaults are unchanged.
- **Search:** blocks a model whose comparison BED is negative, shows which model credits each row, flags fraction sizes above 3.1 Gy, allows a negative minimum change, and offers HYPNO 55 Gy/20 fractions as an acute and as a late limit. Default maximum dose per fraction 5 Gy.
- Optional background growth g, late BED at α/β 2 and 4 Gy, HYHOPE benchmark row, ZIP export with a rerunnable scenario, Sensitivity tab.

Details and evidence: `docs/FIXES.md`.

## Change from version 1

Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days, but hard-codes that regimen's value at a tumor α/β of 10 Gy. At any other α/β the calibration is wrong. For example, at α/β = 5 Gy the DD tumor BED is 73.8 Gy instead of 79.0 Gy. Version 2 recalculates the calibration for any α/β. At α/β = 10 Gy the two versions agree to a relative difference below 1e-8 for tumor BED. Version 2 also marks the DD model as unavailable when no calibration with a repopulation onset after the start of treatment exists, instead of showing a meaningless number.

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
