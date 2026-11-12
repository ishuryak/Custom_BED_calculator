# Version 4.6.1: the page "How the models were fitted (Bayesian)"

2026-11-12. Base revision `52f0d3b` (4.6.0, live at `c80f63b`). Not deployed. Instruction (Igor
Shuryak, 2026-11-12): "check, analyze and implement this review", for
`2026_11_12_bayesian_methods_and_late_toxicity_review.md`.

## 1. What led here

1. Igor Shuryak reported two issues with the live 4.6.0: the Bayesian fitting method is not clearly
   explained anywhere a user can find, and asked whether the same data can be used to refit the
   late-toxicity model as was done for early toxicity.
2. Claude proposed a plan and, on request, audited it
   (`2026_11_12_plan_audit_fitting_page_and_late_refit.md`).
3. An independent review checked that audit against the code, the fits, the live site and the
   primary reports. The review is the governing plan for this change.

## 2. Check of the review

Each claim the implementation depends on was checked against the code or the retained fits before
it was used:

| Review claim | Check | Result |
|---|---|---|
| The tumor brief already says α/β is mostly set by its prior and gives trial baselines | `refit_in_brief.md` | Confirmed, so the audit's "mentions no prior" was too strong |
| The deployed DD model is the calibrated one, conditional on each DI draw, with a soft onset penalty | `fit/stan/dd_cal.stan`, its identity (`onset_min` 7, `onset_scale` 1, λ prior lognormal(log 0.5, 1.5)) | Confirmed |
| 4,000 tumor rows hold 1,000 distinct DI draws, each repeated 4 times | `draws/tumor_di_primary.csv`: 1,000 distinct (α, α/β, λ, Tk, g) tuples | Confirmed |
| The tumor model has arm effects with an estimated common spread, so "hierarchical" cannot be banned by a regex | `di.stan`: `sigma_raw ~ normal(0, 0.5)`, `sigma * z[i]` | Confirmed. The audit's F2 test was dropped |
| Early priors and fixed settings (table in review 3.4) | `step5_model.stan`, `step5_fit.py` PRIOR, U02 diagnostics | Confirmed, and proven by reproducing the U02 run identity (below) |
| The early gate records tree-depth hits but does not use them | `step5_fit.py::gate_from_summary` | Confirmed and stated on the page |
| A lognormal with median 3 Gy and 95% limits 1–10 Gy is impossible | geometric mean of 1 and 10 is 3.162 | Confirmed. The audit's suggestion is withdrawn |
| Late α/β 1.62 → 1.61 Gy is where the default list first loses a dose (29 × 2.30 Gy) | (66.7 × 2.30 − 140) / (75 − 66.7) = 1.6157, and `v4_regimens()` at late α/β 1.5, 1.61, 1.62, 2 and 3 Gy | Confirmed: 65, 68, 69, 69, 69 doses with Fowler, 9 at every value in research mode |
| `test_v45_review.R` asserts six destinations | file, lines 179–186 | Confirmed and updated to seven |
| `docs/`, `fit/`, `tools/` are not deployed | `.rscignore` | Confirmed. The page text is at the app root and the record in `draws/` |

Nothing in the review was found to be wrong. One local discrepancy it raised was explained:
`deployed_revision.txt` on this computer names `faf6817` because the 4.6.0 deployment ran on
another computer. The file is gitignored and written by `tools/deploy.R` at upload, so the README
deployment history (4.6.0 at `c80f63b`) is the record.

## 3. Decisions taken in implementing it

- **One page, contextual entry points kept.** The menu item and the Details button open the page at
  its start. The three "How was the tumor model fitted?" links open it at the tumor section, and
  "How was the early-toxicity model fitted?" at the early section (review 3.1). The two dialogs
  that used to answer those links are removed. The trial-list dialogs stay and are linked from the
  page.
- **No hand-typed fitted numbers on the page.** The narrative (`how_fitted.md`) holds the method.
  The two existing summaries (`refit_in_brief.md`, `early_fit_in_brief.md`) are embedded unchanged as
  the result summaries. Every prior, fixed setting, sampler setting and diagnostic comes from
  `draws/model_methods.json`. The one numerical illustration in the narrative (EQD2 = BED / 1.2 at
  α/β 10 Gy) is a conversion, not a result.
- **The record is proven, not transcribed.** `tools/build_model_methods.py` refuses to write
  unless all of the following hold:
  - the DI and calibrated-DD identities hash to their recorded identity values
  - the current `di.stan` and `dd_cal.stan` have the sha256 the fits recorded
  - the DI Stan data, rebuilt from the data file with the historical anchor and evaluation
    regimens, reproduce the recorded `stan_data_sha256`, which proves every data-passed prior
  - the early U02 run identity reproduces from the current data, settings, priors and Stan model
  - every prior statement the page names occurs in the Stan model of the fit
  - the exported draws were produced by these fits.
- **Early identity environment.** The U02 identity also hashes the cmdstanpy and numpy versions,
  which the early manifest does not record. It reproduces exactly with cmdstan 2.40.0,
  cmdstanpy 1.3.0 and numpy 2.4.3, the versions recorded by the tumor fits. Other numpy versions
  tried (2.4.2, 2.3.5, 2.2.6) do not reproduce it. A sha256 match does not happen by chance.
  These versions are stated in the generator.
- **Fail closed in the app.** `load_model_methods()` compares the record's draw hashes with the
  manifests of the draws the app loads. On a mismatch, or a missing or empty record, the page shows
  a warning in place of the technical tables and keeps the summaries.
- **Rendering.** The page uses `markdown::mark(..., "+latex_math")` with `withMathJax()`, as Model
  and sources does, with its own table and equation styles. Technical details sit in `<details>`
  blocks so a clinician can stop after the summary.
- **Reader checks.** The generic reader-first checker expects a manuscript structure and has
  empty inventories, so it is not used as the acceptance test (review F5). Acceptance is by the page
  tests and the browser walkthrough on desktop and phone.
- **Version.** 4.6.1, a patch: no number, list, draw or export changes.

## 4. Late toxicity

Only step L0 of the review's late plan is implemented: `docs/development/late_toxicity_protocol.md`
fixes the objectives, the endpoint tuples, the rules and the decisions D1–D4 for Igor Shuryak before
anything is extracted. Steps L1–L11 (source inventory, double extraction, analysis sets, evaluator,
identification, model, fit, validation, decision impact, integration) wait for those decisions and
run on a separate branch, as the review recommends. The page states that late toxicity is not
fitted.

## 5. Tests

| Test | What it pins | Seen red |
|---|---|---|
| `tests/python/test_model_methods.py::test_committed_record_is_current` | the committed record equals a rebuild from the sources | yes: one prior string edited in the record |
| `::test_a_changed_prior_bound_or_setting_is_refused` (9 cases) | a changed Tk bound or prior, DD onset penalty, DD identity, DI α/β prior width, early α/β or chemotherapy prior width, early arm-spread prior, or early repair half-time is refused | yes: the early-identity check disabled turned 3 cases red |
| `::test_a_changed_late_assumption_makes_the_record_stale` | the late α/β range on the page follows `R/uncertainty.R` | by construction (`--check` fails) |
| `::test_a_described_prior_missing_from_the_model_is_refused` | the page cannot describe a prior the model does not contain | by construction |
| `::test_record_states_the_active_settings` | the review's table 3.4 values, 45 arms in 21 units, 1,000 DI draws × 4 | not needed (fixture check) |
| `tests/testthat/test_v461_methods_page.R::"every route reaches the one methods page, at its own section"` | menu entry, panel value, five links and their target sections (expected sections written in the test) | yes: red against the 4.6.0 `app.R` (2 errors), and red with the early link mapped to the tumor section |
| `::"the page shows the active priors and settings from the record, with no placeholder left"` | every prior of the record appears, calibrated-DD wording, late not fitted, briefs embedded, download link | yes: red against the 4.6.0 `app.R` |
| `::"the early summary's stated priors agree with the record"` | the early brief and the record agree on all six priors | not needed (fixture check) |
| `::"a record that does not match the deployed draws is refused, and the page says so"` | fail-closed loading and the warning | not needed (fixture check) |
| `::"a placeholder the record cannot fill stops the page instead of showing braces"` | no literal `{{...}}` reaches a reader | not needed (fixture check) |
| `tests/testthat/test_deploy_preflight.R::"the bundle check requires the methods page text and record, present and not empty"` | missing `how_fitted.md`, record or `R/methods_page.R`, and empty or unreadable copies, are refused | yes: fails against the 4.6.0 `check_bundle` (no `PAGE_CONTENT`) |
| `test_v45_review.R`, `test_v46_app_switch.R`, `test_v27_colleague_review.R`, `test_v31_stage_e.R` | seven destinations, routes of the existing links, version 4.6.1 | updated ripple |

## 6. Verification

2026-11-12, working tree on base `52f0d3b` (uncommitted, so the suite marker is not written): `Rscript tools/run_tests.R` gave 6,786 expectations, 0 failed, 0 errors, 0 skipped (4.6.0 baseline: 6,701). `python3 -m pytest -q -m "not sampling" tests/python`: 324 passed, 1 skipped, 2 failed, the 2 failures being `test_early_tox_step5.py::test_run_identity_changes_with_each_defining_input` and `test_early_tox_step6.py::test_cache_without_identity_or_with_another_identity_is_refused`, which call `cmdstan_path()` and fail because this computer has no CmdStan installation (they do not touch the changed files). `python3 tools/build_model_methods.py --check`: current. `stale_constant_audit.py --project . --project-root . --strict`: 0 errors (31 constants). Local browser walkthrough `docs/development/2026_10_28_v4/walkthrough_v4.R`: 73 of 73 (`walkthrough/20261112_123953_local`), including every route to the page, the back link, rendered tables and equations, and the phone view. Every new test was seen red against the unchanged code or a mutation. Not deployed: the full suite must be rerun on the committed revision before `tools/deploy.R`, and deployment needs Igor Shuryak's approval.
