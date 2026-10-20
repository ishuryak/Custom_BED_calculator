# Version 2.7: a colleague's comments, the plan, and its independent review

Date: 2026-10-20. Branch `v2-toxicity-hypofractionation`, starting revision `0c51ae6` (version 2.6.0, deployed 2026-10-19 from `cd5d835`).

Supporting files in `docs/development/2026_10_20_v27/`: `plan_original.md` (the plan as first written), `independent_plan_review.md` (the review of that plan), `make_v26_replay_fixtures.R`, `red_before_change.log`, `mutation_checks.log`, `walkthrough_v27.R` and its run folders under `walkthrough/`.

## 1. The comments

A colleague reviewed the live 2.6 calculator and its code on 2026-10-20. He raised four points:

1. The GitHub repository opens on version 1 (branch `main`). The live 2.6 code is on `v2-toxicity-hypofractionation`.
2. The default candidate search accepts a regimen when either the DI or the DD model supports it and ranks by the more favorable result. The three default regimens therefore look favorable under DI but less favorable under DD. For trial planning, agreement between both models might be the better default.
3. "15 of 15 scenarios" could read as a probability or as 15 independent confirmations. "15 of 15 selected parameter settings" would be plainer.
4. The limitations (exploratory parameters, and no modeling of chemotherapy, patient characteristics, treatment volumes or organ doses) are explained, but a short version near the main results would help.

All four were confirmed against the code at `0c51ae6` (plan, section E1–E4):

- The repository's default branch is `main` (`gh api`), and the app had no source link at all.
- Defaults were `criterion = "either"`, `rank_by = "best"`. At the defaults the search listed 29 × 2.25, 28 × 2.30 and 30 × 2.20 Gy, with ΔDI +1.19, +1.13, +1.12 and ΔDD −2.41, −2.86, −2.05 Gy₁₀.
- The agreement cell read "15 of 15", and the Sensitivity summary read "evaluable scenarios".
- The Results and Model agreement tabs had no limitation line.

## 2. The plan and its independent review

The plan (`plan_original.md`) proposed targeted changes E1–E4 and a 2.7.0 release with unchanged selection rules. An independent review (`independent_plan_review.md`) checked the plan against the code. It reproduced the numerical claims in Python (maximum discrepancy 1.8e-14 Gy) and ran the 2.6 suite (2543 expectations, 0 failed). It found nine problems in the plan. Each was verified here before being acted on.

| Finding | Verified here | Action |
|---|---|---|
| F1. `caveat_strip()` placed "beside `context_strip()`", which is inside `server()`, so the UI could not find it | Yes: `context_strip` is defined at `app.R:653`, inside `server` (line 271). | Defined at top level after `reason_text()`. A mutation that moves it into `server` turns 10 assertions red (`mutation_checks.log`, M6). |
| F2. A `testServer()` run cannot see UI defaults, because its inputs are NULL until set | Yes: inputs come only from `session$setInputs()`. | The defaults are tested in the rendered UI (`<select>` selected options), in the server with an explicit `opening` fixture, and in a fresh browser session. The legacy `set_defaults()` is left as it is. |
| F3. An empty table cannot test the ranking, because "all rows pass" is true for zero rows | Yes. | Separate checks: the selected ranking in the UI, and a −3 Gy fixture whose order differs by ranking (30, 29, 28 conservative vs 29, 28, 30 larger change). |
| F4. Overstatements: "DI or DD lists regimens that only one model supports" (OR is inclusive), "which model is true", "both models support", and "the 15 settings share α" | Yes. "Either" passes regimens meeting the threshold under both, and α differs between presets and C is recalibrated per setting. | The review's replacement text is used in the app message, About and the Sensitivity summary. |
| F5. "Chemotherapy assumed equal in both regimens", "gaps not modeled" and "no uncertainty" were imprecise | Yes. No chemotherapy term exists to be equal, a typed longer overall time does change BED, and uncertainty exists but is not calculated. | Caveat line and About limitations reworded. |
| F6. Missed documentation: About line 67 (default example) and line 84 ("no regimen is supported by both models without a loss"); guide statements that α has no DI effect before onset and that adding a chemotherapy term double-counts it; README history saying version 1 calibrated at the comparison | Yes, all present. The guide's pre-onset claim is false with g > 0 and, slightly, through the smooth threshold. The double-count claim treats the HYPNO cap as containing a chemotherapy term, which it does not. README lines 64 and 95 misdescribed version 1. | All corrected. |
| F7. "Old files replay exactly" was too broad | Yes. Only compatible 2.6 files replay. Older files are refused unless migrated. | Release text says "compatible 2.6 scenario files replay with the same numerical results and selection rules". Tested with four frozen 2.6 files. |
| F8. The copied v2.6 walkthrough hard-codes 2.6, resets to either/best before its "default" checks, and never fails the process | Yes. | New `walkthrough_v27.R`: local or live mode, opening-page checks first, explicit fixtures, nonzero exit on any failure or timeout, cleanup in `on.exit`. Failure path checked with `--self-test-fail` (exit 1). |
| F9. Release order ran final checks before docs were finished, omitted `DESCRIPTION`, and the deploy script takes the first shinyapps.io account | Yes (`tools/deploy.R` uses `acc$name[1]`). | Docs finished before final checks. `DESCRIPTION` is 2.7.0. Deployment is not part of this change (section 5). |

What was considered and not done:
- Changing the defaults of `suggest_regimens()`. Direct callers would change behavior while the app passes explicit values.
- Renaming export fields (`n_scenarios` and others). This would break replay and change scenario IDs.
- Switching the repository's default branch. A README banner on `main` solves discoverability and keeps version 1 where users expect it.
- A −3 Gy default minimum. It would put an allowed loss in front of every user only to fill the opening table.

## 3. Decisions

On 2026-10-20 Igor Shuryak asked for the plan to be merged with the independent review and implemented. Both recommended the same answers to the four open decisions, and those were taken:

1. Repository entry point: a source link on the Model and sources page, plus a README-only banner on `main`. The default branch is unchanged.
2. Opening Candidate screen: both models, smaller change, minimum +1 Gy. It shows an empty table with an explanation.
3. Display wording: "selected parameter settings". Machine-readable `scenario` names are unchanged.
4. Deployment: not done in this change. It waits for Igor's local inspection with `run_app.bat`, then `tools/deploy.R` with the account checked, then the walkthrough in live mode.

## 4. What changed

- `app.R`:
  - top-level `caveat_strip()`, shown on the Results, Candidate and Model agreement tabs
  - `selected = "both"` and `selected = "worst"`
  - a new empty-result message under "both models"
  - agreement help text and the table labels "settings" and "Selected parameter settings with both ≥ allowed change"
  - Sensitivity help and summary wording, and the "Preset" column.
- `about.md`:
  - the source link
  - the default-screen text (unavailable DD, the Tumor and Output bullets, the default example)
  - "settings" in the agreement description
  - the scoped statement on nonnegative changes
  - the chemotherapy and gap limitations.
- `R/bed_models.R` (`APP_VERSION` only) and `DESCRIPTION`: 2.7.0. `SELECTION_VERSION` 2.6, `EQUATION_VERSION` and `SCENARIO_SET_VERSION` are unchanged.
- `README.md`: the release-candidate note, the 2.7 changes, and the corrected version 1 calibration history.
- Tests:
  - new `tests/testthat/test_v27_colleague_review.R`, with fixtures `tests/testthat/fixtures/v26_*.json` frozen from 2.6
  - two existing assertions updated because they pinned display wording that changed on purpose: `test_app_server.R:290` ("scenarios under DI" became "parameter settings under DI") and `test_v25_review_fixes.R:83` ("12 of 15 scenarios" became "12 of 15 settings").

## 5. Evidence

All runs use Windows R 4.5.1 with the release library (`tools/env.R`).

- **Baseline** at `0c51ae6` (clean tree): `tools/run_tests.R` gave 2543 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings.
- **Red before the change:** `red_before_change.log`, the new tests against 2.6 code.
  - Feature tests failed on their own assertions, for example "Expected `selected_value(html, "criterion")` to equal "both"; x[1]: "either"". The failing tests were the UI defaults, caveats, source link, empty-result message, settings wording and version.
  - Compatibility tests passed, as they should: the −3 Gy order, either/best rows, the DD-unavailable message, agreement rows, and the 2.6 replay (117 assertions).
  - Two test-helper bugs were found and fixed on this run. A regex did not span line breaks, and the tab position was read from the nav link instead of the panel.
- **Mutations:** `mutation_checks.log`. Each restores one old behavior, runs the new tests and restores the file:
  - M1 removing the source link fails E1
  - M2 criterion "either" fails the default test
  - M3 ranking "best" fails the default test, independently of the empty table
  - M4 the old count label fails E3
  - M5 removing the agreement caveat fails E4
  - M6 the helper in server scope fails 10 assertions, including the startup test.
- **Full suite** after the change: see `docs/FIXES.md`, version 2.7, for the final count and revision.
- **Browser walkthrough** (`walkthrough_v27.R`, local mode):
  - the run with `--self-test-fail` exited 1, as intended (`walkthrough/20261020_111914_local_selftest_fail`). That run also exposed two timing faults in the script, which were then fixed
  - final local run: 52 of 52 checks, exit 0. The checks cover the opening page before any input, all tabs, the −3 Gy and either/best fixtures with download and replay, the agreement download, the DD-unavailable state, the chemotherapy box (Fowler zone changes, BEDs and the HYPNO limit do not), the Sensitivity wording, 420-px width, reload and the browser console.
