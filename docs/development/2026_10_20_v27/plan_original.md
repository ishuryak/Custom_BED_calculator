---
title: "Plan: BED calculator v2 changes from a colleague's comments"
author: "Igor Shuryak, MD, PhD, is144@cumc.columbia.edu"
date: "October 20, 2026"
---

# Scope

Four comments from a colleague on the live calculator (version 2.6.0, https://ishuryak.shinyapps.io/custom_bed_calculator_v2/). Each was checked against the code on branch `v2-toxicity-hypofractionation` (HEAD `0c51ae6`, working tree clean). Nothing below has been changed yet. Only the locations listed are touched. Equations, parameter defaults, toxicity limits and selection rules are unchanged.

Proposed release: **version 2.7.0**. A default that changes what a fresh page shows is more than a patch, although no rule changes. `SELECTION_VERSION` stays at 2.6, because the selection rules are unchanged and scenario files record the criterion and ranking explicitly. Old scenario files therefore replay exactly.

# Decisions needed (one round, recommendation first)

1. **E1, GitHub.** (a) Add a short "Version 2 is on branch v2-toxicity-hypofractionation" banner to the README on `main`, and a "Source code" link on the app's Model and sources page *(recommended)*. (b) Make `v2-toxicity-hypofractionation` the default branch, plus the same app link. (c) The app link only. Option (a) keeps v1's code and history as they are and adds one README commit to `main`. Option (b) is reversible but changes what every visitor and clone sees first.
2. **E2, default minimum change.** With the two-model default, should the default minimum tumor BED change stay at +1 Gy *(recommended)*? At the defaults the table is then empty, and the message points to the Model agreement tab. The alternative, −3 Gy (the same as the agreement tab), lists the same three regimens as today, now labeled honestly (30 × 2.20, 29 × 2.25, 28 × 2.30 Gy, ΔDD −2.05 to −2.86). It also makes "allowed loss" the first thing a user sees.
3. **E3, wording.** Use "selected parameter settings" in the table and summaries *(recommended)*. Keep the internal and export name "scenario" so that files, column names and scenario IDs do not change.
4. **Deploy.** After review: redeploy with `Rscript tools/deploy.R` *(recommended only after you look at it locally with `run_app.bat`)*.

# E1. The GitHub link opens version 1

**Verified.** `gh api repos/ishuryak/Custom_BED_calculator` gives default branch `main`, which holds v1 (`75d96b8`). The app itself has no GitHub link anywhere: `grep -i github app.R about.md` finds nothing. the colleague probably opened the repository directly.

| Item | Location | Current | Proposed |
|---|---|---|---|
| E1.1 | `about.md`, the `{{VERSIONS}}` paragraph (line 169) | "Version 1 remains available at https://ishuryak.shinyapps.io/custom_bed_calculator/. …" | Prepend: "Source code: https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation (version 2). The repository's default branch holds version 1." |
| E1.2 (option a) | `README.md` on `main`, top | (v1 README, no mention of v2) | One line under the title: "**Version 2** (live at https://ishuryak.shinyapps.io/custom_bed_calculator_v2/) is on branch [`v2-toxicity-hypofractionation`](…/tree/v2-toxicity-hypofractionation). This branch holds version 1, which carries a known DD calibration error at tumor α/β other than 10 Gy." A separate one-commit change on `main`. No v1 code is touched, and the live v1 app is not redeployed. |
| E1.2 (option b) | GitHub settings | default `main` | `gh repo edit ishuryak/Custom_BED_calculator --default-branch v2-toxicity-hypofractionation` |

**Test.** `test_app_server.R::E1 about page links the version 2 source branch`: render the About output and assert that it contains the branch URL. It is red now (no link) and green after the change.

# E2. The default search credits a regimen that only one model supports

**Verified.** The defaults are `criterion = "either"` and `rank_by = "best"` (`app.R:197-202`). Rerunning the default search with `suggest_regimens()` (R 4.5.1, 2026-10-20):

| Setting | Rows |
|---|---|
| Current default (either, larger change, min +1 Gy) | 29 × 2.25 Gy (ΔDI +1.19, ΔDD −2.41), 28 × 2.30 (+1.13, −2.86), 30 × 2.20 (+1.12, −2.05) |
| Both models, smaller change, min +1, 0 or −1 Gy | none |
| Both models, smaller change, min −3 Gy | 30 × 2.20, 29 × 2.25, 28 × 2.30 (same regimens, conservative order) |

So the colleague is right: every default row gains under DI and loses 2–3 Gy₁₀ under DD.

| Item | Location | Current | Proposed |
|---|---|---|---|
| E2.1 | `app.R:197-202`, the `criterion` and `rank_by` inputs | first choices "DI or DD model" and "larger of ΔDI and ΔDD (optimistic)", which are therefore the defaults | `selected = "both"` and `selected = "worst"`. The order of the choices is unchanged. |
| E2.2 | `app.R:711-716`, the `suggest_msg` empty-result text | "...Options: require the change under 'DI or DD model', allow a tumor BED loss…" | When `criterion == "both"`: "No regimen meets the minimum change under both models within the limits. The Model agreement tab shows which regimens both models support when a tumor BED loss is allowed (default −3 Gy against HYPNO). You can also allow a loss here (a negative minimum change), lower the days saved, or choose another limit. 'DI or DD model' lists regimens that only one model supports." Other criteria keep the current text. |
| E2.3 | `app.R:79`, `reason_text(reason, criterion = "either")` | default argument "either" | Unchanged. The argument is only read by callers that pass it. Listed to confirm it was checked. |
| E2.4 | `R/bed_models.R:222`, `suggest_regimens(... criterion = "either", rank_by = "best")` | function defaults | **Unchanged.** Tests and scripts call the function with these defaults, and the app always passes the values explicitly. Changing them would alter every direct caller without changing the app. |
| E2.5 | `about.md:34` | "...and the default "DI or DD" criterion with ranking by the larger gain uses DI alone." | "...and the "DI or DD" criterion with ranking by the larger gain uses DI alone. The default (both models, smaller change) then returns no regimens, with an explanation." |
| E2.6 | `about.md:59` (Tumor bullet) | "...under the chosen model: DI or DD (default), both, DI only, or DD only." | "...under the chosen model: both (default), DI or DD, DI only, or DD only. Requiring both is the default because a regimen that only one model supports depends on which repopulation model is true." |
| E2.7 | `about.md:63` (Output bullet) | "The ranking is by the larger change of the two models (default, optimistic), the smaller change (conservative), or one model's change." | "The ranking is by the smaller change of the two models (default, conservative), the larger change (optimistic), or one model's change." |

If decision 2 is −3 Gy, E2.8 changes `num("min_gain", …, 1, …)` to −3 at `app.R:204`, and the about.md "(default +1 Gy)" in the Tumor bullet to "(default −3 Gy, an allowed loss)".

**Tests** (in a new `tests/testthat/test_v27_colleague_review.R`, plus `test_app_server.R`):

- `E2: a fresh session defaults to both models and the conservative ranking`: start the server with no inputs set and assert `input$criterion == "both"` and `input$rank_by == "worst"`. Red now.
- `E2: the default candidate table lists no regimen that only one model supports`: at the defaults, every listed row has `d_di >= min_gain` and `d_dd >= min_gain`. For +1 Gy, assert zero rows and that the message names the Model agreement tab. Red now (3 rows, each with ΔDD < 0).
- `E2: the empty-result message under both models points to the agreement tab, not back to 'DI or DD'`: the message contains "Model agreement" and does not start its options with "require the change under 'DI or DD model'".
- Mutation check: restore `selected = "either"`, and the first two tests go red.
- Existing tests: `set_defaults()` in `test_app_server.R:11` and `test_v25_review_fixes.R:137` set `criterion = "either"` explicitly, so they keep testing the old path and stay green. I will check each test that relies on the app's own defaults (not `set_defaults()`) and update only those, naming each in the record.

# E3. "15 of 15 scenarios" can read as a probability or 15 confirmations

**Verified.** The agreement table column is headed "Scenarios with both ≥ allowed change", with cells such as "15 of 15" (`app.R:833-842`). The Sensitivity summary reads "%d of %d evaluable scenarios under DI…" (`app.R:876-880`). The help text already says "a count of scenarios is not a chance of benefit" (`app.R:258`), but the table cell itself does not.

| Item | Location | Current | Proposed |
|---|---|---|---|
| E3.1 | `app.R:842`, column name | `Scenarios with both ≥ allowed change` | `Selected parameter settings (of 15) with both ≥ allowed change` |
| E3.2 | `app.R:833-836`, cell text | "15 of 15" / "%d of %d evaluable (%d of %d not evaluable)" | "15 of 15 settings" / "%d of %d evaluable settings (%d of %d not evaluable)" |
| E3.3 | `app.R:876-880`, Sensitivity summary | "...%d of %d evaluable scenarios under DI, %d of %d evaluable scenarios under DD." | "...%d of %d evaluable parameter settings under DI, %d of %d under DD. These are selected settings that share α and the DD calibration, not independent confirmations or a probability." |
| E3.4 | `app.R:238`, agreement help text | "...vary over the 15 scenarios of the Sensitivity tab..." | "...vary over the 15 selected parameter settings of the Sensitivity tab (five presets × three tumor α/β values)..." |
| E3.5 | Column names in the exported CSV files, `n_pass_both` and similar, the scenario set | unchanged | **Unchanged**, so exports, scenario IDs and replay are unaffected. |

**Tests:** `E3: the agreement table says 'settings', not a bare count` (the header contains "parameter settings", and the default cell is "15 of 15 settings"). `E3: the Sensitivity summary states the settings are not independent`. Both are red now. The walkthrough check that greps "15 of 15" still matches.

# E4. Key limitations shown next to the results

**Verified.** The Candidate tab's grey box (`app.R:190-194`) mentions hypothesis generation, uncertainty, safety and organ-at-risk planning, but not chemotherapy, patient or tumor characteristics, or treatment volume. The Results tab has one sentence (`app.R:174-177`). The Model agreement tab has none. The full list is only on the Model and sources page.

| Item | Location | Proposed |
|---|---|---|
| E4.1 | New helper in `app.R`, beside `context_strip()` | `caveat_strip()` returns one short grey line: "**Not modeled:** chemotherapy effects (assumed equal in both regimens), patient and tumor characteristics (HPV status, stage, volume), treatment volumes, organ-at-risk doses, treatment gaps. Default tumor parameters are exploratory, and differences carry no uncertainty interval. Details: Model and sources." |
| E4.2 | Results tab, after `tableOutput("results")` (`app.R:178`) | `caveat_strip()` |
| E4.3 | Candidate tab, directly below the grey box (`app.R:194`) | `caveat_strip()` |
| E4.4 | Model agreement tab, after `agree_table` (`app.R:253`) | `caveat_strip()` |

One shared function means one wording to maintain. The existing grey box and help texts are left as they are.

**Test:** `E4: the caveat line appears on the Results, Candidate and Model agreement tabs` checks the rendered UI for "chemotherapy" and "organ-at-risk" on each tab. Red now on the Results and Agreement tabs.

# Version, records and documents

| Item | Location | Change |
|---|---|---|
| V1 | `R/bed_models.R:11` | `APP_VERSION <- "2.7.0"`. `SELECTION_VERSION`, `EQUATION_VERSION` and `SCENARIO_SET_VERSION` are unchanged. |
| V2 | `README.md:13` | The live line: version 2.7.0, new default rows (none at +1 Gy, or the three at −3 Gy), "15 of 15 settings". |
| V3 | `docs/FIXES.md` | One row per item E1–E4 with all five fields. The first line reads: "**Changes what a fresh page shows: yes.** The default candidate table goes from 3 rows to 0 (or the same 3 at −3 Gy, conservative order). No equation, parameter, limit or selection rule changes, and saved scenarios replay unchanged." |
| V4 | `docs/development/2026_10_20_eric_wang_review.md` | This plan, the colleague's comments, the verification numbers above, your decisions 1–4, and what was ruled out (changing `suggest_regimens()` defaults, renaming export columns). |
| V5 | Guide `BED_calculator_v2_references_and_interface_2026_10_20.md` | Candidate controls table (line 339-340: defaults "both models" and "smaller of ΔDI and ΔDD (conservative)"), the minimum-change row if decision 2 is −3 Gy, the Model agreement outputs ("15 of 15 settings"), and a note on the caveat line. Rebuilt as a tracked and a clean DOCX with a `_pre_` snapshot. Changed table cells must be checked by hand, because the tracked-changes helper does not mark edits inside tables. |
| V6 | Walkthrough | A new `docs/development/2026_10_20_v27/walkthrough_v27.R`, copied from `walkthrough_v26.R`, with check 6 changed from the three default rows to the empty-table message (or the conservative order at −3 Gy), plus checks for the caveat line on three tabs and the "settings" wording. The v26 script is left as the record of 2.6. |

# Order of work and verification

1. Write the new tests (E1–E4) and run them against the current code. Each must fail, and fail on its own assertion.
2. Make the changes E1.1, E2.1-2.7, E3.1-3.4, E4.1-4.4 and V1. For each, restore the old line once and confirm its test goes red.
3. Run the full suite with `Rscript tools/run_tests.R` (baseline 2543 passed, 0 failed), then the browser walkthrough `walkthrough_v27.R` against a local app. Both must pass on the same revision.
4. Update about.md, README, FIXES.md and the development record (V2-V4), then commit on `v2-toxicity-hypofractionation`. After checking `git rev-parse --show-toplevel` is the app's own repository, push.
5. E1.2 on `main` (or the default-branch switch) as a separate commit, per decision 1.
6. Deploy only after you have looked at it locally and approved. Then run the walkthrough against the live URL.
7. Update and rebuild the guide (V5).
8. Draft a short reply to the colleague listing what changed. Not sent.

# What is not changed

- Equations, parameter defaults, presets, toxicity limits, the DD anchor and the 15 settings.
- `suggest_regimens()` function defaults and export column names.
- The live version 1 app and v1's code on `main` (option a adds a README line only).
