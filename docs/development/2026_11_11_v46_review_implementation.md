# Version 4.6 early-toxicity review: what was implemented and why (2026-11-11)

Source: `HNSCC_BED_v46_early_toxicity_review_and_action_plan_2026-11-11.md`, an independent review of the 4.6.0 branch, kept outside the repository with the first review. Igor Shuryak asked on 2026-11-11 to "analyze and check and implement this review", then commit and push. Deployment to the live app was not authorized and has not been done. Evidence rows: `docs/FIXES.md`, section "second independent review of the early-toxicity fit, R1 to R8".

## Findings checked before implementing

Each finding was reproduced or checked against the code before it was changed.

| ID | Reviewer's claim | Check | Verdict |
|---|---|---|---|
| R1 | Banner implies a risk meaning for +5, the chemotherapy term is overstated, and the ascertainment count is confused | Read `app.R`, `early_fit_in_brief.md`, `about.md`, and counted families with any post-RT visit (12) and meeting the coverage rule (2) in `data/early_tox/source/ascertainment.csv` | Confirmed |
| R2 | Trial list hides schedule variants and qualifiers | Compared `draws/early_fit_trials.csv` with `data/early_tox/source/extraction_map.csv` and the papers (CAIR, EORTC 22851, DAHANCA, HYPNO, ORO via Olmi 2003) | Confirmed |
| R3 | Repair uses one calendar for every regimen | Reproduced in R: the same m and T on different calendars gave the same repair term | Confirmed |
| R4 | Exporter accepts non-primary profiles and writes in place | Ran the exporter on U04/U15/U17/U19 manifests | Confirmed |
| R5 | Runtime loads any package whose draws hash matches | Built packages with a failed gate and with shifted medians, and both loaded | Confirmed |
| R6 | E-BFMI and fixed-constant gates are incomplete, and caches are unbound | Recomputed gates for all 114 retained fits from the raw chains in `~/early_tox_chain_archive` | Confirmed. 111 pass, and 3 known binomial fits fail (not used by the app) |
| R7 | A count can be confirmed by a different endpoint's row | Traced the ELAN pharyngeal row | Confirmed |
| R8 | Model identity is missing on other pages, and the Fowler fields are editable but unused | Shiny server probe: setting α/β 50, onset 30 and doubling time 10 in research mode changed no result | Confirmed |

Investigated and dismissed: none of the eight. The review's clinical-plausibility section (shift of the research list to 26 to 30 fractions) and its list of accepted choices were checked against the app and need no change.

## Choices made while implementing

These were open in the review or arose during the work. None changes a fitted number.

1. **Typed overall time that the calendar cannot give (R3 item 3).** The review allowed either an explicit gap pattern or "unavailable with a clear explanation". The app reports early BED as not available and says why under the time field, pointing to the skipped-days field, which already describes gaps exactly. Inventing a gap pattern would put fractions on days nobody chose. Replay refuses a file whose recorded calendar cannot give its recorded time.
2. **Fowler fields in research mode (R8).** Hidden, not disabled, and replaced by a read-only list of the fitted values. A disabled field still looks like a setting in use. Their values are kept, so switching back restores what was entered.
3. **Prior-use record in exports (R8).** Added only to research-mode Uncertainty files. Adding it to Fowler files would change the IDs of every Fowler export without telling a reader anything new, since Fowler mode uses the priors as listed.
4. **Search candidates keep no skipped days.** As before and as the review asked: candidates use the sidebar's frequency and start day without skips. The comparison uses its own calendar (preset comparators carry theirs, for example HYPNO's 66/33 control at 6 a week).

## Defects found during the work

- `compute_all()` already used a local variable `cal` for the DD calibration. A new argument of the same name was overwritten, which surfaced as "$ operator is invalid for atomic vectors". Renamed the argument `early_cal`.
- `suggest_regimens()` built the candidate calendar in Fowler mode too, which needs `default_calendar()` from `R/schedule.R`. Engine-only tests (`test_bed_models.R`, `test_v21_changes.R`) do not load that file and failed with 7 and 3 errors. The calendar is now built only in research mode.
- A mutation run of the exporter (R4) had earlier published dummy files over the real package. `publish_package()` now takes a `root`, tests use temporary folders, and the package was regenerated with byte-identical draws.

## Measurements

- Fixture values (`tests/testthat/fixtures/early_fit_evaluator.json`, independent Python evaluator, unrounded posterior medians): 20 × 2.75 Gy, skipped days 4, 5, 8, day 31: 31.11301. 33 × 2 Gy at 6 a week, day 38: 29.73564. At the app's 4-decimal point values the R engine gives 31.11405 and 29.73699, the reviewer's values to 1e-6.
- Uncertainty, research mode, 20 × 2.75 Gy (skips 4, 5, 8) against 33 × 2 Gy at 6 a week, early limit comparison + 5: toxicity share 0.84325 over 4,000 draws, replayed exactly from the saved file.
- Chain recheck: 114 fits. Output `docs/development/2026_11_10_early_toxicity/chain_recheck_2026_11_11.csv`.

## Still open

- A browser walkthrough of both model states on the final commit (review "local browser checklist").
- Deployment, which needs Igor Shuryak's approval.
- A pure-Python CI job for the fitting-pipeline tests (review step 6). CI runs only the R suite.
