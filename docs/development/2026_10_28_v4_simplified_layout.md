# Version 4.0.0: one clear main page, everything else in a menu (2026-10-28)

Branch `v4`, from `c782ab0` (3.2.1 with Eric Wang's PR #1 merged). Not deployed.

## Request

Igor Shuryak, 2026-10-28: version 3 is too complex for a clinician (many tabs, much text and many numbers that need prior explanation). Make version 4 with the same engine and a simpler layout. The main page lists plausible hypofractionated regimens that meet, against the standard 35 × 2 Gy:

1. DI or DD tumor BED ≥ the standard's DI or DD BED − 5 Gy,
2. early toxicity BED ≤ the standard's early BED + 5 Gy,
3. late toxicity BED ≤ the standard's late BED + 5 Gy.

The other tabs stay reachable through a side menu the user can open.

## What was built

- `R/v4_home.R`: `v4_regimens()` calls the unchanged `suggest_regimens()` with comparison 35 × 2 Gy in 47 days, `criterion = "either"`, `min_gain = -5`, `tol_early = 5`, `tol_late = 5`, `min_days_saved = 1`, fractions 5–34, dose steps of 0.05 Gy from just above 2 Gy to 8 Gy, 5 fractions a week from a Monday. It returns every qualifying grid point (the engine's `feasible` attribute). `v4_summary()` makes one row per number of fractions: dose-per-fraction range, total-dose range, overall time, days saved, which tumor models meet condition 1 (and at which doses), and any published trial arm inside the range (HYPNO 55 Gy/20, HYHOPE 46.5 Gy/15).
- `app.R`: the version 3 controls and tabs are unchanged in content and IDs. The page is now a collapsible menu (closed at start, opened by the **Menu** button or the arrow) and a hidden tab set. The main page is the list. Clicking a row shows each dose with its tumor (DI, DD), early and late BED differences. The version 3 tabs are menu pages. The regimen and comparison inputs moved to the Report card page (two columns), the model parameters to a Model settings page. The banner shows on every page except the main one.
- Wording: every visible "sidebar" and "tab" now names the page where the control lives (app, uncertainty and organ modules, About). `R/export.R` strings are unchanged because they enter exported scenario files.
- About: a "What is new in version 4" section. Version 4.0.0 in `R/bed_models.R` and `DESCRIPTION`.
- Deploy: `tools/deploy.R` targets a new site, `custom_bed_calculator_v4`, and `deploy_preflight()` now also refuses the live version 3 site name, following the version 3 precedent (each major version is its own site).

## Decisions and alternatives

| Question | Chosen | Alternatives considered | Why |
|---|---|---|---|
| "DI or DD" in condition 1 | Either model suffices (`criterion = "either"`). The table says which model(s) and at which doses | Both models required | The request reads "DI or DD". Requiring both empties 16–21 fractions at the defaults, including HYPNO's 55 Gy/20, under the calibrated DD model. The column keeps the stricter reading visible. Confirmed by Igor Shuryak on 2026-10-28 ("you understand right"). |
| One row per regimen or per fraction number | Per fraction number, with the dose range | Every (fractions, dose) pair (≈90 rows at the defaults). One "best" dose per fraction number | A range is what a clinician chooses within. All pairs are one click away. The feasible set is a contiguous dose interval for every fraction number at the defaults and at the refit (checked). |
| Which comparison and parameters | Always 35 × 2 Gy in 47 days. Parameters and DD anchor from Model settings | Follow the Report card comparison | The request names the standard. The parameters stay adjustable so the list can be redone under the refit or another preset. |
| Margins adjustable on the main page | No (fixed at 5 Gy) | Inputs on the main page | The main page must stand alone. The Search page already offers every limit. |
| Search ceiling | 8 Gy per fraction | 5 Gy (the Search page default) | 5 Gy would truncate short courses. At 8 Gy no listed setting reaches the ceiling (late BED + 5 Gy caps 5 fractions near 7.2 Gy). |
| Menu | bslib sidebar, closed at start, plus a visible Menu button | Top navbar with a "More" dropdown | The request asks for a side menu the user opens. The small default arrow alone was easy to miss in the first screenshot, so the button was added. |

## Results at the defaults (Scenario D18, calibrated DD)

Standard 35 × 2 Gy in 47 days: tumor BED 71.51 Gy₁₀ (DI and DD, equal by calibration), early BED 53.11 Gy₁₀, late BED 116.67 Gy₃.

19 schedules, 16–34 fractions. 16–21 and 23 fractions meet condition 1 under the DI model only. 22 and 24–31 meet it under both models at part of their dose range. 32–34 meet it under both models throughout. 20 fractions: 2.65–2.80 Gy (contains HYPNO's 2.75 Gy). 16 fractions (3.20–3.25 Gy) is above the 3.1 Gy benchmark fraction size and marked.

Under the verified-trial refit preset (posterior medians, calibrated DD): 22 schedules, 13–34 fractions, both models at all doses for 18 of them. HYHOPE's 15 × 3.1 Gy is not in the list (15 fractions qualify only at 3.30 Gy).

## Verification

- `tests/testthat/test_v4_home.R` (81 expectations): the list equals an explicit brute-force evaluation of the three conditions over the same grid, at the defaults and under the refit, HYHOPE and version 1 presets. Pinned defaults (rows 16–34, 20-fraction row). Summary wording. Margins widen the list monotonically. Empty list keeps the engine's reason. The server's main page, row detail and independence from the Report card comparison. The UI opens on the main page with every menu link, page and back link, and no visible "sidebar".
- Deliberate breaks, each caught: tumor-margin sign flipped (23 failures), `criterion = "both"` (18), late tolerance 0 (10).
- Full suite and browser walkthrough: see `docs/FIXES.md`, version 4.0.0 row.
- `docs/development/2026_10_28_v4/walkthrough_v4.R`: headless Chrome, 36 checks (main page first, 19 rows, HYPNO marked, menu closed, a real mouse click shows the dose card, every menu page opens without an output error, banner only on detail pages, a setting change reaches the main list, no horizontal scroll at 390 px, no console errors). A synthetic jQuery click does not reach DT's row selection. The walkthrough uses a dispatched mouse event.

## Review of Eric Wang's PR #1 (merged into `fixes_v32_review_2026_10_26` as `c782ab0`)

Findings, kept here because they bear on releasing either 3.2.1 or 4.0.0:

1. **Report-card threshold status (correct fix).** Before, the status said "met" whenever the difference was finite, including when both BEDs were negative. `tumor_threshold_txt()` now applies `tumor_endpoint()`, so a non-evaluable comparison reads "not evaluable". Tested in `test_eric_feedback.R`.
2. **Stale-constant audit broken by the merge (release blocker for 3.2.1).** The 9 lines added near the top of `app.R` moved two registered sites (`verdict_tolerance_gy` at app.R:89, `max_days_input` at app.R:47). `stale_constant_audit.py --strict` on `c782ab0` gives 2 errors, so `.stale_constant_record.json` (2026-10-26) no longer attests to the branch. Fixed on `v4` (sites now app.R:101, 57, 805). The 3.2.1 branch still needs the same three-line fix if it is released.
3. **About text about the `v3` branch is wrong in direction.** It says the v3 branch "has the development code and may differ from the live website". The live site runs 3.2.0 from `review_2026_10_26` (4876ae4). `v3` is at 535b338 (3.1.0), behind it, and 3.2.1 lives on `fixes_v32_review_2026_10_26`. The new deployed-revision link is the reliable pointer.
4. **A default changed without a record.** The Uncertainty page now opens with the tumor criterion "both models" instead of "DI model", so its default joint shares change. The PR states "the equations, draws and search rules stay the same", which is true, but there is no `docs/FIXES.md` row saying a displayed default moved (rule 10, item 5).
5. **Deployed-revision link** points to `/tree/<sha>`. It resolves only if that commit is pushed. `tools/deploy.R` does not require a pushed revision.
6. **PR #2** (README link on `main`) is still open, not merged.
7. Suite on `c782ab0`: 5514 expectations, 0 failed, 0 errors (Windows R 4.5.1, release library). Suite marker written for `c782ab0`.

## Not done

- No deploy (Igor deploys). Committed and pushed to branch `v4` on 2026-10-28 with his approval.
- Live walkthrough of v4 (needs a deploy).
- The user guide (`docs/user_guide`) still describes the version 3 tabs.

**Correction (2026-10-28, review of 4.0.0).** The Verification section above overstated two tests in 4.0.0. "The list matches brute force under other tumor settings" compared row counts only, and "an empty list carries the engine's reason" did not check the reason. Both now check what they claim (version 4.1.0). The brute-force oracle uses the same 0.05 Gy grid as the app, so it checks the selection on that grid, not the exact dose edges. See `2026_10_28_v4_review_and_plan.md`.
