# Interface clarity changes and input-effect audit (version 2.5.1, 2026-10-19)

Requested by Igor Shuryak on 2026-10-19, after using the local copy of version 2.5:

1. The legacy checkbox "Legacy: calibrate at the comparison regimen instead (versions 1 and 2.0)" and the note under it ("DD and DI tumor BEDs are equal at the anchor (35 × 2 Gy in 47 days). DD repopulation there starts on day 23.5. Changing the comparison regimen does not change the DD model") were unclear: what they mean, what they do, and whether the option is needed.
2. "Parameter preset" gave no explanation of what the default means.
3. Make the changes, then review and audit the app again, save the results, and check every menu box of the web app for the intended effect on the outputs.

No result changes. Every number the app computes with its own inputs is unchanged; the selection version stays 2.5.

## 1. The legacy coupling option

**What it did.** The DD model has one constant ($C$, the cell-kill level at which accelerated repopulation starts) that the user does not enter. It is solved so that the DD and DI models give the same tumor BED for one regimen, the anchor (default 35 × 2 Gy in 47 days). Since version 2.1 the anchor is a separate input. Versions 1 and 2.0 instead used the comparison regimen as the anchor, so choosing another comparison (for example HYPNO's control arm) changed the DD model itself: 55 Gy/20 fractions had DD BED 58.02, 64.18 or 70.12 Gy depending on the comparison (`docs/FIXES.md`, version 2.1, item 2). The checkbox restored that behavior.

**Is it needed? No.** The decision to remove it was taken here under the instruction to "make these changes". Reasons:

- It reproduces a behavior that version 2.1 removed as a defect, because a comparison choice should not change a model.
- The same calibration is available without it: enter the comparison regimen in the anchor fields. The test `test_v25_review_fixes.R::version 1 calibration is reproduced by entering the comparison regimen as the anchor` checks that DD equals DI at 20 × 2.75 Gy in 26 days when it is entered as the anchor. One difference: the anchor's overall time is always computed for 5 fractions per week, so a comparison with another schedule (HYPNO's control arm, 6 fractions per week, 38 days) cannot be entered exactly as an anchor. The removed checkbox would have used 38 days.
- Exported scenarios from version 2.5 that recorded the legacy mode (`anchor.mode = "coupled"`) still replay in R, because the replay reads the saved anchor regimen. `R/export.R` is unchanged.

**Alternatives considered.** (a) Keep it with a clearer label: rejected, because a second way to do what the anchor fields already do adds a control that most users should not touch. (b) Move it into a collapsed "reproduce older versions" section: rejected for the same reason. Removal is easy to reverse from git history (`953b284`).

**Code removed.** The checkbox, `coupled()`, `anchor_for()`, `anchor_mode()` and the legacy branches of the anchor note and the agreement note in `app.R`. Every analysis now uses `fixed_anchor()`.

**New anchor note** (sidebar, "DD calibration anchor"), with the defaults:

> The DD model has one constant that is not entered: it is set so that DD and DI give the same tumor BED for this anchor course (35 × 2 Gy in 47 days). The two models therefore agree for the anchor course and differ only for other schedules. With the current tumor values, DD accelerated repopulation in the anchor course starts on day 23.5 (DI: day 28). The anchor is separate from the comparison regimen: changing the comparison changes the differences shown, not the DD model. To reproduce version 1, which calibrated at the comparison regimen, enter that regimen here.

Day 23.5 comes from `dd_calibration()` at the D18 values and is pinned by `test_v25_review_fixes.R::the UI has no legacy coupling switch and the anchor note explains the anchor`.

## 2. Parameter preset

The menu said "Scenario D18 (2026-10-18 proposal, not fitted)" and gave no explanation, and the first entry was named "Calculator defaults (version 1)" although it is not the default. Changes:

- Labels: "Version 1 values (defaults of the published version 1 website)" and "Scenario D18 (default): α 0.35/Gy, doubling time 3 d, onset day 28".
- A note under the menu that changes with the selected preset. For D18: "The default. Scenario D18 was chosen by the developer on 2026-10-18 and is not fitted to outcome data. Its DI values are within published ranges: repopulation at 0.23 per day (doubling time 3.0 days) from day 28, a loss of 0.66 Gy₁₀ of tumor BED per day, at the low end of published estimates. Its DD maximum rate (0.40 per day) is a constructed value, not a published DD estimate. See "Model and sources"." Each other preset gets its source values, and Withers, Hartley and HYHOPE say that they keep the DD maximum rate entered below. All presets end with: "A preset fills the tumor values below; editing a value switches the preset to Custom."

Numbers checked against the code: ln 2 / 0.23 = 3.01 days; 0.23 / 0.35 = 0.657 Gy₁₀ per day; version 1, 0.2 / 0.2 = 1.0 Gy₁₀ per day. The range "0.66–0.84 Gy₁₀/day" and the D18 provenance are in `about.md` and `References/README.md` section G.

## 3. Tests and checks

| Check | Result |
|---|---|
| New or rewritten tests | `test_app_server.R::comparison regimen does not move the DD model; only the anchor does`; `test_v25_review_fixes.R::version 1 calibration is reproduced by entering the comparison regimen as the anchor`; `::the UI has no legacy coupling switch and the anchor note explains the anchor`; `::the preset menu explains the default and every preset`. The old `N6: legacy coupling…` test was removed with the feature. |
| Red without the change | The three tests that check the change failed when run against `app.R` at `953b284`, each on its assertion about the changed behavior. The fourth test (version 1 calibration through the anchor fields) checks behavior that already worked and is a guard. |
| Full suite | `testthat::test_dir("tests/testthat")`, Windows R 4.5.1: **1670 expectations, 0 failed, 0 errors** (baseline at `953b284`: 1651). |
| Mutations | 18 of 18 caught (`docs/development/2026_10_19_review_scripts/mutations_output.txt`). M9 now restores an agreement anchor that follows the comparison; M15–M18 are new: the comparison moves the DD model, the preset note says nothing about the default, the default is not labeled, and the anchor note loses its explanation. |
| Browser walkthrough | `walkthrough.R`; result in `docs/development/2026_10_19_walkthrough/walkthrough_log.txt`. The legacy check now asserts that the switch is gone and that the anchor and preset notes render. **25 of 25 passed.** It ran on the final interface code a minute before the version string was raised to 2.5.1, so its screenshots show "version 2.5"; the suite was rerun after the bump (1670, 0 failed). |
| Input-effect audit | See section 4. |

## 4. Input-effect audit: every menu box

Script: `docs/development/2026_10_19_input_effects/input_effects.R`. Log and per-case table: `input_effects_log.txt` and `input_effects.csv` in the same folder.

Method. A headless Chrome drives the running app. For each case it reloads the page (all defaults), opens every sidebar panel so that no output is suspended, optionally applies base settings, takes a snapshot, changes one menu box through its own widget (number field, selectize menu, radio button, checkbox, slider), and takes a second snapshot. A snapshot holds every output on the four working tabs (each Results row separately), the plot image, and both ZIP downloads split into result rows and the saved scenario file (without its export time). A case passes when every output named as required changed, no output named as must-stay-the-same changed, no output shows an error, and something changed. "Must stay the same" encodes the independence rules: for example the Results regimen must not change the candidate or agreement rows, the comparison regimen must not change the anchor or the agreement (whose default comparison is HYPNO), candidate-only inputs must not change the agreement, the anchor must not change the DI or toxicity rows, and agreement inputs must not change the candidates.

Coverage: all 38 inputs of the interface (the Model tab has none): Results regimen (3), comparison regimen (5), tumor (7 including the preset, each preset tested), anchor (2), early toxicity (5 including chemotherapy, tested under the HYPNO and the Fowler limits), late toxicity (1), candidate tab (11, every choice of each menu), agreement tab (3).

First full run: 42 of 46 cases passed. The four failures were checked one by one and are expectations of the audit script that were wrong, not app defects:

| Case | Observed | Why it is correct |
|---|---|---|
| Fraction range 15–30 | candidate rows unchanged | The default list is 29 × 2.25, 28 × 2.3 and 30 × 2.2 Gy, all inside 15–30 fractions. |
| Highest dose 2.75 Gy | candidate rows unchanged, agreement rows changed | Candidates are 2.2–2.3 Gy; the agreement list (2.80–3.05 Gy) loses its rows. |
| Minimum days saved 3 | candidate rows unchanged | 31 × 2.2 Gy (43 days) and 32 × 2.15 Gy (44 days) gain 1.83 and 1.56 Gy₁₀ under DI and meet the HYPNO acute limit, but their late BED (118.2 and 118.1 Gy₃) exceeds the default late limit of 116.7 Gy₃ (the comparison's). Without the late limit the list has 5 rows. |
| Early limit "comparison regimen" | agreement rows unchanged | The agreement tab compares with HYPNO itself, so "comparison plus 0" is the HYPNO acute limit. |

Those cases were rewritten to values that cut into the default list (fractions 29–30, highest dose 2.25 Gy, minimum days 3 with a 2 Gy late excess). Two cases were added: highest dose 2.75 Gy now requires the agreement rows to change and the candidate rows not to, and re-entering the default minimum days must change nothing. That gives 48 cases.

The second full run passed 48 of 48, but the case that must change nothing showed a changed scenario file. The cause was the audit, not the app: the version was raised to 2.5.1 in `R/bed_models.R` at 13:00:43, during the run, and the scenario's code fingerprint (`code_identity()`) hashes that file at export time. Two direct probes showed identical scenario files before and after re-entering the value, within one session and across three separate page loads. A third run stopped after 25 passing cases on a headless-Chrome page-load timeout; the script now retries a page load up to three times. The final run (13:29, no file edited during it) passed **48 of 48**, and the re-entered default changed nothing. It is the one saved in `input_effects_log.txt`.

Observations from the audit, with no change made:

- Chemotherapy, with the default HYPNO early limit, changes only the Fowler zone text on the Results tab and the saved scenario. This is intended (no BED includes chemotherapy) and the checkbox's help text says so.
- The criterion "DI model" gives the same list as the default "DI or DD", because with the defaults no regimen gains under DD.
- Two defects of the audit script itself were found and fixed during the work: a ZIP hash that mixed result rows with interface context saved in the scenario, and a slider command whose return value Chrome could not pass back.

## 5. Version

`APP_VERSION` 2.5.1 (interface only). `SELECTION_VERSION` stays 2.5, so exports and replays from 2.5.0 remain valid. Not deployed.
