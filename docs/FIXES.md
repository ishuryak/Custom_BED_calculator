# Fixes

## 2026-10-18: version 2 (branch `v2-toxicity-hypofractionation`)

**Changes published numbers (version 1 website):** yes, for any tumor α/β other than 10 Gy. At α/β = 10 Gy (the default) nothing moves.

| Item | Record |
|---|---|
| What was wrong | Version 1 calibrates the DD model so that the DD and DI tumor BEDs are equal for 35 × 2 Gy in 47 days. The calibration hard-codes that regimen's value for α/β = 10 Gy (the constants 84 and 3948). At any other α/β, the DD tumor BED was wrong at every regimen. |
| Proof it was real | With the v1 formula at α/β = 5 Gy (other defaults unchanged), 35 × 2 Gy gives BED_DD = 73.80 Gy, while BED_DI = 79.00 Gy. The model defines them as equal there. At α/β = 15 Gy: 62.14 vs 60.33. |
| Test that pins it | `tests/testthat/test_bed_models.R::DD is calibrated to DI at the reference for every tumor alpha/beta (v1 is not)`. Confirmed failing when the hard-coded value 84 is restored (deliberate-break check, 2026-10-18). |
| Verification | `Rscript -e "testthat::test_dir('tests/testthat')"` gave 1131 expectations, 0 failed, on 2026-10-18 at this branch's first commit. Baseline before version 2: no tests existed. |
| What it moves | BED_DD for every regimen when tumor α/β ≠ 10 Gy. Example: 35 × 2 Gy at α/β = 5 Gy goes from 73.80 to 79.00 Gy. Version 1 values at α/β = 10 Gy are reproduced to a relative difference below 1e-8. |

Other changes in version 2 (new features, review fixes, and decisions) are recorded in `docs/development/2026_10_18_code_review.md`.

## 2026-10-18 (later): records only, no code change

**Changes published numbers:** no. The app code is unchanged. This entry adds the planning record: `docs/development/2026_10_18_handoff.md` (start here), `docs/development/2026_10_18_fix_plan.md` (the external-review verification, work packages and decisions D1–D18), and `docs/development/2026_10_18_dd_parameters/` (DD parameter analysis and grid search). Fixes from that plan will get their own rows when they land.

## 2026-10-18 (night): version 2.1 (plan `docs/development/2026_10_18_synthesis_plan.md`)

**Changes published numbers: yes.** Acute (early) BEDs rise by 0.79 Gy₁₀ for any course longer than 8 days (35 × 2 Gy: 52.31 → 53.11), so the fixed Fowler limit is stricter (default-parameter search under the 59 Gy₁₀ cap: 18 → 17 rows). DD BEDs change whenever the comparison regimen is not 35 × 2 Gy / 47 days (now calibrated at a fixed anchor). The Withers preset changes at tumor α/β ≠ 10 Gy. Default-setting results on the Results tab and the default 7-row search are otherwise unchanged. Nothing is deployed.

Verification for every row: `testthat::test_dir("tests/testthat")` with Windows R 4.5.1 on 2026-10-18 gave **1419 expectations, 0 failed, 0 errors, 0 skipped, 2 warnings** (shiny and DT built under R 4.5.3), at the commit that adds this entry. Baseline at `6f0209c`: 1129 passed plus the same 2 warnings. The 12 new tests in `tests/testthat/test_v21_changes.R` were run on `6f0209c` first: 10 failed (the two that passed guard behavior that must not change). Twelve mutations of the shipped lines (restoring each old line in a temporary copy) were each caught by at least one test. A headless-browser walkthrough of all four tabs, the comparison presets, the legacy option and the presets showed no output errors.

| Item | What was wrong | Proof it was real (code at `6f0209c`) | Test that pins it | What it moves |
|---|---|---|---|---|
| 1. Acute clock | Acute BED used inclusive days with Fowler's onset (day 7), but Fowler's parameters and zone count elapsed days | 35 × 2 Gy gave 52.31 Gy₁₀, Fowler 2003 Table 1 gives 53.1. 22 × 2.75 Gy in 30 days gave 58.92 and passed the 59 Gy₁₀ cap; on Fowler's clock it is 59.71 | `test_v21_changes.R::acute BED follows Fowler's elapsed-day convention (Fowler 2003, 2010)`, red on `6f0209c` | All acute BEDs +0.79 Gy₁₀; Fowler-cap search 18 → 17 rows |
| 2. Fixed DD anchor | The DD model was recalibrated at whatever comparison regimen was chosen, so the comparison changed the model | 55 Gy/20 fx had DD BED 58.02, 64.18 or 70.12 Gy for comparison 35 × 2 Gy, HYPNO control or 55/20 itself | `test_v21_changes.R::the DD anchor is separate from the comparison regimen`; `test_app_server.R::comparison regimen does not move the DD model unless legacy coupling is chosen` | DD values for non-default comparisons. Old behavior kept as a labeled legacy option |
| 3. Presets | Withers preset fixed λ_DI/α at 0.72 (valid only at α/β 10). Presets set only three values | At α/β 5 the reference DI BED was 84.32 Gy instead of 82.04 | `test_v21_changes.R::Withers preset follows the tumor alpha/beta`, `::source presets carry their published values and clocks` | Withers values at α/β ≠ 10. New presets: Hartley 2011, scenario D18. Defaults unchanged |
| 4. Search validity | A negative comparison BED was used in the ranking | λ_DI 1/day at α 0.2: comparison DI BED −11 Gy, and the search listed 26 rows with apparent gains up to +78.04 Gy | `test_v21_changes.R::a negative reference BED blocks that endpoint` | Only invalid parameter settings |
| 4b. Trial-design options | No non-inferiority margin, no HYPNO-anchored limits, no model-membership column, 8 Gy default ceiling | Minimum gain below 0 was rejected by input validation; late limit fixed to the comparison (70/35 allows 116.7 Gy₃ vs HYPNO's 105.4) | `test_v21_changes.R::absolute late cap (HYPNO 55/20) is honored`, `::search rows say which model credits them...`, `test_app_server.R::HYPNO acute and late limits reach the search`, `::a negative minimum change (non-inferiority margin) is accepted` | New options only; default d_max 8 → 5 Gy leaves the default rows unchanged |
| 5. Growth term g | The 2018 paper's background growth could not be included | With the fitted 2018 parameters, 55/20 vs 70/35 ΔDD is +1.41 Gy with g and −0.03 without | `test_v21_changes.R::background growth g subtracts gT/alpha in both models and keeps the calibration` | Nothing at default g = 0 |
| 6. Export | CSV rounded to 3 decimals, without parameters, anchor or versions, so a result could not be reproduced | Download held display columns only | `test_v21_changes.R::an exported scenario reruns to the same rows`; `test_app_server.R::download is a ZIP with full-precision results and a rerunnable scenario` | Download format (ZIP) |
| 7. Display | Late BED at one α/β only; no HYHOPE row; no sensitivity view | 20 × 2.9 Gy is below 35 × 2 Gy at α/β 3 (114.07 vs 116.67) and above it at α/β 2 (142.1 vs 140) | `test_app_server.R::results show late BED at alpha/beta 2 and 4, and the sensitivity tab renders`, `::HYHOPE benchmark row and HYPNO context are shown` | Display only |
| 8. Text and records | HYPNO's 0.65 Gy/day was treated as Gy₁₀ (it is EQD2, 0.78 Gy₁₀); incomplete-repair direction reversed in the analysis notes; Gupta 6.4 Gy attributed to Stuschke; page claimed the Fowler full text was unread | HYPNO p. 389: 66 Gy + 9 × 0.65 = 71.85 Gy in 2-Gy equivalents | Text; `about.md`, `References/README.md`, `docs/development/2026_10_18_dd_parameters/ERRATA.md` | No numbers in the app |

Reasoning, the two reviews and the decisions taken without the PI (A1–A6): `docs/development/2026_10_18_synthesis_plan.md`.

## 2026-10-19: version 2.2 (independent review of 2.1, Fenwick option, CI)

**Changes published numbers: no** for default settings. The search tables change only where a model's comparison BED is negative (its differences now show as NA instead of spurious numbers). Nothing is deployed.

Verification: `testthat::test_dir("tests/testthat")`, Windows R 4.5.1, 2026-10-19: **1462 passed, 0 failed, 0 errors, 0 skipped, plus 2 warnings** (shiny and DT built under R 4.5.3), at the commit that adds this entry. (Earlier rows report the testthat result count, which includes the 2 warnings: 1419 = 1417 passed + 2 warnings.) The new review-fix tests were run against the 2.1 model code (`6225dd4`) first: the three that test changed behavior failed; three others pin 2.1 behavior that was correct but untested. Mutation checks: the seven mutations the reviewer found surviving were rerun, and six are now caught; the seventh (dropping the candidate-positivity rule for DI) cannot be reached with realistic inputs, because a shorter course with larger fractions cannot lose more DI BED than its comparison, so it stays a defensive guard. Four Fenwick mutations (Meade's rounded constant 31.8, option ignored, 50-day domain not enforced, app not passing the option) are each caught.

| Item | What was wrong | Proof it was real (code at `6225dd4`) | Test that pins it | What it moves |
|---|---|---|---|---|
| R1. Unusable model display | A model whose comparison BED is negative was excluded from crediting and ranking, but its differences were still shown and entered `d_min` | λ_DI 1/day, λ_DD 2/day, anchor 20 × 2 Gy / 26 d: comparison DI BED −11 Gy, rows showed ΔBED DI +57.8 to +63.2 Gy | `test_v22_review_fixes.R::an unusable model shows no difference and does not enter the score` | Display and export in invalid settings only |
| R2. Both models negative | Reported as DI only | λ_DI 1, λ_DD 5, anchor 20 × 2 Gy: DI −11, DD −50.7, reason `negative_reference_di` | `::both models negative gives its own reason` | Message only |
| R3. Scenario ID | Hashed interface context (preset label, evaluated regimen), so identical searches got different IDs | Three IDs for one search (e54c3a920b74, 09a68b76c398, f06de2cba5b6) | `::scenario ID ignores interface context but not search settings` | File names of exports |
| R4. Empty export | `scenario.json` did not record why no regimen qualified | `attr(s, "reason")` not written | Rerun test plus `result$empty_reason` in `write_export_zip()` | Export content |
| R5. Withers rule and an empty field | Clearing α/β (NA) while retyping switched the preset to Custom and left λ_DI at the α/β 10 value | Reviewer's trace of `app.R` observer | Code review; browser walkthrough (the observer now returns early on any empty tumor field) | UI only |
| R6. Weak tests | A test could not fail (empty table), and six mutations survived | `any(logical(0))` is FALSE | `test_v22_review_fixes.R` (six tests), `test_v21_changes.R` negative-reference test now asserts the empty result | Tests only |
| F1. Fenwick 2008 boundary | Not available | Fix plan D11 | `test_v22_fenwick.R` (published-formula values 98.574, 74.169, 70.203 Gy₁₀; 20 × 2.9 Gy above, 55/20 below), `test_app_server.R::Fenwick acute limit reaches the search and the Results tab` | New option; default unchanged |
| C1. CI | No automated test run on push | Not applicable | `.github/workflows/tests.yml` (test only, never deploys) | Nothing in the app |
| D1. Text | HYHOPE acute BEDs (51.3, 51.4) are on the protocol's inclusive clock, not the app's; DD smoothing range unstated; archived `dd_floor.R` uses the inclusive acute clock | Review | `about.md`, `ERRATA.md` | Text only |

## 2026-10-19: version 2.3 (new defaults)

**Changes published numbers: yes, all default results.** Decision by Igor Shuryak on 2026-10-19: default tumor parameters become scenario D18, and the default early-toxicity limit becomes the acute BED of HYPNO 55 Gy/20 fractions. With the default regimen (20 × 2.65 Gy, 26 days) against 35 × 2 Gy: tumor BED DI 67.04 vs 71.51 (was 67.04 vs 65.00), DD 60.45 vs 71.51 (was 56.01 vs 65.00). Default search: 3 rows (was 7), led by 29 × 2.25 Gy in 39 days, ΔDI +1.19, ΔDD −2.41, early limit 55.87 Gy₁₀ (was 53.11). Nothing is deployed.

| Item | Record |
|---|---|
| What changed | `default_params()` tumor values α 0.35, α/β 10, λ_DI 0.23, λ_DD 0.40, Tk 28 (were 0.2, 10, 0.2, 0.5, 28). `DEFAULT_EARLY_MODE = "hypno"`, `DEFAULT_PRESET = "d18"`. The interface reads its starting values from `default_params()`. Version 1 values: `legacy_params()` and the "Calculator defaults (version 1)" preset. Tests that pin version 1 results now use `legacy_params()` explicitly. |
| Proof of the previous state | At `a82d8a6`, `default_params()$alpha` was 0.2 and the search tab opened on the comparison-relative acute limit. |
| Tests | `tests/testthat/test_v23_defaults.R` (defaults equal D18; the default search uses the 55.866 Gy₁₀ HYPNO limit and returns 3 rows led by 29 × 2.25 Gy; the rendered interface starts at α 0.35, λ_DI 0.23, λ_DD 0.40, preset D18, HYPNO limit checked). All three failed at `a82d8a6`. |
| Verification | Full suite, Windows R 4.5.1, 2026-10-19: 1474 passed, 0 failed, plus 2 package-version warnings. Headless browser: a fresh page shows preset D18, α 0.35, HYPNO limit selected, and the same 3 rows. |
| What it moves | Every default number on the Results and Candidate tabs. The model page explains the HYPNO limit and its limits, and that D18's DD values conflict with HYPNO's result. |
