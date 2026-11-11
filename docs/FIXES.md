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

## 2026-10-19: version 2.4 (two-model agreement tab)

**Changes published numbers: no.** A new tab only; existing tabs and defaults are unchanged. Nothing is deployed.

| Item | Record |
|---|---|
| What was missing | No view of regimens that both tumor models support against an evidence-based comparison (HYPNO 55/20), or of how robust that is over parameter scenarios. Requested by Igor Shuryak, 2026-10-19, after the analysis in this session. |
| What it does | `agreement_regimens()` in `R/bed_models.R` runs the search with both models required (ranked by the smaller change) at a margin, then adds each model's range over the 15 sensitivity scenarios, the count of scenarios in which both meet the margin, the DD onset day and whether the course ends before the DI onset. Tab "Model agreement" in `app.R`, defaults: HYPNO 55/20, margin −3 Gy₁₀, at least 1 day saved, toxicity limits from the candidate tab. |
| Tests | `tests/testthat/test_v24_agreement.R` (rows 17 × 3.05, 18 × 2.90, 19 × 2.80 Gy with ΔDI −2.460, −2.786, −2.028 and ΔDD −1.590, −1.961, −1.477 Gy₁₀ and DD onsets 14.28, 14.97, 15.42 days from an independent Python implementation; 15 of 15 scenarios; none at margin 0; counts need both models; long courses do not end before the DI onset), `test_app_server.R::model agreement tab lists the regimens both models support against HYPNO`. All failed before the function existed. Four mutations (either model instead of both, count ignoring DD, wrong DI-onset flag, tab ignoring the days input) are each caught. |
| Verification | Full suite, Windows R 4.5.1, 2026-10-19: 1502 passed, 0 failed, plus 2 package-version warnings. Headless browser: on a fresh page the tab shows the three regimens with "15 of 15". |
| What it moves | Nothing on existing tabs. |

## 2026-10-19: version 2.5 (review of 2.4, findings F01–F08 and N1–N6)

**Changes published numbers: no for default settings.** No equation, default or acute policy changed: the default candidate rows (29 × 2.25, 28 × 2.30, 30 × 2.20 Gy; ΔDI +1.1884, +1.1263, +1.1200; ΔDD −2.4080, −2.8648, −2.0495 Gy₁₀) and agreement rows (17 × 3.05, 18 × 2.90, 19 × 2.80 Gy, 15 of 15) are identical to 2.4. What moves: agreement and sensitivity counts and ranges where a scenario has a comparison or regimen BED that is not positive (they were counted; now they are reported as not evaluable); search differences against a non-positive regimen BED (now NA); export format (scenario schema 2; scenario IDs are now SHA-256 of canonical JSON, so 2.4 IDs are not comparable); labels and wording. Nothing is deployed.

Plan, verification of the external review and the decisions A1–A5 (taken without the PI, confirmed by Igor Shuryak on 2026-10-19): `docs/development/2026_10_19_review_and_plan.md`. External review: `docs/development/2026_10_19_external_review_codex.md`.

**Verification** (Windows R 4.5.1, 2026-10-19, revision `725514a`, which is `73b54a5` plus the display part of F01): `testthat::test_dir("tests/testthat")` gave **1651 expectations, 0 failed, 0 errors, 0 skipped** (baseline at `ad16504`: 1504, 0 failed; warnings are package-build-version notes only). The first 26 tests in `tests/testthat/test_v25_review_fixes.R` were run against `ad16504` in a separate worktree: 25 failed, each on the reproduced defect (for example "Enter a number for: Minimum tumor BED change" for F03); the one that passed is a guard on behavior that must not change (fixed anchor ignores the sidebar comparison). The 27th (Results and benchmark display) was added in `725514a`; mutation M14, which restores the old display rule, turns it red. Mutation check: all 14 mutations (`docs/development/2026_10_19_review_scripts/mutations.R`, output in `mutations_output.txt`) turn at least one test red on `725514a`: sensitivity validity bypassed, invalid scenarios counted as passing, wrong agreement comparator, candidate inputs leaking into the agreement, replay ignoring hinge or equation version, key-order-dependent ID, wrong unit label, sidebar anchor in the agreement, unchecked calibration, export requiring the Results regimen, negative candidate credited, Inf onset printed, and the old display rule. Combined smoke on the same revision: the full suite plus a headless-browser walkthrough (`docs/development/2026_10_19_review_scripts/walkthrough.R`, Chrome via chromote 0.5.1): 23 of 23 checks passed, log and screenshots in `docs/development/2026_10_19_walkthrough/`. Probe outputs before and after: `docs/development/2026_10_19_review_scripts/probe_v24.R` (on `ad16504`) and `probe_after_output.txt`.

| Item | What was wrong | Proof it was real (code at `ad16504`) | Test that pins it (red at `ad16504`) | What it moves |
|---|---|---|---|---|
| F01 Validity | Sensitivity and agreement counted scenarios with a negative comparison BED | λ_DI 0, comparison 35 × 2 Gy over 150 days: 32 × 2.10 Gy listed "15 of 15"; version 1 preset at α/β 10 had comparison DI −38 Gy and ΔDI +103.31 Gy, counted as a pass. Now 3 of 3 evaluable, 12 of 15 not evaluable | `test_v25_review_fixes.R::F01: agreement counts only scenarios in which both comparisons are evaluable`, `::F01: sensitivity table hides deltas against a negative comparison BED`, `::F01: endpoint rule covers negative candidate, zero, unavailable and all-invalid`, `::F01: Sensitivity tab and its summary count evaluable scenarios only`, `::F01: the Results and benchmark tables apply the validity rule too` | Counts, ranges and displayed differences in invalid settings only |
| F02, N1, N2 Claims | "Within the parameter uncertainty", "close to HYPNO in tumor effect", "pessimistic bound", "non-inferiority margin" for a BED threshold, "can only lose dose" before the DI onset, "cannot separate regimens 1–2 Gy₁₀", Hartley 15% as a general rule, "does not support escalation" | g 0 → 0.0154 changes ΔDI of 17 × 3.05 vs HYPNO from −2.4602 to −2.3282 Gy₁₀ although both end before day 28; at tumor α/β 2 the same regimen gains +0.2968 Gy at lower acute BED | `::F02: withdrawn claims are gone from the About page and the interface`, `::N1: the D18 HYPNO statement is conditional on a cited dose-response slope`, `::F02: before the DI onset, background growth still rewards a shorter course` | Text only |
| F03, N6 Tab independence | Agreement validated candidate-only inputs; export read the Results regimen; legacy coupling calibrated the agreement at the sidebar comparison | Clearing the candidate minimum gave "Enter a number for: Minimum tumor BED change" on the agreement tab; clearing the Results fractions blocked the candidate download | `::F03: agreement does not depend on candidate-only inputs`, `::F03: the candidate export does not depend on the Results regimen`, `::F03: an empty agreement explains its own reason`, `::N6: legacy coupling calibrates the agreement at its own comparison and says so`; guard `::F03: in fixed-anchor mode the sidebar comparison does not move the agreement` | Legacy-mode agreement values (calibrated at its own comparison); nothing in default mode |
| F04 Agreement export | Only the candidate search could be downloaded; file name said 2.1 | `app.R` had one `downloadHandler`; name `bed_v2.1_...` | `::F04: the agreement analysis downloads and replays in a fresh process` (callr), `::F04: download names carry the app version and the analysis` | New download |
| F05 Replay | Replay accepted any equation/schema version, hinge and tumor clock; ignored the saved calibration; ID changed with key order | Equation 9.9 + hard hinge + elapsed tumor clock replayed to the same 3 rows; reversing `params` changed the ID a98024bb355b → 0f3766cc8351 | `::F05: replay rejects unsupported model versions, hinge and clocks`, `::F05: schema 1 candidate files are migrated explicitly`, `::F05: the scenario ID ignores key order but not scientific inputs` | Export format; old candidate files still rerun (migrated) |
| F06, N3, N5 Units and labels | HYPNO acute cap labeled Gy₁₀ at mucosal α/β 5 (70.99 BED Gy); "Beyond precedent" overstated; Fenwick comment said elapsed days; DD onset printed "Inf"; Results rows had no BED unit | Probe and `app.R:528` at `ad16504` | `::F06: the HYPNO acute cap is labeled with the mucosal alpha/beta in use`, `::F06: Fenwick is not described as an elapsed-day model`, `::F06: the fraction-size flag names the benchmark set literally`, `::N3: an infinite DD onset is shown as 'none'`, `::N5: Results rows carry BED units`, `::F06: the sensitivity record states inherited values` | Labels |
| F07 Records | about.md said 33 trial arms (the 2018 paper used 35 regimens from 16 publications, 7283 patients) and HYHOPE protocol 4.1 (local copy is 4.2); catalogue said Fowler 2003 unread and g omitted; Fowler 2010 rate without its EQD2 unit; DESCRIPTION 2.1.0 | Shuryak 2018 p. 21; HYHOPE protocol signature page; Fowler 2010 text "0.7 Gy EQD 10/2 per day" | `::F07: versions agree and source statements are current` | Text only |
| F08 Release | No dependency lock; launcher tied to `C:\Program Files\R\R-4.5.1`; agreement table's last column cut off at the card edge (found by the walkthrough) | `run_app.bat`; screenshot `03_agreement.png` before the fix | `::Wide tables scroll horizontally instead of being cut off (walkthrough finding)`; launcher dry run through `cmd.exe` (finds R 4.5.1 under Program Files; honors `RSCRIPT`) | `renv.lock`, launcher, layout |

## 2026-10-19: version 2.5.1 (interface clarity, input-effect audit)

**Changes published numbers: no.** No computed value changes; `SELECTION_VERSION` stays 2.5. Record: `docs/development/2026_10_19_interface_clarity_and_input_audit.md`.

| Item | What was wrong | Proof it was real | Test that pins it (red at `953b284`, green after) | What it moves |
|---|---|---|---|---|
| Legacy coupling switch removed | An unexplained checkbox restored a calibration removed as a defect in 2.1; the user could not tell what it did | Reported by Igor Shuryak on 2026-10-19 while using the local 2.5 app | `test_app_server.R::comparison regimen does not move the DD model; only the anchor does`; `test_v25_review_fixes.R::the UI has no legacy coupling switch and the anchor note explains the anchor`; guard `::version 1 calibration is reproduced by entering the comparison regimen as the anchor` (passes on both, by design) | Nothing in default use; legacy-mode values are now obtained by entering the comparison as the anchor |
| Anchor note | "DD and DI tumor BEDs are equal at the anchor … starts on day 23.5" did not say why the anchor exists | Same report | `::the UI has no legacy coupling switch and the anchor note explains the anchor` | Text only |
| Preset menu | The default preset was not explained, and "Calculator defaults (version 1)" was not the default | Same report | `test_v25_review_fixes.R::the preset menu explains the default and every preset` | Text only |

Verification (Windows R 4.5.1, 2026-10-19): full suite **1670 expectations, 0 failed, 0 errors** (baseline 1651 at `953b284`); 18 of 18 mutations caught (`docs/development/2026_10_19_review_scripts/mutations_output.txt`); browser walkthrough 25 of 25 (`docs/development/2026_10_19_walkthrough/walkthrough_log.txt`); input-effect audit of all 38 inputs (`docs/development/2026_10_19_input_effects/input_effects_log.txt`). Live evidence: the walkthrough reads the new notes from the running app, and the input audit shows the preset note changing with every preset.

## 2026-10-19: version 2.6 (external review of 2.5.1, findings R1–R9)

**Changes published numbers: no for default settings.** No equation, default or acute policy changed. The default candidate rows (29 × 2.25, 28 × 2.30, 30 × 2.20 Gy; ΔDI +1.1884, +1.1263, +1.1200; ΔDD −2.4080, −2.8648, −2.0495 Gy₁₀) and agreement rows (17 × 3.05, 18 × 2.90, 19 × 2.80 Gy, 15 of 15) are identical to 2.5.1. What moves: the listed regimen when one model's BED is not positive (R1), agreement ranges when some scenarios are not evaluable (R5), and the replay of files from older versions, which now needs an explicit migration (R3). Selection version 2.5 → 2.6.

Record, verification of each finding and decisions B1–B9 taken without the PI: `docs/development/2026_10_19_v26_review_response.md`. External review: `docs/development/2026_10_19_external_review_v251_codex.md`.

| Item | What was wrong | Proof it was real (code at `451d4c5`) | Test that pins it (red at `451d4c5`) | What it moves |
|---|---|---|---|---|
| R1 Ranking | The conservative ranking scored a regimen by DI alone when its DD BED was not positive. A DD ranking could list a row with no score | α 0.2, λ_DI 0.4, λ_DD 5, 32 fractions, step 0.05 Gy, max 2.5 Gy: listed 32 × 2.50 Gy (DD BED −0.604, score +22) instead of 32 × 2.10 Gy (score +3.312). Step 0.5 with DD ranking: listed 32 × 2.50 Gy with score NA | `test_v26_review_fixes.R::R1: the conservative ranking never lists a regimen whose DD BED is not positive`, `::R1: a DD ranking with no evaluable DD returns an explicit reason, never an NA score`, `::R1: all 16 criterion and ranking pairs list only rankable rows, best per fraction number`, `::R1: the next-dose column names a ranking that becomes non-evaluable`, `::R1: the app explains a search with no rankable regimen` | Listed regimen in such cases; new `rank_eligible`, `fractions_saved`, `early_slack`, `late_slack` columns; `d_min` NA unless both models |
| R5 Agreement ranges | Counts used jointly evaluable scenarios, ranges each model's own | 150-day stress case: 3 evaluable scenarios, DD range −0.214 to +2.389 over 15 | `::R5: both ranges use the scenarios in which both models are evaluable` | That range: +0.988 to +2.389 |
| R2 Replay inputs | Fenwick flag "true", 1, NA or missing read as false; contradictory agreement fields and a wrong onset accepted; no grid bound | A Fenwick search replayed with `"fenwick": "true"` as a different acute policy without error | `::R2: the Fenwick flag must be one true or false value`, `::R2: enumerated fields, calibration and the grid size are validated`, `::R2: agreement files cannot carry search fields that contradict the agreement`, `::R2: source hashes are read from the bytes and flag code changed since loading` | Malformed files refused |
| R3 Replay vs migration | An older file was silently reanalyzed with new rules and called a replay | Version 2.4 file `tests/testthat/fixtures/schema1_v24_r1case.json`: 2.4 engine chose 2.10 Gy, 2.5.1 "replay" chose 2.50 Gy | `::R3: a version 2.4 file is refused by strict replay and reanalyzed only on request`, `::R3: version 2.5 files need migration; unknown and Fenwick-2.1 files are refused` | Older files need `migrate_scenario()` |
| R4 Illustration | About applied γ₅₀ to time-corrected EQD2 without saying so, used a linear approximation and called the DD loss "more than HYPNO allows" | Stuschke 1999 Appendix A (printed p. 119) defines γ₅₀ on total dose; exact logistic DD loss at γ₅₀ 3.2 is 25.0, not 27, points | `::R4: the illustration reproduces the review's independent values`, `::R4: the About page states the illustration's assumptions and no unconditional verdict` | Text only |
| R6 Comparator scope | Sidebar said every difference and limit was relative to it; anchor note said version 1 calibrated at the comparison | `Custom_BED_calculator_2/app.R` hard-codes 47 days and 84 Gy | `::R6: each analysis tab names its own comparison, anchor and numerical limits` | Text and context strips |
| R7 Wording and sources | Default paragraph contradicted its table; Withers note implied Withers gave α; four catalogue entries said "not downloaded" for held PDFs; HYPNO "≤ 1.4 points"; footer said 2.5 | `about.md:66`, `app.R:366`, `References/README.md` A2, A3, B1, C1 | `::R7: default narrative, preset wording, version footer and source records are current`, `::R1/R7: versions agree across the code and DESCRIPTION` | Text only |
| R8 Browser | No walkthrough of 2.5.1; at 420 px the unwrapped title made the page 532 px too wide; the old walkthrough's menu setter could not change select menus | `walkthrough_v26.R` first run: 5 narrow-width failures | Walkthrough check "narrow 420 px: … no page-level horizontal scroll" (5 tabs) | Layout only |
| R9 Environment | Lockfile lacked test packages and could not be restored into an empty Windows library without a compiler | `docs/development/2026_10_19_v26/restore_attempt_cran_old_lock.log`: 48 packages failed to build | Clean restore `clean_restore_empty_library.log`: 68 binaries, 0 source builds; `tools/run_tests.R` in that library | Package versions (B8): shiny 1.14.0, digest 0.6.39, testthat 3.3.2 among others |

**Verification** (Windows R 4.5.1, 2026-10-19, revision `cd5d835`): `Rscript tools/run_tests.R`, release library only, gave **2543 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings** (`docs/development/2026_10_19_v26/suite_release_library.log`; baseline at `451d4c5`: 1670, 0 failed, with package-build warnings). The 19 tests of `tests/testthat/test_v26_review_fixes.R` were run against `451d4c5` in a separate worktree: 18 failed on the reproduced defect (for example "Expected `s$d` to equal 2.1"); the one that passed checks the new illustration script against the review's independent values on unchanged equations. Two assertions (singular "day") were added after that run in `cd5d835`. Mutation check: 14 of 14 restored old lines turn the new tests red (`docs/development/2026_10_19_v26/mutations.py`, `mutations_output.txt`, run on the working tree that became `4253389`/`472542c`; the mutated lines are unchanged in `cd5d835`). Combined smoke on the same revision: browser walkthrough of the review's section 10 procedure, app in the release library, Chrome via chromote 0.5.1: **80 of 80** checks, including a 420-px pass and no browser-console errors (`docs/development/2026_10_19_v26/walkthrough/walkthrough_log.txt`). Four existing tests were changed because they pinned behavior this version changes on purpose (listed in the record, section 4).

## 2026-10-20: version 2.7 (a colleague's comments, E1–E4, and the independent review of the plan, F1–F9)

**Changes what a fresh page shows: yes. Changes any computed number: no.** The opening candidate table goes from 3 rows (29 × 2.25, 28 × 2.30, 30 × 2.20 Gy, each meeting +1 Gy under DI only, ΔDD −2.41, −2.86, −2.05 Gy₁₀) to 0 rows, because the default now requires both models. With the criterion and ranking chosen explicitly, every result is the same as in 2.6: compatible 2.6 scenario files replay to the same numbers (tested on four frozen files). No equation, parameter default, toxicity limit or selection rule changed. Selection version stays 2.6, and app version 2.6.0 → 2.7.0.

Record, verification of each finding, and decisions: `docs/development/2026_10_20_v27_colleague_review.md`.

| Item | What was wrong | Proof it was real (code at `0c51ae6`) | Test that pins it (red at `0c51ae6`) | What it moves |
|---|---|---|---|---|
| E2 Default screen | The opening search credited a regimen that only DI supported, and ranked by the more favorable model | Default search listed 29 × 2.25 Gy (ΔDI +1.19, ΔDD −2.41 Gy₁₀) and two similar rows. UI selected `either` and `best` | `test_v27_colleague_review.R::E2: the page opens with both models, the smaller change and +1 Gy` (red: "x[1]: "either""), `::E2: at the opening settings no candidate qualifies and the message explains it` | Opening table 3 → 0 rows. New empty-result message |
| E3 Counts | "15 of 15" and "evaluable scenarios" could read as a probability or as independent confirmations | Agreement cell "15 of 15", Sensitivity summary "evaluable scenarios" | `::E3: agreement counts are labeled as selected parameter settings`, `::E3: with invalid settings the denominator is the evaluable count, shown with the invalid count`, `::E3: the Sensitivity summary qualifies what the settings are` | Display text only. Export fields unchanged |
| E4 Caveats | Limitations were only on the Model and sources page | No caveat element on the Results or Model agreement tabs | `::E4: one caveat line on each of the Results, Candidate and Model agreement tabs`, `::A: app.R builds its UI` | New line on three tabs |
| E1 Source | The app had no link to its code, and the repository opens on version 1 | `grep -i github app.R about.md`: no hit. Default branch `main` | `::E1: the Model and sources page links the version 2 source branch` | About page. README banner on `main` (separate commit) |
| F4–F6 Text | "DI or DD lists only one-model regimens", "assumed equal chemotherapy", "gaps not modeled", "no regimen is supported by both", and README "versions 1 and 2.0 calibrated at the comparison" were wrong or overstated | `about.md:34, 59, 63, 67, 81, 158, 164`, `README.md:64, 95` | Wording in `::E2: at the opening settings …` and `::E3: the Sensitivity summary …` | Text only |
| Compatibility | Must not move | — | `::D: compatible version 2.6 files replay to the same results` (117 assertions, passes on 2.6 and 2.7), `::E2: the either/larger screen is still available and unchanged`, `::E2: an allowed loss of 3 Gy lists three regimens in conservative order` | Nothing |

**Verification** (Windows R 4.5.1, release library, 2026-10-20):
- **Baseline:** at `0c51ae6` the suite gave 2543 expectations, 0 failed.
- **New tests on the 2.6 code:** they failed on their own assertions (`docs/development/2026_10_20_v27/red_before_change.log`). The compatibility tests passed there, as intended.
- **Mutations:** 6 of 6 restored old behaviors turn the new tests red (`mutation_checks.log`).
- **Full suite after the change:** 2727 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings.
- **Walkthrough:** `walkthrough_v27.R` (local) passed 52 of 52 checks and exits 1 when a check fails (`--self-test-fail`).
- **Changed existing tests:** two assertions were updated because they pinned display wording changed on purpose, `test_app_server.R:290` and `test_v25_review_fixes.R:83`.
- **Revision:** the tested revision is the commit that adds this row. Its SHA, and the rerun of the suite and walkthrough on it, are recorded in the next row.

**Verification on the committed revision** `96ec8ac` (2026-10-20, clean tree, release library):
- **Suite:** `tools/run_tests.R` gave 2727 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings (`docs/development/2026_10_20_v27/suite_96ec8ac.log`).
- **Walkthrough:** `walkthrough_v27.R` in local mode passed 52 of 52 checks with exit 0, with no uncommitted app changes (`docs/development/2026_10_20_v27/walkthrough/20261020_113231_local/`).
- **Not deployed:** the live site is still 2.6.0. Deployment and the live walkthrough wait for local inspection.

**Deployed** 2026-10-20 from `436375a` (code identical to `96ec8ac`) with `Rscript tools/deploy.R`:
- **Account:** the only configured shinyapps.io account, `ishuryak`. The bundle was `about.md`, `app.R`, `DESCRIPTION`, `LICENSE`, `R/bed_models.R`, `R/export.R` and `renv.lock`.
- **First live walkthrough:** 50 of 51 checks. The −3 Gy row check read the table before the live server had re-rendered it. The context strip already showed −3 Gy, and the download replay of the same rows passed (`walkthrough/20261020_115952_live_first_run_timing_fail/`).
- **Script fix:** the row checks now poll for up to 20 s. Rerun live: **51 of 51, exit 0** (`walkthrough/20261020_120129_live/`). The local rerun with the same script was 52 of 52 (`walkthrough/20261020_120255_local/`).

## 2026-10-21: deployment safety (PLAN_v3 section 3), no app change

**Changes published numbers:** no. `app.R`, `R/` and `about.md` are unchanged, so the live 2.7.0 app is unaffected.

| Item | Record |
|---|---|
| What was wrong | (1) The tracked `.Rhistory` held a `rsconnect::setAccountInfo()` call with a full shinyapps.io token and secret, in this public repository since `2b30456` (2024-10-25). (2) `tools/deploy.R` deployed to whichever shinyapps.io account was listed first and did not check for uncommitted changes. (3) `DESCRIPTION` said MIT while `LICENSE` is Apache 2.0. (4) The inactive version 1 deploy workflow was still on this branch. |
| Proof it was real | (1) `git log --oneline -- .Rhistory` lists `2b30456`; the file's line 1 is the `setAccountInfo` call (values not reproduced here). (2) `deploy.R` line 30 at `65e826f`: `account = acc$name[1]`. (3) `DESCRIPTION` line 10 vs `LICENSE` line 2. |
| Test that pins it | `tests/testthat/test_deploy_preflight.R` (10 expectations). Against the previous `deploy.R`, `.rscignore` and `.gitignore`, 5 of them fail on the deploy.R and ignore-file assertions (checked 2026-10-21). |
| Verification | `.Rhistory` untracked (`git rm --cached`), listed in `.gitignore` and `.rscignore`; `.Renviron` added to both. `deploy_preflight()` refuses another app name, the version 1 name, an unconfigured pinned account and uncommitted tracked changes. `DESCRIPTION` now names Apache 2.0, the licence text actually in `LICENSE`. |
| What it moves | Nothing in any result. **Still required, by the account owner:** revoke the exposed token at shinyapps.io (Account > Tokens) and issue a new one, typed in the R console only. Untracking does not remove the old file from the history of this branch or of `main`; a history rewrite is a separate decision. |

## 2026-10-21: version 3.0.0 (branch `v3`, not deployed)

**Changes published numbers:** no live number changes (version 3 is not deployed; versions 1 and 2 are untouched). Within the app, for the same inputs every BED, search, sensitivity and both-model result equals version 2.7. Two defaults moved: the evaluated regimen (20 × 2.65 Gy to HYPNO's 20 × 2.75 Gy) and the both-model screen's comparison and threshold (HYPNO and −3 Gy to the global 35 × 2 Gy and +1 Gy, so it opens empty instead of listing 17, 18 and 19 fractions).

| Item | Record |
|---|---|
| What was wrong | (1) Parameter uncertainty was not quantified, organ-at-risk doses were not modeled, each tab had its own comparison and threshold, and only Mon-Fri courses from a Monday were covered (PLAN_v3 workstreams A-D). (2) Found during development and fixed before commit: `load_draws()` dropped its SHA-256 and ESS attributes. |
| Proof it was real | (2) `attr(load_draws("di"), "sha256")` was `NULL`, and `rerun_uncertainty()` refused every file with "the draws files differ from those recorded". |
| Test that pins it | `tests/testthat/test_v3_app.R`, `test_v3_uncertainty.R`, `test_v3_oar_schedule.R`, `test_deploy_preflight.R`. (2) `test_v3_uncertainty.R::draws files match the manifest, and a changed file is refused` fails when the attribute bug is restored (mutation M1). Twelve deliberate breaks of the new guards are all caught (`docs/development/2026_10_21_v3/mutation_checks.log`). |
| Verification | See the verification line below this table. |
| What it moves | No published number. The defaults listed above. Reasoning and departures from the plan: `docs/development/2026_10_21_v3_build.md`; refit: `docs/development/2026_10_21_v3/fit_report.md`. |

**Verification (version 3.0.0, revision `8fa0b57`, 2026-10-21):** `Rscript tools/run_tests.R` gave 4956 expectations, 0 failed, 0 errors (version 2.7 baseline on `65e826f`: 2727, 0 failed). Version 2 tests were changed only where version 3 changes behavior on purpose (listed in the commit messages of `41f6659` and `8fa0b57`). Deliberate breaks: 12 of 12 caught. Local browser walkthrough on the same revision: 34 of 34 (`docs/development/2026_10_21_v3/walkthrough/20261021_155556_local/`), and the self-test run's deliberate failure gave exit status 1. Not deployed.

## 2026-10-21: first version 3 deployment bundle lacked the draws manifest

**Changes published numbers:** no. The first upload (from `9b29d15`) was replaced within minutes, before any walkthrough or announcement.

| Item | Record |
|---|---|
| What was wrong | rsconnect skips files named `manifest.json`, so `draws/manifest.json` was not uploaded. Every draws load on the server would fail (Uncertainty tab, report card refit column). The bundle also carried local-only files: `external_proofs.yaml` (local paths of licensed PDFs, no PDF content), the audit registries and `.stale_constant_record.json`. |
| Proof it was real | The deploy log's file list for `9b29d15` names `draws/refit_point_values.json` and the four draws CSVs but no manifest, and includes `external_proofs.yaml`. |
| Test that pins it | `tests/testthat/test_deploy_preflight.R::the bundle check refuses a missing draws manifest and local-only files`. `tools/deploy.R` now runs `check_bundle(rsconnect::listDeploymentFiles("."))` before uploading. |
| Verification | The manifest is renamed `draws/draws_manifest.json`; `listDeploymentFiles(".")` lists it and none of the local-only files (checked 2026-10-21). |
| What it moves | Nothing in any result. |

**Deployed** 2026-10-21 from `f951cd3` (suite 4961 expectations, 0 failed, on that revision) as `custom_bed_calculator_v3` with `Rscript tools/deploy.R`, account `ishuryak`. Bundle: `about.md`, `app.R`, `DESCRIPTION`, `LICENSE`, `renv.lock`, `deployed_revision.txt`, `R/` (7 files), `data/` (3 files), `draws/` (6 files). Live walkthrough: 34 of 34 (`docs/development/2026_10_21_v3/walkthrough/20261021_162106_live/`). Versions 1 and 2 not redeployed.

## 2026-10-26: version 3.1.0, the verified-trial refit in the app (review of 2026-10-25, stages A–E)

**Changes published numbers:** not yet. Version 3.1.0 is not deployed, so the live 3.0.0 site is unchanged. Once deployed, the Uncertainty tab, the report card's refit column and the refit preset move. For example, HYPNO 55 Gy/20 against 35 × 2 Gy over the first 50 shipped DI draws gives a mean ΔDI of −2.28 Gy, against +1.75 Gy₁₀ before. Mean-dose organ statuses move (physical comparison, stage A). D18 deterministic results do not move.

| Item | Record |
|---|---|
| What was wrong | (1) The app's draws came from a fit to the 35 arms of the 2018 Table 1 rather than the trial publications. The DD draws were a separate, unidentified fit not paired with the DI draws, and nothing refused draws from a fit with sampling problems. (2) Two data errors in the 10-25 trial file: the KROG 0201 HYPO arm had 8 T2 patients (Moon 2014 Table 1: 9), and RTOG 9003 combined rates from one report with patient numbers from another. (3) Mixed fraction sizes entered as one mean fraction size. (4) BED differences at varying α/β were labeled Gy₁₀. |
| Proof it was real | (1) The 10-21 manifest has no pairing and no gate field, and its free DD fit has 74 divergences. (2) Moon 2014 Table 1 lists T2a 7 + T2b 2. Fu 2000 Tables 1–2 give N 268/263/274/268 with rates 46.0/54.4/47.5/54.5%. (3) Ghoshal 25 × 1.8 + 15 × 1.5 Gy: Σd² 114.75 vs 40 × 1.6875² = 113.9. (4) `paper_restrictions` is the only variant with α/β 10 in every draw. |
| Test that pins it | `tests/python/test_refit_data_and_identity.py` (data, components, identity), `tests/python/test_export_draws.py` (pairing, gate refusal), `tests/python/test_repair_dd_calibrated.py`, `tests/testthat/test_v31_stage_e.R` (manifest, gate refusal, row pairing with anchor equality < 1e-6, units, point values). Each was observed red without its change (on the 10-25 file, the HEAD runner, the 10-21 draws, or by mutating the shipped line) and green with it. Details: `docs/development/2026_10_25_stages_BCDE_refit.md`. |
| Verification | 2026-10-26, base revision 677ce21 plus the working tree: `python3 -m pytest tests/python -q` 38 passed. `Rscript tools/run_tests.R` 5265 expectations, 0 failed, 0 errors (baseline 3.0.0 on `f951cd3`: 4961, 0 failed). Fits: 12 DI fits and 3600 calibrated-DD fits pass the gate, 35 of them after the repair step. |
| What it moves | The uncertainty results and refit preset listed in the first line, organ mean-dose statuses, labels. Not the D18 defaults, the anchor or any deterministic BED. |

### 2026-10-26: deploy bundle limited to committed files (review S6)

| | |
|---|---|
| What was wrong | The deploy check compared tracked files only (`--untracked-files=no`), so an untracked file in the app folder would be uploaded with the bundle while `deployed_revision.txt` named a commit that did not contain it. |
| Proof it was real | Stage A baseline: `.suite_passed` was written and the preflight passed with an untracked file present (`docs/development/2026_10_25_stage_A_app_fixes.md`, S6). `check_bundle()` had no knowledge of git. |
| Test that pins it | `tests/testthat/test_deploy_preflight.R::the bundle check refuses files that are not tracked by git`. Red without the change (`check_bundle` has no `tracked` argument), and red when the shipped line is mutated to check nothing (the error is not thrown); green with it. |
| Verification | 2026-10-26: `Rscript tools/run_tests.R` on `e4ff22f`: 5269 expectations, 0 failed (baseline before the change on `ba1e7a6`: 5267, 0 failed). The deploy of `e4ff22f` ran the check on the real bundle (40 files); its first dry run caught `.pytest_cache/`. `tools/deploy.R` passes `git ls-files` to `check_bundle()`. |
| What it moves | No number. Deploys now stop if an uncommitted file would be uploaded. |

### 2026-10-26: version 3.2.0, response to the review of version 3.1

| | |
|---|---|
| What was wrong | The calibrated-DD fits were accepted without the ESS part of the approved gate (R01), so 748 of 3600 conditional fits behind the shipped and sensitivity draws had bulk ESS below 400. The export did not verify the data or parent-draw bytes it relabeled (R02). Malformed organ mappings could replay with no organ limit (R03). The Results tab showed "not available" for benchmark rows 2-3 (R04). An indeterminate DVH maximum read "within at every α/β" (R05). Pooled DD summaries gave retried fits up to 8 times the weight (R06). Several statements were stale or too categorical (R07). Known source intervals were dropped (R08). The suite marker trusted an unchecked `git` call (R09). Malformed Eclipse rows were silently dropped (R10). Also: default evaluated regimen changed to 18 x 2.9 Gy at Igor's request. |
| Proof it was real | Stored per-fit diagnostics: primary 177/1000, restricted 272/1000, sensitivity sets 26-46 of 200 below bulk ESS 400, minimum 151, all with `gate_pass` true. Code reading confirmed R02, R04, R05, R10 (`app.R` `y[i]` on a scalar, `status_exceeds()` mapping indeterminate to FALSE, regex row filter). The review's probes reproduced R03 and R02 (a changed parent lambda exported with a 6.19 Gy anchor error). |
| Test that pins it | `tests/python/test_gate_and_ess_repair.py`, `tests/python/test_export_draws.py` (identity mutations, ESS-150 fit, negative g, missing variant key, no partial output), `tests/python/test_loso_cache.py`, `tests/python/test_provenance.py`, `tests/python/test_repair_dd_calibrated.py` (retired tool refuses), `tests/testthat/test_v32_review.R` (default, R01 app side, R03, R04, R05, R07, R09, R10, grid bound). Each was observed red against the unchanged code or by mutating the shipped line, then green. Details in `docs/development/2026_10_26_review_v31_fixes.md`. |
| Verification | 2026-10-26: the ESS repair (`fit/run_ess_repair_2026_10_26.sh`, revision f1a81ed) refitted 748 fits, all passing, 0 still failing. Export passed every identity check, anchor residual at most 4.5e-8 Gy. `Rscript tools/run_tests.R` 5379 expectations, 0 failed, 0 errors (baseline 3.1.0 on `e4ff22f`: 5269, 0 failed). `python3 -m pytest tests/python` 70 passed (baseline 38). |
| What it moves | The default evaluated regimen (18 x 2.9 Gy), so every tab's opening values. Display (audit decisions): report card and organ benchmark rows show the α/β-range result and the mean-dose note, the banner flags extrapolated fraction sizes, both Fowler zones are shown, the organ limit column is relabeled. The DD draws (DI draws byte-identical): HYPNO vs 35 x 2 Gy, DD −3.27 (−5.92 to −0.53) Gy in the conditional fits, was −3.3, within Monte Carlo error. DI unchanged (−2.0, full posterior). Sensitivity medians unchanged in substance (DI −2.7 to +0.3 from the full DI posteriors, DD −3.7 to −2.1). Benchmark rows 2-3 on the Results tab. Organ α/β-range text for indeterminate maxima. Uncertainty scenario files of 3.1 no longer replay. |

**Correction (2026-10-26, review of version 3.2).** Two numbers in the 3.2.0 row above are wrong as stated. 752 fits were refitted, not 748: 748 was the count below the bulk-ESS criterion, and 4 more were below the tail criterion only (from `repair_ess.json`). The "anchor residual at most 4.5e-8 Gy" was in log cells, not Gy. In Gy of BED it is 9.8e-7 (primary) and 7.8e-7 (restricted), within the 1e-4 Gy tolerance. The row is kept as written.

### 2026-10-26: version 3.2.1, response to the review of version 3.2

**Changes published numbers:** no. No BED equation, selection rule, draw file or displayed tumor, acute or late value changes. Two manifest residuals change unit (see below). Two displayed texts change: the toxicity share near the required share (one decimal, "meets"/"below") and the Monte Carlo error statement (95% limits "up to about 0.3 Gy", was 0.25).

| | |
|---|---|
| What was wrong | F01: a DVH row with a missing field was recycled and a repeated structure name replaced the first block, so an exceeded organ limit could read "within". F02: the diagnostic gate skipped NaN quantities before the minimum, excluded constant sampled quantities, accepted negative counts (Python) and ±Inf (R). F03: the anchor residual was log cells labelled Gy. F04: an ESS repair could run with a different prior and keep the original identity. F05: a failed second rename left no `draws/`. F06: the guide said the Both-model screen needs 15 of 15 settings. F07: repair count 748 (is 752), selection hash depended on the NumPy integer width, and several wording items. M01: `deploy.R` used unchecked `system("git ...")`. M02: unused DVH fraction input. M03: 88.7% and 90.3% both shown as "90%". M04: `References/README.md` presented 2.6 defaults as current. |
| Proof it was real | F01: on f436e57 the review's parotid file parsed to volumes 100/10/0 (mean 7.5 Gy, "within"), its two "Cord" blocks to the second only. F02: one quantity with NaN bulk ESS beside others at 700 passed; a manifest with ESS `1e999` loaded in R. F03: a pair with log-cell residual 5e-5 (0.0011 Gy) passed the 1e-4 "Gy" gate. F04: prior sd 10 accepted by a folder fitted with sd 3. F05: the f436e57 logic, with the second rename failing, left `draws/` absent. F07: the same 1000 indices hash differently as int32 and int64. |
| Test that pins it | `tests/testthat/test_v321_dvh_parser.R`, `tests/testthat/test_v321_gate.R` and `tests/python/test_v321_gate.py` (shared table `tests/fixtures/gate_cases.csv`), `tests/python/test_v321_anchor_and_publish.py` (F03, F05, F07 hash), `tests/python/test_gate_and_ess_repair.py::test_ess_repair_refuses_a_changed_posterior_target` and `::test_ess_repair_with_prior_sd_10_matches_only_a_prior_10_folder`, `tests/python/test_export_draws.py::test_v321_a_held_export_lock_refuses_and_leaves_draws_untouched`, `tests/testthat/test_v321_review_text.R` (F06 12/15 case, F07 text), `tests/testthat/test_v321_deploy_identity.R`, `tests/testthat/test_v321_ui_items.R`, `tests/testthat/test_v321_golden.R` (review sections 4.3-4.4). Each was observed red on f436e57 (F05 with the old logic placed in `publish()`), then green. Details in `docs/development/2026_10_26_v321_review_fixes.md`. |
| Verification | 2026-10-26: `Rscript tools/run_tests.R` on 2350c2b (all code changes): 5489 expectations, 0 failed, 0 errors, 0 warnings (baseline on f436e57: 5394, 0 failed). `python3 -m pytest tests/python`: 114 passed (baseline 70). Export run for real: draw CSVs byte-identical. Repair target check passed on all ten 2026-10-25c DD folders. App launched locally on the branch and served every tab, with no DVH fraction input. `stale_constant_audit.py --strict`: 0 errors (19 constants). The final suite run on the release commit writes `.suite_passed`. Live browser walkthrough not yet run on this build. |
| What it moves | No result. `draws/draws_manifest.json`: gate policy `gate-2026-10-27`, a `per_quantity_check` note, anchor residuals 3.89e-8 / 4.48e-8 (log cells) become 9.840e-7 / 7.752e-7 Gy. Draw CSVs and `refit_point_values.json` byte-identical. DVH files with ragged rows or repeated names are now refused. The DVH tab no longer asks for the plan's number of fractions. |

### 2026-10-28: version 4.0.0, simplified layout (branch `v4`, from `c782ab0`)

**Changes published numbers:** no. Same engine, parameters, draws and selection rules as 3.2.1, so every version 3 page gives the same numbers for the same inputs. The new main page is a new display of `suggest_regimens()` output. Record: `docs/development/2026_10_28_v4_simplified_layout.md`.

| | |
|---|---|
| What was wrong | Version 3 opened on a report card among nine tabs, with a seven-panel parameter sidebar. A clinician could not get a list of usable hypofractionated schedules without first learning the comparison, threshold, toxicity-limit and ranking controls. Separately, PR #1 moved two registered `app.R` sites, so the strict stale-constant audit failed on `c782ab0`. |
| Proof it was real | Request from Igor Shuryak, 2026-10-28. `stale_constant_audit.py --project . --project-root . --strict` on `c782ab0`: 2 errors (`verdict_tolerance_gy` at app.R:89, `max_days_input` at app.R:47). |
| Test that pins it | `tests/testthat/test_v4_home.R` (81 expectations): list equals a brute-force evaluation of the three conditions at the defaults and under three other presets, pinned default rows (16–34 fractions), summary wording, server main page and row detail, UI opens on the main page with every menu page. Red without the change (`R/v4_home.R` absent). Mutations caught: tumor-margin sign flipped (23 failures), criterion "both" (18), late margin 0 (10). Browser: `docs/development/2026_10_28_v4/walkthrough_v4.R`, 36 checks. |
| Verification | 2026-10-28, Windows R 4.5.1, release library: `Rscript tools/run_tests.R` on the uncommitted v4 tree gave 5597 expectations, 0 failed, 0 errors, 0 warnings (baseline on `c782ab0`: 5514, 0 failed). Eight version 3 expectations were updated for the new wording and version number, and one for the v4 deploy name. `walkthrough_v4.R` (local): 36 of 36 checks. `stale_constant_audit.py --strict`: 0 errors (2 on `c782ab0`). A release-certified suite run needs a commit. |
| What it moves | No BED value. Main page at the defaults: 19 schedules, 16–34 fractions (20 fractions at 2.65–2.80 Gy, containing HYPNO 55 Gy/20). Deploy target is now `custom_bed_calculator_v4`, and the version 3 site name is refused. |

### 2026-10-28: version 4.1.0, review of version 4.0.0 (branch `v4`, from `61ccaed`)

**Changes published numbers:** no BED value changes. The main page shows more: the tumor BED change under each model, refit support and notes, and its wording no longer implies equivalence. Every v3 page keeps its numbers, defaults and text, apart from the link color, a non-sticky banner and menu descriptions. Review, decisions and the applied list: `docs/development/2026_10_28_v4_review_and_plan.md`.

| | |
|---|---|
| What was wrong | Four reviews of 4.0.0 found that the main page read as an equivalence claim ("about as effective"). It also hid that 16–21 fractions pass on the DI model alone while DD estimates losses of 7–11 Gy₁₀. It showed no extrapolation warning (rows 16–20 lie wholly above the 2.55 Gy data limit) and no parameter uncertainty. Its early limit (58.1 Gy₁₀) is above HYPNO's acute BED (55.9) without saying so. The Report card contradicted the list for the same regimen (+1 Gy threshold, HYPNO acute limit). Two tests checked less than recorded, and a missing DD model was not announced. |
| Proof it was real | Re-run on `61ccaed`, defaults: 20 × 2.65 Gy ΔDI −4.47, ΔDD −11.06 Gy₁₀. 20 × 2.75 Gy (HYPNO) −1.39 and −8.67. Early cap 58.11 against HYPNO 55.87 Gy₁₀. Refit support 21 × 2.55 Gy 0.37, 20 × 2.75 Gy 0.95. Report card for 18 × 2.9 Gy: "threshold +1: not met" while the list contains it. The preset test compared `nrow()` only. |
| Test that pins it | `tests/testthat/test_v4_home.R`: "pinned result at the defaults" (DI / DD ranges, notes), "refit support" (0.95 / 0.37 / 1.00, checked against an independent count), "summary wording", "the main page lists the regimens…" (settings disclosure, 0.05 Gy wording, DD notice, table columns, Report card agreement through `v4_report_inputs()`: "−5: met" for 20 × 2.75 Gy, "exceeds" for 20 × 2.85 Gy), and the two corrected tests. Against the 4.0.0 code: 10 failed, 4 errors. Mutations caught: acute flag at "≥" (1 failure), support requiring both models (2). |
| Verification | 2026-10-28, Windows R 4.5.1, release library: `Rscript tools/run_tests.R` 5625 expectations, 0 failed, 0 errors, 0 warnings (baseline `61ccaed`: 5597, 0 failed). `walkthrough_v4.R` (local): 46 of 46, including a real click on "Open on the Report card", the card scrolled into view, pages opening at the top, frozen first column and menu closing at 390 px. `stale_constant_audit.py --strict`: 0 errors (new constant `v4_support_min`). |
| What it moves | No BED value. Main page at the defaults: same 19 rows (16–34 fractions), now with DI / DD changes, refit support (for example 51–99% at 20 fractions) and notes. Version 4.1.0. |

Also in 4.1.0: `tools/deploy.R` titled the shinyapps.io app "version 3". It now reads "version 4". Pinned by `tests/testthat/test_deploy_preflight.R::"deploy.R titles the shinyapps.io app as version 4"`, observed red before the edit and green after. Moves no number.

### 2026-10-30: version 4.2.0, main-page layout requests (branch `v4`, from `8355d3e`)

**Changes published numbers:** no BED value changes. The main list at the defaults drops from 19 rows (16–34 fractions) to 15 rows (30 down to 16); every remaining row keeps its doses and BED changes. Record: `docs/development/2026_10_30_v4_2_main_page.md`.

| | |
|---|---|
| What was wrong | Requests from Igor Shuryak, 2026-10-30: the list included 31–34 fraction schedules that save little time; it was sorted fewest fractions first; seven columns made the font small, three of them (tumor model, refit support, notes) secondary; early and late toxicity changes were shown only after a row click; no overview figure. |
| Proof it was real | Live 4.1.0 at the defaults: 19 rows, first row 16 fractions, no early or late column in the main table, no plot. |
| Test that pins it | `tests/testthat/test_v4_home.R`: "the standard is 35 x 2 Gy…" (search ceiling 30), "pinned result at the defaults" (rows 30:16, early −0.3 to +4.3 and late −16.8 to −8.4 Gy at 20 fractions, no "−0.0"), "summary wording", new "plot data…", and "the main page lists the regimens…" (early/late headers, details hidden until the checkbox is ticked, plot rendered, row 11 = 20 fractions). Against the 4.1.0 code: 17 failed, 2 errors. Browser: `walkthrough_v4.R`, 49 checks. |
| Verification | 2026-10-30, Windows R 4.5.1, release library: `Rscript tools/run_tests.R` 5649 expectations, 0 failed, 0 errors, 0 warnings (baseline `8355d3e`: 5626, 0 failed). `walkthrough_v4.R` (local): 49 of 49. `stale_constant_audit.py --strict --project-root .`: 0 errors (new constant `v4_min_fewer_fractions`). |
| What it moves | No BED value. Main page: 15 rows at the defaults (30 → 16 fractions), new early and late BED change columns, a Details checkbox, and a BED-change plot. Version 4.2.0. |

### 2026-10-30: version 4.2.1, range explanation under Details (branch `v4`, from `ad73a34`)

**Changes published numbers:** no. Display text only.

| | |
|---|---|
| What was wrong | Igor Shuryak, 2026-10-30: the main page did not say how a row's dose range and BED change ranges are derived (asked "this regimen has BED DI between −1 and +1, how was that determined?"). |
| Proof it was real | 4.2.0 footnotes and Details had no description of the 0.05 Gy dose search or of what the two ends of a range mean. |
| Test that pins it | `tests/testthat/test_v4_home.R::range explanation: the doses next to a row fail, for the stated reasons` (hand check of the DI formula at 30 × 2.10 and 2.25 Gy, the reasons 2.05 and 2.30 Gy fail, and that the next dose outside every row fails a condition) and `::the main page lists the regimens…` (explanation shown only with Details). Against `ad73a34`: 1 error, 2 failed. Browser: `walkthrough_v4.R`, 50 checks. |
| Verification | See the 4.2.1 commit message for the suite count on that revision. |
| What it moves | Nothing numeric. Version 4.2.1. |

## 2026-10-31: version 4.3.0, external review of 4.2.1 (branch `review_v4_2026_10_31`)

**Changes published numbers (live 4.2.1):** no. No equation, parameter, draws file or selection rule changed. At the defaults the main list is identical (60 points, 15 rows, same BEDs and refit support, checked in `test_v43_review.R::F01: the default list and its worked example are unchanged`). Text on the main page changes (units, model description, HYPNO counts, transfer scope). Plan, decisions and what was not done: `docs/development/2026_10_31_v4_review_fixes.md`.

| Item | What was wrong | Proof it was real (4.2.1, 917825b) | Test that pins it (red on 917825b, green on 4.3.0) |
|---|---|---|---|
| F01 | Range explanation said every BED rises with dose and that left/right values belong to the lowest/highest dose; its failure helper credited a model whose standard BED was negative | λDD = 2: 30-fraction DD maximum at 2.10 Gy, minimum at 2.25 Gy. Negative-reference fixture printed "2.10 Gy fails ()" | `tests/testthat/test_v43_review.R::F01: a non-monotone DD row names the doses of its extremes`, `::F01: an unusable reference never credits a schedule in the explanation`, `::F01: unavailable DD and a negative candidate are named, not left blank`, `::F01: a dose range with a gap is shown as separate runs` |
| F02 | Lead and criteria fixed at Gy₁₀/Gy₃ and "about 4 Gy"; footnote said the models agree at 35 × 2 Gy for any anchor and mode; a finite negative standard DD BED was not announced | At α/β 2.5 the table showed "Gy (α/β = 2.5 Gy)" while the lead said Gy₁₀. Custom anchor 20 × 2.75: DI 71.51, DD 81.49 at the standard. Direct DD α 0.005: standard DD −284.35, no notice | `::F02: the DD model description names the actual anchor and mode`, `::F02: a finite but negative standard DD BED is unusable`, `::F02/F03: the main page text follows the settings and quotes HYPNO correctly` |
| F03 | HYPNO note: "50.9% of 395 patients" | HYPNO: 190/373 evaluable (395 randomized) | `::F02/F03: ...` (asserts "190 of the 373 evaluable", no "of 395 patients") |
| F04 | "Open on the Report card" claimed every menu page uses the main-page rules | Search keeps criterion both/worst and excludes the DI-only 20 × 2.70 point | `::F04/F05: the transfer states its scope and refuses a dose outside the row` (scope text; no agreement/OAR/Search inputs sent) |
| F05 | Transfer sent any dose the browser supplied | Row 11 with dose 8: `updateNumericInput(d, 8)` sent | `::F05: only a listed point of the current list can be transferred`, `::F04/F05: ...` (dose 8 sends nothing) |
| F06 | One bad chain's NaN E-BFMI could be hidden by `min()`; missing c[]/z[]/delta_eval indices passed; `False` counted as 0 | Constant energy in chain 2 passed the gate, in chain 1 failed | `tests/python/test_v43_review_guards.py::test_one_bad_energy_chain_fails_in_every_position`, `::test_bad_chain_classification_does_not_depend_on_chain_order`, `::test_repair_fit_loop_rejects_a_fit_missing_its_nuisance_vectors`, `::test_bool_and_nonnumeric_diagnostics_are_rejected`, `::test_malformed_sampler_arrays_fail` |
| F07 | ESS repair filled a missing historical anchor from current `C.ANCHOR` | Current anchor changed to 30 × 2/40 was accepted against the 35 × 2/47 identity and reached the compile step | `::test_changed_current_anchor_is_refused_against_a_legacy_identity`, `::test_changed_current_evaluation_regimen_is_refused...`, `::test_fit_draw_calibrates_at_the_checked_anchor...`; all ten real targets still pass (`::test_real_intended_targets_pass_the_target_check`) |
| F08 | No export of the main list | Only the Candidate ZIP existed (Search settings) | `test_v43_review.R::F08: the default main list is exported in full and replays in a fresh R process` and three further F08 tests |
| F09 | Guide was 3.2.1; About said 5 to 34 fractions | About line 14 | Text change: `docs/user_guide/HNSCC_BED_v4_user_guide.md` and `.docx` |
| F10 | CI did not run on pushes to `v4` | `.github/workflows/tests.yml` branches | `test_v43_review.R::F10: CI runs on pushes to the v4 branch` |
| F11 | Where retained fit artifacts live was undocumented | v4 lacks `fit/output_refit_2026_10_25c_main/...` | Documentation: `fit/README.md` |
| F12 | A first failed draws load lasted until restart | Restoring a valid loader still returned NULL | `test_v43_review.R::F12: a failed draw load is retried after a bounded wait; a success is cached` |

Several R tests were red on 917825b because the new helper did not exist yet. The behavioral reds are the server tests (F02/F03, F04/F05) and, for F06/F07, the tests listed with their observed messages in the development record. Mutation checks (F06/F07): restoring the `C.ANCHOR` fallback, a plain `min(ebfmi)`, accepting booleans, `fit_draw` taking T0 from `C.ANCHOR`, and dropping the expected-index check each turned at least one test red.

**Verification (2026-10-31, R 4.5.1, pinned library `library-2026-10-15-R4.5.1`, Python 3 with NumPy 2.4.3):** `Rscript tools/run_tests.R` gave 6119 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings, against 5691/0 on 917825b (baseline before these changes). `python3 -m pytest -q tests/python` gave 202 passed, against 114 passed on 917825b. `stale_constant_audit.py --strict` gave 0 errors after updating four line pointers in `constant_registry.yaml` (one, the HYPNO quote, was already stale on 917825b). Combined run on the final working tree before commit; the release marker is written by a rerun on the committed revision.

## 2026-11-04: version 4.3.1, refit explainer published and linked (branch `review_v4_2026_10_31`)

**Changes published numbers (live 4.3.0):** no. No equation, parameter, draws file or selection rule changed. The site gains a link on the main page and a section at the top of About, both pointing to the refit explainer served from `www/`. The explainer itself is new. Its AICc comparison (`docs/refit_explainer/aicc_check.py`) was computed for it and is not used by the app. Record and audit: `docs/development/2026_11_04_refit_explainer.md`.

| Item | Record |
|---|---|
| What was wrong | The site had no account of how its tumor-model refit was done. The first draft of the explainer had stale numbers: its sensitivity table held calibrated-DD values from before the 2026-10-26 ESS repair. |
| Proof it was real | The prior sd 10 row read −3.3 (−6.1, −0.7), while `fit/output_refit_2026_10_25c_prior_c10/dd_calibrated_onset7_n200/summary.json` gives −3.39 (−6.11, −0.76). Seven of nine rows differed by up to 0.1 Gy₁₀. |
| Test that pins it | `tests/testthat/test_v431_refit_explainer.R::the explainer's sensitivity table matches the fit summaries` failed on all 9 rows with the pre-audit explainer. `::the main page links to the explainer` and `::the About section links to the explainer` failed with the 4.3.0 `app.R` and `about.md`. `::the served page is built from the committed source` failed with the pre-audit source. Each passes with this change. |
| Verification | See the commit message of the 4.3.1 release for the full-suite count on the committed revision. |
| What it moves | Nothing numeric. Version 4.3.1. |

## 2026-11-06: version 4.4.0, refit support always shown and the list of refitted trials (branch `v4-trials-support`, not deployed)

**Changes published numbers (live 4.3.0):** no. No equation, parameter, draws file or selection rule changed. The main table always shows the Refit support column, and Details opens a list of the trials of the Bayesian refit. Record: `docs/development/2026_11_06_v440_trials_and_support.md`.

| Item | Record |
|---|---|
| What was wrong | Requested changes (Igor Shuryak, 2026-11-06), not defects. The Refit support column was hidden unless Details was on, and the site did not list the trials the refit used, their schedules or the drugs given. |
| Proof it was so | On 4.3.1 (`ba9fe51`) the default main table has 7 columns and no "Refit support" header (`test_v4_home.R` asserted that it was absent), and no page or dialog lists the 19 trials. |
| Tests that pin it | `tests/testthat/test_v44_trials_support.R`: `::the Refit support column and its footnote are shown without Details`, `::Details offers the fitted-trial list`, `::the fitted-trial list follows the fitting data named in the draws manifest`, `::a trial list that does not match the fitting data is refused`, `::the trial table gives each trial's reference once, spanning its arms, with links`. With 4.3.1's `app.R` and `R/v4_home.R` they gave 5 failures and 3 errors, and they pass with this change. `tests/python/test_fit_trial_list_builder.py` (5 tests): removing the nimorazole entry or the disagreement check each turned a test red. |
| Verification | 2026-11-06, R 4.5.1, pinned library `library-2026-10-15-R4.5.1`, Windows: `Rscript tools/run_tests.R` on the working tree gave 6188 expectations, 2 failed (the two version pins in `test_v27_colleague_review.R` and `test_v31_stage_e.R`, updated to 4.4.0 and rerun with 0 failures), 0 errors, against 6155/0 on `ba9fe51`. `python3 -m pytest -q tests/python`: 196 passed (191 + 5 new), 11 skipped, and 5 failed for missing `cmdstanpy` and a Windows-R path, the same 5 as on `ba9fe51` in this environment. The certifying run on the committed revision is recorded in the development note. |
| What it moves | Nothing numeric. Version 4.4.0. |

### 2026-11-06, second round on the same branch: footnotes and changeable limits

**Changes published numbers (live 4.3.0):** no at the defaults. With the default limits (5, 5, 5) the list, support and exports are unchanged. Other limits give other lists by design. Record: `docs/development/2026_11_06_v440_trials_and_support.md`, "Second round".

| Item | Record |
|---|---|
| What was wrong | Requested changes: drug columns replaced by footnotes; trials not in the 2018 and 2019 papers marked; the three 5 Gy limits fixed in code. |
| Proof it was so | At `ce41191` the dialog had 9 columns including Chemotherapy and Other drugs, and `V4_MARGINS` was the only source of the limits (`v4_regimens()` and the export used it unconditionally). |
| Tests that pin it | `tests/testthat/test_v44_trials_support.R`: `::the trial table gives each trial's reference once, with footnotes for drugs and new trials`, `::margins: checked, and out-of-range values are refused with a message`, `::the main page follows the limits the user sets, and passes them on`, `::the Report card receives the limits in use`, `::a main-list export with changed limits replays exactly`. Red with `ce41191`'s `app.R`, `R/v4_home.R` and `R/v4_export.R` (9 failures, 4 errors), green with this change. `tests/python/test_fit_trial_list_builder.py` checks the 2018 reference column. |
| Verification | 2026-11-06, Windows, R 4.5.1, pinned library: full suite on the working tree 6222 expectations, 1 failed (`test_v43_review.R`, F08 expectation for an edited tumor margin, updated as described in the record), 0 errors. Affected files rerun with 0 failures. Builder tests 5 passed. Certifying run: see the development record. |
| What it moves | Nothing at the default limits. |

### 2026-11-06, third round: the plot no longer goes blank while a limit is retyped

**Changes published numbers:** no. Record: `docs/development/2026_11_06_v440_trials_and_support.md`, "Third round".

| Item | Record |
|---|---|
| What was wrong | With a limit box empty or out of range (mid-retype, or 25), the plot section went blank and the table kept the old list without explanation. |
| Proof it was real | Headless Chrome against `8da8c18`: after Backspace in the tumor box, `#home_plot` was absent and the plot UI held only an empty error element (68 characters, against about 16,500 with a plot). |
| Test that pins it | `tests/testthat/test_v44_trials_support.R::the main page follows the limits the user sets, and passes them on` (NA, −1 and 25 keep the last valid limits, rows and plot, with a message). Errored with `8da8c18`'s `app.R`, passes with the fix. |
| Verification | Browser probe after the fix: plot present in all 17 states, the note only in the 3 invalid ones. Full suite on the committed revision: see the development record. |
| What it moves | Nothing numeric. |

### 2026-11-06, fourth round: review and audit of the branch

**Changes published numbers:** no. At the default limits the list, support and exports are unchanged (the 4.3.1 export fixture replays exactly). Record with every finding, its evidence and what was left for a decision: `docs/development/2026_11_06_v440_trials_and_support.md`, "Fourth round".

| Item | Record |
|---|---|
| What was wrong | The plot ignored the late limit and dropped rows with a single dose. "−0" for a zero limit. Edited exports replayed. The trial list said no chemotherapy above a note allowing previous chemotherapy, and showed ranges, mixes and whole-arm times as exact values. Stale docs and a failing stale-constant audit. |
| Proof it was real | `app.R` at `b4826f3`: `geom_hline(yintercept = c(-mg$tumor, mg$early))`. Browser at limits 3/8/2: 15- and 14-fraction rows listed with no bar. CAIR 2006 Table 1 gives only whole-arm medians (40, 54 days). `stale_constant_audit.py --strict`: exit 1 with 3 line-site errors. |
| Tests that pin it | `tests/testthat/test_v44_trials_support.R`: `::the dialog's opening text is accurate about drugs, strata and what was fitted`, `::the plot draws one dashed line per limit, in the colour of its bars, and the caption names all three`, `::a box Shiny reports as NULL keeps the limits in use rather than resetting them`, `::replay refuses an edited export and still reproduces a real 4.3.1 export exactly`, `::the trial table gives each trial's reference once, with footnotes for drugs and new trials`, `::the fitted-trial list follows the fitting data named in the draws manifest`. Red with `b4826f3`'s `app.R`, `R/v4_home.R` and `R/v4_export.R`, green here. `::the support column and the failure explanations use the limits in force` and `::the main table has 8 columns, and 10 with Details` also pass on `b4826f3` (coverage of behavior that was already right). |
| Verification | Recorded with the certifying run of this revision in the next section of this file. |
| What it moves | Nothing numeric. |

### 2026-11-06: certifying run of the branch after the review (`34bbe8d`)

`Rscript tools/run_tests.R` on Windows, R 4.5.1, pinned library `library-2026-10-15-R4.5.1`, revision `34bbe8d` (committed, clean tree): 6272 expectations, 0 failed, 0 errors, 0 skipped, 0 warnings, suite marker written for `34bbe8d`. Baselines: 6155/0 on `ba9fe51` (4.3.1) and 6188/0 on `ab9152e` (4.4.0, first round). This run covers the second, third and fourth rounds above, whose earlier counts were on working trees. `python3 -m pytest -q tests/python`: 192 passed, 11 skipped, 4 failed, all four `ModuleNotFoundError: No module named 'cmdstanpy'` (not installed on this machine; the same tests failed the same way on `ba9fe51`). `stale_constant_audit.py --strict`: 0 errors. A deploy from a later revision (this record, or a merge into `v4`) needs the suite rerun on that revision.

## 2026-11-09: version 4.5.0, independent review of 4.4.0 (deployed 2026-11-09 at `faf6817`)

**Changes published numbers (previous live 4.4.0):** yes, at the defaults, because the startup tumor values change from scenario D18 to the trial-fit point scenario (decision D4). Main list: 60 qualifying doses in 15 rows (16 to 30 fractions) before, 69 in 18 rows (13 to 30) after; 58 common, 11 added, 2 removed (19 × 2.75 and 21 × 2.55 Gy); doses meeting the tumor condition under both models 17 before, 65 after. Refit support values are unchanged (they come from the draws). Under D18, chosen explicitly, every number is unchanged. No equation, draws file, fitting dataset or selection rule changed. Review: `HNSCC_BED_v44_independent_review_and_action_plan_2026-11-09.md` (outside the repository). Record: `docs/development/2026_11_09_v45_review_implementation.md`.

| Fix | What was wrong (proof) | Test (red on `3b72f85` or by mutation) | What it moves |
|---|---|---|---|
| Named profiles (4.1) | Refit values with DD rate 0.6 or g 0.02 kept the refit label; D18 after the refit kept the refit's g (review: +0.8102 Gy at 18 × 2.90). | `test_v45_review.R::a named profile label never covers values that differ from it (review 4.1)`; matching on four fields (old rule) gave 4 failures. | Labels only; startup values (with D4). |
| Main snapshot (6.2) | After a transient draw failure the ZIP's points and its scenario could come from different draw states. | `::a Main ZIP agrees with its own scenario after a draw-loading failure recovers (review 6.2)`: 2 failures on `3b72f85`. | Exports only. |
| Exact replay ID (6.3) | Deleting `analysis_id` let an edited file (early margin 5 → 0) recompute 18 points as a "replay". | `::exact replay needs the exported ID; an intentional edit is a new analysis (review 6.3)`: 2 failures, 1 error on `3b72f85`. | Replay API: `main_rerun()` refuses; `main_recompute()` added. |
| Mode messages (5) | Sensitivity called a directly entered DD calibrated; Main claimed a larger DD loss when DD was unusable. | `::Sensitivity never calls a directly entered DD calibrated; Main does not invent a DD loss`. | Text only. |
| Byte-exact checkout (6.1) | With `core.autocrlf=true`, a fresh clone of `origin/v4` matched 0 of 4 draw hashes. | Clone check (this branch: 4 of 4); `.gitattributes`. | Nothing numeric. |
| Bundle and deploy library (6.4) | A bundle without both explainer copies passed; the deploy library was added to the path before it existed. | `test_deploy_preflight.R::the bundle check requires each linked document, present and not empty`, `::the deploy tools library is created before it is put on the library path`. | Release checks only. |
| Six destinations (5) | Ten menu pages, three of them overlapping. | `test_v45_review.R::the menu has six destinations including the main page, and every former task keeps its ids`; walkthrough 54 of 54. | Layout only. |
| EQD2, cautions, count (4.1, 4.2) | Tumor change only in BED; extrapolation and low support only under Details; "15 schedules" counted rows. | `::EQD2 is a display of the same BED change...`, `::the count names dose and fraction combinations...`, `::D18 and the trial-fit point reproduce the review's fixtures (review 4.1)`. | Display only; exports keep BED. |
| Trial evidence (2.2, 2.3) | CAIR note denied period medians that its treatment report gives; Jackson and PMH times lacked their sources; metadata joins accepted a substituted PMID. | `test_v45_trials.R` (8 of 9 red on the old code), `tests/python/test_fit_trial_list_builder.py` (2 red). | Text and table only. |

**Verification.** 2026-11-09, R 4.5.1, pinned library, Windows R from WSL: `Rscript tools/run_tests.R` gave 6485 expectations, 3 failed (the DESCRIPTION version, then fixed), 0 errors; the certifying run on the committed revision is in the development record. Baseline on `3b72f85`: 6269 expectations, 0 failed, 1 error caused by an edit made while that run was in progress (the old `refit_preset()` call). `python3 -B -m pytest -q tests/python`: 208 passed. Strict stale-constant audit: 0 errors, 0 warnings. Local walkthrough: 54 of 54.

**Deployed** 2026-11-09 from `faf6817` as `custom_bed_calculator_v4` with `Rscript tools/deploy.R`, account `ishuryak`, after the shinyapps.io token exposed in the history of `main` (`.Rhistory`) was replaced and the old token deleted; the configured token and secret were checked not to match the exposed ones. Certifying suite on `faf6817` (committed, clean tree, Windows R 4.5.1, pinned library): 6485 expectations, 0 failed, 0 errors, marker written. Bundle: 55 files. Live walkthrough: 54 of 54 (`docs/development/2026_10_28_v4/walkthrough/20261109_173609_live/`). A first live run gave 53 of 54: the phone check for the frozen first column ran 4 s after the page opened and the live table took about 5 s. A separate probe found 18 rows and 20 frozen cells; the walkthrough now waits for the table (up to 60 s) instead of sleeping. It also no longer prints a failure note beside a passing served-file check. Rollback point: 4.4.0 at `62a0a67`.

## 2026-11-10: early-toxicity plan step 2, historical label tables (local, not pushed, not deployed)

**Changes published numbers:** no. The app does not read `data/early_tox/` yet, and the main list is unchanged (69 doses in 18 rows).

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| Historical tables from the Pardo-Montero 2021 supporting file (`tools/early_tox/build_historical.py` writes `data/early_tox/{schedules,fractions,labels,reports}.csv`) | The exploratory refits of 2026-11-10 included the 8 first-part-of-gap records that Pardo-Montero excludes, and keyed schedules on (n, total dose, days), which merges 4 pairs with different timing (plan R1, R2). | `tests/python/test_early_tox_historical.py`: treating partial records as complete gave 2 failures (`::test_raw_labels_and_selection_fixtures`, `::test_partial_courses_never_enter_a_likelihood`); a coarse key gave 1 (`::test_full_vector_identity_separates_the_false_merges`, 92 against 96 groups); keeping disputed labels in the primary selection gave 1; accepting an unknown label gave 1 (`::test_unknown_label_is_refused`). | Data tables only. Selections: primary 103 records (13 events), sensitivity 107 (15). |

**Verification.** 2026-11-10, revision recorded in the commit that adds this row: `python3 -m pytest -q tests/python/test_early_tox_historical.py` gave 5 passed. Rebuilding from the spreadsheet after the mutation checks reproduced the committed tables byte for byte. Record: `docs/development/2026_11_10_early_toxicity_step1_target_definition.md` (step 1) and the step 2 record.

### 2026-11-10, second part: randomized-arm tables (local, not pushed, not deployed)

**Changes published numbers:** no. The app does not read `data/early_tox/`.

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| Arm-level mucositis counts from 36 trial reports, each extracted twice independently (`tools/early_tox/build_arms.py`, `data/early_tox/{arms,outcomes,arm_schedules,arm_fractions,evidence_families}.csv`) | No arm-level outcome data existed: the early-toxicity parameters rest on Fowler 2003 and historical labels only (plan section 2). | `tests/python/test_early_tox_arms.py`: accepting one extraction gave 1 failure (`::test_a_count_found_in_only_one_extraction_is_refused`); Saturday as a W5 day gave 9; a gap off by one gave 2 (`::test_schedule_patterns_reproduce_the_protocols[RTOG9003_AFXS-...]`); ignoring the same-day interval gave 2; ignoring the screening log gave 1 (`::test_duplicate_reports_add_no_counts_even_when_both_extractions_agree`); accepting an unknown pattern token gave 1. | Data tables only: 45 outcomes, 43 arms, 15 candidate families. |

**Verification.** 2026-11-10: `python3 -m pytest -q tests/python` gave 232 passed (213 before); `stale_constant_audit.py --project . --project-root . --strict` gave 0 errors. Rebuilding after the mutation checks reproduced the committed tables byte for byte. Record: `docs/development/2026_11_10_early_toxicity_step2_data.md`.

### 2026-11-10, third part: early-toxicity schedule evaluator (local, not pushed, not deployed)

**Changes published numbers:** no. The app does not use `tools/early_tox/evaluator.py`.

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| One repair-aware evaluator for explicit fraction schedules (plan R8) | The exploratory fitter used repair with a 3.2 h half-time while `exploratory_after_audit/implic.py` used complete repair, so fit and prediction disagreed (plan R8). | `tests/python/test_early_tox_evaluator.py`: an inclusive clock gave 1 failure (`::test_fowler_no_repair_fixtures_from_the_step2_arm_schedules`); an undoubled repair cross-term gave 2; dropping the onset hinge gave 1 (`::test_onset_cancels_from_a_contrast_only_when_both_courses_pass_it`); half-time used as a rate gave 2; no ordering check gave 1. | Nothing in the app. Fixtures: 53.10544, 55.86597, 55.43495 (from the step 2 arm schedules), 84.12426 and +0.28402 Gy10. |

**Verification.** 2026-11-10: `python3 -m pytest -q tests/python` gave 239 passed (232 before).

### 2026-11-10, fourth part: early-toxicity steps 3 (mapping) and 4 (identification) (local, not pushed, not deployed)

**Changes published numbers:** no. Nothing in the app reads `data/early_tox/` or the step 4 outputs.

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| Endpoint mapping M1-M4 (`tools/early_tox/build_mapping.py`, `data/early_tox/endpoint_mapping.csv`) and identification run (`tools/early_tox/step4_identification.py`, outputs in `docs/development/2026_11_10_early_toxicity/step4/`) | The plan forbids fitting before identification is checked: randomized contrasts were claimed to identify Tk in plan version 1, which is false when both arms pass onset (plan R3). Confirmed here: changing Tk from 0 to 16 days changes no primary contrast (largest change 1e-10). | `tests/python/test_early_tox_step3_step4.py`: Tandon window changed, CAIR severe confluent in the primary, onset hinge removed, likelihood sign flipped, alpha/beta division removed; each gave 1 failure. | Nothing in the app. Recorded: slope and time effect identified, alpha/beta weakly informed, Tk not identified up to 16 days. |

**Verification.** 2026-11-10: `python3 -m pytest -q tests/python` gave 245 passed (239 before). Step 4 run on `e2de68e`, seed 20261110, 200 replicates, all converged. `sim_provenance.py validate sim_provenance.yaml`: registry OK. Record: `docs/development/2026_11_10_early_toxicity_step4_identification.md`.

### 2026-11-11: early-toxicity step 5, first refit (local, not pushed, not deployed)

**Changes published numbers:** no. Nothing in the app reads the step 5 outputs. The Fowler 2003 limit stays the default.

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| Arm-level refit (`tools/early_tox/step5_model.stan`, `step5_fit.py`) with C1 chemotherapy term, C2 Tk scenarios, C3 overdispersion; SBC (`step5_sbc.py`); historical classifier (`step5_historical.py`) | Fowler's fixed parameters were never fitted to randomized-arm toxicity counts (plan goal). Fitted: the binomial model fails its predictive check (p = 0.0025); the overdispersed candidate gives r 1.51 (1.00-2.41) BED/day against Fowler's 0.79, and an odds ratio of 2.27 (1.50-3.49) for 20 x 2.75 Gy against 70/35, where the Fowler benchmark gives 1.22. | `tests/python/test_early_tox_step5.py` (8) and `test_early_tox_step5_historical.py` (6); 16 mutations each gave at least 1 failure (list in the step 5 record, section 6). | Research results only. Step 6 validation decides whether anything reaches the app. |
| Fowler acute alpha/beta registered (`constant_registry.yaml::fowler_acute_ab`); `evaluator.py` FOWLER_2003 split one key per line, same values | The value 10.0 sat unregistered in two files, so a change in one would not be caught. | `stale_constant_audit.py --strict` gave 1 error with the evaluator value changed to 11.0, 0 with it restored. | Nothing. |

**Verification.** 2026-11-11, revision `f47cd4d` for all fits: `python3 -m pytest -q tests/python` gave 259 passed (245 before). All 27 primary fits pass the sampler gate (no divergences, R-hat < 1.01, ESS > 400). SBC with 500 replicates: uniform ranks for all 7 parameters (smallest p = 0.063). `stale_constant_audit.py --project . --project-root . --strict`: 0 errors. `sim_provenance.py validate`: registry OK. Record: `docs/development/2026_11_10_early_toxicity_step5_fit.md`.

### 2026-11-11, later: early-toxicity data expanded (10 trials) and step 6 validation (not deployed)

**Changes published numbers:** no. The app does not read these outputs. Fowler 2003 stays the default.

| Change | What was wrong (proof) | Test (red by mutation) | What it moves |
|---|---|---|---|
| Builder accepts a count only from mucositis rows (`build_arms.is_mucositis_endpoint`) | ELAN-RT's pharyngeal grade 3 count 11/91 in extraction E1 equals the mucositis count, and was accepted as "both extractions". | `test_early_tox_arms.py::test_only_mucositis_rows_confirm_a_count` and `::test_every_count_is_in_both_independent_extractions` (4 failures with the filter removed) | Nothing for the 45 earlier outcomes (byte-identical); O068 now cites the mucositis row |
| One-extraction counts allowed only for sensitivity arms with status `adjudicated_single` | ELAN-RT's mucositis row is hidden by a print banner; one extractor recovered it | `::test_adjudicated_single_extraction_is_refused_for_a_candidate_arm` (red without the role check) | ELAN-RT enters as sensitivity only |
| 10 trials added (7 primary, 3 sensitivity) from 12 new references, double extraction C-E | Step 5 rested on 15 trials; the time effect relied mainly on two arms with chemotherapy differences | Updated expectations in `test_early_tox_step3_step4.py` and `test_early_tox_step5.py`; new chemotherapy codes red when mutated | Research results: r 1.55 (1.07-2.45), previously 1.51 |
| Step 6 validation (`tools/early_tox/step6_validate.py`) | The plan requires held-out validation before any default change | `tests/python/test_early_tox_step6.py`, 7 tests, each red by a targeted mutation | Research results: the candidate passes the numerical rule on held-out trials but HYPNO and ELAN-RT favor Fowler; recommendation is not to promote |

**Verification.** 2026-11-11: `python3 -m pytest -q tests/python` gave 268 passed (261 before step 6). Step 5 expanded on `439d466`: 34 of 34 primary fits pass the gate, SBC uniform (smallest p = 0.28). Step 6 on `e80aa12`: 46 of 46 fits pass the gate; rescoring with new seeds reproduced the totals to within 0.1. `stale_constant_audit.py --strict`: 0 errors. Record: `docs/development/2026_11_11_early_toxicity_step6_validation.md`.

