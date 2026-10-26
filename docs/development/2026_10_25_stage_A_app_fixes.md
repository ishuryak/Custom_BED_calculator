# Stage A: app fixes from the 2026-10-25 review

Plan: `docs/development/PLAN_codex_review_2026_10_25.md`, section 4, stage A (approved by Igor Shuryak on 2026-10-25). Review: `HNSCC_BED_v3_comprehensive_review_and_action_plan_2026-10-25.md` (sections 4, 5, 8. Appendix A S1-S9. Appendix B M3, M4. Appendix D OAR-01 to OAR-05. Appendix F D1-D3).

## Outcome

Every stage A finding is fixed and pinned by a test that was observed red on the unchanged code (revision 677ce21, working tree) and green after the fix. Each fix was also reverted once, line by line, and its test turned red again (two mutation batches, below). The full suite ends with 0 failures and 0 errors (counts at the end).

Not changed: the equations, the D18 defaults, the fixed anchor, the both-model criterion, the +1 Gy threshold, the conservative ranking rule (apart from M3's invalid-draw handling), `fit/`, `data/refit_*`, `draws/`, `References*`, the draws file format and `load_draws()`. APP_VERSION stays 3.0.0, EQUATION_VERSION 3.0, SELECTION_VERSION 2.6, SCHEMA_VERSION 3.

Versions bumped because meaning changed:
- `UNCERTAINTY_VERSION` 3.0 -> 3.1. Uncertainty 3.0 files are refused with a message: they did not record the applied organ limits or the fixed organ α/β, so they cannot be replayed exactly.
- `OAR_VERSION` 3.0 -> 3.1. Deterministic exports now record `search$oar$oar_version`. A schema-3 file whose search applied organ limits under the old rules is refused with a message ("Run the search again in the app"). Files without organ limits replay exactly as before (the frozen v2.6 fixtures still pass).

Golden fixtures (test file `tests/testthat/test_v31_stage_a.R`):
1. HYPNO 20 x 2.75 Gy/26 d against itself, cord 100%, QUANTEC, limits applied: `share_tox` 0 for prior and fixed toxicity, and the downloaded file replays as 0 in every column (before: replay gave 1).
2. D18 against 35 x 2 Gy/47 d, threshold -3, cord 100% applied: zero deterministic candidates. Clearing the cord α/β (or the other-organ α/β) now blocks the candidate search, the both-model screen and the uncertainty tab with "Organ limits are applied ..., but an organ input is invalid" (before: three unchecked rows).
3. Nonempty: threshold -3 lists 30/29/28 without organ limits and with a non-binding cord (50%) applied.

## Per finding

Test paths are in `tests/testthat/test_v31_stage_a.R` unless named. "Red" quotes the failure on the unchanged code.

| ID | What was wrong (symptom) | Reproduction | Test | Red then green | Moves a displayed number? |
|---|---|---|---|---|---|
| S9 | Clearing an applied organ's α/β, or entering an organ with no verified constraint, turned the check into "no check" | D18 vs 35x2/47, -3 Gy, cord 100%: 0 rows. Cord α/β cleared: 3 rows. Mandible 200% (QUANTEC): rows identical to box off, uncertainty share_tox 1 | `S9: clearing an applied organ's alpha/beta blocks...`. `S9: an entered organ with no verified constraint...`. `S9: the checkers refuse...` | Red: "Expected `suggestions()` to throw a error" (x4), unc-summary error "Enter the cord α/β." without the organ reason. Green. Mutation (policy returns NULL on error): red | Yes: such inputs now block with a reason instead of showing rows. Report card shows "Organs at risk: not assessed: an organ input is invalid" and lists entered unassessed organs |
| S1 | After the search, changing settings left old rows on screen and in the ZIP, under a JSON of the new settings | HYPNO, -3, 15-19 fx, d max 3.5. Change threshold to +20: table and CSV unchanged, JSON +20 | `S1: changing any effective input...` (10 input classes: threshold, refit, tox source, rule, required share, d max, fraction range, late mode, late tolerance, organ dose) | Red: 40 of 50 expectations ("listed" still shown, search CSV still exported, `sc$search` present). Green. Mutation (ignore key): red | No numbers move. Stale results are hidden, a warning says the settings changed, and the ZIP holds only the current summary |
| S2 | Uncertainty export dropped the applied organ limits, so replay turned share_tox 0 into 1 | Golden fixture 1 | `S2 golden: HYPNO against itself...` | Red: `sc$oar$applied` FALSE. Replay share_tox 1 vs 0. Share_joint 0.998 vs 0. Green (all columns, prior and fixed). Mutation: red | Replayed values now equal the exported ones |
| S8 | Fixed-toxicity uncertainty used hidden organ α/β 2/3, not the tab's controls | Cord 83% (45.65 Gy/20 fx), α/β 0.5: tab "exceeds", uncertainty share_tox 1 | `S8: fixed-toxicity uncertainty uses the organ alpha/beta...`. `S8: one fixed draw gives the deterministic organ check...` | Red: share_tox 1 (expected 0), replay 1. Green: 0 at α/β 0.5, 1 at 2. Mutation: red | Yes, fixed-toxicity share_tox with organ limits applied and nondefault organ α/β |
| S3 | Deterministic README omitted `R/schedule.R` and `R/oar.R` | Six-day scenario in a clean R process: "could not find function schedule_time" | `S3: the deterministic README replays six-day, Thursday-start and organ scenarios in a clean R process` (callr, emitted README lines) | Red: exactly that error. Green. Mutation: red | No |
| S4 | Uncertainty replay accepted wrong versions, prior, RNG, seed, counts, enums | 15 single-field mutations accepted | `S4: replay refuses incompatible or malformed uncertainty metadata, and old 3.0 files`. `S4: the search part ... replays ... exactly` (new `rerun_uncertainty(sc, part = "search")`) | Red: 14 mutations accepted, 3.0 accepted. Search replay did not exist. Green: all refused with "Cannot rerun". Search replays every column. Mutation (drop prior check): red | No |
| S5 | Uncertainty candidates used entered skipped days. Deterministic candidates did not | 20 fx with day 2 skipped: deterministic T 26, uncertainty T 29 | `S5: both searches build candidates on the regular calendar...` | Red: candidate T with gaps (and sidebar text absent). Green. Mutation: red | Yes, uncertainty candidate T when skipped days are entered. The evaluated regimen keeps its gap (T 29) |
| S7 | `schedule_time()` coerced 2.5 to 2 and accepted NA, negatives, 0, 401, text | Direct API calls | `S7: the direct calendar API rejects...` | Red: 7 inputs accepted. Green. Mutation: red | No |
| OAR-01 | Source physical mean limits were compared as EQD2 | Parotid 25 Gy physical mean in 35 fx: EQD2 18.57 "within" the 20 Gy goal | `OAR-01: a parotid physical mean of 25 Gy...`. `OAR-01: the search checks use the same physical mean comparison` | Red: status "within", limit 26.6 (converted). Green: "exceeds the source limit", limit 20, EQD2 shown as illustration. Mutation: red | Yes: every mean-dose organ status and every search with a mean-dose organ applied. Maxima unchanged. Cord 2-Gy anchor kept |
| OAR-02 | HYHOPE PRVs not offered. Strict "<" read as "≤". Unassessed mandatory metrics silent | Cord 37.0 Gy/15 fx passed. No cord + 5 mm row | `OAR-02: HYHOPE PRVs are organs of their own...`. `OAR-02: HYHOPE's strict '<' limit...`. `OAR-02: unassessed mandatory metrics are listed...` | Red: no PRV rows. 37.0 "within 5%". No coverage text. Green. Mutation (operator ignored): red | Yes at exact boundaries. Two new organ rows |
| OAR-03 | HYHOPE "Max" could be read as a defined metric | Source metric not stated | `OAR-03: HYHOPE 'Max' keeps its unspecified volume definition...` | Red. Green. No D0.03cc anywhere in R/ or the CSV | No |
| OAR-04 | Upper end of a source range presented as the source's limit. Cochlea 35 Gy alternative hidden | Parotid 17-22 Gy flag | `OAR-04: the upper end of a source range is labeled...` | Red. Green | No (labels) |
| OAR-05 | "HYHOPE protocol 4.1.4.2" looked like a version. 12-fraction level not flagged | CSV source_short | `OAR-05: HYHOPE rows name protocol version 4.2...` | Red. Green | No (labels) |
| D1 | DVH Dmax was the top bin's midpoint | Edges 0/40/50 Gy, volumes 100/1/0, cord 20 fx: Dmax 45, "within 5%" | `D1: the DVH maximum is an interval...` | Red: "within 5% of the converted constraint". No bounds. Green: "indeterminate at this DVH resolution", bounds 40-50. Mutation: red | Yes, DVH maximum statuses. The table shows "40.0 to 50.0 (DVH bin)" |
| D2 | NaN-free but invalid curves accepted | 12 synthetic malformed files | `D2: malformed DVHs are refused...` | Red: 12 accepted. Green. Mutation: red | Malformed files are refused |
| D3 | Quoted, BOM or padded identifier headers passed. UI said "never stored" | Synthetic `"Patient ID"` CSV accepted | `D3: quoted, BOM-prefixed and padded identifier headers are refused` | Red: 4 accepted, "temporary file on the server" absent. Green. Mutation: red | No |
| M3 | Conservative score dropped invalid draws | 1 valid +10 draw and 99 invalid: min_q10 +10 | `M3: one valid +10 draw...`. `M3: the 10th percentile is evaluable only while...` (0, 1, 9 % evaluable. 10, 90, 99, 100 % unavailable. No flooring) | Red: +10. N_valid_both absent. 10/90/99 % gave 10. Green. Mutation: red | Yes: scores fall when some draws are not evaluable. At the default search (35x2/47, +1, 90%), all 26 choices are unchanged and screen scores fall by 0.05 to 0.57 Gy |
| M4 | The 1,000-draw screen's choice was not re-checked on all draws | Required share 100%: 8 x 4.5 Gy/10 d and 11 x 3.8 Gy/15 d listed at 99.95% | `M4: at a required toxicity share of 100%...`. `M4 and wording: at the default 90%...` | Red: no role column, failed choices listed. Green: the failed choices are marked "failed the all-draw check" and the next-ranked dose is listed. None fail at 90%. Mutation: red | Yes at 100%: those two fraction counts list the next-ranked dose |
| Sec. 5 | "eligible" read as tumor eligibility. 5 x 5 Gy looked like it met +1 | Default search message | `M4 and wording: ... 5 x 5 Gy is not shown as meeting +1` | Red: old wording, no tumor column. Green: the review's wording, a "Tumor threshold met (share of draws): DI / DD / both" column, 5 x 5 Gy shows both < 5% | No numbers. New columns |

Existing tests changed because the behavior they pinned was the defect (each change is commented in the test):
- `test_v3_app.R::organ rows reach the report card...`: "Parotid ..., mean EQD2" -> "mean physical dose" (OAR-01).
- `test_v3_oar_schedule.R::cGy and Gy, % and cc are read...`: midpoint `max_dose` 30 -> bin bounds 29.99-30 (D1. The review asked for this line to be revisited).
- `test_v3_oar_schedule.R::organs without a verified QUANTEC constraint...`: an entered unsupported organ now reads "not assessed: no verified constraint in this library" (S9).

## Mutation checks (each fix reverted, its test red)

Batch A (one run, each reversion hits its own test): S7 (skip days coerced again), OAR-01 (mean compared as EQD2), D1 (interval override off), D2 (CSV validator off), D3 (identifier columns ignored), M3 (old NA-dropping quantile). The first batch-A S9 reversion was malformed (it made validation unconditional) and was redone in batch B.
Batch B: S1 (ignore run key), S2 (organ block dropped), S8 (default organ α/β), S3 (old README line), S4 (prior check off), S5 (skipped days in candidates), M4 (no all-draw check), OAR-02 (operator ignored), S9 (policy returns NULL on invalid input). Every listed test failed. Files restored byte for byte from a backup (checked with `cmp`).

## Decisions the plan did not settle

1. **S5 policy.** Candidates in both searches use the regular calendar (fractions per week and first day). Skipped days apply to the evaluated regimen only. This keeps the deterministic selection (2.6) unchanged. Stated under the skip-days input and in the uncertainty search help. Recorded as `search$candidate_calendar`.
2. **S1.** Invalidation rather than a saved context: results carry the scenario ID of the settings they were run with. Any difference hides them and drops them from the ZIP.
3. **OAR-01 when fractionations differ.** Per decision 12, a mean-dose limit is compared on physical dose with the source fractionation stated, also when the evaluated regimen has another fraction count (for HYHOPE, the nearest table). For a regimen much more hypofractionated than the source, a physical comparison can understate the biological effect. This is flagged in every such row. It needs Igor's confirmation.
4. **M4.** Both options: the failing choice stays visible as "failed the all-draw check", and the next-ranked dose is evaluated.
5. **Deterministic replay of organ searches.** Old files with an organ block are refused (organ rules 3.0 vs 3.1), instead of bumping SELECTION_VERSION, which would have broken replay of all files without organs.
6. **DVH validation thresholds.** A % DVH must start within 0.5 points of 100%. The curve must end at zero volume (≤ 1e-6 of the total). Doses above 300 Gy are refused as a likely unit error. Blank identifier columns are allowed. Any value in one is refused.
7. **Operators.** HYHOPE rows "<". QUANTEC rows from the quote: cord, brainstem, larynx and cochlea "≤". Parotid and submandibular "<" ("less than ≈20 Gy", "<35 Gy").

## Full suite

Command (repository root, Windows R 4.5.1, release library): `Rscript tools/run_tests.R`

- Before (revision 677ce21, unchanged tree): expectations 4961, failed 0, errors 0, skipped 0, warnings 0 (5 min 16 s).
- After (same tree plus these changes, 2026-10-25): expectations 5230, failed 0, errors 0, skipped 0, warnings 0 (8 min 56 s). The new file adds 269 expectations.

## Deferred, and why

- `about.md` lines 83, 125 and 209 still describe statuses only as "converted constraint" and say uploads are "never stored". Plan stage E owns the Model and sources page wording.
- S6 (release Git identity): stage F. The baseline run wrote `.suite_passed` although an untracked file was present (`--untracked-files=no`), as the review describes.
- Larynx V50 ≤ 27% cannot be checked from one entered value. It is listed as unassessed. cc DVHs are not checked against a structure total.
- An optional prespecified minimum tumor-pass rule (review section 5) was not added. The tumor-pass shares are shown instead.
- M2 (units), M8 (MCSE wording), M1 and the draws themselves: stage E.
