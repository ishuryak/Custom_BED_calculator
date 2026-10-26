# Plan: response to the v3.1 review of 2026-10-26

Review: `../HNSCC_BED_v31_review_and_action_plan_2026-10-26.md` (project root), evidence archive beside it.
Reviewed revision: 535b338 (deployed e4ff22f). Plan written 2026-10-26. Nothing has been edited yet.

## 1. Verification of the findings before planning

Each finding was checked against the code or data before being accepted. All ten hold.

| ID | Checked how | Result |
|---|---|---|
| R01 | Read the gate in `fit/fit_dd_calibrated.py:153-154` and `fit/repair_dd_calibrated.py:71-72`; counted `min_ess_bulk` and `min_ess_tail` in every conditional-DD file | **Confirmed.** The conditional-DD gate checks divergences, tree depth, R-hat and E-BFMI, not ESS. Below bulk ESS 400: primary 177/1000 (tail 15), restricted 272/1000 (tail 24), and in the eight sensitivity sets 26-46 of each 200 (total 299). **748 conditional fits in all.** Minimum bulk ESS 151. |
| R01 (claim) | Read README line 59, my 2026-10-26 report to Igor | The README says ESS applies to DI, but "Every fit ... passes" reads as the full gate. My report to Igor ("every fit passes the convergence gate you set") was wrong for ESS. |
| R02 | Read `fit/export_draws.py:48-96, 125-129` | **Confirmed.** Export checks only the parent identity string, not `parent_di_draws_sha256` or the data hash, compares only α and α/β of each pair, and hashes whichever data file has the expected name. Holdout cache not yet read; the reproduction will be a Phase 0 test. |
| R03 | Evidence fixtures named in the review | Accepted on the reviewer's reproduction; becomes a Phase 0 red test before any fix. |
| R04 | Read `app.R:575-581, 672-673` | **Confirmed.** `tumor_delta_txt` indexes the scalar reference `y[i]`, which is NA after row 1. Live visible defect. |
| R05 | Read `R/oar.R:292, 403-421` | **Confirmed.** `status_exceeds()` maps "indeterminate" to FALSE and `flip_text()` reads FALSE as "within". |
| R06 | Development record already discloses it (equal-weight check) | Confirmed; the record calls it negligible, which is a magnitude, not a fix. |
| R07 | Text collected below | Confirmed for every item. |
| R08 | Builder lines cited | Accepted; Phase 0 test confirms `value_ci` is dropped. |
| R09 | `tools/run_tests.R:10-14` | Accepted; the marker is written from `git rev-parse` without checking its exit status. |
| R10 | Read `R/oar.R:222-223` | **Confirmed.** Rows are filtered by a numeric regex before validation, so a malformed interior row disappears. |

Nothing was dismissed. One review statement is now out of date: the branch is pushed (535b338 on origin), so the in-app source link can resolve.

## 2. Proposed release structure

- **3.1.1 (runtime fixes and wording, same draws):** R03, R04, R05, R07, R09, R10, the replay grid bound, and the correction of the gate claim. No number from the draws changes, so saved 3.1 scenarios keep replaying exactly.
- **3.2 (new draw generation):** R01, R02, R06 and the targeted refits, then R08 provenance. New draw version; 3.1 draw identities archived so 3.1 scenarios replay or are refused explicitly.
- **Later (separate scientific work package, not in this plan):** the review's Phase 4 sensitivity questions (estimator-aware uncertainty, prescription mixtures, population exclusions, incomplete repair, Monte Carlo precision study, external validation) and the optional tumor-share filter.

Every change follows the change gate: a test observed red without the change and green with it, a mutation check, evidence that the fix runs in the app, one combined smoke and the full suite before any run or deploy. Each lands with a `docs/FIXES.md` row.

## 3. Items for 3.1.1

### P1. R04 benchmark deltas (app.R:575-581)
- Current: `tumor_delta_txt <- function(x, y) { vapply(seq_along(x), function(i) { if (is.na(x[i]) || is.na(y[i])) ...`
- Proposed: require `length(y) == 1` (or `length(x)`), otherwise stop; use the scalar for every row.
- Test: at defaults the three benchmark rows read HYPNO test −1.39/−8.67, HYPNO control +1.11/−2.87, HYHOPE −10.60/−15.14 and equal Compare regimens; mismatched length errors. Walkthrough gains a check on all three rows.

### P2. R05 indeterminate through the α/β sweep (R/oar.R:403-421)
- Proposed: classify each grid point as within / exceeds / indeterminate / not assessed; `flip_text()` says "within at every α/β" only when every point is definitively within, otherwise names the indeterminate range. CSV export carries the same state.
- Test: the reviewer's cord fixture (edges 0/40/50 Gy, volumes 100/1/0, 20 fractions) reads indeterminate throughout; plus definite-below, definite-above, mixed and boundary-equality cases.

### P3. R10 Eclipse parser (R/oar.R:222-223)
- Proposed: once the data header is found, every following non-blank line up to the next section must parse as finite numbers, otherwise stop with the line number.
- Test: "0 100 / 40 NaN / 50 0", Inf, missing volume and a text token each fail; existing valid Gy/cGy/%/cc fixtures still parse with the same row count.

### P4. R03 one organ-mapping validator
- Proposed: one function used by deterministic replay, uncertainty replay, direct API and search: named list, unique nonblank supported keys, finite scalar values, valid library/version/mode; an active policy with unusable entries is refused rather than treated as "no organ restriction". Disabled/empty policies keep working.
- Test: the reviewer's `accepted_unnamed_oar.json` and `accepted_duplicated_oar.json` are refused in both replay paths; blank, unknown, non-scalar, nonfinite keys refused; well-formed cord-only and multi-organ scenarios replay exactly; historical no-organ scenarios replay.

### P5. Replay grid bound (R/export.R:567-586)
- Proposed: share the deterministic search's grid-size calculation and refuse an uncertainty replay above the same documented candidate limit before allocation.
- Test: below / at / above the bound; normal exported searches unaffected.

### P6. R09 suite marker fails closed (tools/run_tests.R)
- Proposed: check the exit status and format (40 hex characters) of `git rev-parse HEAD`, refuse the marker when tracked files are dirty, and print "numerical check only, not release-certified" when Git is unavailable.
- Test: failing git, empty output, dirty tree and failed tests write no marker; a clean passing run writes the HEAD revision.

### P7. Correct the gate claim now (README:59, About, Model and sources)
- Current README: "Every fit behind the shipped draws passes 0 divergences, R-hat ≤ 1.01, bulk and tail ESS ≥ 400 (DI) and E-BFMI ≥ 0.3, and the app refuses draws whose manifest does not say so."
- Proposed: "The DI fits pass 0 divergences, R-hat ≤ 1.01, bulk and tail ESS ≥ 400 and E-BFMI ≥ 0.3. The conditional DD fits pass the same checks except ESS: in version 3.1, 177 of 1000 primary and 272 of 1000 restricted conditional fits have bulk ESS below 400 (minimum 151). These are refitted in version 3.2."
- The same correction goes into the stages B-E record ("all fits pass the gate") and the About refit section, and the manifest's DD `gate_pass` is renamed to state its narrower policy (P8 of 3.2 replaces it).

### P8. R07 wording (exact replacements)

| Location | Current | Proposed |
|---|---|---|
| app.R:857-861 (search empty message) | "The Model agreement tab explores a separate comparison: by default HYPNO 55 Gy/20 fractions, an allowed change of −3 Gy under each model, and at least 1 day saved. A limit defined relative to the comparison also uses that tab's comparison." | "The Both-model screen uses the global comparison and threshold by default; you can override them there. An allowed BED loss changes the hypothesis being screened." |
| about.md:107 (cohort) | "...randomized trials of definitive radiotherapy alone in head and neck squamous cell carcinoma other than nasopharyngeal cancer..." | "...randomized trials of definitive fractionation, mainly in head and neck squamous cell carcinoma without cytotoxic chemotherapy. Exceptions within arms: RTOG 83-13 included some nasopharyngeal and sinus cases, and DAHANCA included nasopharynx and gave nimorazole to nonglottic patients; arm-level results cannot remove them." |
| about.md:114 (chemotherapy sensitivity) | "trials with the same chemotherapy in both arms" | "trials with a common chemotherapy policy in both arms (planned or delivered exposure can still differ, e.g. four vs five weekly courses in HYPNO)" |
| about.md:116 (α) | "α is about 0.04 per Gy because it is a population-level slope that includes patient-to-patient variation." | "α is about 0.04 per Gy. It is an effective population-level slope given the endpoint, likelihood, priors, clock and arm-level data, not a measurement of intrinsic cellular radiosensitivity. Patient-to-patient variation may contribute, but the fit does not isolate the cause." |
| about.md:103 and 143 (TCP) | "...needs a dose-response relation this calculator..." / "The calculator has no dose-response relation, so it cannot say what these deficits mean for tumor control." | "D18 is not fitted to outcomes. The Uncertainty tab gives a model-conditional tumor control illustration from the refitted α and an entered reference control rate; it does not establish clinical efficacy or non-inferiority." |
| about.md:124 and mod_oar.R:29 (QUANTEC) | "...cochlea mean 45 Gy, all in 2-Gy fractions)" / "QUANTEC 2010 (2-Gy fractions, converted)" | "...cochlea mean 45 Gy). Maximum-dose values are EQD2 and converted; mean-dose values are physical means from conventionally fractionated courses and compared physically." / "QUANTEC 2010 (conventional fractionation: maxima converted, means physical)" |
| about.md:122 and organ table flag (mean scope) | "...this comparison can understate the effect of larger fractions, and the row is flagged." | "...the row is flagged: physical source comparison only; the effect of changed fractionation is not validated and may understate toxicity of hypofractionation. For HYHOPE, a mismatch with the selected table's fraction count is flagged." |
| R/uncertainty.R:119-120 (acute onset prior) | "Mucosal onset (elapsed days): Fowler 2003 day 7; Pardo-Montero 2021 refit near day 19 (as reported in the v3 plan review, to be re-checked)." | "Mucosal onset (elapsed days), analyst sensitivity prior Uniform(7, 19): Fowler 2003 day 7; Pardo-Montero 2021 refitted Tk 19.44 inclusive days (about 18.4 elapsed) jointly with other parameters. Neither is a direct biological range." |
| about.md (agreement count) | "92% of the 55 regimen pairs" | Generated from the regimen list (11 regimens, 55 pairs) and the text names the 11. |
| about.md (predictive coverage) | "45 of 45 arms inside their 90% predictive interval" | "...an in-sample check using fitted intercepts and arm effects, not prospective validation." |

Ripple: the matching About-page tests (`test_v3_uncertainty.R` patterns) change with the text; listed per item when applied.

Small items folded into 3.1.1: rename `loss_per_day_after_onset` in the DD point values to `lambda_over_alpha` (it omits the kill-probability factor); the Stan comment "keeps every onset" becomes "discourages onsets below".

## 4. Items for 3.2

### Q1. R01 one diagnostic policy
- One versioned policy object (thresholds, quantities checked, missing or nonfinite diagnostics fail) used by DI fit, conditional DD fit, repair, export and the app's manifest check.
- Recompute acceptance from the saved per-fit diagnostics, never from a stored `gate_pass`.
- Test: a synthetic fit with 0 divergences, R-hat 1.001 and ESS 150 fails; missing ESS fails; repair and export give the same verdict.

### Q2. Targeted refits (the costly step)
- Refit only the 748 conditional fits below ESS 400 (449 main, 299 sensitivity), as a new repair generation with its own identity referencing the original fits; originals kept. The DI fits are not refitted.
- Estimated cost: bulk ESS minimum 151 needs about 3 times the draws, so attempts at 4x and 8x chain length. From the earlier runs (1000 fits in about 1.5-2 hours at 1x, 35 repairs in 16 minutes at 8x), roughly 5-10 hours. Checkpointed per fit, as before.

### Q3. R02 identity chain and atomic export
- Export checks data bytes, parent draw bytes, parent α, α/β, g, λ, Tk per pair, draw counts and finiteness; recomputes anchor equality for every pair (≤ 1e-4 Gy after CSV, actual residual recorded); writes into a staging folder and swaps the whole set only when everything passes. Manifest hash comes from the fitted identity, not the file on disk.
- Holdout cache gets a full identity and refuses mismatches.
- Test: each mutation the reviewer lists is refused before any output; the 6.19-Gy broken-anchor probe never exports.

### Q4. R06 equal outer weights
- Every pooled DD summary weights each outer DI draw 1/1000 regardless of chain length; one machine-readable summary feeds the record and the About numbers.
- Test: duplicating all rows of one outer fit leaves summaries unchanged.

### Q5. Re-export, compare and release as 3.2
- New draw version; 3.1 draws archived by hash. Old vs new distributions and candidate lists compared and any material shift explained. The HYPNO numbers quoted in Pardee are rechecked against 3.2.

### Q6. R08 provenance sidecar
- A provenance table linked to the dataset carries the eight known source CIs, estimator, prescription categories (DAHANCA 62 Gy/31 for well-differentiated T1 glottic, weights unknown) and population exceptions. The numerical fitting file is unchanged, so no refit follows from it.

## 5. Decisions for Igor

1. **R01 policy:** refit the 748 conditional fits to meet ESS ≥ 400, as approved in decision 9 (recommended), or adopt a different documented gate for conditional fits? Recommended: refit. The reviewer accepts either if documented, but a changed gate after the fact is weaker.
2. **Two releases:** ship 3.1.1 (runtime fixes and corrected wording, same draws) first, then 3.2 with the refits? Recommended: yes. R04 is a visible defect on the live site and the gate claim is wrong today.
3. **R08 as a sidecar** rather than new columns in the fitting CSV? Recommended: sidecar, so it is clear nothing was refitted.
4. **Defer** the optional tumor-share filter and the review's Phase 4 scientific sensitivities to a separate work package? Recommended: defer.
5. **Deploy 3.1.1** after its checks pass, on the same guarded path (v3 only, rollback e4ff22f), and push? Recommended: yes.
