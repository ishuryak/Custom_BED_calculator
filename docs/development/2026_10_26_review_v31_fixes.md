# Version 3.2: response to the 2026-10-26 review of version 3.1

Review: `HNSCC_BED_v31_review_and_action_plan_2026-10-26.md` (project root, outside this repository).
Plan: `PLAN_codex_review_2026_10_26.md`. Branch `review_2026_10_26`, from 535b338.

## Decisions (Igor Shuryak, 2026-10-26)

- "Make this change and the other changes based on the review, then audit." Taken as approval of
  the plan's recommendations. The conditional DD fits that fail ESS are refitted rather than the
  gate changed. Provenance is kept in a table beside the fitting file. The optional tumor-share
  filter and the review's Phase 4 scientific analyses are deferred.
- The default evaluated regimen changes from 20 x 2.75 Gy (HYPNO test arm) to 18 x 2.9 Gy
  (52.2 Gy in 24 days). Igor first wrote "2.65 x 20". Asked, he confirmed he meant the 2.75 x 20
  default on the report card. The comparison (35 x 2 Gy) and DD anchor are unchanged. (Version 2
  used 20 x 2.65 Gy as its default, which is likely where that figure came from.)
- One release (3.2.0) instead of 3.1.1 followed by 3.2, because deployment was not requested.
- Not deployed and not pushed. Both wait for Igor.

## Verification of the findings before any change

All ten findings were confirmed against the code or data (plan, section 1). The most serious: the
calibrated-DD loop never checked ESS, so 748 of 3600 conditional fits were below bulk ESS 400
(primary 177 of 1000, restricted 272 of 1000, sensitivity sets 26-46 of each 200, minimum 151).
None had a divergence or R-hat above 1.01. The status report to Igor on 2026-10-26 ("every fit
passes the gate you set") was wrong on ESS. A correction note heads the 2026-10-25 record.

One review item was not applied as written. The Stan comment in `dd_cal.stan` ("keeps every
onset") is unchanged, because its bytes are bound into every conditional fit's identity and the
export verifies them. The description is corrected in `fit_dd_calibrated.py` and here: the
penalty discourages onsets before day 7 and does not keep them above it (0.5% of DD draws still
have a fitted-arm onset before day 7).

## What changed, by finding

| ID | Change | Test (each observed red without the change) |
|---|---|---|
| R01 | `fit/gate.py` is one fail-closed gate, ESS included, used by the DI fit, the DD loop, the new repair, the export and, in R, `draws_gate_failures()`. The app refuses draws older than version 3.2. `fit/repair_dd_ess.py` refits only failing fits (4x, 8x, 16x chains, seeds offset from every earlier attempt), keeps the originals in `repaired_originals_ess/` and records `repair_ess.json`. The first repair script is retired and refuses to run. | `tests/python/test_gate_and_ess_repair.py`, `tests/python/test_repair_dd_calibrated.py`, `test_v32_review.R` "R01" |
| R02 | The export verifies data bytes, parent draw bytes, parent link, variant, outer selection, Stan sources, all five parent values per pair, finiteness and ranges, and the anchor equality of every exported pair after CSV round trip. It stages the whole set and swaps `draws/` in one step. The manifest data hash comes from the verified identity. The holdout cache carries an identity. | `tests/python/test_export_draws.py`, `tests/python/test_loso_cache.py` |
| R03 | `validate_oar_percents()` is one validator for both replay paths and both direct checkers. | `test_v32_review.R` "R03" (9 malformed forms, both replay paths, JSON round trip) |
| R04 | `tumor_delta_txt()` takes a reference of length 1 or length(x). | `test_v32_review.R` "R04", walkthrough check on all three rows |
| R05 | A four-way status class runs through the α/β sweep. | `test_v32_review.R` "R05" |
| R06 | DD summaries weight each outer DI draw 1/N with weighted step-CDF quantiles, exactly invariant to duplicated rows. | `test_gate_and_ess_repair.py::test_equal_outer_weights_ignore_chain_length` |
| R07 | Wording replacements listed in the plan. | `test_v32_review.R` "R07" |
| R08 | `fit/build_provenance.py` and `data/*_provenance.csv`. The builder refuses ambiguous matches and absent quote checks. All 8 fitting files rebuild byte-identical. | `tests/python/test_provenance.py` |
| R09 | `release_revision()`: the suite marker needs a resolvable, committed, unchanged revision. | `test_v32_review.R` "R09" |
| R10 | Eclipse data rows are checked one by one, and a malformed row is refused with its line number. The table header must name both dose and volume, so metadata lines such as "Max Dose [cGy]:" are not taken for it. | `test_v32_review.R` "R10" |
| Grid | `grid_size()` is shared by both replay checks. Uncertainty replay is bounded at 60,000 regimens, above the app's largest grid (99 x 598 = 59,202). | `test_v32_review.R` grid block |
| Default | 18 x 2.9 Gy. | `test_v32_review.R` default block |

Mutation checks of the shipped lines, each of which turned its test red: no uncertainty grid
bound, the old 20,000 bound, a marker ignoring a dirty tree, an app gate ignoring bulk ESS, no
minimum draws version, a validator allowing duplicates or treating NaN as blank, unchecked
Eclipse rows, the old Eclipse header match, a gate without bulk ESS, a gate passing missing
values, unweighted rows, a repair trusting the stored `gate_pass`, a repair skipping identity, an
export without the anchor, parent-draw, variant or g checks, the old holdout cache and the old
builder.

## Independent code audit

A subagent audited `git diff 535b338..daa94e2` read-only. It confirmed three defects, all fixed
before the export:

1. The app still loaded the 3.1 draws (DD minimum bulk ESS 151) and would label results 3.2,
   because the new check applied only from draws version 3.2. Fixed: draws older than 3.2, or
   without a version, are refused.
2. The replay cap of 20,000 regimens refused grids the app can produce, because the fraction
   slider spans 1 to (comparison fractions - 1) and the comparison may have up to 100 fractions.
   Fixed: 60,000, derived from the app's input limits.
3. The first repair script still applied the narrow gate. Fixed: retired.

It also suspected that real Eclipse exports with "Min/Max Dose [cGy]:" metadata lines would be
misparsed (fixed), that the export's variant check defaulted open when the key was missing
(fixed) and that g, Tk and c/α were not range-checked (fixed). It reported that NaN organ doses
were read as blank (fixed). Checks that found nothing: key agreement between exporter and app
gate, weighting of every pooled DD statistic, file layouts across all three fit generations,
seed separation, the organ validator on UI and JSON inputs, Eclipse line numbers, the suite
marker, the holdout cache and every call site of the shared helpers.

## Costs and measurements

- ESS repair: 09:10 to 10:47 on revision f1a81ed, five folders in parallel (20 of 22 cores).
  748 fits refitted, every one passing at 4x, 8x or 16x chain length, none failing at 16x.
  A smoke run on a copy repaired two fits at 4x in about 10 s each.
- Export: anchor residual at most 4.5e-8 Gy. DI draw files byte-identical to 3.1. Both DD files
  regenerated.
- Full R suite 5379 expectations, 0 failed (3.1.0 baseline 5269). Python 70 passed (baseline 38).

## Results

Full conditional fits, each outer DI draw weighted equally. HYPNO 55 Gy in 20 fractions against
35 x 2 Gy in 47 days, median tumor BED change in Gy at each draw's α/β (80% range):

| Set | DI | DD | λ_DD median | Order agreement |
|---|---|---|---|---|
| main, primary | −1.95 (−4.42 to +1.19) | −3.27 (−5.92 to −0.53) | 0.405 | 0.919 |
| main, paper restrictions | −2.00 (−4.00 to +0.41) | −2.91 (−4.74 to −0.66) | 0.361 | 0.940 |
| ARTSCAN added | −2.54 | −3.66 | 0.379 | 0.926 |
| common chemotherapy policy added | −1.16 | −2.39 | 0.427 | 0.912 |
| effective N halved | −1.83 | −3.11 | 0.417 | 0.917 |
| Kaplan-Meier rates only | −2.00 | −3.35 | 0.413 | 0.915 |
| three later trials added | −1.48 | −2.90 | 0.434 | 0.919 |
| nominal clock | +0.14 | −2.13 | 0.481 | 0.882 |
| intercept prior sd 10 | −2.09 | −3.39 | 0.412 | 0.914 |
| RTOG 9003 updated report | −1.88 | −3.32 | 0.429 | 0.916 |

Version 3.1 reported −2.0 and −3.3 for the main set, sensitivity DI medians −2.7 to +0.3 and DD
−3.7 to −2.1. The main results moved by less than 0.1 Gy, the sensitivity medians by up to
0.2 Gy, partly from the change to equal weighting. The planned (nominal) clock still moves the DI
result to about zero, so the sign of the DI conclusion depends on that timing choice.

## Deferred (not in this release)

The optional tumor-share requirement in the uncertainty search, and the review's Phase 4
questions: estimator-aware uncertainty, prescription-mixture and exclusion sensitivity,
population and systemic-modifier exclusions, incomplete repair in BID/TID trials, a nested Monte
Carlo precision study and validation against evidence not used in the fit.
