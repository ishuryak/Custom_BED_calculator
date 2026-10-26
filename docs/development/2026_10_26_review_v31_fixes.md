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
- Full R suite on b31db54: 5379 expectations, 0 failed, suite marker written for that revision (3.1.0 baseline 5269). Python 70 passed (baseline 38). Local browser walkthrough 35 of 35, including all three benchmark rows (`2026_10_21_v3/walkthrough/`).

## Results

HYPNO 55 Gy in 20 fractions against 35 x 2 Gy in 47 days, median tumor BED change in Gy at each
draw's α/β (80% range). DI from the full DI posterior (8000 draws in the main set). DD from the
conditional fits (1000 DI draws in the main sets, 200 in the sensitivity sets), each outer draw
weighted equally.

| Set | DI | DD | λ_DD median | Order agreement |
|---|---|---|---|---|
| main, primary | −2.02 (−4.31 to +0.97) | −3.27 (−5.92 to −0.53) | 0.405 | 0.919 |
| main, paper restrictions | −2.01 (−3.87 to +0.37) | −2.91 (−4.74 to −0.66) | 0.361 | 0.940 |
| ARTSCAN added | −2.69 | −3.66 | 0.379 | 0.926 |
| common chemotherapy policy added | −1.03 | −2.39 (to +0.20) | 0.427 | 0.912 |
| effective N halved | −1.89 | −3.11 | 0.417 | 0.917 |
| Kaplan-Meier rates only | −1.88 | −3.35 | 0.413 | 0.915 |
| three later trials added | −1.66 | −2.90 | 0.434 | 0.919 |
| nominal clock | +0.26 | −2.13 (to +1.03) | 0.481 | 0.882 |
| intercept prior sd 10 | −1.93 | −3.39 | 0.412 | 0.914 |
| RTOG 9003 updated report | −1.89 | −3.32 | 0.429 | 0.916 |

The DI results are unchanged from version 3.1, because the DI fits and draws are unchanged. The DD
main median moved from −3.25 to −3.27, within its Monte Carlo standard error. Monte Carlo
standard errors of the app's 4000-draw summaries, by resampling the 1000 DI draws: about 0.1 Gy
for medians, up to 0.25 Gy for 95% limits, about 1.5 percentage points for threshold shares. The
DI 80% range includes zero in every analysis. The planned (nominal) clock moves the DI median to
about zero.

## Statistical and clinical audits (2026-10-26, read-only subagents)

Statistician. Confirmed and fixed in the text: the DI numbers I had first written for 3.2 came
from the 1000-draw thinned subset, not the full posterior, and the claimed "sensitivity medians
moved by up to 0.2 Gy" was Monte Carlo noise (the 3.1 figures were right and are restored).
Background growth g is essentially its prior (posterior 0.005-0.020 per day, prior uniform
0.0046-0.021), now stated. The predictive check uses a new arm-level effect, not the fitted one,
and 45 of 45 inside 90% intervals shows the check has little power, now stated. "The sign of the
DI result depends on the clock" overstated, because the DI 80% range spans zero in every
analysis, now reworded. Monte Carlo errors are now quantified. Checked and correct: the repair is
not a source of selection bias. ESS failures were correlated with lower λ_DD (Mann-Whitney
p = 0.001), so dropping them would have biased DD (−3.32 instead of −3.27), but refitting them
kept every outer draw and moved per-fit medians by about 0.01 Gy. Equal outer weighting matches
the modular target. All 1000 main conditional fits needed more than the base chain length, which
the 2026-10-25 record understated. Not changed (noted for Igor): the order agreement shares
structure between the models, the DD rate interval is largely prior and penalty at its lower end,
and the uncertainty search re-checks its pick on draws that share the same DI draws.

Radiation oncologist and radiobiologist. Confirmed and fixed: four places still described the
cohort as "radiotherapy alone", and the approved sentence on population exceptions had not been
added, now both done, with the twice-daily arms (no incomplete-repair correction) and the single
arm above 2.25 Gy per fraction stated. The claim that "every table flags" extrapolation was
false, now corrected to the four tables that do. Reported for Igor's decision, not changed: the
report card and organ benchmark table show organ status without the α/β-range result or mean-dose
note; the "Limit, physical in this regimen" column label for mean limits; the extrapolation flag
is absent from the banner and the Results tab; the default's help text gives no scope; the
Results tab shows the radiotherapy-alone Fowler zone although the default acute limit is framed
for chemoradiation. Checked and sound: no stale text naming HYPNO as the default, calendar, clock
conventions, DD calibration description, α and α/β statements, HYPNO and HYHOPE descriptions,
QUANTEC values, Fowler and Fenwick scope, no recommendation language, organ status logic.

## Audit decisions (Igor Shuryak, 2026-10-26: "as recommended")

All five display items from the clinical audit were approved and applied, each with a block in
`test_v32_review.R` ("audit 1" to "audit 5") observed red against db4bd3a and green after:

1. Report card organ rows add the α/β-range result for maxima and the cross-fractionation note
   for physical means, and the coverage sentence ("selected constraint metrics, not full protocol
   compliance") appears under the card. The organ benchmark table gains the same column. Example:
   cord 44 Gy in 18 fractions reads "within 5%" at α/β 2 Gy but now also "exceeds for α/β 0.87-1.2 Gy".
2. The banner flags a fraction size outside the tumor data (the default 2.9 Gy is flagged on every tab).
3. The default's help text says it is an example schedule, not a published or trial-tested
   regimen and not a recommendation.
4. The Results tab shows both Fowler zones (59-61 Gy₁₀ radiotherapy alone, 54-58 Gy₁₀ with
   chemotherapy) and names the one the search uses.
5. The organ table's limit column is "Limit, physical (Gy)", each value marked "(source
   fractionation)" for means or "(converted to this regimen)" for maxima.

Deploy of 3.2.0 to `custom_bed_calculator_v3` (rollback e4ff22f) and push of this branch also approved.

## Deferred (not in this release)

The optional tumor-share requirement in the uncertainty search, and the review's Phase 4
questions: estimator-aware uncertainty, prescription-mixture and exclusion sensitivity,
population and systemic-modifier exclusions, incomplete repair in BID/TID trials, a nested Monte
Carlo precision study and validation against evidence not used in the fit.
