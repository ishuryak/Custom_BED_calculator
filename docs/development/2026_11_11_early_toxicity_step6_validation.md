# Early-toxicity plan: expanded data (10 trials added) and step 6 validation


> **Correction, 2026-11-11 (independent review).** Several statements in this record are superseded.
> See `2026_11_11_early_toxicity_review_response.md`, section 4:
> - The evaluable totals are 4,599 (15 families) and 6,042 (22 families), not 4,500 and 6,256.
> - HYPNO's ranking depends on its planned 4 vs 5 cisplatin courses.
> - The conditional coverage figures were computed without conditioning on the reference arm.
> - The benchmark is a Fowler-parameter statistical model, not the app's Fowler limit.
>
> The numbers below are kept as the record of what was run.

Written 2026-11-11. These are research results only: the app is unchanged, and Fowler 2003 stays the
default. Nothing is deployed. The decision in section 5 is Igor Shuryak's.

## 1. New references and what they added

On 2026-11-11 Igor supplied 12 PDFs, found by a PubMed search (`References/README.md`, section H),
and Claude downloaded 3 open-access ones. Every report went through the step 2 double extraction:
two independent extractors per batch (batches C, D and E, files `2026_11_10_early_toxicity/extraction/out_{C,D,E}{1,2}.json`),
following the unchanged `extraction/BRIEF.md`.

| Report | Decision | Count used (grade ≥3 or confluent, per arm) |
|---|---|---|
| Mortensen 2012, DAHANCA 6&7 | Included (DAHANCA 7 part, by treatment received, since DAHANCA 6 glottic is context C2) | 187/383 (5 a week) vs 269/398 (6 a week), confluent |
| Ghosh-Laskar 2016, Tata 3-arm | Included | 12/57 RT, 23/65 RT + cisplatin, 20/64 6 a week |
| Grover 2021 | Included (CTCAE v4.03, so group G2) | 23/51 (70/35) vs 20/50 (66/30), cisplatin in both |
| Gupta M 2015 | Included | 9/54 (5 a week) vs 18/55 (6 a week), radiotherapy alone |
| Gupta P 2020 | Included | 8/30 conventional + cisplatin, 10/30 concomitant boost, 15/30 boost + cisplatin |
| Tallari 2017 | Included | 9/48 (2 Gy daily) vs 14/47 (1.1 Gy twice daily), same chemotherapy |
| Singh 2024 | Included | 19/40 (5 a week) vs 26/41 (6 a week), same cisplatin |
| Sanchiz 1990 | Sensitivity only: about 11% nasopharynx, internally inconsistent tables, WHO scale | 14/277, 11/282, 31/300 |
| Dwivedi 2023 | Sensitivity only: alternate allocation, not randomized | 5/34 vs 11/34 |
| Ortholan 2026, ELAN-RT | Sensitivity only: context C6 (frail elderly, split course). See below | 18/93 (70/35) vs 11/91 (55/20 split) |
| Bentzen 2001, CHART | Not derivable: 75% vs 44% with no per-arm denominators | none |
| Poulsen 1999, TROG 91.01 | Not derivable: 90% vs 65% matches more than one count | none |
| Marcial 1987, RTOG 79-13 | Not derivable: 13% vs 23% "severe", no denominator, and the extractors disagreed (22/94 vs 21/92) | none |
| Miszczyk 2014, CHA-CHA | Excluded: weekly median Dische scores only | none |
| Maciejewski 1990 | Excluded: not randomized | none |

The data now hold 22 primary trial families, 51 arms and 6,256 patients. Before, there were 15
families, 35 arms and 4,500 patients.

**ELAN-RT.** The supplied PDF was printed through a proxy, and the banner hides the mucositis row of
Table 3. Extractor E2 recovered the row from the page's vector layer, and Claude confirmed it: the
grades sum to the any-grade totals (71 and 79), and grade ≥2 matches the text (44 and 54).
- Extractor E1 could not read the row.
- The count enters with status `adjudicated_single`. The builder allows that only for a
  sensitivity-only arm with a note.
- A clean publisher copy would turn it into a normal double extraction.

**Builder defect found and fixed.** The builder accepted a count when both extractions held the same
numbers, whatever their endpoint. E1's ELAN-RT pharyngeal grade 3 count (11/91) equals the mucositis
count, so it passed as "both extractions".
- Now only mucositis rows confirm a count (`build_arms.is_mucositis_endpoint`).
- Rebuilding left all 45 earlier outcomes byte-identical.
- The tests are `test_early_tox_arms.py::test_only_mucositis_rows_confirm_a_count` and
  `::test_adjudicated_single_extraction_is_refused_for_a_candidate_arm`. Both were red when their
  line was mutated.

**Chemotherapy coding (C1), new arms:**
- Tata RT + cisplatin: +1.
- Gupta P concomitant boost alone: −1.
- Sanchiz + 5-FU: +1.
- Gupta P boost + cisplatin: 0. It receives 3 × 75 mg/m² as in its reference arm, on different days
  (`CHEMO_SAME_DOSE`).

## 2. Step 4 and step 5 on the expanded data

- **Step 4** (`2026_11_10_early_toxicity/step4_expanded/`, 200 replicates):
  - Design rank 3, condition number 2.1.
  - Expected posterior SD relative to the prior: log r 0.21, log α/β 0.67 (was 0.69).
  - α/β is still weakly informed.
  - Tk still cancels from every contrast (largest effect 1e-10).
- **Step 5** (`2026_11_10_early_toxicity/step5_expanded/`, revision `439d466`, same settings as
  the first fit):
  - All 34 overdispersed fits pass the sampler gate.
  - The binomial model again fails its predictive check (U02 p = 0.012), so the overdispersed set
    stays primary.
  - SBC (500 replicates): uniform ranks for every parameter, smallest p = 0.28.

| Unit | Arms | r (BED/day) | α/β (Gy) | γ (chemo) |
|---|---:|---|---|---|
| U02 primary | 51 | 1.55 (1.07–2.45) | 12.9 (7.6–21.8) | 0.55 (0.26–0.84) |
| U04 Tk 21 | 51 | 1.64 (1.11–2.66) | 13.0 (7.6–21.9) | 0.57 (0.28–0.87) |
| U07 α/β prior SD doubled | 51 | 1.51 (1.04–2.35) | 20.1 (8.4–51.9) | 0.56 (0.26–0.85) |
| U08 chemotherapy-differing trials dropped | 32 | 1.20 (0.69–2.04) | 12.0 (7.1–20.9) | n/a |
| U10 W1 windows only | 32 | 1.55 (1.00–2.62) | 12.0 (7.1–20.0) | 0.52 (0.12–0.88) |
| U11 plus sensitivity trials | 64 | 1.55 (1.10–2.28) | 13.4 (8.3–22.2) | 0.63 (0.38–0.89) |
| U14 no chemotherapy term | 51 | 1.52 (1.00–2.52) | 12.4 (7.5–20.9) | n/a |

Leaving out one family at a time, r ranges 1.37–1.65. The historical classifier still prefers
Fowler: leave-one-source-out log score −25.90 for Fowler against −30.29 for the step 5 point.

## 3. Step 6 method

`tools/early_tox/step6_validate.py`, revision `e80aa12`.

**Fits.**
- Benchmark U01: Fowler's α/β and r fixed.
- Candidate U02.
- Both use the step 5 primary settings, so the two models differ only in α/β and r.
- For each of the 22 primary families, both models are refitted without it: 44 fits, all passing
  the sampler gate.
- Leakage control: no prior in either model is estimated from any trial in the data. A committed
  test checks that no training set contains its held-out family, and the run asserts it for every
  fit.

**Tasks** (plan step 6). The score is the log predictive density of the held-out family:
- New trial: the family's baseline is drawn from the fitted population.
- Conditional contrast: p(other arms | reference arm), which is the app's use, a schedule against
  70/35.

**Settings.**
- Monte Carlo: 100 inner draws per posterior draw, with common random numbers for the two models.
- Bayesian bootstrap over families: Dirichlet(1,…,1), 20,000 replicates, seeds 20261113 and
  20261114.
- External checks: models fitted on all 22 families predict HYPNO (the R6 benchmark, never fitted)
  and the 6 sensitivity trials.

**Stability.** Rescoring with new Monte Carlo seeds and 200 inner draws reproduced every total to
within 0.1 and the bootstrap bounds to within 0.01.

## 4. Step 6 results

**Held-out primary trials (22 families):**

| Task | Total log score, Fowler | Total, candidate | Sum of differences | Candidate better | Mean difference per family (95% Bayesian bootstrap) |
|---|---:|---:|---:|---:|---|
| New trial | −195.0 | −188.8 | +6.2 | 18 of 22 | +0.28 (+0.11 to +0.55), both seeds |
| Conditional contrast | −107.5 | −101.4 | +6.1 | 18 of 22 | +0.28 (+0.10 to +0.55), both seeds |

- **Influence:** the leave-one-family mean stays positive whichever family is removed
  (0.18–0.31).
- **Largest gains:** Vienna (+2.3), RTOG 9003 (+1.1) and GORTEC 99-02 (+0.8).
- **Fowler better:** PMH, Grover, EORTC 22791 and EORTC 22851 (−0.2 to −0.3).
- **Calibration:** 23 of 29 held-out contrasts fall inside the 90% predictive interval for both
  models.

**External checks** (conditional contrast, fitted on all 22 families):

| Trial | What it tests | Difference (candidate − Fowler) | Observed log OR | Candidate prediction (90%) |
|---|---|---:|---:|---|
| **HYPNO** | 55 Gy/20 vs 66/33 at 6 a week | **−0.68 (Fowler better)** | −0.15 | +0.42 (−0.21 to +1.07) |
| **ELAN-RT** | 55/20 split vs 70/35 | **−0.42 (Fowler better)** | −0.54 | −0.16 (−0.85 to +0.48) |
| IAEA-ACC | 6 vs 5 a week | +0.20 | +0.76 | +0.50 |
| Poddar | 6 vs 5 a week | +0.14 | +0.55 | +0.50 |
| Pinto | hyperfractionation | +0.07 | +0.51 | +0.15 |
| Dwivedi | 64/25 vs 70/35 (not randomized) | +0.13 | +0.96 | +0.83 |
| Sanchiz | twice daily, + 5-FU | +0.98 | −0.26 and +0.75 | +0.10 and +0.56 |

## 5. Interpretation and the decision for Igor

**What passes.** The candidate meets the plan's numerical promotion convention on held-out primary
trials:
- the lower 95% bootstrap bound is above zero in both seeds, and
- no single family reverses the sign.

**What does not.** The plan's convention also requires "no unresolved endpoint or context mismatch".
- The candidate's advantage comes from **accelerated** trials: more fractions per week or per day,
  at 2 Gy or less.
- The app's decision region is **moderate hypofractionation**: fewer, larger fractions.
- The two randomized tests of that region, HYPNO and ELAN-RT, both favor Fowler. The candidate
  over-predicts the toxicity of the shorter arm in both.
- That is a context mismatch in the region where the app makes its recommendations.

**Possible explanations (inferences, not tested).**
- In accelerated trials, overall time may stand in for fractions per day, interfraction interval
  and field technique.
- A single regrowth rate fitted mostly to acceleration may not carry over to hypofractionation.

**Recommendation.** Do not promote the candidate to the default. Keep Fowler as the default, and offer
the candidate only as a labeled research scenario ("fitted to accelerated trials, over-predicted
hypofractionation toxicity in HYPNO and ELAN-RT").

**Alternatives.**
- (a) Promote, since the numerical rule passes. Not recommended, because the app's region
  contradicts it.
- (b) Neither model in the app until a model that separates acceleration from hypofractionation is
  fitted. This would need another plan step.

## 6. Tests and verification

- **Unit tests.** `tests/python/test_early_tox_step6.py` has 7 tests. Each was turned red by a
  targeted mutation:
  - the binomial constant
  - the linear predictor dropped
  - the reference-arm marginal
  - the bootstrap scaled
  - training leakage
  - the G2 offset dropped
  - the HYPNO reference arm swapped
- **Updated expectations.** `test_early_tox_arms.py`, `test_early_tox_step3_step4.py` and
  `test_early_tox_step5.py` now expect the expanded data. The new chemotherapy codes and the
  same-dose rule were each confirmed red when mutated.
- **Suite:** `python3 -m pytest -q tests/python` gave 268 passed.
- **Stale-constant audit:** `stale_constant_audit.py --project . --project-root . --strict` gave
  0 errors.
