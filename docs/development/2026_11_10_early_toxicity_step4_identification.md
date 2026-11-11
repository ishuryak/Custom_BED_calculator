# Early-toxicity plan, steps 3 and 4: mapping applied, identification results

Written 2026-11-10. Step 3 applies decisions M1-M4, which Igor Shuryak approved on 2026-11-10
("proceed autonomously as recommended"). Step 4 uses the trial DESIGN only. No posterior for the
real outcomes has been computed, the app is unchanged, and step 5 (the first refit) has not
started.

## 1. Step 3: mapping applied

`tools/early_tox/build_mapping.py` writes `data/early_tox/endpoint_mapping.csv`.

| Set | Outcomes | Families |
|---|---:|---:|
| Primary (G1 and G2, windows W1 and W2) | 35 | 15 |
| Sensitivity: W1 windows only | 23 | 9 |
| Secondary: CAIR severe confluent | 2 | 1 |
| Sensitivity trials: IAEA-ACC, Pinto, Poddar | 6 | 3 |
| Benchmark: HYPNO | 2 | 1 |

**M1 check.** The NCI CTCAE v4.03 quick reference (downloaded from evs.nci.nih.gov on 2026-11-10,
local copy `NCI_CTCAE_v4.03_2010_QuickReference.pdf`) defines oral mucositis grade 3 as "Severe
pain; interfering with oral intake". It is symptom-based, as assumed. The trials cite v4.0. The
check used v4.03, which may differ in wording.

**What M1 can and cannot do.** The scale-family offset sits in the baseline. Every G2 family uses
one scale for all of its arms, so within-trial contrasts cancel the offset. The offset is therefore
informed only through the distribution of baselines across trials. It matters for absolute-risk
prediction, including the HYPNO benchmark, and not for the schedule parameters.

## 2. Step 4: identification

Run: `python3 tools/early_tox/step4_identification.py --reps 200`.
- Revision `e2de68e`, seed 20261110, one thread, 622 s.
- Python 3.12.3, NumPy 2.4.3, SciPy 1.17.1.
- Outputs are in `2026_11_10_early_toxicity/step4/`.
- The simulation parameters are registered in `sim_provenance.yaml` (ids `et_rec_*` and
  `et_prior_*`), and the registry validates.

Set-up:
- **Model:** the plan's arm-only candidate, logit p = η_trial + κ·(E − E_reference).
- **Exposure E:** from `tools/early_tox/evaluator.py`, with a 3.2 h repair half-time and the elapsed
  clock.
- **Reference arm:** the conventional or control arm of each family.
- **CAIR:** E is the patient-weighted mean of its two dose-per-fraction eras.

### Per-parameter results

| Parameter | Result | Evidence |
|---|---|---|
| Design | Separable | Across 20 contrasts the dose, quadratic and time differences have rank 3 and condition number 2.2. Dose and time differences correlate at 0.78, because accelerated arms are shorter and lower in dose |
| Slope κ | **Identified** | Recovery at κ = 0.1: 90% of estimates within 0.086-0.113 and 95% Laplace intervals cover 96%. Expected posterior SD of log κ is 0.09 (0.18 if L* = 0.5) |
| r (time effect) | **Identified** | Expected posterior SD of log r is 0.22 of the prior SD (0.40 if L* = 0.5, 0.13 if L* = 2). Recovery: 90% of estimates within 0.65-0.96 around 0.79, coverage 95.5% |
| α/β | **Weakly informed** | Expected posterior SD of log α/β is 0.69 of the prior SD (0.87 if L* = 0.5, 0.54 if L* = 2). Recovery: 90% of estimates within 7.5-12.7 around 10, coverage 97.5%. The posterior will depend visibly on the prior, so prior sensitivity must be reported (D6) |
| Tk (onset) | **Not identified for Tk up to 16 days** | Moving Tk from 0 to 16 days changes no contrast at all (largest change 1e-10). Above 16 days the contrasts change only through the shortest arms (Vienna V-CHART, 16 elapsed days, and GORTEC very accelerated, 23). By plan R3, Tk is fixed within named scenarios and no fitted Tk interval is reported |

### Prior predictive (D3 priors, 2,000 draws)

| Quantity | 5% | 50% | 95% |
|---|---:|---:|---:|
| Largest within-trial risk difference | 0.6 points | 14.3 points | 65.3 points |
| Median within-trial risk difference | 0.1 | 3.0 | 21.1 |
| Mean arm risk | 3% | 57% | 98% |
| Odds ratio per 10 BED units | 1.05 | 1.92 | 6.90 |

The priors cover the observed range of within-trial differences: for example, about 50 points in
Vienna, 26 in CAIR and 21 in RTOG 9003. The baseline prior is wide: it allows arm risks from 3% to
98%.

### What step 4 did NOT test

The recovery simulations assume the model is right. They do not test four things:
- chemotherapy differences between arms (section 3)
- scale offsets
- arm-level overdispersion
- misspecified repair or onset

Each is a step 5 question.

## 3. Decisions before step 5 (recommendation first)

- **C1. Chemotherapy differences within trials.** Five primary families compare arms that also
  differ in planned chemotherapy:
  - ORO 93-01 adds carboplatin/5-FU at an identical schedule: 8/54 against 24/50, an observed odds
    ratio of about 5.3 with zero BED difference.
  - GORTEC 99-02 gives 3 cycles, 2 cycles and none.
  - RTOG 0129 gives 3 against 2 cisplatin cycles.
  - Vienna adds mitomycin C to one arm.
  - Bartelink uses 6 against 10 mg/m² of daily cisplatin.

  Without a term for this, the slope κ absorbs or is distorted by chemotherapy effects.
  *Recommended:* add one within-trial main effect, "arm has concurrent chemotherapy that its
  comparator lacks or has less of", with a normal(0, 1) prior on the log-odds scale, and run a
  sensitivity analysis that drops the five families. Alternatives: drop the chemotherapy-differing
  arms from the primary fit, or keep the model as planned.
- **C2. Onset scenarios.** Tk = 7 days (Fowler reference), plus 14 and 21 days as named scenarios.
  *Recommended.* Between 7 and 14 the fit is identical by construction. Only the 21-day scenario
  can differ.
- **C3. Overdispersion.** Fit without an arm-level random effect, and add one, half-normal(0.5),
  if posterior predictive checks show more spread than binomial. *Recommended.*
- **Gate items still open:** your spot-check of the step 2 counts (list in the step 3 proposal,
  section 4), and your agreement to start step 5, which the plan requires.

## 4. Tests

`tests/python/test_early_tox_step3_step4.py` has 6 tests. Each was turned red by a targeted
mutation:
- the Tandon window changed
- CAIR severe confluent put in the primary
- the onset hinge removed
- the likelihood sign flipped
- the α/β division removed

The ordering was tests, then a 3-replicate smoke run, then the tests again, then the 200-replicate
run.
