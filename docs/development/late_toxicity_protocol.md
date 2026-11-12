# Late-toxicity analysis protocol (step L0)

Version 1, 2026-11-12. Status: **draft for Igor Shuryak's decisions D1–D4 below**. Records only:
no data have been extracted, no model has been fitted, and the app is unchanged (late BED uses
α/β = 3 Gy). Source: section 4–7 of `2026_11_12_bayesian_methods_and_late_toxicity_review.md`,
which this protocol implements as step L0. Later steps (L1–L11) follow the review's table in
section 7 and its stopping rules. The protocol is fixed before any late outcome is extracted, so the
choice of endpoint cannot depend on which one gives a convenient α/β.

## 1. Two objectives, kept separate

- **A. Fractionation sensitivity of a defined late tissue endpoint.** The aim is an estimate of α/β,
  and possibly repair half-time, for one tissue endpoint, labeled as that endpoint (for example
  "subcutaneous fibrosis"), never as "late toxicity" in general.
- **B. A clinically important late outcome for comparing schedules.** Severe late mucosal injury,
  swallowing dysfunction, feeding-tube dependence, necrosis and any-site severe toxicity are
  different outcomes. A model of one is not a general late-complication model.

## 2. Endpoint tuples

Every candidate outcome is specified completely before extraction:

| Field | Meaning |
|---|---|
| Tissue or composite | the organ or tissue, or the named composite |
| Grading | scale and version (RTOG/EORTC, LENT-SOMA, CTCAE, trial-specific) and threshold |
| Time origin and horizon | start of radiotherapy, end of radiotherapy or randomization, with a fixed horizon or none |
| Estimator | cumulative incidence (death competing), net risk (Kaplan–Meier, death censored), prevalence among those assessed, or worst grade ever recorded |
| Population | randomized, treated, assessed survivors, locoregionally controlled |
| Handling | death, recurrence, salvage surgery, missing assessments |
| Dose definition | prescription dose, or a tissue dose where reported (hot spots) |

Candidate tuples, from the review's extraction leads (to be confirmed against the full texts in
`References/late_toxicity_review_2026_11_12/` during L1–L2, not from this list):

1. Severe late mucous-membrane toxicity, RTOG/EORTC grade 3 or worse, worst grade recorded from
   90 days after the start of treatment (RTOG 9003 first report, Fu 2000).
2. Any-site grade 3–5 late toxicity, prevalence at 5 years among those assessed (RTOG 9003 final
   report, Beitler 2014). A distinct estimand from tuple 1 in the same trial family, never pooled
   with it. The PDF of this report could not be downloaded, so nothing from it enters the data until
   a full-text PDF is in `References/`.
3. Severe subcutaneous or soft-tissue fibrosis, peak incidence (counts) or 5-year Kaplan–Meier net
   risk (DAHANCA 6&7, Mortensen 2012), with counts and actuarial estimates kept as two observation
   models.
4. Tissue-specific late endpoints of CHART (laryngeal edema, subcutaneous fibrosis and telangiectasia,
   Bentzen 1999), actuarial, for the repair analysis only. The repair half-times derived in that paper
   are not used as priors if CHART outcomes enter the likelihood.
5. Late effects of RTOG 83-13 (Fu 1995), with death treated as a competing event as reported.

## 3. Rules fixed now

- **Primary-endpoint rule.** The primary analysis set is the endpoint tuple with the most compatible
  within-trial dose-per-fraction contrasts, judged on design information only (number of trial
  families, denominator completeness, follow-up, dose quality) after L3 and before any model is fitted.
  If no tuple supports a defensible fit, the result is an evidence synthesis and the app keeps its
  assumed α/β, now with explicit alternative scenarios.
- **No broadening for sample size.** Tuples are not merged into "mucosal or any site" to gain arms.
- **One record per cohort, endpoint, horizon and estimator**, unless a dependence model is specified.
  Later reports of the same trial are not independent trials.
- **Likelihood follows the estimator** (review section 6.2). A binomial likelihood is used only for
  genuine integer counts with a valid denominator, and an approximate likelihood for a fixed-horizon
  estimate with its standard error. An exact count is never made by rounding a percentage times a
  guessed denominator.
- **Exposure.** One tested evaluator for fitting and prediction, with incomplete repair
  $Q(h) = \sum_i d_i^2 + 2\sum_{i<j} d_i d_j \exp[-\ln 2\,(t_j - t_i)/h]$ and late BED $D + Q(h)/a$.
  Reported, minimum, reconstructed and assumed inter-fraction intervals are flagged separately.
  The acute repair half-time (3.2 hours) is not carried over.
- **Priors.** No prior uses an estimate derived from trials whose outcomes enter the likelihood
  (Stuschke and Thames 1999 overlaps EORTC 22791, and Bentzen 1999 overlaps CHART). Overlap is tracked
  at the trial-family level. Priors are chosen after identification (L5) from prior predictive
  checks, and every lognormal is stated with its actual quantiles.
- **Identification gate (L5).** Recovery simulations on the actual design and denominators decide
  whether α/β, repair half-time and slope can be separated. If they cannot, named repair scenarios
  are fixed or the joint fit stops. A converged sampler is not evidence of identification.
- **Integration.** A feasible fit enters the app as a research option at most, as the early fit did.
  Its +5 Gy margin keeps no clinical risk meaning unless a validated risk model gives it one.

## 4. Decision relevance, already known

At the default settings, changing late α/β alone from 3 Gy down to 1.62 Gy leaves the default list
of 69 doses unchanged, and the research-mode list of 9 doses is unchanged at every positive α/β
(review section 5, recomputed with the app engine). The value of this work lies mainly in replacing
the assumed late range on the Uncertainty page, in repair uncertainty, and in an endpoint-specific
interpretation of the late allowance, not in moving the default list.

## 5. Decisions for Igor Shuryak

- **D1.** Approve objectives A and B as separate analyses. Recommended: yes.
- **D2.** Approve the candidate tuples 1–5 as the extraction scope for L1–L2, adding others only from
  already-held reports with usable content. Recommended: yes.
- **D3.** Approve the primary-endpoint rule of section 3. Recommended: yes, including the
  evidence-synthesis outcome when no tuple qualifies.
- **D4.** Supply or obtain a full-text PDF of Beitler 2014 (RTOG 9003 final report) before tuple 2
  is extracted. Recommended: yes, through the library, because the open-access copy could not be
  downloaded automatically.

Next step after D1–D4: L1 (source-family inventory, `data/late_tox/sources.csv` and
`trial_families.csv`) on a separate branch from `v4`.
