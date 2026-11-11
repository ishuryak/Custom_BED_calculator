# Plan: early (acute mucosal) toxicity model with data-based uncertainty

Version 2, 2026-11-10. Status: PLAN ONLY. No refit has been run under this version, and nothing in
the app has changed. Applies to version 4.5.0 (revision `faf6817`, live at
https://ishuryak.shinyapps.io/custom_bed_calculator_v4/).

Version history and inputs:

- Version 1 (same date) is kept as `2026_11_10_early_toxicity_plan_pre_2026_11_10.md`. The diff is
  `2026_11_10_early_toxicity_plan_diff_2026_11_10.patch`.
- Inputs to version 2:
  - the audit `2026_11_10_early_toxicity_plan_audit.md` (findings F1–F12)
  - the independent review `2026_11_10_early_toxicity_independent_review.md` (corrections R1–R10,
    implementation steps 1–8, acceptance tests 1–18)
  - the Igor Shuryak instruction of 2026-11-10: "fix the plan, then commit and push to github, do
    not do the refit yet".
- Every numeric claim of the review that could be recomputed was recomputed independently and
  agrees to the printed digits. The script is
  `2026_11_10_early_toxicity/verification/verify_review_claims.py`. Section 10 lists what was
  and was not rechecked.

## 1. Why this matters

Under the default settings, the early-toxicity limit is the binding constraint. The defaults are
the trial-fit tumor profile, a reference of 70 Gy in 35 fractions over 47 inclusive days, and
margins of tumor −5, early +5 and late +5. Under them, the early limit sets the highest allowed
dose per fraction from 5 to 29 fractions. At 30 fractions the early and late ceilings both round
to 2.25 Gy on the 0.05-Gy grid (continuous values 2.26762 and 2.29693 Gy, so early is still
slightly tighter). The current engine lists 69 qualifying doses in 18 fraction-count rows
(`exploratory/binding.R`).

| Fractions | Inclusive days | Tumor needs d ≥ | Early allows d ≤ | Late allows d ≤ |
|---|---|---|---|---|
| 12 | 16 | 3.90 | 3.85 | 4.20 |
| 13 | 17 | 3.65 | 3.65 | 4.00 |
| 15 | 19 | 3.30 | 3.30 | 3.65 |
| 16 | 22 | 3.15 | 3.25 | 3.50 |
| 20 | 26 | 2.65 | 2.80 | 3.00 |
| 25 | 33 | 2.30 | 2.45 | 2.60 |
| 30 | 40 | 2.05 | 2.25 | 2.25 |

**Scope of this finding.** These are conditional model results. At 12 fractions the two bounds
differ by one grid step, so the 12/13 boundary is fragile. It must not be read as a biological
claim that 12 fractions are unacceptable and 13 are acceptable. Report continuous ceilings beside
grid-rounded ones near this transition. "Early toxicity binds under the calculator defaults" is
credible. "These data establish a clinical lower limit on course length" is not established.

## 2. What the app uses now, and how good the evidence is

**Model.** Early BED = n·d·(1 + d/(α/β)) − (λ/α)·max(0, T_elapsed − Tk), with the Fowler et al.
2003 values: α/β 10 Gy, α 0.35 Gy⁻¹, Tk 7 days, Tp 2.5 days. These give
r = λ/α = ln 2/(α·Tp) = 0.79217 Gy₁₀ per day. The code replaces max(0, ·) with a softplus of
sharpness 3 per day, which is at most about 0.2 Gy₁₀ away from the hinge at the kink. The main
page uses fixed values and the limit "reference early BED + 5 Gy₁₀". The reference is 53.1 Gy₁₀
and HYPNO 55 Gy/20 is 55.9 Gy₁₀. Toxicity uncertainty appears only on the Uncertainty tab.

**Evidence by component.** These are qualitative judgments of applicability and identification,
not formal certainty grades.

| Component | Source | Assessment | Main weakness |
|---|---|---|---|
| Fowler parameter set as a whole | Fowler 2003 | Reasonable reference scenario | Limited precision and transportability. The 59–61 Gy₁₀ zone belongs to that parameterization and to historical radiotherapy-alone schedules |
| α/β 10 Gy | Denham 1995: primary abstract supports about 7–10 Gy (61 palliative patients, 30–42 Gy in 10 fractions) | Order of magnitude supported | The often-quoted 9.3 (5.8–17.9) comes from secondary slides. It is not used for any prior until the full text is read |
| Tk 7 days (elapsed) | Fowler 2003, with biological support including human mucosal observations | Weakly identified | Schedule classification hardly constrains it. Pardo-Montero's 19.44-day breakpoint is a classification optimum, not a measured onset |
| r 0.79 Gy₁₀/day | Fowler 2003 (Tp 2.5 days as a selected datum) | Time effect clearly real, value imprecise | Bentzen 2001 (CHART): 0.80 (0.7–1.1) Gy/day in 2-Gy equivalents, averaged over days 12–45. It assumes recovery throughout and omits incomplete repair, and the authors estimate that repair with a 3.2 h half-time could lower it by about 25% |
| Tolerance margin | App policy "+5 Gy₁₀" | Policy, not evidence | No established patient-risk meaning |
| 13–20 fractions of 2.6–3.7 Gy | F31 (16 × 3.28 Gy tolerable, 16 × 3.44 Gy disputed, leaning intolerable), K18 (20 × 2.55 Gy) | Sparse | A few historical schedules lie in or next to the region. Modern transportability is not established |
| Clinical BED signal | Hartley 2011: 7 studies, 8 comparisons, trend P = 0.06 (with α = 0.30) | Some support | Not enough for an absolute-risk model |
| Concurrent chemotherapy | Meade 2014 offsets, HYPNO | Weak | A fixed offset does not rule out a schedule interaction |

**Uncertainty tab priors** (`R/uncertainty.R::toxicity_prior`):
- The priors are α/β lognormal (95% 5–20), α uniform 0.30–0.40, Tp uniform 2–4 days and Tk uniform
  7–19 days, all independent and labeled as assumptions.
- Only r = ln 2/(α·Tp) and Tk enter the model, so drawing α and Tp separately obscures the
  assumption.
- The ranges are not fitted, and they ignore the strong correlation between Tk and r in the data.

**Clinical benchmark.** HYPNO (Bentzen 2026) randomized 55 Gy/20 against 66 Gy/33 given six
fractions a week. Grade ≥ 3 mucositis was 190/373 against 208/380 evaluable patients. The
unadjusted difference, recomputed from these counts, is −3.80 percentage points (Wald 95% CI
−10.92 to +3.33). This is a descriptive calculation, not a trial-reported noninferiority test.
Some features limit how far it generalizes:
- Planned cisplatin cycles differed by arm.
- Acute mucositis was not a coprimary endpoint.
- The technique was not exclusively IMRT.

The Fowler values give no-repair early BEDs of:
- 53.10544 Gy₁₀ for 70/35 (46 elapsed days)
- 55.86597 Gy₁₀ for 55/20 (25 days)
- 55.43495 Gy₁₀ for the 66/33 comparator (37 days).

The randomized contrast is therefore +0.43 Gy₁₀. The +2.76 Gy₁₀ contrast against 70/35 does not
apply, because HYPNO did not randomize against 70/35.

## 3. Exploratory label-only refit: provisional, superseded

Two exploratory refits of the Pardo-Montero 2021 schedule labels were run on 2026-11-10 (scripts
in `exploratory/` and `exploratory_after_audit/`). Their posterior tables are **not carried into
this version**, for these reasons:

1. **Partial courses.** Both runs included the 8 partial-course records that Pardo-Montero
   excludes ("eight additional tolerable schedules that correspond to the first part of gap
   treatments that have not been used in the analysis"). They are spreadsheet rows 541–576
   (labels K64, K87 ×4, K88, F15, K65).
2. **No archive.** No posterior draws or run logs were archived, so the numbers cannot be
   regenerated exactly.
3. **Patient risk.** The label-only model predicts schedule labels, not patient risk.

The only conclusions kept are qualitative:
- Labels alone do not separate Tk from r.
- The fitted ceilings tend to fall below the app's ceilings.
- The label model's prediction of a positive BED difference for HYPNO's arms does not by itself
  show disagreement with HYPNO.

On the last point: whether a BED difference contradicts HYPNO depends on a patient-risk slope that
the label model does not estimate. With baseline 208/380, the same +3.9 Gy₁₀ difference gives
+1.92 percentage points at a slope scale of 50 Gy₁₀ per logit, but +9.37 points at a scale of 10.
These are illustrative values, not fitted ones.

The audit's statement that the label-only refit "contradicts" HYPNO is therefore withdrawn. It
remains true that the label-only refit must not replace the Fowler values.

**No refit is run under this version.** The next fit happens only at step 5 of section 6, after
the data and provenance gates pass.

## 4. Goal and success criteria

Goal: keep the Fowler model as a frozen reference, build provenance-checked data, and fit models
that the data can identify. Uncertainty is shown in separate layers (fitted, assumed, structural,
transportability). Absolute patient risk is shown only when a clinical context and baseline are
declared. Done when:

1. Every early-toxicity quantity used by the app traces to a fitted posterior, a cited source or a
   labeled analyst assumption (evidence ledger, constant registry).
2. The predictive comparison in step 6 is complete, whether conclusive or inconclusive. An
   inconclusive result blocks default promotion, not completion.
3. Each main-page row can show the share of joint draws meeting the early limit. Candidate and
   reference use the same draw, and failure states are kept explicit (section 5, R9 and R10).
4. A "Fowler 2003 (version 4.5)" preset reproduces the 4.5.0 list exactly (69 doses in 18 rows).

## 5. Corrections adopted from the review (R1–R10)

Each item gives what changes and its acceptance test. The test numbers refer to section 7.

- **R1. Partial courses.** Add a source-backed `course_scope` field (`complete`, `partial`,
  `unresolved`). Keep all raw records, exclude `partial` from the complete-course likelihood, and
  report each exclusion. Do not infer partial status from a low total dose alone. Selection
  fixtures, all verified 2026-11-10:

  | Selection | Records | Intolerable labels |
  |---|---:|---:|
  | Entire spreadsheet | 116 | 13 plus 2 disputed leaning intolerable |
  | Unresolved `Disp` dropped | 115 | 15 |
  | All disputed dropped, partial retained | 111 | 13 |
  | Complete, unresolved `Disp` dropped, leanings kept | 107 | 15 |
  | Complete, all disputed dropped | 103 | 13 |

  Analyses:
  - Primary historical analysis: 103 records.
  - Sensitivity analysis: 107 records.
  - The published 108-schedule selection is reproduced separately, and the 107-record analysis is
    not called a replication of it.

  Tests 1–2.
- **R2. Schedule identity and evidence identity are different things.** The audit's
  `(n, total dose, T)` key is replaced (F6 replaced). In the 115-record set, that key gives 92
  groups, while the full fraction-dose and relative-time vectors give 96. The coarse key wrongly
  merges four groups: K34/K17 (76.8 Gy/64), K73/K72 (70.4 Gy/64), K46 twice (70.2 Gy/40) and
  K88/F15 (38.4 Gy/24). Both keys give 13 event-bearing groups.
  - Create `report_id`, `cohort_id`, `trial_id`, `arm_id`, `endpoint_id` and `schedule_id`.
  - Repeated publications of one cohort and endpoint count once.
  - Independent cohorts on one schedule stay separate.
  - Conflicting reports on one cohort are adjudicated or labeled uncertain.
  - The group count is recomputed after removing partial courses. 96 is not enforced as a final
    count.

  Test 3.
- **R3. Onset is not identified by most randomized contrasts.** The version 1 sentence "The
  within-trial contrast identifies α/β, Tk and r" is wrong and is replaced by: *Within-trial
  contrasts inform identifiable combinations of fractionation and time effects. The onset
  parameter is not identified by contrasts for which both schedules end after onset. Its
  estimation requires suitable duration contrasts, longitudinal observations, or explicitly
  external assumptions.* With E = D + Q/a − r·max(0, T − Tk), the difference
  ΔE = ΔD + ΔQ/a − r·ΔT contains no Tk once both arms end after Tk, and a free trial intercept
  absorbs r·Tk.
  - Before any complex fit: evaluate the contrast design matrix and prediction Jacobian, profile
    Tk, and run recovery simulations on the real trial design.
  - If Tk is not identified, fix it within named scenarios and report no fitted Tk interval.

  Test 8.
- **R4. Tk prior.** Tk = 7 days stays in the Fowler reference (F4 revised). The audit's lognormal
  (median 8, about 4–16 days) becomes a labeled **biological sensitivity scenario**, not a fitted
  or literature-derived prior. A broad scenario uses one bound, the same in text and code (the
  version 1 text said uniform(0, 35) and the code used (0, 40), so this is reconciled to one
  value). No prior is chosen because it makes HYPNO pass.
- **R5. Schedule intolerance is not patient mucositis risk.** The two likelihoods are fitted
  separately by default (F8 and F10 extended).
  - Borrowing between them is allowed only when the source crosswalk shows no overlapping
    evidence and a no-borrowing analysis shows how far the labels move predictions.
  - Absolute risk needs a declared reference population, endpoint, assessment window, planned
    systemic therapy and baseline. Without them the app shows relative BED and model support only,
    with no Gy-to-risk translation table.

  Test 13.
- **R6. HYPNO is an external compatibility benchmark, not a pass/fail gate** (F5 replaced).
  - Removed: the 90% interval-overlap rule (excessively uncertain predictions pass it) and the
    "55/20 passes 70/35 + 5 in at least 50% of draws" rule (a policy condition against a
    comparator HYPNO did not randomize).
  - Kept: comparison of predicted counts or risk differences with the observed 190/373 against
    208/380, using HYPNO's actual arms, with baseline and sampling uncertainty.
  - A material disagreement is investigated, never removed by tuning a prior.
  - HYPNO has already shaped model criticism, so it is not called an untouched holdout.

  Test 12.
- **R7. Chemotherapy is part of the randomized strategy** (F9 revised). Extract planned and
  delivered chemotherapy separately.
  - The primary analysis targets the assigned strategy, with planned therapy in its definition.
  - Delivered dose, interruptions and actual overall time are described, not adjusted for
    automatically, because toxicity can affect them.
  - An arm is not rejected because delivered cisplatin is missing. The missingness is recorded.
- **R8. One schedule evaluator for fitting and prediction.** The exploratory fitter used repair
  with a 3.2 h half-time, while `implic.py` used complete repair. The new evaluator:
  - uses exact fraction times in hours, with the calendar-day convention stored separately
  - uses the repair term Q = Σdᵢ² + 2Σ_{i<j} dᵢdⱼ·exp(−μ(tⱼ − tᵢ)), μ = ln 2/h, with h = 3.2 h
    for reproduction and 2–4 h as a sensitivity range, not a fitted distribution
  - keeps Denham's repair estimates out of a pooled prior, since they come from a different
    irradiation design.

  The audit's planned test "below 0.2 Gy₁₀ for once-daily schedules" is not universal (F12c
  revised): 70/35 gives 84.12426 against 84.00000 Gy₁₀, but a stress fixture of five 8-Gy
  fractions 24 h apart adds 0.28402 Gy₁₀. When α/β varies, results are labeled "BED at the
  sampled α/β", and a CHART-type rate c in 2-Gy equivalents converts as r = c·(1 + 2 Gy/a). A
  separate legacy evaluator reproduces the current app. Tests 5–7.
- **R9. What an uncertainty share means.** For a fixed candidate the support is
  S_E = (1/M)·Σ_b I{E_cand,b − E_ref,b ≤ m}, with candidate and reference evaluated in the same
  joint draw b. It is support for a model inequality, not a probability of clinical safety.
  - A fixed-dose support is distinct from the probability that some dose in a row qualifies.
  - Tumor and early fits are paired by randomized index with a fixed seed, as an explicitly
    assumed product distribution, and the dependence-free bounds
    max(0, S_T + S_E − 1) ≤ S_joint ≤ min(S_T, S_E) are reported.
  - The draw export includes everything the prediction needs (baseline or context and risk slope
    when risk is shown), not only biological parameters.

  Tests 9, 14.
- **R10. Failure states and reproducible outputs.**
  - No numeric sentinels. The 2.00 Gy placeholder in `exploratory_after_audit/implic.py`
    lines 6–9 is not reused.
  - Store a no-feasible-dose flag, its probability, and conditional quantiles among feasible
    draws.
  - Reconcile the reference-validity rule: the deterministic search rejects a nonpositive
    reference BED (`R/bed_models.R:263`), while the uncertainty path is reported to check the
    candidate. Keep separate flags for computational failure, invalid model domain and "no
    qualifying dose".
  - Archive chains, warmup, draws, diagnostics, data hashes, package versions and configuration.
  - Use median predictions or a representative joint draw, never a point model built from
    marginal medians.

  Tests 11, 15.

## 6. Implementation sequence

Each step has a completion gate, and the next step starts only when it passes.

**Step 1. Freeze scope and define the clinical question (Claude, 2 hours).**
- Keep the plan, audit, review and both exploratory folders unchanged, and keep the Fowler
  evaluation as a frozen fixture.
- Write a structured target definition: disease site, definitive or postoperative treatment,
  planned systemic therapy, mucosal endpoint, assessment window, scoring system and reference
  regimen.
- Treat early glottic radiotherapy, broad-field definitive chemoradiotherapy and reduced-volume
  favorable-risk oropharyngeal treatment as separate contexts.
- Gate: every planned prediction has a stated context and comparator.

**Step 2. Data and provenance tables (Claude extracts, Igor spot-checks, 3–5 days).**
- Create these files under `data/early_tox/`:
  - `reports.csv`
  - `schedules.csv` (with `course_scope`, intended or delivered status, source clock, normalized
    elapsed clock)
  - `fractions.csv` (times in hours)
  - `labels.csv` (original and adjudicated labels, inclusion rule and reason)
  - `arms.csv` (assigned regimen, context, technique, planned and delivered chemotherapy,
    dose-volume information, missingness)
  - `outcomes.csv` (exact grading definition and version, threshold, assessment start, end and
    frequency, evaluable n, events, missing outcomes)
  - `evidence_families.csv`
- Rules:
  - No patient event counts are inferred from schedule labels.
  - A percentage is not rounded into a count unless the count is unique.
  - A missing outcome is not a zero.
- Extraction:
  - Every historical record is read against its source.
  - New arm counts are extracted twice independently and reconciled.
  - A search and screening log is kept (date, query, inclusion reason, missing full text).
  - The 19-trial tumor list is only a starting index. KROG 0201, PMH and Yamazaki are already in
    it (F7).
  - Candidate additional data sets, all to be verified before use: the eight trials listed by
    Lester in HYPNO's references, Madhava (55/20 with carboplatin), Sanghera 2007, Jegannathen
    and Bratman phase 2.
  - Every value comes from a downloaded full text, and missing PDFs are reported by name.
- Gate: the R1 fixtures reproduce, every denominator is explicit, the primary analysis has no
  unresolved duplicate evidence, and no evidence family crosses a validation split.

**Step 3. Harmonize endpoints (Claude, 1 day, Igor approves the mapping).**
- These stay separate endpoints:
  - prolonged intolerability
  - maximum grade ≥ 3 mucositis
  - duration of confluent mucositis
  - feeding tube use
  - unplanned breaks
  - dose-limiting toxicity
- Record the exact CTCAE, RTOG, WHO or other scale version, and write a mapping with a primary
  compatible group, secondary groups and exclusions.
- The observation window must be able to catch peaks after a short course. Otherwise the study
  is a limited-ascertainment sensitivity analysis or a qualitative check.
- HYHOPE's persistent-toxicity dose-limiting endpoint is not treated as maximum grade 3
  mucositis, and its abandoned 44.4/12 level is not an observed intolerable cohort.
- Gate: every pooled event has a defensible endpoint mapping and window.

**Step 4. Identification and prior predictive work, before any outcome-based fit (Claude,
1 day).**
- Contrast design matrix, Jacobian, Tk profile and recovery simulations on the real design (R3).
- Prior predictive checks for the clinical benchmark priors below.
- Gate: identification results recorded per parameter.

**Step 5. Fit, in this order (Claude, 1–2 days, the first refit run of this plan).**
1. **Clinical benchmark** (one compatible endpoint and context group): y_ij ~ Binomial(N_ij, p_ij),
   with logit p_ij = η_j + κ·(E_ij − E_j,ref).
   - E is computed from the fixed Fowler set.
   - η_j is the reference-arm log odds of trial family j: η_j ~ Normal(m₀, τ₀²),
     m₀ ~ Normal(0, 2.5²), τ₀ ~ HalfNormal(1).
   - The slope prior is set on L* = κ·10 BED units: HalfNormal(1), with HalfNormal(0.5) and
     HalfNormal(2) as sensitivities.
   - These are analyst-proposed regularization choices, recorded before fitting.
   - Arm-level overdispersion is assessed.
2. **Arm-only candidate:** the same baseline and slope structure, with only the biological
   combinations supported by step 4 allowed to vary.
   - Exposure is centered at each trial's reference arm.
   - Tk is fixed within named scenarios if it is not identified.
   - Starting sensitivity priors are α/β lognormal(log 10, 0.354) and r lognormal(log 0.8, 0.5),
     labeled analyst-chosen.
3. **Separate historical classifier:** 103 complete-course records (primary) and 107
   (sensitivity), with the source and cohort structure.
   - Compare Fowler with an alternative time-boundary model.
   - It is not used to predict patient risk.

General rules:
- Bentzen 2001 is a consistency check, not an independent validation and not a prior, when CHART
  outcomes are in the data (F3 retained and extended).
- No chemotherapy-by-time or technique-by-dose interaction unless the contrast design can inform
  it.
- Long fits follow the analysis-script convention: seeds, bounded workers, thread limits,
  checkpoints, logs and a resumable orchestrator.

Gate: no divergences, rank-normalized R-hat < 1.01, and bulk and tail ESS > 400 for parameters
and decision quantities. Simulation-based calibration is registered in `sim_provenance.yaml`.

**Step 6. Validation without leakage (Claude, 1 day).**
- Hold out whole evidence families, including their labels, derived summaries and any priors
  estimated from them (CHART and Bentzen's derived estimate included).
- Report two tasks separately: new-trial prediction, and conditional contrast prediction (only the
  held-out comparator arm informs its baseline).
- Report paired family-level log-score differences between candidate and benchmark, calibration,
  observed against predicted contrasts and interval widths.
- Coverage is descriptive only. The 80–97% band is dropped (F11 remedy revised).
- Promotion convention:
  - The lower 95% bound of the paired advantage, from a Bayesian bootstrap of family weights
    (Dirichlet(1,…,1), recorded seed, repeated with a second seed), must be above zero, with no
    unresolved endpoint or context mismatch.
  - A leave-one-family influence analysis is required.
  - If the comparison is unstable, it is called inconclusive: Fowler stays the default and the
    candidate becomes a labeled research scenario.
- HYPNO is the R6 benchmark, and HYHOPE and Jacinto are contextual checks with their own volume,
  therapy and endpoint definitions.
- Gate: no leakage, the same task for all models, and favorable and unfavorable benchmarks both
  reported.

**Step 7. Decision display (Claude, 2 days, Igor decides policy).**
- Keep "+5 Gy₁₀" as an explicit user-selected margin. When α/β varies, state whether the margin
  is a numerical BED increment or defined in EQD2 and converted per draw.
- A patient-risk rule is added only with a calibrated reference context and a clinically chosen
  tolerance δ. Neither δ = 5 points nor an 80% support threshold is an empirical safety cutoff.
  Without δ the risk-difference distribution is shown without accept or reject labels.
- Export joint draws with source hashes, model identity, clock and repair conventions and
  endpoint. Mismatched manifests are refused.
- Then: `toxicity_prior()` uses r directly (α and Tp are no longer drawn separately), the
  constant registry is updated (`stale_constant_audit.py --strict`), and the About page, user
  guide and README are updated in the same pass.
- Gate: the legacy preset reproduces 4.5.0 exactly.

**Step 8. Acceptance suite, release (Claude, 1 day, Igor approves deployment).** Tests 1–18 of
section 7.
- Each change gets one committed test, seen red without the change and green with it.
- Then the combined smoke run and the full suite (`tools/run_tests.R`), followed by the local
  walkthrough.
- The `docs/FIXES.md` rows state on their first line whether the main list moves.
- `tools/deploy.R` runs only after Igor approves.

## 7. Acceptance tests (adopted from the review, section 4 step 8)

1. Raw label counts and the R1 selection counts reproduce. Unknown labels fail explicitly, and
   original labels are kept.
2. The 8 partial-course records never enter the complete-course likelihood.
3. Fraction times separate the four false merges of R2. Independent cohorts stay distinct, and
   duplicate publications do not.
4. No evidence-family identifier appears in both training and held-out data, including indirect
   prior sources.
5. Dose and time masks agree, fraction counts and totals reconcile, times are ordered, and unit
   conversions are explicit.
6. Fowler no-repair values reproduce 53.10544, 55.86597 and 55.43495 Gy₁₀, and the section 1 grid
   fixtures are kept.
7. The 84.12426 Gy₁₀ repair fixture reproduces. Infinite separation tends to complete repair, and
   coincident fractions act as one combined dose.
8. Exact Tk cancellation for the hinge when both courses pass onset, and no cancellation when one
   course ends before onset. The softplus approximation tolerance is documented, with its
   sharpness in inverse days.
9. Candidate equal to reference gives zero difference in every draw and support 1 for m ≥ 0.
10. Exposure at fixed n increases with dose per fraction under positive parameters. Monotonicity
    in time is checked per model, not assumed.
11. No-feasible-dose states are kept without numeric sentinels.
12. HYPNO fractions and the descriptive interval reproduce. A validation call that substitutes
    70/35 for the randomized comparator is refused.
13. Patient-risk output is refused without endpoint, context and baseline draws. A label
    probability never enters that function.
14. Risk and support summaries reproduce in an independent implementation. The joint-support
    bounds hold and posterior correlations are kept.
15. Sampler diagnostics as in step 5.
16. Syntax checks, end-to-end smoke fit and frozen validation report.
17. Every app page uses the intended evaluator and draw manifest. Labels, user guide and source
    descriptions agree with the endpoint and units.
18. The existing app suite and a local walkthrough pass. Deployment is a separate approval.

**Stopping rules.**
- A data or provenance defect blocks fitting.
- Nonidentification blocks an empirical claim for that parameter but allows a labeled sensitivity
  model.
- An inconclusive predictive comparison blocks default promotion, not completion.
- Missing clinical context blocks patient-risk output.
- When all checks pass, the work proceeds without another broad review. Reopen only for a failed
  gate or new evidence.

## 8. Uncertainty layers shown to the user

| Layer | Treatment | Not claimed |
|---|---|---|
| Parameter estimation | Joint posterior draws, candidate and reference in the same draw | Independent uniform ranges are fitted uncertainty |
| Weak identification | Prior sensitivity, profiles, recovery results | A narrow posterior proves identification |
| Between-study prediction | New-context baseline uncertainty | Pooled-mean uncertainty equals a new clinic's |
| Structural | Named onset, repair, time-function and borrowing variants, reported as a sensitivity envelope | An arbitrary mixture has a calibrated 95% interval |
| Exposure | Recorded timing assumptions, intended against delivered, mucosal-dose limits | Prescription dose is the dose to all relevant mucosa |
| Endpoint | Separate scales and windows, or justified mapping sensitivities | Trial intercepts make different outcomes identical |
| Patient variability | Patient-risk model when suitable data exist | Parameter draws represent individual patients |

The first implementation includes the Fowler reference, a broad biological sensitivity, an
early-onset sensitivity and repair variants. "Extrapolation" is flagged separately from support
when data do not cover the fraction size, duration, volume or systemic therapy. Supported
schedules are never labeled clinically safe.

## 9. Decisions for Igor (recommendation first)

Settled by this revision, unless Igor objects:
- Fowler stays the default until step 6 shows a conclusive advantage.
- The live app is unchanged until then.
- Bentzen 2001 is a consistency check only.
- HYPNO is a benchmark, not a gate.
- The two likelihoods are fitted separately by default.

**Status 2026-11-10: D1–D5 decided as recommended below, D6 replaced, and D7 added** (Igor
Shuryak: "all as recommended"). Denham 1995 could not be obtained, so D6 now reads: cite Denham's
interval as reported by Kaanders 1999 (secondary), keep Fowler's 10 Gy as the center, and check in
step 4 whether the hyperfractionation trials inform α/β. D7: follow the gates, with no exploratory
pilot fit. Record: `2026_11_10_early_toxicity_step1_target_definition.md`. The original wording
of D1–D6 follows unchanged.

Open, needed before the step shown:

- **D1 (before step 3). Primary endpoint group.** Maximum grade ≥ 3 acute mucositis with the scale
  version recorded and mapped. Recommended. Alternative: confluent mucositis only (fewer trials).
- **D2 (before step 1). Primary clinical context.** Broad-field definitive chemoradiotherapy for
  locally advanced HNSCC. Recommended, because it matches the app's reference and HYPNO.
- **D3 (before step 5). Prior scales.** Adopt the review's analyst-proposed baseline and slope
  priors (section 6 step 5) as the starting specification. Recommended.
- **D4 (before step 7). Margin definition when α/β varies.** Keep +5 as a numerical BED increment
  at the sampled α/β (the simplest continuation of 4.5.0). Recommended, with the EQD2 version as a
  sensitivity.
- **D5 (before step 7). Clinical tolerance δ for a risk rule.** No default. If none is chosen, show
  the risk-difference distribution without labels. Recommended.
- **D6 (before step 2). Denham 1995 full text.** Obtain it. Until then, no prior uses its interval.
  Recommended.

## 10. Verification record for this revision

Recomputed on 2026-11-10 with `2026_11_10_early_toxicity/verification/verify_review_claims.py`
(argument: the Pardo-Montero supporting `.xls`), all matching the review:
- the 116 records, with partial-course rows 541–576 and their labels
- the five selection counts
- 92 coarse against 96 full-vector groups, the four false merges, and 13 event-bearing groups
- the Fowler BEDs 53.10544, 55.86597 and 55.43495, and the difference 0.43102
- the repair fixture 84.12426 (schedule F15 in the spreadsheet) and the stress fixture 0.28402
- the HYPNO −3.80 (−10.92 to +3.33)
- the slope counterexample +1.92 and +9.37 points.

Also confirmed:
- The Pardo-Montero text excludes the 8 first-part-of-gap schedules.
- `exploratory_after_audit/implic.py` uses the 2.00 Gy sentinel.
- `R/bed_models.R:263` rejects a nonpositive reference BED.

Not rechecked here:
- the review's claim about the uncertainty path's candidate check
- its clinical source summaries for Hartley, Denham, HYHOPE and Jacinto (taken as reported, and
  to be read in full at step 2)
- the binding table, which was reproduced earlier by `binding.R` and by the review.
