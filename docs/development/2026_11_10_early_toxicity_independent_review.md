# Early toxicity BED evidence and implementation review

Igor Shuryak, MD, PhD — is144@cumc.columbia.edu

November 10, 2026

Internal development review for `HNSCC_BED_calculator/Custom_BED_calculator_v45`

## 1 Decision and practical meaning

**The objective is worthwhile, but the plan needs substantial revision before fitting a replacement model. Several recommendations in the existing audit also need revision.** Early mucosal biologically effective dose (BED) is the binding constraint in the calculator under the stated defaults. This does not establish that it is the principal clinical limitation of every shorter regimen, or that the displayed fraction-dose ceiling is accurately known.

The evidence supports a clinically important effect of treatment duration on acute mucosal reactions. It supports the Fowler parameter set as a useful historical reference. It is substantially weaker for a universal onset time, a precise proliferation rate under different repair assumptions, and an absolute mucositis-risk function for contemporary hypofractionated chemoradiotherapy. The current “reference BED + 5” margin is a calculator policy choice without an established patient-risk interpretation.

The recommended development path is to retain the reproducible Fowler reference, correct and document the historical data, and fit patient-level event probabilities using aggregate trial-arm counts with compatible endpoints. Keep historical schedule tolerability separate initially. Add uncertainty in parameters, model form, and clinical context. Do not replace the current values with the corrected label-only posterior or the audit's narrower onset prior.

This review provides one implementation specification and conditional stopping rules. Successful implementation does not require repeated broad reviews. A failed check should trigger the specific correction listed below. Clinical risk tolerances and deployment remain explicit decisions, rather than quantities that can be estimated from these data.

### What was independently reproduced

The binding calculation was rerun against the local application functions. At the specified trial-fit tumor profile, 70 Gy in 35 fractions over 47 inclusive days, and tumor/early/late margins of −5/+5/+5, the early constraint gives the lower toxicity ceiling from 5 through 29 fractions. At 30 fractions, early and late ceilings both round to 2.25 Gy on the 0.05-Gy search grid. These are conditional model results.

| Fractions | Inclusive days | Minimum tumor-compatible dose per fraction | Early ceiling | Late ceiling |
|---|---:|---:|---:|---:|
| 12 | 16 | 3.90 Gy | 3.85 Gy | 4.20 Gy |
| 13 | 17 | 3.65 Gy | 3.65 Gy | 4.00 Gy |
| 15 | 19 | 3.30 Gy | 3.30 Gy | 3.65 Gy |
| 16 | 22 | 3.15 Gy | 3.25 Gy | 3.50 Gy |
| 20 | 26 | 2.65 Gy | 2.80 Gy | 3.00 Gy |
| 25 | 33 | 2.30 Gy | 2.45 Gy | 2.60 Gy |
| 30 | 40 | 2.05 Gy | 2.25 Gy | 2.25 Gy |

The apparent boundary between 12 and 13 fractions is particularly fragile: the 12-fraction bounds differ by one search-grid increment. It must not become a biological claim that 12 fractions are unacceptable and 13 are acceptable. Report continuous ceilings alongside grid-rounded values when evaluating this transition. At 30 fractions, the continuous early and late ceilings are 2.26762 and 2.29693 Gy, respectively. Early remains tighter, although both round down to 2.25 Gy on the search grid. The current main-page engine reproduces 69 qualifying doses across 18 fraction-count rows for this configuration.

## 2 How strong is the evidence

The assessments below are qualitative judgments about applicability and identification, not formal evidence-certainty grades.

### Fowler parameters and the tolerance boundary

Fowler and colleagues proposed α/β = 10 Gy, α = 0.35 Gy⁻¹, onset Tₖ = 7 days, and proliferation doubling time Tₚ = 2.5 days. The resulting rate is r = ln(2)/(αTₚ) = 0.79217 BED units/day. The 7-day choice has biological support, including human mucosal observations, but was weakly identified by schedule classification. The doubling time was a selected modeling datum. The reported 59–61 Gy₁₀ tolerance zone applies to that parameterization and historical radiotherapy-alone schedules. It is neither a patient-risk confidence interval nor a threshold transferable unchanged to another model. [1]

**Assessment:** reasonable reference scenario, limited precision, and limited transportability. Preserve the entire parameter set as a named legacy model. Do not mix its threshold with a newly fitted onset, repair model, or proliferation rate.

### Fractionation sensitivity

Denham and colleagues studied acute oropharyngeal reactions in 61 palliative patients receiving 30–42 Gy in 10 fractions across a wide range of dose rates. The primary abstract supports α/β values around 7–10 Gy and explicitly questions transport of the repair estimates to high-dose-rate multifraction schedules. [2] The often-repeated 9.3 Gy estimate with limits 5.8–17.9 should not be treated as a modern peak grade ≥3 mucositis prior without checking its precise endpoint and statistical model in the full paper.

**Assessment:** useful support for the order of magnitude, weaker support for a single precise distribution in the intended application. The current lognormal distribution centered at 10 Gy with an approximate 5–20 Gy interval is a defensible broad sensitivity assumption, but is not the published sampling distribution of this study. The primary Denham PDF was not available in the inspected local reference collection. Its precise interval must remain excluded from quantitative prior construction until that source is obtained.

### Proliferation rate and repair

The continuous hyperfractionated accelerated radiotherapy (CHART) trial analysis estimated an effective recovery of 0.80 Gy/day in equivalent 2-Gy fractions, with 95% confidence limits of 0.7–1.1. However, this is an average over the interval between the trial's 12-day and 45-day schedules. It assumes recovery throughout that interval and does not identify onset. The main estimate omitted incomplete repair. The authors estimated that allowing a 3.2-hour repair half-time could reduce it by roughly 25%, and described that repair half-time as an informal estimate. [3]

**Assessment:** substantial evidence that time matters, but weak justification for an independent narrow prior on r in a model with freely varying onset and explicit repair. Existing audit F3 is correct on this point. Its “validation only” recommendation also needs a qualification: if CHART outcomes enter fitting, the CHART-derived rate is a consistency check, not independent validation.

### Historical schedule labels

Fenwick's analysis and the later Pardo-Montero study provide valuable schedule-level information about prolonged or intolerable reactions. They do not provide a representative denominator for individual patient risk. The Pardo-Montero spreadsheet contains 116 entries, but the paper explicitly excludes eight partial-course observations from its 108-schedule analysis. Its fitted Fowler-type breakpoint of 19.44 days is a classification optimum under that analysis, not a measured human onset distribution. [4,5]

**Assessment:** useful for checking historical boundaries and alternative time functions. Limited for calibrating patient probabilities or resolving modern dose-volume and chemotherapy effects. A good sampler cannot repair an inappropriate observation unit or a mismatched endpoint.

### Clinical checks of the BED signal

Hartley and colleagues examined seven studies with eight comparisons. The reported association between acute mucosal BED differences and grade 3 mucositis differences was a trend at P = 0.06. Their acute parameterization used α = 0.30 Gy⁻¹, rather than Fowler's 0.35 Gy⁻¹. [6]

**Assessment:** some empirical support for a dose-time signal, insufficient evidence for a calibrated absolute-risk model. The modest evidence base is a reason to compare simple models carefully, not a reason to force a more complex biological fit.

### HYPNO and contemporary treatment

HYPNO provides the strongest directly relevant randomized comparison: 55 Gy in 20 fractions versus 66 Gy in 33 fractions delivered six fractions per week. Grade ≥3 mucositis occurred in 190/373 versus 208/380 evaluable patients. The independently recomputed unadjusted difference is −3.80 percentage points, with an approximate Wald 95% confidence interval of −10.92 to +3.33 points. This is a descriptive calculation from published evaluable counts, not a trial-reported acute-mucositis noninferiority test. [7]

The trial compared assigned treatment strategies, including different numbers of planned weekly cisplatin cycles. Acute mucositis was not its coprimary noninferiority endpoint. Its technique mix was not exclusively intensity-modulated radiation therapy (IMRT). Therefore, HYPNO is informative about this specific strategy comparison, not a universal guarantee for 55/20 or an equivalence test against 70/35.

Meade's chemotherapy-adjusted tolerance calculations are useful historical comparisons, but a fixed chemotherapy BED offset does not establish absence of schedule interactions. [8] HYHOPE supports feasibility of 46.5 Gy in 15 fractions in a selected human papillomavirus-associated oropharyngeal population with reduced elective volumes. Its abandoned 44.4/12 level should not be entered as an observed intolerable HYHOPE cohort. [9] Jacinto's 55/20 study also used lower elective doses than the gross-tumor prescription. [10]

**Assessment:** these studies improve applicability checks, but pooling them as if they used the same mucosal exposure and outcome definition would create false precision.

## 3 Corrections required before a replacement fit

### R1 Remove partial-course observations from the complete-course likelihood

**Location:** both `exploratory/parse.py` and `exploratory_after_audit/parse.py`, and corrected `fit.py` lines 6–7. The parser reads all spreadsheet records. The corrected fitter removes the unresolved `Disp` category but does not remove partial courses.

The final eight spreadsheet schedules, numbered 109–116 in extraction order, are first portions of gap treatments. Their source labels are K64, K87, K87, K87, K87, K88, F15, and K65. Their reference rows in the original spreadsheet are 541, 546, 551, 556, 561, 566, 571, and 576. Pardo-Montero excludes these from its complete-course analysis [5]. Treating them as independent complete tolerated regimens adds inappropriate negative examples, particularly at short durations.

Implement a source-backed `course_scope` field with values `complete`, `partial`, or `unresolved`. Preserve all raw records. Exclude `partial` from the complete-course likelihood and report each exclusion. Do not infer partial status only from a low total dose. If eventually used in a longitudinal model, link a partial observation to its parent course and observed time window.

The independently verified count fixtures are:

| Selection before publication or cohort reconciliation | Records | Intolerable labels |
|---|---:|---:|
| Entire spreadsheet | 116 | 13 unequivocal plus 2 disputed leaning intolerable |
| Current corrected script, unresolved `Disp` dropped | 115 | 15 |
| All disputed categories dropped, partial courses retained | 111 | 13 |
| Complete courses, unresolved `Disp` dropped, leanings retained | 107 | 15 |
| Complete courses, all disputed categories dropped | 103 | 13 |

These are parser fixtures, not final independent sample sizes. For a conservative primary historical-label analysis, use the 103 unambiguous complete-course records before cohort reconciliation. Use the 107-record leaning rule as a sensitivity analysis. Separately reproduce the published 108-schedule selection and adjudications if comparing directly with the published classifier. Do not call the 107-record analysis an exact replication of that publication.

### R2 Separate publication duplication from schedule identity

**Location:** plan Phase A1 and audit F6. The proposed `(n, total dose, T)` uniqueness rule is insufficient and its stated independence interpretation is incorrect.

In the current 115-record set, that coarse key produces 92 schedule groups. A key containing the full fraction-dose vector and relative fraction-time vector produces 96 groups. Four coarse groups merge genuinely different timing patterns: K34/K17 at 76.8 Gy in 64 fractions, K73/K72 at 70.4 Gy in 64 fractions, two K46 records at 70.2 Gy in 40 fractions, and K88/F15 at 38.4 Gy in 24 fractions. Consequently, the audit's approximate count of 11 unique event-bearing schedules is not supported by its key. The full-vector check gives 13 event-bearing groups in that same 115-record set.

Create separate `report_id`, `cohort_id`, `trial_id`, `arm_id`, `endpoint_id`, and `schedule_id`. Repeated publications about the same cohort and endpoint must not contribute independent observations. Independent cohorts using the same schedule may legitimately remain separate. Conflicting reports about the same cohort require adjudication or explicit uncertain labeling. A shared-schedule random effect does not by itself fix duplicate publications.

**Acceptance:** the source crosswalk accounts for every retained observation and every exclusion. No training/validation split separates reports from the same evidence family. Recompute group counts after removing partial courses. Do not enforce 96 as the final cohort count.

### R3 Correct the claim that randomized contrasts identify onset

**Location:** plan Phase B, lines 138–146, especially “The within-trial contrast identifies α/β, Tk and r.”

Let a = α/β in Gy, D be total dose in Gy, Q be the quadratic dose term in Gy² including repair, T be elapsed days, and r be BED units/day. The proposed effective exposure is

$$E=D+Q/a-r\max(0,T-T_k).$$

When both schedules end after Tₖ, their difference is

$$\Delta E=\Delta D+\Delta Q/a-r\Delta T.$$

Tₖ cancels exactly. An unrestricted trial intercept also absorbs the common rTₖ contribution to the absolute exposure. Thus many conventional and accelerated trial contrasts cannot identify onset. Strong partial pooling or a narrow prior can make the posterior look precise without supplying the missing randomized information.

**Replacement wording:** “Within-trial contrasts inform identifiable combinations of fractionation and time effects. The onset parameter is not identified by contrasts for which both schedules end after onset. Its estimation requires suitable duration contrasts, longitudinal observations, or explicitly external assumptions.”

**Acceptance:** before fitting the complex model, evaluate the contrast design matrix and prediction Jacobian, profile the onset parameter, and perform recovery simulations using the actual trial design. Demonstrate separately whether a, r, and Tₖ are identified. If onset is not identified, fix it within named structural scenarios and do not report a data-estimated onset interval.

### R4 Do not promote the proposed narrow onset distribution to evidence

**Location:** audit F4 and its implementation step 5.

The suggested lognormal median of 8 days and approximate 4–16-day interval is an analyst-selected distribution. Human observations can support early regeneration without estimating that distribution for the effective clinical breakpoint. Mouse recovery experiments do not supply a human confidence interval. Biological proliferation onset and the best breakpoint of an aggregate exposure model need not coincide.

Keep Tₖ = 7 days as part of the Fowler reference. The proposed 8-day distribution may be a labeled biological sensitivity scenario. A broad onset scenario may also be used, with one explicitly chosen bound across text and code. Neither should be described as a fitted human prior. Do not select a prior because it makes HYPNO pass.

### R5 Separate schedule intolerance from patient mucositis risk

**Location:** plan Sections 3–5, both likelihoods, and audit F8/F10.

A probability that a curated review labels a schedule intolerable is not a patient's probability of grade ≥3 mucositis. Different slopes for the two likelihoods do not establish that they share identical biological parameters. Endpoint differences, scoring, ascertainment, and treatment volume can alter the apparent exposure relationship.

Use compatible trial-arm patient counts for the primary risk model. Fit the historical-label model separately. Permit borrowing only after the source crosswalk excludes overlapping evidence and a no-borrowing analysis shows how much the labels move predictions. A separate historical model is the default, not a failure of the development plan.

For risk predictions, define the reference population, endpoint, assessment window, planned systemic therapy, and baseline risk. An unknown trial intercept is not a usable baseline for a generic calculator user. If these cannot be specified, display relative BED and model support only. Do not display absolute risk or a universal Gy-to-risk translation table.

### R6 Replace the HYPNO gate

**Location:** audit F5 and implementation step 6.

Remove the requirements that a predicted interval merely overlap the observed interval and that 55/20 pass the 70/35-plus-5 rule in at least half the draws. The latter is an arbitrary calculator-policy condition against a comparator that HYPNO did not randomize. The former can be satisfied by predictions that are excessively uncertain.

With no repair and the Fowler parameters, independently calculated early BED values are 53.10544 Gy₁₀ for 70/35, 55.86597 Gy₁₀ for 55/20, and 55.43495 Gy₁₀ for the six-per-week 66/33 comparator. These use 46, 25, and 37 elapsed days, respectively. The relevant randomized difference is +0.43102 Gy₁₀, not the +2.76053 Gy₁₀ difference against 70/35.

A positive predicted BED difference alone does not establish disagreement with HYPNO. Its interpretation requires a calibrated patient-risk slope. As a mathematical counterexample, at baseline risk 208/380, a +3.9 BED-unit difference produces about +1.92 percentage points with slope scale 50, but +9.37 points with slope scale 10. These are illustrative assumptions, not fitted clinical values. They show why the label-only posterior cannot be declared incompatible with patient outcomes from the sign of its BED difference alone.

HYPNO has already influenced model criticism and prior proposals. Designate it an external compatibility benchmark, not an untouched confirmatory holdout. Compare predictions with its actual randomized treatment strategies and endpoints. Reserve a genuinely new dataset for future independent confirmation when one becomes available.

### R7 Preserve the randomized treatment strategy when handling chemotherapy

**Location:** audit F9 and plan decision D4.

Extract planned and delivered chemotherapy separately. In HYPNO, the planned numbers of weekly cycles differ by arm, and omission can follow preceding toxicity [7]. Delivered dose is therefore potentially affected by treatment assignment and toxicity. Adjusting for it as an ordinary baseline covariate can remove part of the treatment effect or introduce selection bias.

The primary randomized analysis should target the assigned treatment strategy. Use planned systemic therapy as part of that definition. Describe delivered dose and interruptions, but do not automatically adjust for them. An analysis intended to isolate radiation's direct effect requires a different causal estimand and richer data than aggregate arm averages.

Do not reject an otherwise usable randomized arm merely because delivered cisplatin dose is missing. Record missingness and its consequence. The same caution applies to actual overall treatment time if toxicity caused treatment breaks: a prolonged delivered course is not an unconfounded experiment in protection from mucositis.

### R8 Make repair and dose units coherent across fitting and prediction

**Location:** corrected `fit.py` lines 8–11 and `implic.py` line 4, plan Phase E/F.

The exploratory fit includes repair with a fixed 3.2-hour half-time, but the implication script uses a complete-repair dose expression. The discrepancy is small for ordinary once-daily treatment, but still breaks strict parity. With the spreadsheet's 70/35 schedule and 3.2-hour half-time, the independently calculated cell-kill BED is 84.12426 Gy₁₀, compared with 84.00000 under complete repair. The audit's proposed “below 0.2 Gy₁₀” once-daily test is not universal. A mathematical stress fixture of five 8-Gy fractions exactly 24 hours apart adds 0.28402 Gy₁₀ at the same half-time. This is a numerical test case, not a recommended clinical regimen.

Use one schedule evaluator for fitting and prediction, plus a deliberately separate legacy evaluator that reproduces the current app. Include exact fraction times in hours. Store the calendar-day convention separately. Do not accidentally subtract one day twice or silently replace fractional-day timing with integer duration.

For instantaneous fractions at times tᵢ, with doses dᵢ, the repair term is

$$Q=\sum_i d_i^2+2\sum_{i<j}d_i d_j\exp[-\mu(t_j-t_i)],\qquad \mu=\ln(2)/h.$$

Here h is repair half-time in hours and μ is the repair rate in h⁻¹. This expression is the model used in the corrected exploratory fitter. Retain 3.2 hours for reproduction and examine source-supported alternative repair assumptions, including the 2–4-hour range discussed by Bentzen [3]. Treat that range as a sensitivity range, not a fitted probability distribution. Denham's shorter repair estimates arise from a different irradiation design [2]. Do not average them into a single pooled prior without a transport argument.

If α/β varies, label the result BED at the sampled α/β, rather than Gy₁₀. For a recovery rate c expressed in equivalent 2-Gy fractions per day, the conversion is

$$r=c(1+2\,\mathrm{Gy}/a).$$

The factor 1.2 applies only at a = 10 Gy. This conversion alone does not fix the onset and repair dependencies of the CHART estimate.

### R9 Specify what an uncertainty percentage means

**Location:** plan Phase D/E and audit Section 4.

For each fixed candidate regimen, use the same joint parameter draw for both candidate and reference. For M draws and a declared BED margin m, compute

$$S_E=\frac{1}{M}\sum_{b=1}^{M} I\{E_{\mathrm{candidate},b}-E_{\mathrm{reference},b}\le m\}.$$

This is support for a model inequality from equally weighted posterior draws, conditional on the specified model and priors. It is not the probability of clinical safety. If draws receive explicit weights, replace the unweighted average by the sum of weighted indicators divided by the sum of weights, and declare the weights. Distinguish a fixed-dose support value from the probability that some dose exists within a fraction-number row. They answer different questions.

When patient risk is calibrated, one possible conditional transformation is

$$p_{\mathrm{candidate},b}=\operatorname{logit}^{-1}\{\operatorname{logit}(p_{\mathrm{reference},b})+\kappa_b\Delta E_b\},$$

where κ is the positive risk slope per BED unit and the reference probability is for the declared clinical context. Export their joint draws and context parameters. The original proposed export of only biological parameters and label thresholds is insufficient for this calculation.

Pairing separate tumor and early-toxicity draws approximates an explicitly assumed product distribution. Separate fits do not establish biological independence. Use randomized pairing with a fixed seed, preserve dependence within each fit, and report the assumption. For two criteria with supports Sₜ and Sₑ, the dependence-free limits are

$$\max(0,S_T+S_E-1)\le S_{\mathrm{joint}}\le\min(S_T,S_E).$$

If conclusions change substantially across plausible dependence assumptions, do not present one combined percentage as settled evidence. Keep late-toxicity and organ constraints visible as separate modeled conditions.

### R10 Preserve failure states and reproducible outputs

**Location:** corrected `implic.py` lines 6–9 and 19–20, both exploratory fitting scripts, and plan Phase F.

The corrected implication script represents no qualifying dose above its grid floor as 2.00 Gy before calculating quantiles. This is a display sentinel, not an observed dose limit. It can contaminate quantiles near the transition between infeasible and feasible draws. Store an explicit no-feasible-dose flag, the probability of that state, and conditional dose quantiles among feasible draws. Alternatively, report quantiles of an explicitly defined ordered distribution with a separate infeasible category. Do not assign it a fabricated numeric dose.

Also reconcile reference-validity rules across deterministic and uncertainty calculations. The current deterministic search rejects a nonpositive reference BED, while an uncertainty path checks the candidate. Preserve a separate invalid-model-domain flag and report its frequency. Never silently discard such draws and describe the renormalized support as unconditional support. Distinguish a computational failure, a model-domain restriction, and a valid draw with no qualifying dose.

Preserve chain identity, warmup, posterior draws, diagnostics, data hashes, package versions, and model configuration. Use rank-normalized convergence diagnostics and bulk/tail effective sample sizes for the actual reported derived quantities. Good diagnostics establish computational adequacy, not biological identification. Do not construct a point model by combining marginal medians when the posterior is strongly correlated. Use median predictions or a representative joint draw.

The reviewed exploratory directories contain scripts but no archived posterior draws or run logs. The binding values, parser counts, repair fixture, and HYPNO arithmetic were reproduced here. The audit's posterior medians and credible intervals were not independently regenerated. They should remain provisional and cannot justify an app change, especially because the complete-course selection still needs correction. Archive the future corrected analysis under a new identity. Do not describe a fit with different inclusion rules as reproduction of the missing original posterior.

## 4 Implementation sequence with completion criteria

### Step 1 Freeze scope and define the clinical question

Preserve the two existing documents and both exploratory directories. Create a new implementation branch or directory after approval of the changes. Retain the current Fowler evaluation as a frozen compatibility fixture.

The intended output is a comparison of named schedules for a declared clinical setting. Write a structured definition containing disease site, definitive versus postoperative treatment, planned systemic therapy, mucosal endpoint, assessment window, scoring system, and reference regimen. Separate early glottic radiotherapy, broad-field definitive chemoradiotherapy, and reduced-volume favorable-risk oropharyngeal treatment when their exposures are not comparable.

**Completion:** every planned prediction has a stated target context and comparator. Generic use without a target baseline produces relative model quantities, not absolute risk. Journal-specific submission formatting is not relevant to this internal implementation review.

### Step 2 Build the data and provenance tables

Create the following files under `data/early_tox/`:

| File | Required contents |
|---|---|
| `reports.csv` | Source identifier, citation, digital object identifier (DOI) or source link, local PDF path, page/table locator, cohort and trial links, extraction status |
| `schedules.csv` | Schedule identifier, complete/partial status, intended versus delivered status, source clock, normalized elapsed clock, dose and timing certainty |
| `fractions.csv` | Schedule identifier, fraction index, relative time in hours, dose in Gy, dose level or target region where available |
| `labels.csv` | Cohort/report/schedule links, original label, adjudicated label, uncertainty category, inclusion rule and reason |
| `arms.csv` | Trial/cohort/arm links, assigned regimen, treatment site/context, technique, planned chemotherapy, reported delivered chemotherapy, dose-volume information and missingness |
| `outcomes.csv` | Arm and endpoint identifiers, exact grading definition, threshold, assessment start/end/frequency, number evaluable, event count, missing outcome count or status |
| `evidence_families.csv` | All reports, labels, derived estimates and priors sharing trial participants or source outcomes |

Do not infer patient event counts from a schedule label. Do not round a published percentage into an exact count unless the denominator and rounding yield a unique justified value. Preserve bounds or use the percentage as qualitative information when they do not. A missing outcome is not a zero event.

Read every historical record against its source mapping. Use two independent extractions for new claim-bearing arm counts and reconcile all differences. The plan's check of ten random schedule rows is insufficient for the exclusions and conflicting labels now identified. Complete a search and screening log for new trials, including date, query, inclusion reason and missing full text. The tumor-refit trial list is a starting index, not an acute-toxicity eligibility list.

**Completion:** all raw records accounted for, the selection fixtures in R1 reproduced, all count denominators explicit, no unresolved duplicate evidence in the primary analysis, and source-family overlap eliminated across validation splits. Missing optional dosimetry or delivered-drug information is recorded rather than used to silently discard studies.

### Step 3 Harmonize endpoints before combining counts

Keep prolonged intolerability, maximum grade ≥3 mucositis, duration of confluent mucositis, feeding-tube outcomes, unplanned breaks, and dose-limiting toxicity as separate endpoints. Record the precise version of the Common Terminology Criteria for Adverse Events (CTCAE), Radiation Therapy Oncology Group (RTOG), World Health Organization (WHO), or other scale.

Do not assume “grade 3” means the same observation under each system. Create a written mapping with a primary compatible group, secondary groups, and exclusions. Require observation windows capable of detecting the endpoint after short-course treatment. Where post-treatment ascertainment is unclear, retain the study as a limited-ascertainment sensitivity analysis or qualitative check. A trial intercept cannot repair systematically missed peaks.

**Completion:** every pooled event has a defensible endpoint mapping and time window. HYHOPE's persistent-toxicity dose-limiting endpoint is not substituted for maximum grade 3 mucositis. The four-week and longer HYPNO arms are compared using the published collection window, with its limitations stated.

### Step 4 Start with models that the data can identify

Fit these in sequence rather than immediately fitting one shared biological model:

1. **Clinical benchmark:** compatible randomized arm counts with trial baselines and a simple calibrated exposure score computed from the fixed Fowler parameter set. Calibrate the risk slope using training data. A fixed BED rule without a risk slope is not a fair probability-prediction comparator.
2. **Arm-only candidate:** use the same baseline and risk-calibration structure, allowing only biological combinations supported by the contrast design to vary. Center exposure at the trial's reference arm. If onset is unsupported, fix it within prespecified scenarios. Examine collinearity among total-dose, quadratic-dose, and duration differences before interpreting ratios of coefficients as a or r.
3. **Separate historical classifier:** fit complete-course labels with the source/cohort structure and disputed-label sensitivities. Assess Fowler and an alternative time-boundary model here. Do not force a label-boundary model to predict patient risk.

Primary estimation of radiation-specific parameters uses contrasts with comparable planned systemic therapy and mucosal exposure. Trials in which these factors change with fractionation remain strategy-level benchmarks unless an identifiable model explicitly accounts for the differences. Randomization identifies the assigned-strategy contrast. It does not separate radiation biology from other components deliberately changed by that strategy.

For the initial clinical benchmark within one compatible endpoint/context group, use the following explicit model:

$$y_{ij}\sim\operatorname{Binomial}(N_{ij},p_{ij}),\qquad
\operatorname{logit}(p_{ij})=\eta_j+\kappa(E_{ij}-E_{j,\mathrm{ref}}).$$

Here j indexes trial families, i indexes arms, y is the event count, N is the evaluable denominator, and ηⱼ is the log odds for the reference arm. Shared ηⱼ handles dependence of comparisons within a multiarm trial. The model assumes binomial sampling conditional on the trial baseline and exposure. Assess remaining arm-level overdispersion rather than assuming that this hierarchy eliminates it.

For a reproducible development specification, proposed proper baseline priors are

$$\eta_j\sim\operatorname{Normal}(m_0,\tau_0^2),\qquad
m_0\sim\operatorname{Normal}(0,2.5^2),\qquad
\tau_0\sim\operatorname{HalfNormal}(1).$$

The Normal distribution is written using mean and variance. HalfNormal denotes the absolute value of a zero-mean Normal with the stated standard deviation. All baseline quantities are dimensionless. These numerical priors are analyst-proposed regularization choices, not literature estimates. Record them as proposals in the implementation specification, not as previously approved user inputs. Examine doubling and halving the heterogeneity scale as sensitivity analyses.

Parameterize the risk-slope prior through a declared exposure contrast E* = 10 BED units and its log odds ratio L* = κE*. A proposed development prior is HalfNormal(1) for L*, with HalfNormal(0.5) and HalfNormal(2) as prespecified sensitivity cases. This makes the assumption interpretable on an outcome scale. The reference scale 10 is a numerical modeling choice, not a tolerance margin. Alternatively, an elicited positive 95th-percentile bound Lmax gives a HalfNormal scale Lmax/1.959964. Record the bound, its rationale, and its elicitor before fitting. Do not select it from the observed benefit of a model.

Use the plan's broad α/β prior and the corrected exploratory broad r prior only as explicitly analyst-chosen starting sensitivity assumptions. Their definitions are `lognormal(log(10), 0.354)` and `lognormal(log(0.8), 0.5)` in the corresponding units. They are not literature-derived uncertainty distributions. Compare with the frozen Fowler values. Record all changes in the prior specification before outcome-based model comparison.

The risk-slope and baseline priors require a separate prior predictive check. Show the induced reference risks and changes in risk across the observed exposure contrasts and the intended prediction range. Record any rejected prior and its scientific rationale before outcome-based model comparison. If prediction or default-promotion conclusions reverse across the prespecified prior scales, retain the model's sensitivity designation. A prior can regularize weak information, but cannot turn a nonidentified parameter into an empirical estimate.

Do not add chemotherapy-by-time, technique-by-dose, or multiple random slopes simply because they are biologically plausible. Add an interaction only when the contrast design can inform it. Otherwise use explicit context-specific models or sensitivity scenarios and label the unresolved interaction.

**Completion:** benchmark and candidate use the same endpoint and calibration task. Identification results are recorded. Prior predictive checks are plausible for the stated context. Any shared-biological model is optional and must show predictions with and without historical-label borrowing.

### Step 5 Represent realistic uncertainty in distinct layers

| Uncertainty | Required treatment | What must not be claimed |
|---|---|---|
| Parameter estimation | Joint posterior draws for identifiable parameters, with candidate and reference evaluated in the same draw | Independent uniform ranges are fitted uncertainty |
| Weak identification | Prior sensitivity, profiles, and parameter-recovery results | Narrow posterior alone proves identification |
| Between-study prediction | New-context baseline uncertainty and supported heterogeneity | Uncertainty in a pooled mean equals uncertainty in a new clinic |
| Structural assumptions | Named onset, repair, time-function, and borrowing variants | An arbitrary weighted mixture has a calibrated 95% interval |
| Exposure uncertainty | Recorded missing timing assumptions, intended/delivered differences, and mucosal-dose limitations | Prescription dose is the dose to all relevant mucosa |
| Endpoint uncertainty | Separate scales/windows or justified mapping sensitivities | Trial intercepts make incompatible outcomes identical |
| Patient variability | Patient-risk outcome model when suitable data exist | Parameter draws represent individual patients |

For unsupported model weights, report variant-specific intervals and the range across variants. Call this a sensitivity envelope. Use data-estimated predictive weights only if the validation design supports estimating them, and report uncertainty in those weights. More Monte Carlo draws do not resolve structural uncertainty.

For the first implementation, include the Fowler reference, a broad biological-parameter sensitivity, an early-onset sensitivity, and repair variants. Add an alternative time function only within its stated domain. When clinical data do not cover the intended fraction size, duration, mucosal volume, or systemic therapy, show “extrapolation” separately from model support. A wide interval should not be relabeled a reliable probability merely because it is wide.

**Completion:** the user can identify which variation is fitted, which is an analyst assumption, and which is unresolved transportability. The main table reports dose-specific support and an applicability indicator. It does not identify supported schedules as clinically safe.

### Step 6 Validate the prediction task without leakage

Hold out whole trial evidence families. Exclude their labels, derivative summaries, and any priors estimated from them. Apply the same rule to CHART and Bentzen's derived estimate. Freeze one candidate before outer validation. If candidate selection uses outcomes, perform it inside each outer training fold. All fitting and prior tuning driven by outcomes must respect this partition.

Define two distinct evaluations if both are useful:

- **New-trial prediction:** integrate over the distribution of an unseen trial baseline and predict its arm counts jointly.
- **Conditional contrast prediction:** use only the held-out comparator arm to inform its baseline, then predict the experimental arm. Do not use the experimental outcomes to estimate that baseline. Label this conditional task explicitly.

Compare benchmark and candidate using the same task. Report paired trial-level predictive log-score differences, uncertainty in that comparison, observed-versus-predicted contrasts, calibration, and interval widths. Report predictive coverage with its sampling interval as a diagnostic. Do not use the proposed 80–97% coverage band as a pass/fail rule for a small trial collection.

For default promotion, a conservative proposed decision convention is that the candidate's paired predictive advantage has a lower 95% uncertainty bound above zero, with no unresolved material endpoint or context mismatch. One reproducible summary is a Bayesian bootstrap of the paired held-out trial-family scores: draw Dirichlet(1,...,1) family weights using a recorded seed, calculate the weighted mean score difference, and take its 2.5th and 97.5th percentiles. Use enough reweightings for stable reported quantiles and repeat with an independent seed. This interval describes uncertainty in family weighting conditional on the fitted validation predictions. It does not capture all dependence from overlapping training sets or guarantee calibrated performance in a new clinical context.

Accompany that interval with leave-one-family influence analysis and the prespecified prior/model sensitivities. If the collection is too small or influential for a stable comparison, call it inconclusive. This is a development convention, not a biological threshold. Retain the simpler reference as default and make the alternative available as a labeled research scenario when the evidence is inconclusive. Do not manufacture a precise promotion probability from an unreliable score-uncertainty calculation.

Use HYPNO as the explicit compatibility benchmark in R6. Compare predicted counts or risk differences with the observed data, including baseline and sampling uncertainty. Investigate a material disagreement rather than changing the prior until it disappears. Use HYHOPE and Jacinto as contextual checks with their own volume, systemic-therapy, and endpoint definitions.

**Completion:** no duplicated evidence leaks into a held-out family, all model comparisons are on the same task, and both favorable and unfavorable benchmarks are retained in the report. Simulation-based calibration checks the computation under its own model. Separate misspecification experiments assess whether plausible endpoint, volume, or repair shifts cause poor prediction.

### Step 7 Implement the decision display after validation

Keep the legacy +5 BED rule as an explicit user-selected model margin during development. If α/β varies, specify whether the margin remains a numerical BED increment or is defined in equivalent 2-Gy dose and converted within each draw. These are different policies. Do not silently change between them.

Add a patient-risk difference rule only after specifying a calibrated reference context and a clinically chosen tolerance δ. Neither δ = 5 percentage points nor an 80% support threshold is established by the literature reviewed here. They may be selected policy settings, but must not be presented as empirical safety cutoffs. If no clinical tolerance is chosen, display the risk-difference distribution without automatic accept/reject labeling.

Export all parameters required for the selected prediction task, including baseline/context and risk slope. Preserve joint draws, source hashes, model identity, clock convention, repair convention, and applicable endpoint. Reject mismatched manifests rather than silently mixing model versions.

**Completion:** every displayed percentage is reproducible independently, the reference and candidate share each draw, unavailable baseline information prevents absolute-risk output, and the legacy preset reproduces the original calculation exactly.

### Step 8 Run a focused acceptance suite and stop when it passes

The following tests cover scientific and implementation failure modes. They should accompany the implementation rather than trigger another broad review.

1. Reproduce all raw label counts and the complete-course selection counts in R1. Unknown labels fail explicitly. Preserve original labels.
2. Confirm the eight partial-course records never enter the complete-course likelihood. Verify complete-source replication separately from conservative-label analysis.
3. Confirm fraction times distinguish the four false-merge examples in R2. Independently repeated cohorts remain distinct, duplicate publications do not.
4. Confirm no evidence-family identifier appears in both training and held-out data, including indirect prior sources.
5. Confirm dose/time column masks agree, fraction counts and total doses reconcile, times are ordered, and Gy/hour/day conversions are explicit.
6. Reproduce the no-repair Fowler values 53.10544, 55.86597, and 55.43495 Gy₁₀ for the declared schedules. Preserve the grid fixtures in Section 1.
7. Reproduce the 84.12426-Gy₁₀ repair fixture for the source 70/35 schedule. Infinite separation tends to complete repair, and a same-time pair agrees with a single combined dose in the quadratic term.
8. Demonstrate exact onset cancellation for the hinge evaluator when both courses are beyond onset. Demonstrate noncancellation when one course lies before onset and the other after it. Finite-sharpness softplus does not cancel exactly. Document its approximation tolerance, verify convergence toward the hinge, and include an explicit sharpness with units of inverse days.
9. Verify candidate equals reference gives zero exposure difference in every draw and support of one for a nonnegative margin. This catches accidental independent resampling.
10. Verify fixed-fraction-dose exposure increases with dose under the fitted positive parameters. Do not assume monotonicity in overall time for every structural model without checking it.
11. Preserve no-feasible-dose states without numeric sentinels. Report both their probability and the clearly labeled conditional distribution of feasible ceilings.
12. Reproduce HYPNO's event fractions and the descriptive interval in Section 2. Refuse a validation call that silently substitutes 70/35 for its randomized comparator.
13. Refuse patient-risk output without endpoint, reference context, and baseline draws. A schedule-label probability must never enter that function as patient risk.
14. Reproduce risk and support summaries in an independent implementation. Check joint-support bounds and preserve posterior correlations.
15. Require no unresolved sampler divergences, rank-normalized R-hat below 1.01, and bulk/tail effective sample sizes above 400 for fitted parameters and decision-relevant quantities. These are the plan's computational conventions, not evidence of clinical validity.
16. Run syntax checks, an end-to-end smoke fit, and the frozen validation report. Long fits need recorded seeds, bounded parallel workers, thread limits, logs, memory monitoring, checkpoints, and a resumable orchestrator.
17. Verify every app page uses the intended evaluator and draw manifest. Confirm labels, user guide, and source description agree with the implemented endpoint and units.
18. Run the existing application test suite and a local walkthrough. Deployment remains a separate action requiring approval.

**Stopping rules:** a data/provenance defect blocks fitting. Nonidentification blocks an empirical claim about that parameter, but allows a labeled sensitivity model. An inconclusive predictive comparison blocks default promotion, not completion of the review. Missing clinical context blocks patient-risk output. Passing all checks permits the reviewed implementation to proceed without reopening unrelated project components.

## 5 Disposition of the existing audit findings

| Existing finding | Disposition | Required refinement |
|---|---|---|
| F1 disputed-label coding | Retain | Correct coding is verified. Add complete-course selection and explicit primary/sensitivity label policies |
| F2 short-course data coverage | Retain | Sparse coverage exists, including the two F31 schedules. It does not establish modern transportability |
| F3 CHART prior mismatch | Retain and extend | Preserve repair/onset dependence, sampled α/β conversion, and whole-family leakage rules |
| F4 narrow biological onset prior | Revise | Keep as an elicited sensitivity only. Do not call its limits an evidence-derived interval |
| F5 HYPNO gate | Replace | Correct comparator and patient-risk scale, remove arbitrary pass share and interval-overlap rule |
| F6 duplication | Replace its key and interpretation | Full fraction timing defines schedule identity. Cohort/report identity defines independent evidence |
| F7 overlap | Retain and extend | Apply source-family exclusions to every validation fold and to derived priors |
| F8 risk-scale confusion | Retain and extend | Require baseline/context draws and separate the two endpoints |
| F9 chemotherapy | Revise | Extract planned and delivered treatment separately. Do not mandate naive adjustment for delivered dose |
| F10 endpoint mismatch | Retain and extend | Separate slopes alone do not solve it. Require explicit grading and observation-window mapping |
| F11 validation thresholds | Retain criticism, revise remedy | Fair calibrated benchmarks, family-level predictive comparison, and uncertainty in model ranking |
| F12 minor corrections | Retain with qualifications | Document clock and softplus conventions, repair parity, source access, and correct HYHOPE interpretation |

Do not apply the existing audit as one approved block. In particular, its F4/F5/F6/F9 prescriptions would leave consequential problems in the revised plan.

## 6 Scope and verification record

This review examined both named development documents, all six exploratory Python scripts and the binding R script, relevant application calculations, the source spreadsheet, and the principal local clinical PDFs. Separate methods, code, and clinical reviews were reconciled against direct source reads and numerical checks. This section is an internal verification record, not text for a manuscript or the public calculator.

Methods verdicts are: conceptual validity **partial**, implementation correctness **partial**, parameterization **needs revision**, reproducibility **partial**, fit for purpose **partial**, and strategic alignment **partial**. The objective and use of joint uncertainty are strengths. The main weaknesses are the observation unit, incomplete-course inclusion, unidentified onset, endpoint transport, and the proposed clinical-validation gate.

The surprising-result check gives two different conclusions. “Early toxicity binds under the specified calculator defaults” is **credible as a conditional calculation**. “These data establish a clinical lower limit on course length” is **not established**. The dependency chain is reference schedule and clock → tumor and toxicity parameters → modeled margins → grid feasibility → displayed regimen list. Uncertainty in the toxicity model and margins propagates directly to the reported shortest course.

| Required review route | Status and outcome |
|---|---|
| Methods audit | RUN. Identification, estimands, uncertainty, and validation reviewed. Revision required |
| Pipeline audit | RUN, scoped to exploratory calculations. Canonical runner reports 12 missing input/output references to six absent intermediates. Binding and parser/repair fixtures rerun, syntax checked. Posterior archives absent |
| Paper-versus-code audit | RUN. Formula, selection, prior, clock, repair, and prediction inconsistencies identified |
| Radiation-oncology review | RUN. Clinical applicability, endpoints, volumes, chemotherapy and source assumptions reviewed |
| Surprising-result skeptic | RUN. Conditional binding reproduced, clinical extrapolation not established |
| Simulation audit | RUN. Existing registry schema passes. Scan of corrected fitter gives 37 advisory literals, including unregistered biological/prior constants |
| Subagent scan | RUN. Independent domain findings consolidated and key claims verified |
| Manuscript audit | N/A. Internal development plan and review, not a manuscript |
| Tough grant review | N/A. No grant or proposal |
| Biomedical data audit | N/A. Aggregate schedules and trial-arm data, no patient-level dataset reviewed |
| Release evidence audit | N/A. No submission or application release requested. Source and numeric checks performed for this review |

The registry's successful schema validation does not establish that its distributions are empirically justified or that the new exploratory fit is registered. Its entries include analyst assumptions pending review. The 37 scan advisories include numerical implementation constants as well as substantive parameters and must not be interpreted as 37 scientific errors. No strict simulation run trace was available for these exploratory fits.

No new full posterior fit, systematic extraction of all trial outcomes, patient-level validation, or deployment was performed. The corrected posterior tables therefore remain provisional. No existing plan, audit, application source, or live deployment was changed by this review.

### Equation checks

Every displayed equation in this review was checked. The effective-exposure expression has dose units and reduces to cell-kill BED before onset. Its contrast removes Tₖ exactly when both schedules are beyond onset. The repair term has units of Gy², reduces to the sum of squared fractions at complete repair, and becomes squared total dose for coincident fractions. The rate conversion follows the definition of equivalent 2-Gy dose. Support is a bounded average of indicators. The conditional risk transformation returns the reference risk at zero exposure difference. The joint-support limits are the two-event probability bounds and contain the independent product. The arm model has integer counts and probabilities between zero and one, with a dimensionless linear predictor. The baseline hierarchy uses proper distributions and positive heterogeneity. The HalfNormal percentile conversion uses the 97.5th percentile of a standard Normal. These are mathematical checks, not validation of the biological assumptions.

### Local reference proofs and remaining access gap

Existing full-text PDFs were reused rather than duplicated. Except for Bentzen 2001, the proof files listed below are in `C:\Users\Igor\Desktop\Office\HNSCC_BED_calculator\Custom_BED_calculator_v2\References\`:

- [1] `Fowler_2003_acute_mucosa_tolerance.pdf`, especially PDF pages 2–4.
- [2] Primary Denham 1995 PDF not located. The primary abstract supports the limited claims made here. The exact parameter interval must not be used for prior construction until the full text is obtained.
- [3] `C:\Users\Igor\Desktop\Office\HNSCC_BED_calculator\References_refit_trials\eric_2026_10_22\Part 3\Bentzen_SM_2001_11439207_Radiotherapy_related_early_morbidity_in_head_and_neck_cancer_quantitat.pdf`, especially PDF pages 9–10, Sections 4.2–4.3.
- [4] `Fenwick_2008_mucosal_tolerance.pdf`, including the schedule tables and boundary definition.
- [5] `PardoMontero_2021_mucosal_tolerance_model.pdf`, PDF page 3 Section 2.D, parameter tables, and discussion. Data: `PardoMontero_2021_SupportingInfo_108_schedules.xls`.
- [6] `Hartley_2011_mucosal_parameters.pdf`, PDF pages 1–2.
- [7] `Bentzen_2026_HYPNO.pdf`, PDF pages 5 and 8 for treatment, ascertainment and acute outcomes.
- [8] `Meade_2014_chemo_mucosal_tolerance.pdf`.
- [9] `Moon_2025_HYHOPE_outcomes.pdf`, PDF pages 2–4 and Table 2. The related protocol is `Moon_2024_HYHOPE_protocol_NCT04580446.pdf`. The local outcomes PDF and current publisher record differ in pagination, so the bibliography uses current publisher pagination and these locators use PDF page numbers.
- [10] `Jacinto_2018_hypofractionated_55Gy20fx_cisplatin.pdf`, PDF pages 2–3 for dose levels and grading.

## Verified References

1. Fowler JF, Harari PM, Leborgne F, Leborgne JH. Acute radiation reactions in oral and pharyngeal mucosa: tolerable levels in altered fractionation schedules. *Radiotherapy and Oncology*. 2003;69(2):161–168. [https://doi.org/10.1016/S0167-8140(03)00231-7](https://doi.org/10.1016/S0167-8140(03)00231-7).

2. Denham JW, Hamilton CS, Simpson SA, O'Brien MY, Ostwald PM, Kron T, Dear KB. Acute reaction parameters for human oropharyngeal mucosa. *Radiotherapy and Oncology*. 1995;35(2):129–137. [https://doi.org/10.1016/0167-8140(95)01545-R](https://doi.org/10.1016/0167-8140(95)01545-R). [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/7569021/).

3. Bentzen SM, Saunders MI, Dische S, Bond SJ. Radiotherapy-related early morbidity in head and neck cancer: quantitative clinical radiobiology as deduced from the CHART trial. *Radiotherapy and Oncology*. 2001;60(2):123–135. [https://doi.org/10.1016/S0167-8140(01)00358-9](https://doi.org/10.1016/S0167-8140(01)00358-9).

4. Fenwick JD, Lawrence GP, Malik Z, Nahum AE, Mayles WPM. Early Mucosal Reactions During and After Head-and-Neck Radiotherapy: Dependence of Treatment Tolerance on Radiation Dose and Schedule Duration. *International Journal of Radiation Oncology Biology Physics*. 2008;71(2):625–634. [https://doi.org/10.1016/j.ijrobp.2008.02.008](https://doi.org/10.1016/j.ijrobp.2008.02.008).

5. Pardo-Montero J, Parga-Pazos M, Fenwick JD. Classification of tolerable/intolerable mucosal toxicity of head-and-neck radiotherapy schedules with a biomathematical model of cell dynamics. *Medical Physics*. 2021;48(7):4075–4084. [https://doi.org/10.1002/mp.14834](https://doi.org/10.1002/mp.14834).

6. Hartley A, Sanghera P, Kazi W, Mehanna H, McConkey C, Glaholm J, Fowler J. Correlation of currently used radiobiological parameters with local control and acute and late mucosal toxicity in randomised studies of altered fractionation for locally advanced head and neck cancer. *Clinical Oncology*. 2011;23(1):29–33. [https://doi.org/10.1016/j.clon.2010.08.007](https://doi.org/10.1016/j.clon.2010.08.007).

7. Bentzen SM, Gupta T, Jacinto AA, et al. Hypofractionated Versus Normofractionated Accelerated Radiation Therapy With or Without Cisplatin for Locally Advanced Head and Neck Squamous Cell Carcinoma (HYPNO): A Randomized, Open-Label, Phase 3, Noninferiority Trial. *International Journal of Radiation Oncology Biology Physics*. 2026;125(2):387–398. [https://doi.org/10.1016/j.ijrobp.2025.12.056](https://doi.org/10.1016/j.ijrobp.2025.12.056).

8. Meade S, Sanghera P, Glaholm J, Hartley A. Models of acute mucosal tolerance to radiotherapy alone applied to synchronous chemoradiation schedules in head and neck cancer. *Tumour Biology*. 2014;35(3):2017–2023. [https://doi.org/10.1007/s13277-013-1267-0](https://doi.org/10.1007/s13277-013-1267-0).

9. Moon DH, Avkshtol V, Hughes RS, et al. HYHOPE: A Phase 1 Study of Deintensified Hypofractionated Radiation Therapy for Human Papillomavirus-Associated Oropharynx Cancer. *International Journal of Radiation Oncology Biology Physics*. 2025;123(3):694–702. [https://doi.org/10.1016/j.ijrobp.2025.05.009](https://doi.org/10.1016/j.ijrobp.2025.05.009).

10. Jacinto AA, Batalha Filho ES, Viana LS, et al. Feasibility of concomitant cisplatin with hypofractionated radiotherapy for locally advanced head and neck squamous cell carcinoma. *BMC Cancer*. 2018;18:1026. [https://doi.org/10.1186/s12885-018-4893-5](https://doi.org/10.1186/s12885-018-4893-5).
