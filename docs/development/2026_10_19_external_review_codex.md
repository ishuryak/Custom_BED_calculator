# HNSCC BED calculator v2.4 review and implementation plan

19 October 2026 | Independent mathematical, software, source and clinical review

**Reviewed project:** `C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator`

**Active application:** `Custom_BED_calculator_v2`, commit `ad16504284ed4955b211e981b7d48cc0e31e0bcf`, app version 2.4.0, equation version 2.2. The working tree was clean when reviewed. File and line references below refer to this snapshot. Earlier reviews are historical evidence; this document evaluates the implementation that actually exists now.

## 1 Overall assessment

The core DI and DD calculations are now internally consistent with the calculator's declared model, and the principal mathematical corrections from the previous review have been implemented. The Fowler elapsed-day conversion, fixed DD anchor, dynamic Withers conversion, background growth term, negative-reference protection in the main search, and Fenwick boundary all pass the checks described below. The default D18/HYPNO search and the three new agreement rows reproduce independently. There is no evidence here of a remaining general arithmetic error in the DI or DD BED equations over the tested domain.

The application still needs corrections before its new agreement output should be used as a documented trial-design screen. The most important remaining computational defect is that sensitivity and agreement counts can credit scenarios with negative reference BED, bypassing protections already present in the main search. The most important clinical defect is overinterpretation: a selected BED loss is described as within parameter uncertainty or close in tumor effect, and D18 is described as contradicted by HYPNO and as a pessimistic bound. Those conclusions do not follow from the current calculations or the cited trial.

Preserve the explicitly selected D18 tumor defaults and HYPNO acute reference policy. The synthesis plan records the PI's 19 October decision to use them. This review does not propose reversing that decision. Instead, make their analyst-selected status and interpretation exact, complete the validity and export contracts, and finish the verification gates below.

**Priority order:** repair scenario validity; correct clinical interpretation and units; separate tab-specific calculation inputs; export and replay agreement analyses faithfully; reconcile sources and release metadata; then complete the exact-release checks. No refit, new TCP model, DD onset floor, or switch of default policy is required to complete these software corrections.

The proposed changes have not been applied to the calculator. This is an implementation plan with reproducible evidence and acceptance criteria, not a declaration of clinical validation.

## 2 Scope and strength of verification

The review covers `app.R`, `R/bed_models.R`, `R/export.R`, all seven test files, application and source documentation, the recent fix record, the synthesis plan and D18 plausibility records. It distinguishes current code from archived analysis scripts. The local source collection contains 19 PDFs and one supporting XLS workbook; the source register later in this document identifies their roles. The original 2018 and 2019 model publications were also used to check the tumor equations and fitted-model claims.

- **Existing tests:** 72 test cases, 1,502 passed assertions, zero failures, errors or skips. Two within-suite warnings report Shiny and DT built under R 4.5.3 while this run used R 4.5.1; testthat's package-build warning occurred separately at startup. Explicit-English test runs had no locale warning. Other runtime/oracle startup messages were recorded separately.
- **Independent equation oracle:** 800 parameter/schedule cases, with DD calibration solved by numerical bisection rather than the app's inverse formula. There were 532 calibratable cases, 220 infeasible cases and 48 zero-accelerated-repopulation cases; no status or numerical mismatches, maximum BED discrepancy approximately 2.04e-10 Gy.
- **Independent search enumeration:** eight settings covering all four acute policies, legacy and D18 parameters, agreement, nonzero background growth and an unusable DI endpoint. No membership or ranking mismatches; maximum numerical difference approximately 3.6e-13 Gy.
- **Additional server and replay probes:** reproduced invalid sensitivity credit, unrelated-input dependencies, mislabeled acute units, missing agreement export state and permissive scenario replay. Passing the existing suite does not cover these new counterexamples.
- **Runtime:** the current application started and listened successfully on a local port. A fresh browser walkthrough could not be completed because the browser automation runtime failed before page loading with a kernel-assets path error. The current review therefore verifies server/reactive behavior, not browser layout, MathJax rendering or download interaction. Earlier browser checks recorded in `docs/FIXES.md` are useful historical evidence but are not a substitute for the final release walkthrough.

Tests ran on a copied snapshot so application source was not modified. A successful arithmetic test establishes implementation consistency, not clinical accuracy. The finite oracle sample is strong regression evidence but is not a proof for every value allowed by the broad UI ranges.

## 3 Changes that should be considered closed

| Earlier issue | Current implementation and evidence | Status |
| --- | --- | --- |
| Fowler acute clock was one day too long | `compute_all()` converts inclusive time to elapsed time once; `bed_early()` still accepts model-clock time. 70/35 now gives 53.10544 Gy10. | Corrected; preserve single conversion. |
| Changing comparator changed the DD model | A separate fixed 35 x 2 Gy/47-day anchor is now the default; legacy coupling is explicit. | Corrected in normal mode. Clarify agreement's coupled-anchor semantics. |
| Withers preset assumed tumor alpha/beta 10 | Lambda follows 0.6 x alpha x (1 + 2/r); Hartley's elapsed onset is converted to inclusive day 23. | Corrected arithmetic. Some source-attribution wording remains too broad. |
| Negative reference BED credited attractive gains | Main search excludes the unusable endpoint, hides its delta and records reasons; both-negative state is distinct. | Corrected in main search; not propagated to sensitivity/agreement. |
| Background growth omitted | Both tumor models subtract gT/alpha, including the zero-AR branch; shared g cancels in anchor calibration. | Corrected; explanation of pre-onset shortening still needs updating. |
| Rounded CSV lacked provenance | Candidate ZIP includes full-precision results, all qualifying rows, parameters, anchor, clocks, version and empty-search reason. | Substantially corrected; replay guards and agreement export remain open. |
| Late BED shown at one alpha/beta | Results include late sensitivity at 2 and 4 Gy in addition to the selected value. | Corrected display; still a prescription-dose surrogate. |
| HYPNO compared with the wrong control | Test and six-fractions/week control presets and contextual benchmark text are present. | Corrected. Keep nominal and actual timing distinct. |
| Fenwick absent or liable to use the wrong equation | Published 69.5, 32.2, minus 3.5 boundary, inclusive time, cell-kill BED and 50-day limit are implemented. | Correct for the declared once-daily approximation. |
| No automated test workflow | A test-only workflow exists and does not deploy. | Implemented; clean dependency restoration remains unfinished. |
| Defaults and agreement lacked explicit decisions | D18/HYPNO default decision and agreement feature are recorded and tested. | Authorized changes, not defects to undo. |

The new default behavior differs materially from the old seven-row legacy search. That is a documented change of scientific scenario and acute policy, rather than an arithmetic regression. Do not restore old expected outputs by changing the new defaults.

## 4 Mathematical specification and source mapping

### 4.1 Notation and clocks

Let n be the integer fraction count, d the uniform dose per fraction in Gy, r = alpha/beta in Gy, alpha in Gy^-1, and T the inclusive first-to-last treatment day. The current generated calendar starts on Monday, treats once daily Monday to Friday, and has no holidays or interruptions:

```text
T(n) = 7 floor((n - 1)/5) + ((n - 1) mod 5) + 1
K(n,d,r) = n d (1 + d/r)
EQD2(B,r) = B / (1 + 2/r)
s(x) = log(1 + exp(k x))/k, with k = 3 day^-1
```

The formula gives T = 47, 26 and 19 inclusive days for 35, 20 and 15 fractions. Fowler acute calculations use T_e = T - 1; Fenwick uses inclusive T directly. The built-in HYPNO control is a specifically defined six-fractions/week schedule, 33 x 2 Gy in 38 nominal inclusive days. Do not apply the five-day automatic calendar to that preset.

The softened hinge is an explicit implementation variant of the sharp source threshold. It has s(0) = log(2)/3 = 0.2310490602 days and a small positive tail before onset. For finite x, its excess over max(0,x) is at most log(2)/3. The code uses stable `log1p`, `expm1` and a large-argument branch. A switch to a hard hinge would be an equation change requiring recalibration and new fixtures, not a cosmetic cleanup.

### 4.2 Dose-independent tumor model

```text
B_DI = K - (lambda_DI/alpha) s(T - Tk) - (g/alpha) T
```

This follows the log-survival representation of Shuryak et al. 2018, equation 1, with the declared soft hinge. Lambda and g are rates in day^-1, not Gy/day. The post-onset accelerated BED cost is lambda_DI/alpha; its EQD2/day value is that quantity divided by 1 + 2/r. The whole-course background term remains present before accelerated repopulation begins.

Withers' approximately 0.6 Gy/day compensation in 2-Gy fractions corresponds to 0.6(1 + 2/r) BED Gy/day. It is 0.72 Gy10/day at r = 10 and 0.84 BED Gy/day at r = 5. The current conversion is correct. Likewise, HYPNO's design discussion uses 0.65 Gy/day in EQD2 terms: 66 + 9 x 0.65 = 71.85 Gy; it corresponds to 0.78 Gy10/day, not 0.65 Gy10/day.

### 4.3 Dose-dependent tumor model and calibration

```text
p_kill = 1 - exp(-alpha K/T)
q = C/alpha
t_onset = q T/K
B_DD = K - (lambda_DD/alpha) p_kill s(T - t_onset) - (g/alpha) T
```

The cell-kill response, onset law and rate law follow Shuryak et al. 2018 equations 2-5. C is dimensionless log cell kill, q has units Gy, and qT/K is in days. `p_kill` is a daily average response based on the whole schedule, not a simulation of individual daily fractions or weekend cell counts.

For the fixed anchor with K0 and T0, equality of the DI and DD BEDs requires:

```text
y = lambda_DI s(T0 - Tk) / (lambda_DD p_kill(K0,T0,alpha))
x0 = s_inverse(y) = log(expm1(k y))/k
t0 = T0 - x0
q = K0 t0/T0
```

The app implements this derivation correctly for positive rates and a finite positive t0. The same alpha and g are used by both models, so g cancels in the equality. Zero lambda_DI intentionally disables accelerated repopulation in both coupled models, using infinite onset and B = K - gT/alpha. Positive DI with zero DD rate, or no positive feasible onset, produces an unavailable DD model. Do not replace unavailable values with zero or force an onset floor to make a parameter set work.

This is a calibrated DD variant. The 2018 study fitted DI and DD separately to outcome data, with its own DD alpha, maximum rate and C. The app's shared-alpha, equality-at-anchor construction is not that fitted model. Reproducing equality at the anchor is a mathematical identity imposed by the calculator, not independent biological evidence. The updated About page correctly explains this distinction.

The fitted 2018 DD comparisons quoted in the current documentation also reproduce: approximately +3.1687 Gy10 for the GORTEC comparison, +2.0607 for TROG and +0.5375 for HYPNO with the stated fitted construction. These are model calculations, not proof that either the fitted model or current calibrated variant is clinically correct for a new population.

### 4.4 Acute mucosal model and policies

```text
B_acute = n d (1 + d/r_e) - (lambda_e/alpha_e) s(T_e - Tk_e)
lambda_e = log(2)/Tp_e
Default: r_e = 10 Gy, alpha_e = 0.35 Gy^-1,
         Tp_e = 2.5 days, Tk_e = 7 elapsed days.
```

These are the Fowler 2003 assumptions, also used in the HYHOPE rationale. The current app gives 53.10543995 Gy10 for 70/35 and 55.86597229 for nominal 55/20. The latter is the default HYPNO reference policy, computed using the same acute parameters as a candidate. If those parameters change, its numerical value changes; this is an internally matched model comparison, not a trial-estimated toxicity threshold.

Four different policies are present and should remain distinct: reference-relative acute BED, HYPNO-reference acute BED, the lower Fowler zone edge, and Fenwick's cell-kill boundary. Fowler's 59-61 Gy10 RT-alone range and possible 3-5 Gy10 chemotherapy decrement produce the displayed 54-58 range and a conservative lower-edge policy of 54. They do not create a toxicity probability. The app appropriately blocks the fixed Fowler policy when the specified Fowler acute parameters are changed.

The chemotherapy checkbox changes Fowler's displayed zone/policy only. It does not add tumor radiosensitization, predict systemic benefit, or account for chemotherapy in the HYPNO, relative or Fenwick policies. Preserve this explicit separation.

### 4.5 Fenwick boundary

```text
B_limit(T) = 69.5 (T/32.2)/sin(T/32.2) - 3.5 Gy10
0 < T <= 50 inclusive days; angle in radians
Compare with K(n,d,10), without a second repopulation subtraction.
```

The implementation uses the modified boundary in Fenwick 2008, equation 10, printed page 630, rather than Meade's rounded time constant. The original cell-kill framework accounts for incomplete repair; the present once-daily grid uses the source's complete-repair approximation for intervals of at least 24 hours. This is reasonable for the current schedule generator. It must not silently become a model for twice-daily treatment, short intervals, mixed fraction sizes or multiple phases.

The source boundary is an empirical discriminator developed from historical RT-alone schedules. It is not a modern chemoradiation risk model. Its domain enforcement at 50 days and independence from custom Fowler parameters are correct. The comment at `R/bed_models.R:68-71` incorrectly says Fenwick counts elapsed days; correct that comment without changing the working Fenwick calculation.

### 4.6 Late BED, search and ranking

```text
B_late = n d (1 + d/r_l)
Default r_l = 3 Gy; no accelerated-repopulation term.
```

This is correct LQ arithmetic at the prescribed tumor dose. It is not an organ dose-volume calculation. A late-BED comparison can reverse with r_l: 20 x 2.9 Gy is below 70/35 at r_l = 3 (114.0667 versus 116.6667 Gy3), but above it at r_l = 2 (142.1 versus 140 Gy2). Showing several r_l values is useful; replacing actual OAR doses with prescription dose remains a substantive limitation.

The main search correctly enforces its discrete dose grid, higher dose per fraction than the comparator, shorter duration, selected tumor criterion and acute/late policy. The next-dose-step column evaluates the actual next grid point. Ranking by the best model is explicitly optimistic; ranking by the smaller change requires both models to be usable. One displayed dose per fraction count is a selected point from a finite grid, not a continuous mathematical optimum.

The agreement tab currently chooses the best current-parameter point for each fraction count and then evaluates those selected points across 15 scenarios. It is not an exhaustive optimization over scenario robustness. The 15 scenarios are five presets crossed with three alpha/beta values, not independent model validations or samples from a probability distribution. DI-only presets inherit the current DD maximum rate; all inherit current g and the anchor convention. These dependencies must be visible in the analysis record.

## 5 Numerical results and their clinical interpretation

All figures below use the declared soft hinge, fixed 35/2/47 anchor and g = 0 unless stated otherwise. Acute BED uses elapsed time. They are acceptance fixtures, not proposed prescriptions.

| Scenario and schedule | DI BED | DD BED | Acute BED | Late BED3 or note |
| --- | --- | --- | --- | --- |
| Legacy 35 x 2 Gy / 47 d | 65.00000000 | 65.00000000 | 53.10543995 | 116.66666667 |
| Legacy 20 x 2.75 Gy / 26 d | 70.12417477 | 58.02143210 | 55.86597229 | 105.41666667 |
| D18 35 x 2 Gy / 47 d | 71.51428571 | 71.51428571 | 53.10543995 | Anchor equality |
| D18 20 x 2.75 Gy / 26 d | 70.12445771 | 62.84730123 | 55.86597229 | HYPNO test schedule |
| D18 33 x 2 Gy / 38 d | 72.62857143 | 68.64033914 | 55.43495381 | HYPNO control schedule |
| D18 15 x 3.1 Gy / 19 d | 60.91500000 | 56.37052984 | 52.20114973 | HYHOPE schedule arithmetic |

D18 gives HYPNO test-minus-control differences of -2.50411372 Gy10 under DI and -5.79303791 under DD. These values alone do not establish contradiction with HYPNO's noninferiority result. A BED-to-LRC relationship, population and treatment context, uncertainty and the trial's actual statistical margin would be required for that assessment.

With the current candidate-search defaults, three rows pass the +1 Gy tumor criterion under DI only: 29 x 2.25 Gy/39 days, 28 x 2.30/38 and 30 x 2.20/40. Their DI gains are approximately +1.1884, +1.1263 and +1.1200 Gy10; DD changes are -2.4080, -2.8648 and -2.0495. Thus the default list is a list of model-discordant candidates, not a both-model efficacy gain. The membership labels correctly reveal this.

At D18 with the same +1 Gy/7-day criteria and a 0.05-Gy grid up to 5 Gy, the acute policies produce zero displayed rows for the reference-relative policy, three for HYPNO, 11 for Fowler 59 and 12 for Fenwick. The corresponding all-qualifying counts are zero, three, 14 and 17. Differences reflect explicit policies; they are not evidence that a more permissive policy is more accurate.

The default agreement analysis compares with HYPNO 55/20, permits a BED loss up to 3 Gy and requires at least one day saved:

| Schedule | Delta DI | Delta DD | DD onset | Selected scenarios meeting both margins |
| --- | --- | --- | --- | --- |
| 17 x 3.05 Gy / 23 d | -2.46020777 | -1.58951637 | 14.28035 d | 15/15 |
| 18 x 2.90 Gy / 24 d | -2.78645905 | -1.96142511 | 14.97343 d | 15/15 |
| 19 x 2.80 Gy / 25 d | -2.02848474 | -1.47714451 | 15.42370 d | 15/15 |

These three numerical rows and their 15/15 counts are correct at the default settings. They save one to three nominal calendar days and all have lower modeled tumor BED than HYPNO under both current models. They satisfy a user-chosen BED-loss budget; they have not been shown clinically equivalent, noninferior or safe. The two models share alpha and are linked by calibration, and the sensitivity sets share assumptions. Agreement is conditional consistency within this model family.

The all-scenario result does not resolve uncertainty in lambda_DD, anchor choice, clinical population, systemic therapy, actual timing, dose distribution or the translation from BED to outcomes. In particular, the app's agreement evidence does not justify describing a 1.5-3 Gy loss as negligible or within a quantified uncertainty interval.

## 6 Findings and exact corrective work

### F01 High priority Apply one validity rule to every scenario

**Files:** `R/bed_models.R:401-446` (`sensitivity_table`, `agreement_regimens`); `app.R:695-705` (`sens_summary`).

**Defect.** The main search requires finite, positive reference BED for the endpoint and positive candidate BED. Sensitivity simply subtracts BEDs; agreement only checks that DD delta is not NA before counting success. A finite negative reference can therefore produce a very large favorable delta and be counted as supporting evidence. Calibration success is not the same as usable reference/candidate BED.

**Reproduction.** Start with D18, set lambda_DI = 0 and g = 0, use a custom 35 x 2 Gy comparator over 150 inclusive days, keep the fixed 35/2/47 anchor, use the HYPNO acute cap, and set agreement comparator to the sidebar with margin -3. The app lists 32 x 2.10 Gy/44 days and 33 x 2.05/45 with 15/15 scenarios. For the first row, only three scenarios have positive DI and DD reference and candidate BEDs. Twelve have a negative DI reference. The legacy r = 10 scenario has reference DI = -38 Gy and candidate DI = 65.312 Gy, yielding an apparent +103.312 Gy delta and a counted pass.

The 150-day example is deliberately a domain stress test within accepted input limits, not a plausible trial schedule. Its purpose is to prove that the new screen can violate the explicit validity contract. Adding an arbitrary shorter maximum time would hide this instance without fixing the duplicated logic.

**Implementation.** Add a shared endpoint evaluator that returns raw reference and candidate BED, `valid_reference`, `valid_candidate`, an enumerated reason, and a decision delta only when usable. Use it in the main search, sensitivity and agreement. Treat nonfinite or nonpositive values as ineligible under the existing policy; keep raw values separately for diagnostics. Retain distinct calibration-unavailable and reference/candidate-domain reasons. Do not silently clip a BED to zero.

For each row return `n_total`, `n_valid_both`, `n_pass_both`, `n_invalid` and, if useful, per-model counts. Display, for example, "3 of 3 evaluable scenarios meet the margin; 12 of 15 ineligible", not an unqualified 3/3 or 15/15. Exclude ineligible deltas from favorable ranges. If none are evaluable, return NA with a reason, not Inf/-Inf or zero. Whether a regimen must pass all valid or all prespecified scenarios is a separate, explicit policy; invalid scenarios must never be counted as passing.

**Verify.** Add a regression for the exact counterexample: 15 total, 3 valid in both, 12 ineligible, 3 valid passes for 32 x 2.10. Cover negative candidate BED, zero BED, unavailable DD, both endpoints invalid, and all-scenarios-invalid ranges. The three default agreement rows must remain unchanged at 15 valid passes out of 15. Run the same cases through the standalone Sensitivity tab and exported scenario detail, not just a helper function.

### F02 High priority Remove unsupported outcome and uncertainty claims

**Files:** `about.md:10-12, 58-60, 66-82, 94-99, 115, 145`; `app.R:140, 656-663`; relevant source notes and fix-record descriptions.

**Defect.** Several statements convert model arithmetic into unsupported clinical conclusions: a few Gy are "within" parameter uncertainty; agreement schedules are close to HYPNO in tumor effect; D18's DD deficit is not supported by HYPNO; D18 is a pessimistic bound; the acute model cannot distinguish exactly 1-2 Gy10; and a negative BED threshold is described as a noninferiority margin. These claims conflict with the app's otherwise correct statement that it lacks an outcome dose-response relation.

**Replace with these meanings.**

- D18 is an analyst-selected scenario informed partly by DI literature, not a fitted DD parameter set or established upper/lower bound. It predicts the stated HYPNO BED deficit; compatibility with observed LRC cannot be determined from the BED difference alone.
- The agreement rows meet a specified modeled BED-loss allowance under both models and the selected scenarios. No uncertainty interval or probability of noninferiority has been calculated.
- Use "allowed BED loss" or "BED screening margin" for the numerical control. A clinical NI margin concerns a defined outcome at a defined time and needs an independently justified statistical design.
- HYHOPE illustrates that similar nominal acute BED does not ensure similar toxicity. It does not establish a measured 1-2 Gy resolution limit for this model.
- Hartley's observed relation between an approximately 15% reduction in late mucosal BED and sparing is an aggregate finding in the analyzed trials, not a universal minimum required for a clinical benefit.

**Pre-onset explanation.** Replace the categorical claim that a course ending before Tk gets no DI benefit from shortening with a conditional explanation: accelerated DI repopulation is approximately absent sufficiently before Tk when g = 0; the soft hinge has a tail, and positive g still gives a whole-course time benefit. For 17 x 3.05 versus HYPNO under D18, increasing g from 0 to 0.0154 changes delta DI from -2.46020777 to -2.32820777 Gy, an improvement of 0.132 Gy = g x 3/0.35 even though both courses end before day 28. At g = 0.1 the improvement is 0.85714286 Gy. DD shortening can also change daily cell-kill response and onset; it does not invariably create a gain.

Different tumor and mucosal fractionation sensitivities also invalidate the blanket claim that this comparison can only lose tumor BED under an acute cap. With D18 except tumor r = 2 and g = 0, the same 17 x 3.05 schedule has delta DI = +0.29679223 Gy while its acute BED is 55.78172691 Gy10, below HYPNO's 55.86597229. This is a model counterexample within the input range, not evidence that r = 2 is an appropriate default for HNSCC.

**Verify.** Inspect all visible strings and active documentation for each claim, not just the About page. A test should confirm the g-dependent explanation and unavailability messages. Retain the exact default numerical examples. A clinical text checklist should reject unsupported equivalence, quantified uncertainty, guaranteed bounds and a direct Gy-to-percentage-point translation.

### F03 Medium priority Separate analysis-specific inputs and anchors

**Files:** `app.R:489-510, 600-608, 636-653`.

**Defect.** Agreement calls `search_args()` and only afterwards overwrites the candidate fraction range, gain and days. As a result, clearing candidate `min_gain` or `min_days`, or an invalid candidate fraction range, blocks the independent agreement calculation. Separately, clearing the evaluated regimen's fraction count leaves a valid three-row candidate search but blocks its export because the scenario builder reads that unrelated regimen for context.

**Implementation.** Split configuration into shared dose-grid/toxicity policy inputs, candidate-specific search inputs, and agreement-specific inputs. Each calculation validates only the values that determine its result. Compute relative limits against that analysis's comparator. Nonessential context must be optional, not a dependency capable of invalidating a valid search export.

Make the effective comparator and calibration anchor explicit on each tab. In legacy-coupled mode, agreement currently receives the sidebar-derived anchor even when its own comparator is HYPNO. Either preserve this as a clearly named "anchor coupled to sidebar comparator" mode, or define coupling per analysis and record that deliberate behavior change. Do not let the phrase "coupled to the comparison" conceal two different comparisons. Fixed-anchor mode should remain the default: comparator changes must preserve the anchor and absolute model predictions, while changing the analysis's own comparator intentionally changes relative deltas.

Reuse a reason-to-message formatter with the agreement calculation's actual endpoints. An unavailable agreement result currently sends the user to the candidate tab for an explanation, although that tab may use a different comparator or permit DI alone. The agreement tab should explain its own failure directly.

**Verify.** Start at defaults; clear candidate gain and days one at a time. Agreement must retain its three rows. Clear the Results regimen fraction count; candidate export must still reproduce the valid search, with context absent or explicitly invalid. Change sidebar comparator while agreement remains HYPNO in fixed-anchor mode: absolute and relative agreement values must not change. Test the chosen coupled-mode contract separately and display its actual anchor.

### F04 Medium priority Export the analysis the user actually reviewed

**Files:** `app.R:600-619, 636-693`; `R/export.R`; new agreement export tests.

**Gap.** The only download serializes the candidate search. At defaults it records comparison 35/2/47, minimum gain +1 and seven days saved. The agreement analysis instead uses HYPNO 20/2.75/26, margin -3, one day saved and 15 sensitivity scenarios. Agreement inputs and detailed scenario results are absent. A candidate ZIP is correct for its own tab but cannot reproduce the new agreement conclusion.

**Implementation.** Add an explicit analysis type and separate agreement download. Save the effective comparison, fixed/coupled anchor and coupling target, current parameters, g, clocks, hinge, grid, acute/late policies and numerical caps, criterion, loss allowance, days saved, scenario definitions and all inherited values. Include the exact selected rows, all qualifying grid rows, and a long-form candidate-by-scenario table containing raw BEDs, usable deltas, statuses and reasons from F01. Preserve the one-row-per-fraction ranking rule as a declared selection rule.

Store selected-scenario min/max as deterministic ranges, not confidence intervals. If the same numerical BED allowance is applied at r = 7, 10 and 15, say so. It corresponds to different EQD2 allowances: -3 BED Gy is approximately -2.333, -2.500 and -2.647 EQD2 Gy respectively. Alternatively offer an explicitly defined fixed EQD2 allowance and convert it per scenario; do not change the interpretation silently. Neither choice supplies an outcome NI margin.

Add equation, selection-policy and scenario-set versions, source identifiers, code revision/hash and dependency provenance. Do not rely solely on `EQUATION_VERSION` to identify changes to filtering or ranking that leave the formulas unchanged. Human-readable README text should identify which analysis the ZIP reproduces and explain unavailable endpoints. Update the filename from hard-coded `bed_v2.1` to the actual version.

**Verify.** Export and rerun both default candidate and agreement analyses in a new process. Match all row identities, all qualifying rows, full-precision deltas, statuses, the 15 scenario definitions and counts within 1e-8 Gy. Include empty, infeasible-DD, invalid-sensitivity, nonzero-g, custom acute alpha/beta, coupled-anchor and each acute-policy cases. No stale sidebar comparison may replace the agreement comparator.

### F05 Medium priority Enforce the saved model during replay

**Files:** `R/export.R:19-63`, especially `rerun_scenario()`.

**Defect.** Replay accepts fabricated schema and equation versions, a hard hinge, an elapsed tumor clock and a changed saved calibration without warning; it returns the same current-code results. The unsupported schema/equation/hinge/tumor-clock metadata changes the scenario identifier while being ignored by calculation. Saved calibration is also ignored, but does not change the identifier because it is excluded as derived data. Current same-version exports replay correctly, but the contract becomes silently wrong when models evolve or a scenario file is edited.

**Implementation.** Validate required fields, types, lengths, finite numeric values, ranges, enumerated clocks and policy names at public calculation/replay boundaries. Accept only supported schema/equation/selection/scenario-set versions or route through an explicit migration. Verify the requested hinge type and k against the implementation, or implement dispatch to the actual requested engine. Recompute calibration from authoritative inputs and compare the saved derived calibration with tolerance; report inconsistency rather than using it blindly or ignoring the mismatch.

Canonicalize named object keys recursively for the scientific identifier while preserving ordered arrays such as grids. The current ID changes when parameter key order is reversed even though the calculation is identical. Choose a stable UTF-8 byte representation and newline convention; test across clean processes. Exclude presentation-only context, as the recent fix already does, but include every scientific selection setting. A checksum is provenance, not a scientific validation certificate.

**Verify.** Each unsupported schema, equation, hinge and clock must produce a specific error without results. Missing, nonfinite, malformed or out-of-domain required inputs must fail before arithmetic. A known supported scenario must still reproduce. Reordered object keys must preserve identity; changes to actual model, anchor, grid, policy or scenario definitions must change it. A corrupted stored calibration must be detected.

### F06 Medium priority Correct units and scenario provenance at the point of use

**Files:** `app.R:528`, `R/bed_models.R:68-71, 309-350, 399-416`, `about.md:93-99`, `References/README.md`.

**Reproduced unit defect.** Set acute alpha/beta to 5 and retain the HYPNO policy. The actual cap is 70.99097 BED Gy at alpha/beta 5, but the search note labels it 71.0 Gy10. Use the existing `early_unit()` helper in this branch; retain Gy10 only for the fixed-alpha/beta Fowler/Fenwick definitions. Include alpha/beta in table/export units wherever a numerical BED might be mistaken for Gy10.

**Provenance defect.** "Each preset sets only the values its source gives" is inaccurate. Withers supplies compensation and approximate onset information, not this app's chosen alpha = 0.2. D18's lambda_DD = 0.4 is an analyst construction informed by DI-rate plausibility, not a direct published DD estimate. A value lying in an individually plausible range does not establish a valid joint parameter distribution or transferability between model meanings.

**Implementation.** For every preset parameter store value/rule, units, source locator, source clock, transformation, and origin: direct source estimate, derived conversion, chosen representative value, inherited current value, or analyst scenario. In sensitivity show/export the inherited DD rate, g and anchor explicitly. Fix the erroneous Fenwick clock comment. Clarify that a missing model means unavailable, not a zero delta. Where only one model is usable, name its score rather than calling it the smaller of two available changes.

**Verify.** Acute r = 5 labels and exported values must agree; default r = 10 labels remain Gy10. Source-origin records must correctly identify Withers alpha, representative Tk, Hartley conversion, HYHOPE protocol values and D18 DD assumptions. Check that changing inherited g or lambda_DD changes the scenario record and identifier.

The 3.1-Gy flag should be named literally, for example "Above included benchmark fraction-size range". It describes the selected benchmark collection, whose highest fraction size comes from favorable-risk HYHOPE. A false flag is not evidence that the whole regimen has clinical precedent in the proposed population, with its proposed systemic therapy and volumes.

### F07 Medium priority Make the current source catalogue and release record consistent

**Files:** `References/README.md`, `about.md` reference 6, `DESCRIPTION`, `README.md`, `docs/FIXES.md`, model/export file headers and archived-record pointers.

Some catalogue entries still say full texts were not seen even though later sections use the downloaded full papers. The catalogue's Fowler 0.7 Gy/day entry must preserve EQD2 units; its E1 statement that g is omitted is now obsolete for current code. The local HYHOPE protocol is version 4.2 while the About bibliography identifies version 4.1. `DESCRIPTION` and some source/export headers still say 2.1 despite app version 2.4.0.

Correct `about.md:36` from 33 trial arms to the 2018 paper's stated 35 treatment regimens from 16 publications, totaling 7,283 patients (printed p21, Data sets, Table 1). If a different subset is intended, identify its exclusions explicitly. The fitted parameters and rounded GORTEC/TROG/HYPNO calculations on that same paragraph otherwise reproduce.

Create one current source ledger with DOI/study ID, exact local file and version, pages/equations/tables, evidence type, extracted parameter or endpoint, transformation and limitations. Link current explanatory text to that ledger. Keep old analyses frozen with visible snapshot/date and errata; do not rewrite old outputs as if they came from the current engine. Any new D18 plausibility run should call the current engine and serialize its full construction, including rate versus maximum-rate distinctions and acute clock.

**Verify.** No active source record has contradictory access/version/unit status. Confirm all 19 local PDFs and the workbook are accounted for. A single version source should generate UI, package metadata and download names. The fix record must distinguish formula corrections, filtering corrections, new default scenarios and changed clinical wording.

### F08 Release requirement Finish reproducibility and interface verification

**Files:** dependency manifest/lock, launcher, CI, test suite and release checklist.

The existing test-only CI is a real improvement. Dependency versions are not locked, and the launcher is tied to a local R installation. Add `renv.lock` or an equivalent reproducible restore specification with a documented R version. Make the launcher discover R or accept an explicit configured path. Keep deployment separate from testing.

Add regression tests that first reproduce F01-F06 on this snapshot and then fail for the intended reason before the repair. Use independently calculated fixture values and meaningful nonempty cases; avoid assertions that pass on empty tables. Add a small mutation set: bypass validity in sensitivity, use the wrong agreement comparator, ignore hinge/version during replay, leak candidate gain into agreement, and restore the wrong unit label. The tests should catch each mutation.

Complete a fresh browser walkthrough on the exact final artifact: all tabs and equations; changing and clearing inputs; presets and Withers rule; fixed/coupled anchor; negative/unavailable endpoints; nonzero g; each acute policy and custom acute alpha/beta; empty and populated searches; agreement scenario counts; horizontal table access; both downloads and their contents. This review's browser tooling failure is why this remains an explicit gate.

**Verify.** A clean dependency restore and all applicable tests pass on the final revision; assertion failures, skips and package warnings are separately reported. Independent equation/search fixtures still agree. The walkthrough record identifies commit, browser and checked behaviors. No deployment is part of this plan; before any separately authorized deployment inspect the deployment-file list and target, excluding local references, historical grids and review evidence from the app bundle.

## 7 Clinical use appropriate to HNSCC trial exploration

### 7.1 What the current benchmark evidence supports

**HYPNO.** The randomized phase III trial compared 55 Gy/20 fractions at five fractions/week with 66 Gy/33 at six/week, not 70 Gy/35. It randomized 395 and 397 patients. Three-year LRC was 50.7% versus 51.2%; the failure hazard ratio was 1.098 with 90% CI 0.93-1.30. Three-year grade 3+ late effects were 18.8% versus 20.2%, HR 0.926 with 90% CI 0.68-1.25. The prespecified NI margin was 10 percentage points. These are trial-level clinical results, not an observed BED-to-LRC calibration. The current app reproduces the main reported numbers correctly. [S01]

Grade 3+ acute mucositis was 190/373 (50.9%) versus 208/380 (54.7%), using acute-evaluable denominators. The 50.9% is pooled across the hypofractionated arm, including RT-alone patients, not a chemotherapy-only risk estimate. Cisplatin was prescribed to 597/792; the paper also summarizes approximately 75% receiving it. The trial used weekly 35 mg/m2, with planned four versus five administrations and differing cumulative exposure. HPV status was not tested; population, tobacco exposure, techniques, target volumes and systemic therapy affect transferability. The nominal 26/38-day schedules are not patient-level delivered times. [S01, printed pages 389, 393-397 and Table 3]

The selected HYPNO cap is useful as an explicit modeled reference with a known clinical burden. Equal radiation-only acute BED under one model does not imply the same mucositis probability, duration, feeding-tube dependence, dysphagia or treatment completion. Retain this limitation alongside the policy; do not hide the observed burden behind a pass/fail label.

**HYHOPE.** This was a phase I deintensification study in favorable-risk HPV-associated oropharyngeal cancer with reduced elective volumes and planned three weekly cisplatin doses of 40 mg/m2. It was designed against 60 Gy/30, not as equivalence to 70/35 in unselected HNSCC. Twenty-four evaluable patients received the selected 46.5 Gy/15-fraction level; grade 3 acute toxicity of any type occurred in 7/24 and oral mucositis in 5/24. No late grade 3+ events were reported at median 26.6-month follow-up. Two-year locoregional recurrence was 13.8%, 95% CI 3.2-32%, a cumulative-incidence endpoint rather than a randomized comparative LRC estimate. Grading used CTCAE v5. DLT meant RT-attributed grade 4 during or within three months, or persistent grade 3 mucositis at three months; it was distinct from any acute grade 3 event. [S13-S14]

The protocol's withdrawal of the 12-fraction level after concerns in a parallel postoperative study remains relevant to the limits of acute BED screening. It is not a calibrated estimate of model error or a basis for assigning a universal numeric toxicity-resolution threshold. Protocol assumptions, protocol amendments and the outcomes paper must remain separate records.

**Other 55/20 cohorts.** Jacinto studied 20 patients with cisplatin, grade 3 mucositis 40%, dermatitis 30%, and median actual duration 29 days (27-34). Miller reports 119 referred, 113 completing, and grade 3 mucositis in 72/113 at treatment end, 26% four weeks later, predominantly in older patients without concurrent chemotherapy and using RTOG grading. Madhava's carboplatin escalation study had 19 patients, grade 3 mucositis in 16 and prolonged grade 3 for at least four weeks in three. [S07, S10, S12]

These series provide context for burden and feasibility; their rates are not interchangeable estimates of a causal chemotherapy effect. Preserve grading system, assessment time, denominator, selection and delivery details. In particular, the uncorrected LQ BED10 of 55/20 is 70.125 Gy10, whereas uncorrected EQD2 is 58.4375 Gy: neither can be relabeled 70 physical Gy. Time correction is a separate, model-dependent comparison.

### 7.2 A concrete trial-exploration workflow

1. **Define the intended population and treatment strategy.** Record site, HPV/p16 status where relevant, stage/risk group, definitive versus postoperative intent, prior RT, target/elective volumes and systemic therapy. Do not transfer a favorable-risk deintensification benchmark to unselected locally advanced disease solely because fraction sizes overlap. The current calculator is not specified for reirradiation, SBRT, mixed-dose phases, twice-daily treatment or spatial dose summation.
2. **Choose a clinical comparator and declare the anchor separately.** A trial comparator should match the proposed clinical question; the numerical DD anchor is an assumption of this model. Show both in every analysis and export. Compare HYPNO and HYHOPE with their own trial context before discussing a new regimen.
3. **State the efficacy question in two layers.** Specify a model-space BED allowance for screening, and separately the intended clinical endpoint/timepoint and clinically acceptable effect difference. The app can calculate the first. A statistical trial design or a validated outcome model is required for the second. Do not convert the selected -3 Gy into an assumed percentage-point margin.
4. **Evaluate both models and perturb consequential assumptions.** Retain discordance, empty searches and invalid scenarios. Vary DI loss/onset, DD maximum rate and anchor, alpha/beta, g and acute-policy assumptions in explicitly named scenarios. Inspect whether a candidate depends on one attractive model, a particular cap or a particular calendar. The current 15 scenarios are a starting panel, not a complete joint uncertainty analysis.
5. **Use actual planning and delivery information for toxicity assessment.** Prescription-dose late BED does not represent cord, brainstem, constrictors, larynx, oral cavity, mandible or salivary-gland dose-volume exposure. Obtain protocol-specific OAR plans and assessment endpoints. Account for the proposed systemic-treatment schedule and likely interruptions rather than assuming shorter RT preserves the same chemotherapy exposure and adherence.
6. **Create an auditable hypothesis package.** Save the exact calculation, full settings, validity statuses, all selected scenario values, source context and planned clinical endpoints. Any phase I dose ladder, escalation rule, toxicity stopping rule or NI design should be independently specified and reviewed by the clinical/statistical team. The calculator can support that package; its pass/fail fields do not supply those rules.

These are requirements for using the app's outputs responsibly in trial development, not requirements to implement an entire treatment-planning system before correcting this release. A compact clinical-context form and export fields can capture the key distinctions without introducing unvalidated risk calculations.

### 7.3 Scientific extensions that remain separate

Refitting DI and DD to outcome data is a separate scientific analysis. It requires reconciled trial-arm data, compatible endpoints and timepoints, appropriate chemotherapy/population covariates, an explicit likelihood/TCP relation, comparable model flexibility, parameter uncertainty and evaluation outside the fitting data. Counting favorable signs, plausible parameter combinations or passing schedules cannot replace it. Do not label the current shared-alpha calibration as a reproduction of the fitted 2018 model.

Adding actual fraction times and incomplete repair, an OAR module, a chemotherapy effect, or an NTCP model would likewise need a separate specification. Meade's cross-trial chemotherapy estimates should not be turned into an automatic cisplatin or carboplatin dose conversion. Pardo-Montero's refitted mucosal framework and workbook are potential source-reproduction material, not an already validated alternative engine in this app.

The archived D18 grid scripts are now explicitly marked as historical and not directly rerunnable. Preserve that mitigation. For a future reproducibility refresh, freeze the old code/parameters for historical reconstruction and write a separate current-engine pipeline with explicit input/output paths. Some archived scripts mix a fixed lambda_DI = 0.2 grid with current default parameters, or write to hard-coded external output locations; simply changing their source path to v2.4 is not a valid update.

## 8 Publication and evidence register

The register covers all 19 local PDFs and the XLS workbook, including supporting and educational material. Relevant full-text readings from the previous pass remain applicable to unchanged files; this pass rechecked the passages supporting current numerical inputs, new interpretation and clinical claims. It does not imply that every cited but absent historical paper or unavailable supplement was newly obtained. Page locators below distinguish local PDF pages from printed pages when necessary.

| ID and local file | Use in this app or review | Source locator and limitation |
| --- | --- | --- |
| S01 Bentzen_2026_HYPNO.pdf | Trial schedule/comparator; chosen acute/late model reference; clinical outcomes | 12 PDF pages. Printed p389/PDF3 design units; pp393-394 delivery/acute events; p395/PDF9 model versus trial distinction; Table 3 p396/PDF10. Not an acute NTCP threshold. |
| S02 Fenwick_2008_mucosal_tolerance.pdf | Implemented optional boundary | 10 pages. Eq 2 p626 for cell kill/repair; Eq 10 p630 for boundary, Eq 11 and discussion for daily repair approximation. Historical RT-alone classifier. |
| S03 Fowler_2003_acute_mucosa_tolerance.pdf | Acute equation/defaults, 59-61 zone, tentative chemotherapy decrement | 8 pages. BED equation and Tables 1-2; elapsed-day examples. Empirical tolerance range, not a probability. |
| S04 Fowler_2010_21_years_BED.pdf | Source fixtures and later context | 15 pages. Tables 2, 4, 5; printed p562 EQD2/day narrative. Its later 59-63 zone should not be silently merged with the 2003 range. |
| S05 Gupta_2020_hypofractionated_HYPNO_perspective.pdf | Design rationale and historical series | 6 pages. Printed p2/PDF3 unpublished alpha/beta 6.4 overview; Table 1. Secondary description is not Stuschke's estimate; use Jacinto's original toxicity ordering. |
| S06 Hartley_2011_mucosal_parameters.pdf | DI preset and aggregate toxicity context | 5 pages. Table 1 p31/PDF3 and Figs 1-5. Elapsed onset conversion is appropriate; approximately 15% late observation is not universal. |
| S07 Jacinto_2018_hypofractionated_55Gy20fx_cisplatin.pdf | Clinical feasibility and delivered-time context | 9 pages. Methods/results and toxicity table, PDF2-5. Twenty patients; actual median 29 days, not ideal 26. |
| S08 Joiner_2023_LQ_fractionation_UofT_lecture.pdf | Educational endpoint-specific alpha/beta context | 28 slides. Slide 16/Table 10.1. Secondary teaching material, not primary validation of every quoted estimate. |
| S09 Kirkpatrick_2010_QUANTEC_spinal_cord.pdf | Explains why one late alpha/beta is insufficient | 8 pages. Table 1 printed S44/PDF3 and cord fit discussion. Cord alpha/beta about 0.87 with uncertainty; no cord NTCP module is implemented. |
| S10 Madhava_2006_carboplatin_hypofractionation.pdf | Toxicity and persistent-event context | 5 pages. Methods/results PDF2-4. Separate 16/19 any grade 3 from 3/19 prolonged; not a causal drug increment. |
| S11 Meade_2014_chemo_mucosal_tolerance.pdf | Multiple-model caution and future chemotherapy context | 7 pages. Methods/Table 4 PDF5; discussion PDF5-7. Phase I caution when models predict intolerability cannot be inverted to guarantee safety below limits. |
| S12 Miller_2026_hypofractionation_cohort.pdf | Contemporary retrospective burden and denominators | 8 pages. Methods PDF2; toxicity table PDF4. 72/113 at RT end, RTOG grading; 119 referred. Do not equate BED10 with physical dose. |
| S13 Moon_2024_HYHOPE_protocol_NCT04580446.pdf | Preset parameters, comparator, design and planning context | 47 pages. Amendment 4.2 dated 2024-01-30 on PDF1 and version 4.2 on PDF2; some body headers still 4.1. PDF15/printed p10 section 4.1.2 has BED assumptions. Cite the actual amendment. |
| S14 Moon_2025_HYHOPE_outcomes.pdf | Clinical phase I schedule and applicability | 9 pages. PDF2 eligibility/DLT; PDF3-4 delivery/toxicity; outcomes/discussion for LRR. Local printed pagination differs from indexed pages; use DOI plus PDF/table locator. |
| S15 PardoMontero_2021_mucosal_tolerance_model.pdf | Alternative mucosal-model context | 10 pages. Table I/PDF5 parameters; Table II cross-validation; Appendix A repair/equations. Its refit and performance are not those of the current Fowler defaults. |
| S16 PardoMontero_2021_SupportingInfo_108_schedules.xls | Supporting schedule dataset, currently not an app engine | Worksheet dataset2 contains 116 schedule blocks: 108 analyzed and eight extra partial courses. Has fraction times/doses and disputed tolerability labels. Requires reconciliation; not independent external validation. |
| S17 Stuschke_1999_HN_alpha_beta.pdf | Tumor and late fractionation sensitivity context | 9 pages. Table 1, p117/PDF5 tumor fits and 1.1-2.55 Gy range; PDF6 late slope-dependent fits. Extrapolation at larger fractions; not a jointly fitted app parameter set. |
| S18 vanLeeuwen_2018_alpha_beta_tumours.pdf | Heterogeneity of tumor alpha/beta estimates | 11 pages. Results PDF3, Fig 1/PDF4; 39 HNSCC estimates. No universal pooled value; negative fitted estimates are not admissible physical inputs. |
| S19 vanLeeuwen_2018_alpha_beta_tumours_AdditionalFile1.pdf | Underlying study/model details | 15 pages. Table S4/PDF7-8 and S6/PDF15. Supporting evidence, not patient-level calculator validation. |
| S20 Withers_1988_accelerated_repopulation.pdf | Direct compensation-rule source | 16 pages. Abstract p131 approximate 4 +/- 1 weeks; p140/PDF10 compensation/doubling assumptions. Alpha = 0.2 is an app completion assumption. |

The original tumor-model sources are S21 Shuryak 2018, equations 1-6, model-fitting section and trial table; and S22 Shuryak 2019, modeled fractionation optimization. Their equations and fitted-model context are necessary even though the supplied local References folder is primarily the later clinical/radiobiological collection. The current calculator does not implement the papers' complete fitted TCP/NTCP framework.

Keep unresolved historical source discrepancies out of model-performance claims. In particular, the EORTC/Horiot outcome discrepancy in earlier plausibility comparisons needs reconciliation to the original trial endpoint/timepoint before being used to score DI against DD. Also retain absent-supplement status for HYPNO Appendix E1/Figure E1, HYHOPE outcome Tables E1/E2 and Miller Supplemental Table 1 unless those exact files have since been acquired. They are not needed for the core corrections above, but should be obtained before importing detailed constraints or supplement-specific endpoints.

### Source identifiers

- S01 Bentzen et al., HYPNO, 2026: https://doi.org/10.1016/j.ijrobp.2025.12.056 (PMID 41711624).
- S02 Fenwick et al., 2008: https://doi.org/10.1016/j.ijrobp.2008.02.008 (PMID 18474317).
- S03 Fowler et al., 2003: https://doi.org/10.1016/S0167-8140(03)00231-7 (PMID 14643953).
- S04 Fowler, 2010: https://doi.org/10.1259/bjr/31372149 (PMID 20603408).
- S05 Gupta et al., 2020: https://doi.org/10.1016/j.oraloncology.2020.105045 (PMID 33091846).
- S06 Hartley et al., 2011: https://doi.org/10.1016/j.clon.2010.08.007 (PMID 20829006).
- S07 Jacinto et al., 2018: https://doi.org/10.1186/s12885-018-4893-5 (PMID 30352576).
- S08 Joiner lecture: https://radonc.utoronto.ca/media/2690/download (educational material; no established DOI).
- S09 Kirkpatrick et al., 2010: https://doi.org/10.1016/j.ijrobp.2009.04.095 (PMID 20171517).
- S10 Madhava et al., 2006: https://doi.org/10.1016/j.clon.2005.08.014 (PMID 16477924).
- S11 Meade et al., 2014: https://doi.org/10.1007/s13277-013-1267-0 (PMID 24104502).
- S12 Miller et al., 2026: https://doi.org/10.1016/j.clon.2026.104332 (PMID 42759386).
- S13 HYHOPE protocol NCT04580446: https://cdn.clinicaltrials.gov/large-docs/46/NCT04580446/Prot_SAP_000.pdf (local amendment version 4.2).
- S14 Moon et al., HYHOPE outcomes, 2025: https://doi.org/10.1016/j.ijrobp.2025.05.009 (PMID 40383163).
- S15-S16 Pardo-Montero et al., 2021 and supporting workbook: https://doi.org/10.1002/mp.14834 (PMID 33704792).
- S17 Stuschke and Thames, 1999: https://doi.org/10.1016/S0167-8140(99)00042-0 (PMID 10435801).
- S18-S19 van Leeuwen et al., 2018 and supplement: https://doi.org/10.1186/s13014-018-1040-z (PMID 29769103).
- S20 Withers et al., 1988: https://doi.org/10.3109/02841868809090333 (PMID 3390344).
- S21 Shuryak, Hall and Brenner, 2018: https://doi.org/10.1016/j.radonc.2018.02.015 (PMID 29534828).
- S22 Shuryak, Hall and Brenner, 2019: https://doi.org/10.1016/j.ijrobp.2019.02.025.

## 9 Ordered implementation plan with acceptance gates

The following order is designed to complete the identified corrections without another broad review. The detailed reproduction and repair instructions are in F01-F08; this sequence defines dependencies and when each work item is finished.

### Step 0 Preserve the reviewed state

Record the commit, source hashes, package versions and current test results. Keep this document, the evidence bundle and the existing defaults as baseline. Create new tests for the defects using the current implementation; observe that they fail because of the reproduced behavior rather than setup or package errors. Do not execute the archived grid scripts in place.

**Gate:** a clean snapshot reproduces the reported default rows and F01/F03/F05/F06 counterexamples. No unintended source edits or output overwrites occur.

### Step 1 Centralize endpoint validity and correct clinical claims

Implement F01 before any new aggregation/export work. Add the shared validity result and use it in main search, sensitivity, agreement and display. In parallel, apply the precise clinical text changes in F02 and preset/unit corrections in F06; these do not require refitting the equations or changing the user-selected defaults.

**Gate:** the 150-day counterexample shows 3 valid passes and 12 ineligible scenarios, not 15 confirmations; all-invalid cases have no favorable range; default arithmetic remains unchanged. No interface text asserts clinical NI, a quantified uncertainty range or a proven pessimistic bound from BED alone. The custom acute r = 5 cap is correctly labeled.

### Step 2 Separate calculation configurations

Implement F03 using independent candidate, agreement and evaluated-regimen configuration objects plus shared toxicity/grid settings. Decide and document the coupled-anchor target once; display the actual anchor. Remove optional Results context from candidate-export validation. Use one reason formatter with analysis-specific context.

**Gate:** independent tabs survive unrelated blank/invalid inputs; real shared inputs still affect both correctly. Fixed-anchor agreement is invariant to a sidebar comparator change when HYPNO is selected. Coupled-mode behavior is explicit and tested. Unavailable DD is explained on the agreement tab itself.

### Step 3 Define and implement the analysis export schema

Specify the schema before adding fields ad hoc. Include analysis type, effective scientific configuration, scenario definitions, per-scenario validity, selection policy and code/runtime identity. Implement F04 and F05 together so a new file can be both faithfully written and strictly replayed. Preserve a tested migration or explicit rejection path for prior candidate-only files.

**Gate:** fresh-process replay reproduces candidate and agreement rows, all qualifying rows and full scenario detail; unsupported specifications fail explicitly. Identity is stable under object-key reordering and independent of irrelevant UI labels. A change to an effective scientific input changes the analysis identity. Invalid optional context cannot block a valid analysis download.

### Step 4 Reconcile current sources and version metadata

Implement F07 using the source register above. Correct the HYHOPE protocol amendment, source-access status, Fowler units, g status, preset origins and current package/download versions. Update `docs/FIXES.md` with before/after behavior and the regression proving each fix. Preserve old tables and errata as historical records with their original model snapshot.

**Gate:** every current numeric input or clinical claim has an identifiable source or explicit analyst/inherited origin. No current note contradicts the running model. All visible version identifiers agree while schema/equation/selection versions retain separate meanings.

### Step 5 Complete exact-release verification

Implement F08: reproducible restore, tests, fresh-process replay, independent equations/grids and browser walkthrough. Review the final diff to ensure that no default, cap, hinge or onset law changed accidentally while fixing presentation/validation. Document expected output changes: valid default rows should stay the same; invalid-scenario counts, unsupported replay and erroneous unit labels should change.

**Gate:** all required checks pass on the same final revision, with warnings explained. Browser controls, rendered equations, large tables, plots and both download paths are checked. Source and exported model identity match the tested artifact. This is the completion criterion for the corrective release; increasing the test count alone is not.

### Step 6 Optional improvements after correction

Add a concise clinical-context record and a richer deterministic sensitivity panel, including DD rate and anchor perturbations. If a robust-search mode is desired, evaluate all feasible grid points across the prespecified scenarios before selecting one point per fraction count. Keep the existing current-parameter selection mode explicitly named; no omitted superior-robustness alternative was established for the default case, which has only three feasible points.

Create a portable archived-analysis reconstruction only if those historical analyses will be reused. Separate it from live-app correction. New scientific modules such as fitted TCP/NTCP, repair-aware schedules, chemotherapy effects or actual OAR dose handling require their own specification and validation; they are not hidden prerequisites for closing F01-F08.

## 10 Verification fixtures and evidence use

The accompanying evidence archive contains portable R/Python review scripts, their outputs, source identity and test logs. It excludes the full app snapshot and publications. It provides reproducibility of this review; it is not a replacement for adding permanent regression tests to the repository. Existing probes intentionally report the current defects, so a future fix should change their diagnostic output in the directions specified above.

Use the evidence README for the exact script invocations and the app-root argument. Run them from a scratch directory, with the application root pointing to the version being checked. The scripts must not be run by substituting current code into an archived scientific analysis without declaring that model change.

Recommended numerical tolerances are 1e-8 Gy for the printed application fixtures and 1e-9 to 1e-8 Gy for independent finite-domain calculations, adjusted only with a numerical reason. Published one-decimal source examples should be checked at their rounding precision. Compare row identity, validity and membership separately from floating-point closeness.

| Check | Required expected result | Why it matters |
| --- | --- | --- |
| Fowler 70/35, inclusive T47 | Acute 53.10543995 Gy10 | Prevents restoring the one-day clock defect. |
| Fowler 22 x 2.75, inclusive T30 | Acute 59.70979946; fails cap 59 | Tests a decision changed by the clock fix. |
| Fowler 2010 source examples | 60 x 1.3, elapsed T39: 62.79061740; T42: 60.41411278 | Independent source-equation fixtures; not authorization for a multidaily search. |
| Fenwick boundary | T47: 98.57419; T26: 74.16888; T19: 70.20307 Gy10 | Pins constants, shift, radians and inclusive clock. T > 50 unavailable. |
| Withers at tumor r5 | Legacy alpha0.2, lambda0.168, Tk28: reference DI82.04 Gy | Pins EQD2-to-BED conversion. |
| Anchor equality | Legacy tumor r5: DI = DD = 79.00 Gy at 35/2/47 | Pins nondefault-alpha/beta calibration. |
| Background growth | g0.0154, alpha0.2 lowers T47 BED by3.619 and T26 by2.002 in both models | Pins whole-course g and zero-AR branch behavior. |
| Late-model reversal | 20 x 2.9 versus70/35: lower at r3, higher at r2 | Prevents overgeneralizing one tissue parameter. |
| Default agreement | Exact three rows and deltas in section5; 15 valid scenarios each | Protects correct new feature arithmetic. |
| Invalid sensitivity counterexample | 32 x 2.10/44: 3 valid passes,12 ineligible of15 | Prevents false agreement from negative reference BED. |
| Acute units at r_e5 | HYPNO cap70.99097229 BED Gy at r_e5 | Prevents relabeling as Gy10. |
| Replay incompatibility | Unsupported schema/equation/hinge/clock rejected | Prevents silent scientific reinterpretation. |

The negative-reference policy is an explicit model-domain safeguard. Positive BED alone is not sufficient evidence of clinical plausibility. Very large fractions, unusual times or out-of-population use can still produce finite positive calculations while remaining extrapolations. Keep those two kinds of limitation distinct.

### Snapshot identity

```text
Commit: ad16504284ed4955b211e981b7d48cc0e31e0bcf
App version: 2.4.0
Equation version: 2.2

SHA256 app.R
3448906B720165159E4FCC8E5DADE8DBEDA8424E31400495093C5B6D78C30C44
SHA256 R/bed_models.R
BB3B29695381C552C37CC6DF50DD78455EC6E43B475719BDEAA70E1CDB69299C
SHA256 R/export.R
FC2EFBC237FB8C4E010A146868460CE1D5588C8E6D3E7F3E22A04575DB7B3EB7
SHA256 about.md
F1AC0A592AB09A7221EC3673AFF6919CD48D6DF195B549B7B6B9FF22AE7EC5B8
```

The existing candidate export implementation also passed eight additional ZIP round trips: current defaults, zero accelerated repopulation, unavailable DD with DI permitted, nonzero g, an empty result, Fowler, Fenwick and HYPNO acute/late caps. Those successful checks should be retained. They do not negate the separate replay-metadata and missing-agreement-export findings.

This plan is complete for the defects documented at the reviewed snapshot. Closing the specified acceptance gates is sufficient to verify those repairs without repeating the broad audit. A materially changed equation, parameter-fitting method, supported treatment domain or clinical outcome claim would require targeted scientific validation of that new scope.
