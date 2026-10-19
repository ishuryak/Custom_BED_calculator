# HNSCC BED calculator v2.5.1 review and action plan

19 October 2026 | Mathematics, source use, clinical interpretation and all application tabs

**Project:** `C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator`

**Reviewed application:** `Custom_BED_calculator_v2`, commit `451d4c511df3c5370ad0afc7623078a48a9fd761`, app 2.5.1, equation 2.2, selection 2.5, scenario set 1. File and line references refer to this snapshot. This review supersedes earlier reviews for current findings; the earlier documents remain useful history. No calculator code was changed by this review.

## 1 Main conclusions

The DI, DD, acute mucosal, Fenwick and late BED equations remain correctly implemented for the model the app declares. Independent calculations reproduce the current defaults and the main search grids. The recent fixes substantially improve endpoint validity, separation of tab inputs, exports, provenance and clinical wording. The app is useful for exploring how assumptions affect proposed schedules, provided its outputs are treated as model screens and the defects below are corrected.

The most important new computational defect is in ranking. A candidate whose DD BED is non-positive can now lose its DD delta, retain a positive DI delta, and win the supposedly conservative “smaller change” ranking. The code correctly rejects that endpoint for qualification, but does not consistently reject it for ranking. This is a regression introduced by the new validity handling; it does not affect the three default candidate rows or the default agreement rows.

The most important scientific interpretation issue is the new Stuschke/HYPNO illustration. Its BED differences and source slope numbers are correctly transcribed, but applying a normalized dose-response slope directly to repopulation-subtracted EQD2 requires an additional assumption or a coordinate transformation. The illustration should remain conditional, as the PI requested, with that assumption and the exact curve made explicit. It does not establish that D18 DD is disproved or that D18 DI is validated.

Further fixes are needed in replay validation, historical migration, agreement range definitions, comparator/anchor explanations, and source records. The ordered plan in section 9 gives implementation and acceptance criteria without requiring another review to discover the intended behavior. Preserve the authorized D18 defaults, HYPNO acute policy, fixed anchor, BED-based agreement allowance and selected 15-scenario set. None of the findings requires silently replacing those choices or adding a new outcome model.

**Release status:** suitable for exploratory investigation with the identified limitations; not yet ready to serve as a fully verified, reproducible trial-design screen. Correct R1-R7, complete R9's environment work and pass R8's browser checks before freezing a scientific result from the next version. A passing calculation or scenario count is not clinical validation.

## 2 Review scope and verification limits

The review covers all active mathematical functions, the entire UI/server implementation, export/replay code, tests, current source documentation, recent fix records and the 19 October development plan. Separate mathematical, clinical and software audits were reconciled. The local collection has 19 PDFs and one supporting XLS workbook. Prior complete readings were retained, and passages supporting current or newly changed claims were rechecked. The source inventory in section 8 distinguishes implementation sources, contextual papers and unavailable supplements; this is not a claim that every historical citation in every paper has been independently obtained.

The unchanged core was checked with an independent numerical calibration solver and independent grid enumeration. Additional probes targeted the new validity paths, malformed replay fields, old-schema replay, every selector option and numeric-input boundaries. Tests used copied snapshots; source hashes tie the results to the active project. The app starts locally.

**Current browser coverage is incomplete.** The in-app browser timed out connecting to the local server. Windows Computer Use then stopped because it could not determine the browser URL reliably enough to enforce its safeguards. No further UI automation was attempted. Server/reactive tests and code inspection cover the menu behavior below, but do not establish clickability, visual clipping, MathJax rendering, responsive layout or actual browser download behavior on this commit. The repository's earlier 23-check browser record is historical evidence for an earlier commit, not a current walkthrough. Section 10 lists the exact remaining checks.

**Test environment qualification:** the unmodified full suite encountered five source-hashing errors because Windows R's readability precheck rejected files that could in fact be opened. A process-local diagnostic workaround used actual binary reads for the digest operation without editing application source. Both raw and supplemental results are retained. A clean-host unmodified run is still required; this environment issue must not be reported as either a proven numerical defect or an unconditional passing release suite.

## 3 Recent fixes that are confirmed

| Earlier concern | Current behavior | Assessment |
| --- | --- | --- |
| Invalid sensitivity scenarios counted as successes | Shared endpoint validity excludes non-positive reference or candidate BED and unavailable DD. The 150-day stress case now reports 3 of 3 evaluable, with 12 of 15 not evaluable. | Corrected counts; range definition still needs R5. |
| Results and benchmark deltas bypassed validity | Invalid tumor differences are labeled not evaluable instead of showing a misleading positive gain. | Corrected. Raw values remain available. |
| Agreement depended on unrelated candidate inputs | Shared toxicity settings are separated from candidate-only gain, days and fraction range. | Corrected. Preserve independence. |
| Export required a valid Results regimen | Evaluated-regimen context is optional, so an unrelated invalid Results input does not block a search export. | Corrected. |
| Agreement could not be exported | It has its own ZIP with scenario, listed and qualifying rows, and candidate-by-scenario detail. | Implemented and round-trip tested. |
| Replay ignored model versions and key order | Schema 2 checks equation/selection/scenario versions, clocks and hinge; canonical key ordering stabilizes IDs. | Major improvement; R2-R3 identify remaining gaps. |
| Fixed Gy10 labels after changing acute alpha/beta | Main acute labels follow the selected alpha/beta; Results name the BED basis. | Corrected in primary outputs; optional consistency improvements remain. |
| Infinite DD onset displayed as Inf | No accelerated repopulation is displayed as “none.” | Corrected. |
| Clinical equivalence inferred from small BED differences | Active text distinguishes a BED allowance from clinical noninferiority and disclaims a confidence interpretation for scenario counts. | Substantially corrected; R4 and R7 remain. |
| Presets described as wholly published parameter sets | New tables distinguish direct, converted, chosen, inherited and analyst values. | Corrected main ledger; short Withers menu wording remains ambiguous. |
| Dataset count, protocol version and optional growth | The 2018 dataset is 35 regimens, 16 publications, 7283 patients; HYHOPE protocol is version 4.2; g is documented. | Corrected. |
| Comparator automatically changed DD calibration | Current UI has a fixed independent anchor; the coupling control was deliberately removed in 2.5.1. | Appropriate authorized design. Do not restore the checkbox as a purported fix. |

The 2.5 plan's main clinical and validity decisions should be retained. Its old sentences calling D18 DI “consistent” with HYPNO and the 0.4-Gy acute ordering reversal a demonstrated model “resolution limit” should be annotated if that plan remains current guidance: neither is established by the active calculations. Historical implementation records should not be silently rewritten as though the old checks were run on the new code.

## 4 Equations and derivations

### 4.1 Common definitions and clocks

Let n be integer fraction count, d the uniform prescription dose per fraction, r the relevant alpha/beta ratio, alpha the tumor radiosensitivity, and T the inclusive first-to-last treatment day. Rates lambda and g are per day; alpha is per Gy; BED and EQD2 are in Gy on their stated model scales.

```text
T(n) = 7 floor((n - 1)/5) + ((n - 1) mod 5) + 1
K(n,d,r) = n d (1 + d/r)
EQD2(B,r) = B / (1 + 2/r)
s(x) = log(1 + exp(k x))/k, k = 3 per day
```

The generated calendar assumes Monday start and once-daily Monday-Friday delivery. It gives 47, 26 and 19 days for 35, 20 and 15 fractions. Fowler acute calculations use elapsed time T_e = T - 1, converted once. Fenwick uses inclusive T. The HYPNO control is a specifically supplied six-fractions/week schedule, 33 x 2 Gy over 38 nominal inclusive days. It is not generated by the ordinary five-day calendar, which would give 45 days.

The soft threshold is an explicit variant of the sharp threshold in the source equations. It leaves a small positive pre-onset tail and has s(0) = log(2)/3 = 0.231049 days. Its excess over max(0,x) is at most that amount. Stable logarithm/exponential branches are used. Changing to a hard hinge would change equations and calibration and needs new versioned fixtures. Manually increasing T represents a longer aggregate course; the model does not simulate the position of gaps, weekends or individual fractions.

### 4.2 DI tumor BED

```text
B_DI = K - (lambda_DI/alpha) s(T - Tk) - (g/alpha) T
```

This is the log-survival form of the DI model in Shuryak et al. 2018, equation 1, with the declared soft hinge. The post-onset accelerated BED loss is lambda_DI/alpha per day. The EQD2/day loss divides that quantity by 1 + 2/r. Background growth applies to the whole course, including before accelerated repopulation starts.

Withers' approximately 0.6 Gy/day compensation in 2-Gy fractions converts to 0.6(1 + 2/r) BED Gy/day: 0.72 at r=10 and 0.84 at r=5. The current dynamic conversion is correct. Withers supplies the compensation rule and approximate four-week lag, not a jointly estimated alpha=0.2 parameter set. Similarly, HYPNO's design rationale uses 0.65 Gy EQD2/day, or 0.78 Gy10/day; these are not interchangeable units. Hartley's day 22 elapsed becomes day 23 inclusive. The HYHOPE protocol uses its own inclusive convention and its tumor onset of 21 does not need the same conversion.

The complete tumor preset contract is:

| Preset | Values filled by the preset | Source status |
| --- | --- | --- |
| D18 | alpha=0.35, r=10, lambda_DI=0.23, lambda_DD=0.40, Tk=28 | Analyst scenario; DI values fall within cited ranges but the joint set and DD rate are not fitted estimates. Default g=0 is a chosen omission of background growth. |
| Version 1 values | alpha=0.2, r=10, lambda_DI=0.2, lambda_DD=0.5, Tk=28 | Historical calculator choices, not a published jointly fitted HNSCC set. |
| Withers rule | alpha=0.2; current r; lambda_DI=alpha x 0.6(1+2/r); Tk=28 | Compensation rule from source; alpha chosen; onset representative of about 4 +/- 1 weeks. |
| Hartley | alpha=0.3, r=10, lambda_DI=ln2/3, Tk=23 inclusive | Source parameters with elapsed-to-inclusive onset conversion. |
| HYHOPE protocol | alpha=0.35, r=10, lambda_DI=ln2/3, Tk=21 inclusive | Protocol's modeling assumptions, not new fitted estimates from its 24 outcomes. |

Withers, Hartley and HYHOPE retain the entered DD maximum rate. Every preset retains entered g; selecting a named preset therefore does not necessarily restore the entire default scenario. Sensitivity overrides r to 7,10,15 and recalibrates DD at the same anchor for each row. D18's DI loss is 0.23/0.35=0.657143 BED Gy/day, or 0.547619 EQD2 Gy/day at r=10. Its quoted three-day doubling time is rounded: ln2/0.23=3.0137 days. The selection of D18 as a default was confirmed on 19 October, following its construction on 18 October.

### 4.3 DD tumor BED and fixed-anchor calibration

```text
p = 1 - exp(-alpha K/T)
q = C/alpha
t_onset = q T/K
B_DD = K - (lambda_DD/alpha) p s(T - t_onset) - (g/alpha) T
```

The onset and response laws follow Shuryak et al. 2018, equations 2-5. C is dimensionless log cell kill; q is in Gy; qT/K is in days. The response p is an average daily cell-kill response derived from the whole course, not a daily biological simulation.

For an anchor with K0,T0, matching accelerated losses between models gives:

```text
y = lambda_DI s(T0 - Tk) / (lambda_DD p0)
x0 = s_inverse(y) = log(expm1(k y))/k
t0 = T0 - x0
q = K0 t0/T0
```

The implemented inverse and units are correct. The common g term cancels. A finite positive t0 is required by the declared calibration domain. Positive DI with no positive DD rate, or no positive feasible anchor onset, makes DD unavailable. Zero DI deliberately disables accelerated repopulation in both models, using infinite onset and B=K-gT/alpha. Do not substitute zero for unavailable DD, force an onset floor, or tune an equation to conceal an infeasible calibration.

The app uses a **calibrated variant**, with shared alpha and equality imposed at an anchor. The 2018 paper fitted DI and DD independently to outcomes, with separate fitted DD parameters. The 2019 paper addresses additional fitted-model/survival implications; it does not validate this app's shared-alpha anchor construction. The current distinction is appropriate. Anchor equality is an identity, not biological corroboration. The separately quoted fitted-2018 contrasts reproduce approximately +3.1687 Gy10 for GORTEC, +2.0607 for TROG and +0.5375 for HYPNO under that fitted construction; they must not be presented as D18 outputs.

### 4.4 Acute mucosal BED and four alternative policies

```text
B_acute = n d (1 + d/r_e) - (lambda_e/alpha_e) s(T_e - Tk_e)
lambda_e = log(2)/Tp_e
Defaults: r_e=10, alpha_e=0.35, Tp_e=2.5 days, Tk_e=7 elapsed days
```

These Fowler 2003 assumptions are implemented correctly. The four selectable acute policies have different meanings:

| Policy | Implemented rule | Correct interpretation |
| --- | --- | --- |
| Relative to comparison | Candidate acute BED no greater than this analysis's comparison plus entered excess | Matched model comparison; assumes a clinically comparable treatment context. |
| HYPNO | Candidate acute BED no greater than 20 x 2.75 Gy/26 days, recomputed with current acute parameters | Trial-referenced model cap, not a fitted probability or guaranteed tolerance threshold. |
| Fowler | Lower edge 59 Gy10 without chemotherapy or 54 Gy10 with it | Chosen conservative edge of historical zones; requires the original acute parameters. |
| Fenwick | Cell-kill BED at r=10 no greater than the time-dependent boundary | Separate empirical RT-alone boundary, with no Fowler repopulation subtraction. |

Fowler 2003's RT-alone zone is 59-61 Gy10; the suggested chemotherapy decrement is possibly 3-5 Gy10. The displayed 54-58 range is constructed from those edges, and the search uses 54. Fowler 2010 later describes 59-63. These are not estimated probabilities. Blocking the fixed Fowler cap after changing its acute parameters is correct. The chemotherapy checkbox changes the Fowler zone/policy only; it changes no tumor, acute or late BED and has no effect on the default HYPNO cap.

### 4.5 Fenwick boundary

```text
x = T/32.2
L(T) = 69.5 x/sin(x) - 3.5 Gy10, for 0 < T <= 50 inclusive days
Candidate cell-kill BED = n d (1 + d/10)
```

The constants, radians, subtraction, time domain and alpha/beta=10 are correct for the app's declared once-daily use of Fenwick 2008 equations 10-11. Neither the Fowler clock conversion nor its repopulation term belongs here. The original analysis includes repair/delivery assumptions not fully modeled by this scalar app; do not extend it to twice-daily or interrupted schedules by changing only T. Historical schedules falling below a boundary do not prove safety for a new chemoradiation regimen.

### 4.6 Late effects, validity and search logic

Late BED is n d(1+d/r_l), without repopulation; its EQD2 conversion is correct. The conventional r_l=3 and sensitivity at 2 and 4 are useful prescription-dose summaries, not estimates of spinal cord, swallowing, laryngeal or mandibular injury. Organ doses, volume effects, prior RT and endpoint-specific biology are absent.

The shared endpoint policy uses a model delta only when the model is available and both reference and candidate BED are finite and positive. This is a declared interpretability/domain policy, not a statement that a negative model expression cannot be calculated. Raw values and reasons should remain available. The new policy is correctly applied to qualification and scenario counts but incompletely to ranking (R1).

The grid correctly searches multiples of the dose step strictly above the comparison dose, shorter aggregate calendars, chosen toxicity limits and minimum tumor change. It lists one highest-ranked dose per fraction count and exports all qualifying rows. “Both” qualification means both models under the current parameter set, not all 15 alternative scenarios. Agreement applies the same two-model criterion, then summarizes sensitivity afterward. Its -3 BED-Gy allowance is applied at each scenario's own tumor alpha/beta: equivalent allowances are about -2.33, -2.50 and -2.65 EQD2-Gy at r=7,10,15. This authorized choice is coherent if labeled; it is not a fixed clinical margin or uncertainty interval.

## 5 New and remaining findings

### R1 High priority Correct candidate ranking when a requested endpoint is invalid

**Files:** `R/bed_models.R:238-281`, especially candidate filtering at `:268-270` and scores at `:276-281`; related next-dose explanation later in `suggest_regimens`.

With the new validity rule, an invalid DD candidate gets an NA delta. `pmin(d_di,d_dd,na.rm=TRUE)` then treats the remaining DI delta as the conservative score. The early reference-level guards do not detect candidate-level invalidity. Ranking by DD can also list a row with an NA score when it qualifies under DI.

**Minimal UI-range reproduction:** alpha=0.2, tumor r=10, lambda_DI=0.4, lambda_DD=5, Tk=28, g=0; reference and anchor 35 x 2 Gy/47 days; default acute/late parameters; both toxicity policies relative to reference with +30 allowances; criterion either, rank smaller change, minimum gain 0, days saved 1; n=32 only; step 0.05, maximum d=2.5.

The app selects 32 x 2.50 Gy/44 days: DI BED=68, DD BED=-0.60436688, DD delta=NA, score=+22 from DI alone. The valid two-model optimum is 32 x 2.10: DI BED=49.312, DD BED=54.76422319, conservative score=+3.312. With a 0.5 dose step and DD ranking, the only 32 x 2.50 candidate is still listed with score=NA. These broad exploratory inputs are allowed by the UI; the defect need not occur at the default to be real.

**Fix:** keep qualification and ranking eligibility explicit. Smaller-change ranking requires both endpoint deltas to be finite. DI and DD ranking require the selected endpoint. Larger-change ranking may use an available endpoint if that documented fallback is intended. Apply eligibility before selecting one row per fraction count. Compute a two-model minimum only when both deltas are valid; otherwise use NA or a distinctly named available-model statistic. Decide whether the all-qualifying export includes unrankable rows, then expose a separate `rank_eligible` flag if it does. Do not silently alter what “qualifying” means between CSV and display.

Use the same eligibility in next-dose explanations, distinguishing failure of the tumor qualification rule from an unavailable ranking endpoint. If no rankable row remains, return an explicit reason. Increment the selection version because output membership/ranking changes; the BED equation version can remain unchanged.

**Acceptance:** the minimal case selects 2.10 for smaller-change ranking; the DD-only ranking case returns no rankable row, never a row with NA score. Exercise all 16 criterion/rank combinations with invalid DI only, invalid DD only, unavailable DD, both invalid, and fully valid endpoints. Preserve default candidate/agreement outputs, and test the next-dose message when the next step first becomes non-evaluable.

### R2 Medium priority Reject malformed replay fields instead of changing their meaning

**Files:** `R/export.R:validate_scenario`, particularly the `fenwick` conversion around `:187` and agreement/search validation.

`isTRUE(z$fenwick)` silently maps the string “true,” numeric 1, NULL and NA to FALSE. A genuine Fenwick search that produces 12 rows can therefore replay as a different acute policy producing zero rows, without a validation error. Generated schema-2 exports round-trip, but the public replay contract also needs to reject malformed or manually edited files.

Agreement files redundantly carry search criterion, rank and minimum gain as well as the agreement margin. Contradictory values are accepted and hashed, then ignored by the agreement runner. This produces different IDs for inputs that do not describe the calculation actually performed and undermines the claim that validated saved settings are authoritative.

**Fix:** validate `fenwick` as one nonmissing logical value; require it in schema 2. Default missing historical fields only in an explicit schema-1 migration. Validate all enumerated strings and numeric shapes before coercion. For agreement, either remove redundant fields in a versioned schema or enforce criterion=both, rank=worst and min_gain=agreement margin. Document which descriptive context fields do not affect computation. Add a practical upper bound/deduplication rule for replay fraction grids and a total grid-size check so edited files cannot cause uncontrolled enumeration.

Apply the same authoritative-versus-derived distinction to saved calibration and margin fields: corrupted `onset_anchor=999` is currently accepted if C/alpha agrees, and derived EQD2 margins need not describe the computation unless checked or recomputed. Validate every field presented as verified provenance, or explicitly regenerate and label it as derived. Capture source identity when the model is loaded, or detect disk changes before export, if live editing during a session is supported; hashing files only at export time can otherwise identify code different from that loaded in memory.

**Acceptance:** TRUE and FALSE preserve the intended policy; each string/numeric/null/vector invalid form yields a named actionable error before computation. Every contradictory agreement field is rejected or normalized by an explicit migration that records the original. Reordered object keys preserve IDs; changed scientific settings change IDs. Existing valid candidate and agreement round trips remain unchanged.

### R3 Medium priority Separate historical replay from migration to current selection rules

**Files:** `R/export.R:validate_scenario` schema-1 branch and `rerun_scenario`; `README.md` replay description; `about.md:62`.

A genuine v2.4 schema-1 fixture for the R1 setting selects 2.10 Gy in its original engine. Current migration silently assigns selection version 2.5 and selects 2.50 Gy, with only a `migrated_from='1'` attribute. The ranking defect makes the numerical change conspicuous, but the contract problem remains after that defect is fixed: historical schemas did not specify all current validity/selection semantics. “Replay” and “reanalyze with new rules” are different operations.

**Fix:** provide separate explicit behavior: strict replay requires a supported recorded engine/selection version; migration creates a new scenario with original schema/app/equation/selection provenance, target versions and a visible statement that results may change. If no historical engine is maintained, reject strict historical replay and give a command for deliberate migration. Do not pretend that a schema conversion proves equality of results. Preserve the original file and ID, create a new migrated ID, and retain a small change report when previous results are present. State exactly which equation-2.1 non-Fenwick compatibility path is supported.

**Acceptance:** freeze actual historical schema-1 fixtures and their old results. A strict call either reproduces those results or refuses with an explicit version message; migration never masquerades as exact replay. Use the 2.10-versus-2.50 fixture as a regression even after R1 is repaired, and include older invalid-scenario behavior. Invalid future versions, hinges and clocks continue to be rejected.

### R4 Medium to high priority Derive the conditional dose-response illustration on a declared coordinate

**Files:** `about.md:94,99`; `References/README.md:267`; development plan section N1.

The source numbers are real: the full Stuschke paper's logit estimate is gamma50=3.2 (95% CI 1.6-4.8), and it discusses RTOG 83-13 at 0.8 +/- 0.4. The abstract's 3.1 (1.5-4.7) is a different reported fit, not proof that the full-text value was misquoted. Current D18 HYPNO contrasts are also correct: -2.5041 Gy10 DI and -5.7930 Gy10 DD. Their percentages use the respective time-corrected control BEDs 72.6286 and 68.6403, giving 3.4% and 8.4%; EQD2 conversion cancels from each ratio.

However, Stuschke Appendix A, printed page 119, equations 1-3 defines gamma50 as D times dP/dD at TCD50, with D physical total dose or isoeffective conventional-fractionation dose. It does not define the slope on the calculator's repopulation-subtracted EQD2 coordinate. Normalized slopes are not invariant under an additive time correction. For E=f(D), the transformation at the corresponding midpoint is:

```text
gamma_E = [E50 / (D50 * f_prime(D50))] * gamma_D
```

For DI at fixed time and 2-Gy fractions, E is physical dose minus a fixed time penalty, so f_prime=1 but E50 differs from D50. For DD, the repopulation loss also changes with dose, so its derivative is required. Different-time comparisons add further assumptions about transporting the dose-response relation. There is no uniquely correct replacement probability derivable from this app alone.

The HYPNO HR conversion is conditional too. With proportional hazards and fixed control probability 0.512, the upper HR 1.30 gives `0.512 - 0.512^1.30 = 0.093157`, or 9.32 percentage points. This is not an independently estimated confidence bound for the absolute difference. The current gamma-times-fractional-dose calculation is a local linear approximation; predicting 27 points is large enough to show the curve rather than assume local linearity.

Under the **extra hypothetical assumption** that the same gamma applies directly on effective EQD2, and an exact logit curve fitted through control probability 0.512, gamma 1.6-3.2 gives about 13.36-25.01 points loss for DD and 5.55-10.95 for DI; the current linear calculation gives 13.50-27.01 and 5.52-11.03. These figures illustrate approximation effects; they do not fix the unvalidated scale transfer. At the central RTOG gamma=0.8 both losses are within the illustrative PH bound, but at its quoted upper value 1.2 the DD loss is about 10.13 points, exceeding 9.32. “Both compatible” therefore must not describe that whole quoted range.

**Fix preserving the approved illustration:** either explicitly assume gamma transfers unchanged to effective EQD2, or derive a stated physical/isoeffective-dose mapping and its Jacobian. State baseline endpoint, timepoint, curve, dose coordinate, source population, chemotherapy transfer and HR assumptions next to the result. Use exact sigmoid values alongside or instead of the linear approximation. Replace “more than HYPNO's upper limit allows” with a conditional comparison under the chosen assumptions. Do not call D18 DI validated or D18 DD a bound. Record this scientific derivation in a short source-linked script, not only narrative.

**Acceptance:** independent fixtures reproduce raw BEDs, denominators, PH conversion, selected curve and coordinate transformation. Include central and interval slope values. Link Stuschke pages 114, 115/Table 1, 117 and 119, and HYPNO pages 391-392 and Table 3. Preserve the D18 defaults and the PI's decision to retain a conditional illustration.

### R5 Medium priority Make agreement ranges use the documented scenario population

**Files:** `R/bed_models.R:scenario_summary`, `about.md:75`, agreement/export descriptions.

Counts correctly define an evaluable scenario as valid under both models. Ranges instead use each model's own validity mask. In the 150-day stress case there are only three jointly evaluable scenarios, but the DD range -0.21437 to +2.38883 uses all 15 DD-valid scenarios. The DD range among the three jointly evaluable scenarios is +0.98809 to +2.38883. The existing phrase “ranges over evaluable scenarios” makes these different denominators appear identical.

**Fix:** the least ambiguous implementation of the recorded A3 decision is to calculate both ranges using the joint-validity mask. Alternatively, retain per-model ranges but label them explicitly, show n_valid_di and n_valid_dd, and distinguish them from the joint pass denominator in UI and CSV. Do not change pass counts or discard rows merely because some scenarios are non-evaluable; the PI explicitly chose to retain such rows and disclose invalid counts.

**Acceptance:** the stress case has either the joint-valid DD range +0.98809 to +2.38883, or a clearly labeled DD-only range with denominator 15 beside the jointly evaluable denominator 3. All-invalid produces NA ranges. Default 15/15 rows stay unchanged. Scenario detail must let a reader reconstruct every displayed range and count exactly.

### R6 Medium priority Correct comparator scope and historical anchor instructions

**Files:** `app.R:123-124,382-395`; `about.md:32`; `README.md:42,73`; v2.5.1 interface audit.

The sidebar says all differences and toxicity limits are relative to its comparison. Agreement defaults to its own HYPNO comparison, while HYPNO, Fowler and Fenwick acute limits are fixed/reference-specific policies. The agreement note itself is accurate, but the global instruction is not. A user can otherwise interpret a delta or late limit against the wrong clinical schedule.

The current text also says versions 1 and 2.0 calibrated at the selected comparison. Inspection of both local v1 sources shows a hard-coded 35 x 2 Gy/47-day anchor and no arbitrary comparator coupling; v1 also had the alpha/beta calibration error. The comparison-coupled history applies to v2.0, not the original v1. This corrects the historical characterization repeated in the previous external review as well as current app text.

Even for v2.0-style reconstruction, entering n and d in the current anchor cannot reproduce every comparator: 33 x 2 becomes 45 days, not HYPNO control's 38 days; an interrupted comparator's T cannot be represented either. Matching n,d alone is insufficient. Reintroducing the removed coupling switch is not necessary.

**Fix:** name the analyses affected by the sidebar comparison and distinguish relative from fixed policies. Show each analysis's actual comparator, anchor time and effective numerical toxicity limits. Correct v1/v2.0 history. Either add an advanced explicit anchor-time control with an explicitly supported calendar policy, or narrow reproduction instructions to schedules the UI can actually represent. A 38-day six-day anchor must not be blocked by a five-day validation rule if that advanced capability is chosen. Saved scenarios already carry anchor T; preserve that exact value in replay.

**Acceptance:** with sidebar HYPNO control and default agreement HYPNO test, each tab names its real comparator and limits. Changing only the sidebar comparison leaves fixed-anchor DD BED for a fixed schedule unchanged. A 33-fraction anchor visibly states 45 days unless an explicit supported 38-day anchor was entered. Historical reproduction examples identify the exact app/equation version and limitations.

### R7 Medium priority Finish active wording and source provenance cleanup

**Files:** `about.md:66,96-112`; `app.R:366-367`; `References/README.md` original entries and sections F/G; About version footer.

The default-search paragraph first reports DI gains, then says shorter courses at HYPNO's acute BED lose tumor BED under both models. Its own displayed candidates refute that general statement. Separate “DI-only gains versus 70/35” from “both-model agreement versus HYPNO under a chosen loss allowance.” A negative allowance is not required for every possible hypofractionation exploration.

The short Withers menu note still reads as though alpha=0.2 came from Withers. State that the paper supplies the compensation rule/onset and the app chooses alpha to complete the preset. The claim that editing any value switches to Custom is too broad: inherited DD/g values do not define every named DI preset. Say which values the preset controls.

Several original catalogue entries still say full text was not downloaded despite the appended inventory: Fowler 2010, Hartley, Stuschke, HYPNO and the section containing Meade. Reconcile these with one current availability/read-status field while retaining acquisition history as history. Replace the HYPNO LRC shorthand with its exact -0.5-point test-minus-control estimate. Mark Meade's 4.4-Gy10 cisplatin estimate as study/model-specific rather than a universal measured drug conversion. The “every numeric input” ledger should include chosen search settings or narrow its stated scope. Make the About footer show 2.5.1 through the same version source as other outputs.

**Acceptance:** menus, About, reference ledger and exported readme assign the same source/chosen/inherited status. No current “not downloaded” entry contradicts a present PDF. Missing supplements remain marked missing. Default narrative matches actual table signs and comparators. These changes need no alteration of default numerical results.

### R8 Verification gate Finish current browser and clean-environment checks

The app's mathematical and server evidence is substantial, but current browser behavior is unverified for the reasons in section 2. This is a verification gap, not evidence of a visual defect. Keep it open until the exact post-fix version passes the walkthrough in section 10. Separately run the unmodified suite and both downloads on a clean supported R environment. The diagnostic hashing workaround is evidence about this sandbox, not a release change to ship without investigation.

Acceptance requires a version/commit-stamped browser record, working downloads with reproducible CSVs, readable equations and tables, and zero unexplained failures in the unmodified clean-host suite. A historical screenshot or test record is insufficient for changed controls.

### R9 Medium priority Make dependency restoration an exercised release path

**Files:** `renv.lock`, `DESCRIPTION`, `run_app.bat`, `.github/workflows/tests.yml`.

The new lockfile and more portable launcher improve the earlier release setup. However, the lockfile contains 55 packages and omits testthat and callr, which are declared test dependencies. The project has no renv activation/bootstrap files, the launcher runs Shiny without explicitly selecting a restored project library, and CI uses setup-r-dependencies rather than restoring the lockfile. There is no evidence yet that the locked environment can reproduce the release tests from a clean machine. This is partial closure of the previous release finding, not a reason to remove the lockfile.

**Fix:** choose one documented workflow for restoring and activating the project library. Include test dependencies in its reproducible specification or a separately pinned test environment. Make the launcher/CI use that environment explicitly and record R/package versions. A dependency installer based on DESCRIPTION can remain useful, but should not be described as verifying the lockfile. See the official renv and r-lib action documentation linked in section 8.

**Acceptance:** on a clean supported Windows environment, restore from the declared specification, start the app, run the unmodified suite and generate/replay both exports without the review harness workaround. Record exact commit, R version, dependency versions and exit results. The current review did not perform a network-dependent clean restore.

## 6 Review of every menu and tab

The following matrix covers every current control group and every choice within it. “Verified” here means code, numerical or server/reactive verification; browser-only behavior remains subject to section 10. The detailed probe CSV in the evidence package records numeric boundary/cleared-input checks and selector combinations.

### 6.1 Sidebar controls

| Menu and controls | Intended and observed behavior | Clinical usefulness and action |
| --- | --- | --- |
| Regimen to evaluate: n, d, automatic time, manual T | Drives Results and Sensitivity, not candidate generation. Automatic time uses five weekdays from Monday; manual time cannot be shorter than that calendar. Values, range and integer n are checked on the server. | Useful single-schedule inspection. State that extra T is aggregate prolongation, not gap-placement modeling. Add a candidate-row inspection workflow later. |
| Comparison: Custom, HYPNO test, HYPNO control | Custom n,d and auto/manual T drive candidate and Results comparisons; presets supply 20/2.75/26 or 33/2/38. Candidate fraction range resets below its comparison n. | Correct actual trial arms are available. R6 fixes the overbroad help. The six-day control is a supplied benchmark exception, not general support for six-day candidate schedules. |
| Tumor preset: Version 1, Withers, Hartley, HYHOPE, D18, Custom | Each named preset sets its declared values; DI-only presets retain sidebar DD rate and all retain g. Withers responds to tumor alpha/beta. Custom leaves current values editable. | Useful contrasting hypotheses; none is a joint modern-HNSCC fit. R7 clarifies inherited fields and chosen alpha. Preset selection/update rendering still needs live browser confirmation. |
| Tumor numeric fields: alpha, r, lambda_DI, lambda_DD, Tk, g | Change the corresponding equations and DD calibration; rate and alpha units are consistent. Zero DI disables accelerated repopulation in both; infeasible DD is disclosed. | Wide limits support stress exploration, not clinically plausible ranges. Display implied doubling time/BED loss and source range separately from allowed input range. R1 is essential for invalid candidate states. |
| DD calibration anchor: n, d | Independent of comparisons; T is always generated from a five-day calendar. Recalibrates C/alpha for each parameter set. | Useful way to expose an important model assumption. No manual anchor T exists; R6 narrows reproduction claims or adds an advanced explicit time. |
| Acute mucosa: r_e, alpha_e, Tp_e, Tk_e | Derives lambda_e=ln2/Tp_e; uses elapsed clock once. A custom parameter changes candidate and HYPNO reference cap consistently; it blocks fixed Fowler policy. | Appropriate mechanistic sensitivity. Show that changing these is a new model assumption, not a measured patient parameter. |
| Concurrent chemotherapy checkbox | Lowers the displayed Fowler zone and the Fowler search cap only. Does not change any BED, the HYPNO cap, Fenwick or late model. | Current explanatory text is clear. Label could be more direct: “Apply Fowler chemotherapy decrement.” Keep actual systemic regimen/cycles in clinical context. |
| Late-tissue r_l | Changes late BED/EQD2 and matched late caps; does not need alpha. | Useful endpoint sensitivity, but no organ-specific risk inference. Current prescription-dose warning is appropriate. |

### 6.2 Results tab

| Item | Behavior and assessment | Required or useful change |
| --- | --- | --- |
| DD and negative-BED messages | Calibration failures and negative model values are disclosed; invalid tumor deltas/verdicts are not credited. | Preserve raw values and reasons. For exactly zero BED the row's not-evaluable message remains important even though the general negative warning uses a strict negative test. |
| Results table and verbal comparison | Correct T, total dose, DI/DD BED, acute/late BED and EQD2; late sensitivity at 2 and 4 in addition to selected r. | Keep dynamic units. Verbal “higher/lower” means model value, not probability of harm/benefit. |
| Fowler zone and Fenwick statement | Display independently applicable policies; changed Fowler parameters are flagged. Fenwick uses r=10 cell kill and a time domain. | Useful cross-check, provided the two scales are not conflated. Show policy names and numerical margin to each boundary. |
| Published regimen table | HYPNO test/control and HYHOPE are recomputed under current app parameters, with deltas versus sidebar comparison. | Label explicitly “recomputed with current parameters.” Add benchmark population, volumes, chemotherapy, original design comparator and source drawer. |
| HYPNO/HYHOPE clinical note | Correctly distinguishes HYPNO's actual control, endpoints and HYHOPE's selected population/60 Gy comparator. | Useful context. An optional 60/30 comparator can help inspect HYHOPE but must not imply it is a standard for all HNSCC. |

### 6.3 Candidate regimens tab

| Item and all choices | Behavior and assessment | Required or useful change |
| --- | --- | --- |
| Tumor criterion: either, both, DI, DD | Qualifications use valid endpoint changes and entered minimum. “Either” admits disagreement by design; “both” is current-parameter agreement. | Preserve the visible credited-by column. Avoid describing either-model candidates as jointly supported. |
| Ranking: larger, smaller, DI, DD | Selects one dose for each n after qualification. Larger is optimistic; smaller is intended conservative. | R1 must enforce per-candidate endpoint availability and finite scores. All 16 criterion/rank combinations need edge-case tests. |
| Fraction range slider | Restricts candidate n below comparison n; changes when the comparison changes. | Browser verify slider bounds/reset. Explain that agreement searches its own full shorter range, independent of this slider. |
| Minimum gain and saved days | Defaults +1 BED Gy and 7 calendar days; negative gain permits modeled loss. | Add visits saved. The input is a model-screening allowance, not a clinical NI margin. |
| Dose step: 0.05, 0.1, 0.2, 0.25, 0.5; maximum dose | Enumerates absolute multiples strictly above reference d, with default max 5. | Correct. The grid's fine numerical resolution does not establish equally fine biological discrimination. Very high allowed d needs explicit extrapolation context. |
| Acute policy: comparison plus excess, Fowler, HYPNO, Fenwick | Each branch applies its own rule; excess matters only for comparison mode. Custom acute parameters block Fowler; Fenwick uses its own fixed r=10. | Correct source separation. Show effective cap and clock in the tab, especially after changing settings elsewhere. |
| Late policy: comparison plus excess, HYPNO | Relative cap follows the analysis comparison; fixed HYPNO cap follows current r_l. | Correct. Do not let sidebar wording imply HYPNO cap is 70/35-based. |
| Search/cap/empty messages | Describe actual grid and limits, with explicit unavailable or no-solution states. | R1 adds unrankable-state reason. Empty-result suggestions should not encourage relaxing a clinical limit without saying it changes the hypothesis. |
| Candidate table, sort/paging/scroll, fraction-size flag | One dose per n; delta signs, credited-by and next-dose reasons are informative. >3.1 means above these displayed benchmarks, not beyond all clinical precedent. | R1 fixes next-dose ranking validity. Browser sorting, paging and width remain unverified. Add row inspection and numerical toxicity slack later. |
| Plot | Plots DI/DD deltas for the listed dose at each n, omitting unavailable values. | Explain it is not a fixed-dose response curve: d changes by point and model-specific availability may differ. Browser legend/axes still need inspection. |
| ZIP download | Full-precision rows, qualifying grid, scenario, provenance and reasons; independent of valid Results input. | R2-R3 strengthen replay. Server exports passed probes; actual browser save is still an open acceptance item. |

### 6.4 Model agreement tab

| Item | Behavior and assessment | Required or useful change |
| --- | --- | --- |
| Comparison: HYPNO test or sidebar | Has its own effective comparator; fixed anchor remains shared. | Correct. Always display the actual comparator and numerical limits. |
| Allowed change and days saved | Defaults -3 BED Gy and 1 day; independent of candidate gain/days/range. | Correct and clinically interpretable as a screen. Retain explicit separation from NI endpoints. |
| Shared dose grid and toxicity | Uses candidate-tab step, maximum dose and acute/late policy; relative caps are computed against agreement's own comparison. | Works but settings are physically on another tab. Repeat their values or move shared controls into a common analysis panel. |
| Agreement table | Both current-model changes must qualify; sensitivity counts exclude non-evaluable scenarios; onset is “none” for zero AR. | R5 aligns ranges with their denominators. The tab is not a robust optimizer across all scenarios. |
| Scenario counts and BED allowance | Fifteen selected preset/r combinations, with inherited DD/g and shared anchor; counts are not probabilities. | Keep invalid counts visible. Optional current custom scenario should be a separate contextual row, not silently a sixteenth denominator. |
| Agreement ZIP | Contains listed rows, all qualifying rows and details for every listed row/scenario. | Useful new audit trail. R2 validates redundant fields; R3 separates migration. State that scenario_detail covers listed rows rather than every feasible grid row. |

### 6.5 Sensitivity tab

The tab correctly evaluates the **sidebar regimen**, not the currently highlighted candidate. It varies five presets across r=7,10,15; reconstructs/recalibrates DD with each scenario; preserves inherited g/DD where appropriate; and shows valid-only deltas and explicit failure reasons. Its summary threshold comes from candidate minimum gain, which is another cross-tab dependency that should be visible locally. Its toxicity table uses current acute parameters and late r=2,3,4, independent of tumor scenarios.

This is useful deterministic sensitivity, not parameter uncertainty propagation. Add a short context strip identifying the regimen, comparison, anchor and summary threshold. An optional CSV/ZIP export for standalone sensitivity would improve auditability; currently agreement exports detail for its own listed rows, which is not a substitute for every arbitrary sidebar sensitivity analysis. Do not count the scenario proportion as a probability or “15 independent confirmations.”

### 6.6 Model and sources tab

All equations and most provenance are accessible in the explanatory page, with source links and limitations. The scientific definitions in section 4 are correct; R4-R7 identify the remaining interpretation/history issues. MathJax, links, long-table readability and narrow-screen behavior need the browser checks in section 10. Add anchors to equation, parameter, benchmark and limitation sections so a user can reach the relevant explanation from a control. Display application/equation/selection versions together in an audit footer, using a single version source.

## 7 Clinical numbers and plausibility

### 7.1 Reproducible numerical anchors

All values below use current D18 tumor parameters, default Fowler acute parameters, r_l=3, fixed 35 x 2/47-day anchor and nominal schedule times. Rounding is for this report; evidence CSVs retain full precision.

| Schedule | Tumor DI and DD BED Gy10 | Acute BED Gy10 and late BED Gy3 |
| --- | --- | --- |
| 35 x 2 Gy, 47 days | 71.514286; 71.514286 | 53.105440; 116.666667 |
| HYPNO test 20 x 2.75 Gy, 26 days | 70.124458; 62.847301 | 55.865972; 105.416667 |
| HYPNO control 33 x 2 Gy, 38 days | 72.628571; 68.640339 | 55.434954; 110.000000 |
| HYHOPE 15 x 3.1 Gy, 19 days | 60.915000; 56.370530 | 52.201150; 94.550000 |
| 18 x 3 Gy, 24 days | 70.199999; 63.142893 | 57.525309; 108.000000 |

The 18 x 3 schedule from the 2019 optimization analysis is not an optimum of this app: under current D18 it loses 1.3143 DI and 8.3714 DD BED Gy10 versus 70/35 and exceeds the selected HYPNO acute cap by 1.6593. The 2019 paper used fitted tumor parameters, TCP/LNTCP objectives and different acute assumptions; its predicted outcomes are not observed clinical evidence and its Monte Carlo analysis is not the app's 15-scenario exercise.

The current candidate defaults produce three DI-only gains versus 70/35: 29 x 2.25/39 days (+1.1884 DI, -2.4080 DD), 28 x 2.30/38 (+1.1263, -2.8648), and 30 x 2.20/40 (+1.1200, -2.0495). Their acute BEDs 55.17-55.45 are below the selected HYPNO cap but above the 70/35 acute value. They are model-discordant hypotheses; claiming simultaneous preserved tumor effect would contradict the DD results.

The default agreement comparator is HYPNO 55/20. The three rows are:

| Candidate | DI and DD change in Gy10 | Practical interpretation |
| --- | --- | --- |
| 17 x 3.05 Gy, 23 days | -2.4602; -1.5895 | 3 fewer visits and 3 calendar days; acute BED 55.7817, just 0.0842 below the model cap. |
| 18 x 2.90 Gy, 24 days | -2.7865; -1.9614 | 2 fewer visits and 2 days; acute 54.6633. |
| 19 x 2.80 Gy, 25 days | -2.0285; -1.4771 | 1 fewer visit and 1 day; acute 54.6291. |

All 15 selected scenarios are evaluable and within the -3 allowance for these default rows. That establishes the chosen arithmetic screen, not efficacy equivalence. The logistical benefit is modest and must be weighed against uncertainty, systemic therapy delivery and normal-tissue planning. Differences of hundredths of a modeled acute BED Gy do not justify fine clinical ranking. The appropriate next clinical question is whether fewer visits offer enough benefit to motivate further evidence/planning, not which row is clinically “best” because its score is largest.

### 7.2 Clinical source values and what they support

HYPNO randomized 792 patients to 55/20 at five fractions/week versus 66/33 at six/week; 75% received weekly cisplatin 35 mg/m2, and HPV was not tested. Three-year locoregional control was 50.7% versus 51.2%, HR 1.098 (90% CI 0.93-1.30); grade 3+ late effects were 18.8% versus 20.2%, HR 0.926 (0.68-1.25). Both co-primary endpoints met the trial's 10-percentage-point noninferiority framework. Grade 3+ acute mucositis was 190/373 (50.9%) versus 208/380 (54.7%). These numbers and denominators are correctly used in the active app. The trial supports its tested comparison in its population; it does not supply a generic acute safety threshold or decide between this app's DI and DD assumptions. [HYPNO primary paper](https://doi.org/10.1016/j.ijrobp.2025.12.056).

HYHOPE was a phase 1 study of 24 favorable-risk HPV-positive oropharynx patients, reduced elective volumes and three weekly cisplatin doses of 40 mg/m2, designed against 60 Gy/30 fractions. The 46.5/15 level had grade 3 oral mucositis in 5/24, any acute grade 3 toxicity in 7/24, and two-year locoregional recurrence 13.8% (95% CI 3.2-32%); absence of late grade 3+ at median 26.6-month follow-up is not proof of long-term safety. The app now keeps this selected-population context. The 12-fraction level was dropped after signals in a parallel postoperative study, rather than a direct randomized comparison proving the acute BED formula wrong. Protocol acute calculations use inclusive time: about 51.3 for 12 fractions and 51.4 for 15; the latter becomes 52.2 on this app's elapsed clock. Preserve the explicit clock distinction. [Outcomes](https://doi.org/10.1016/j.ijrobp.2025.05.009), [protocol record](https://clinicaltrials.gov/study/NCT04580446).

The additional acute series are correctly contextualized: Jacinto 20 patients with cisplatin, 40% grade 3 mucositis and median actual T=29 days; Madhava 16/19 grade 3 with carboplatin, three lasting at least four weeks; Miller 72/113 completers at treatment end, mainly older patients treated without chemotherapy, with a different four-week assessment. These are not exchangeable cohorts from which to subtract rates and estimate a chemotherapy effect. [Jacinto](https://doi.org/10.1186/s12885-018-4893-5), [Madhava](https://doi.org/10.1016/j.clon.2005.08.014), [Miller](https://doi.org/10.1016/j.clon.2026.104332).

Tumor r=10 is a conventional exploratory choice with wide uncertainty. Stuschke's 10.5 Gy (6.5-29) estimate comes from hyperfractionation trials, not validation at every d permitted by the UI. The 6.4-Gy estimate mentioned in the HYPNO rationale is from an unpublished overview; it must not be promoted to an independently checked primary fit. The late 3.3-4.0 Gy estimate and broad limits concern a particular late endpoint and slope assumption. The cord's substantially lower estimate illustrates why one prescription-level late model cannot represent all OARs. Hartley's approximately 15% late-BED observation is an aggregate association in those trials, not a universal minimum improvement threshold. [Stuschke](https://doi.org/10.1016/S0167-8140(99)00042-0), [Hartley](https://doi.org/10.1016/j.clon.2010.08.007), [QUANTEC cord](https://doi.org/10.1016/j.ijrobp.2009.04.095).

### 7.3 Appropriate role in trial exploration

Use the app to define and compare explicit hypotheses for a specified population: site, stage, HPV status, gross/elective dose volumes, systemic regimen and actual calendar. Freeze a candidate scenario and compare it with its appropriate clinical reference. Inspect why DI and DD differ, how the anchor affects that difference, and whether source-supported parameter alternatives change membership. Then evaluate dose-volume feasibility and clinical precedent outside the scalar prescription model.

A protocol still needs endpoint-specific efficacy and toxicity assumptions, a defensible clinical noninferiority or phase 1 objective, statistical operating characteristics, treatment/chemotherapy delivery details, OAR constraints, follow-up and stopping rules. The calculator does not supply these. Adding a TCP/NTCP module would require a separate specification, source dataset, fitting/validation strategy and population transport assumptions; it should not be implemented by attaching one gamma to existing BED values and calling the result validated.

## 8 Source inventory and provenance

All files below are in the active application's `References` collection. Page counts and SHA-256 hashes are in `clinical/publication_inventory.json` in the evidence package. “Context” means read and assessed for scope, not implemented as an additional model. Current main-paper claims have been distinguished from unavailable supplementary information.

| Local source | Pages or content | Role and key locator |
| --- | --- | --- |
| Bentzen 2026 HYPNO | 12 pages | Trial benchmark; design p389, analysis pp391-392, acute toxicity pp393-394, Table 3 p396. |
| Fenwick 2008 mucosal tolerance | 10 pages | Implemented boundary; equations 10-11 p630, repair formulation p626. |
| Fowler 2003 acute mucosa | 8 pages | Implemented acute defaults/zone; Table 1 elapsed clock and conclusions. |
| Fowler 2010 BED review | 15 pages | Gray-zone and parameter context; schedule table/Table 4, pp561-562; distinguish EQD2/day. |
| Gupta 2020 HYPNO perspective | 6 pages | Design rationale/context; unpublished 6.4-Gy estimate and historical series table. |
| Hartley 2011 parameters | 5 pages | DI preset Table 1 p31; endpoint-specific aggregate correlations. |
| Jacinto 2018 55/20 cisplatin | 9 pages | Feasibility and actual duration/acute events; not randomized validation. |
| Joiner 2023 lecture | 28 slides | Secondary teaching source; slide 16/Table 10.1; does not verify each cited original study. |
| Kirkpatrick 2010 QUANTEC cord | 8 pages | Endpoint-specific OAR limitation, not implemented cord risk model. |
| Madhava 2006 carboplatin | 5 pages | Small dose-escalation cohort and mucositis duration. |
| Meade 2014 chemo tolerance | 7 pages | Context/Table 4, pp2021-2023; study-specific drug-effect estimates, not active automatic correction. |
| Miller 2026 cohort | 8 pages | Retrospective population and assessment-time-specific toxicity. |
| Moon 2024 HYHOPE protocol | 47 pages | Version 4.2 title/signatures; PDF p15/printed p10 dose/model table. Some older body headers persist in the source. |
| Moon 2025 HYHOPE outcomes | 9 pages | Phase 1 population, DLT/results, recurrence endpoint and follow-up. |
| Pardo-Montero 2021 | 10 pages | Refit and cell-dynamics validation context, Tables I/II and Appendix A. The refit is not silently substituted for the active Fowler set. |
| Pardo-Montero supporting XLS | 116 schedule blocks | 108 analyzed schedules plus 8 partial-course blocks; not all represent scalar schedules directly supported by the app. |
| Stuschke 1999 | 9 pages | Tumor/late estimates and slope; pp114,115/Table 1,117 and119/Appendix A. R4 addresses its newly expanded use. |
| van Leeuwen 2018 review | 11 pages | Heterogeneous tumor parameter literature; HN estimates are not one universal fit. |
| van Leeuwen Additional File 1 | 15 pages | Study/model endpoint detail, Table S4 PDF pp7-8 and S6 p15. |
| Withers 1988 | 16 pages | Compensation rule and lag; abstract/p140; does not supply the chosen alpha=0.2. |

Additional original sources used for the tumor derivation are Shuryak et al. 2018, pp21-22/equations 1-5 and fitting section, and the 2019 optimization paper's methods, Figure 1 and Table 1. These are available outside this References collection and were checked in the mathematical audit. [2018 model](https://doi.org/10.1016/j.radonc.2018.02.015), [2019 optimization](https://doi.org/10.1016/j.ijrobp.2019.02.025).

The following are still unavailable in the local collection: HYPNO Appendix E1/Figure E1; HYHOPE outcomes E1/E2; Miller Supplemental Table 1; full primary texts for several catalogue-only historical references, including Bentzen 2001, Denham 1995, Fowler 2007/2008/2009, Fenwick 2012, Strigari 2012, Rezvani 1991, Maciejewski 1986 and Sanghera 2007. Secondary/model-dataset discrepancies for historical EORTC outcomes need the original endpoint definitions before being promoted to validation fixtures. Do not block the current code corrections on obtaining every historical paper, but obtain any source before using its unverified numbers as new model inputs or validation targets.

Further primary links: [Fowler 2003](https://doi.org/10.1016/S0167-8140(03)00231-7), [Fowler 2010](https://doi.org/10.1259/bjr/31372149), [Fenwick 2008](https://doi.org/10.1016/j.ijrobp.2008.02.008), [Withers 1988](https://doi.org/10.3109/02841868809090333), [Pardo-Montero 2021](https://doi.org/10.1002/mp.14834), [Meade 2014](https://doi.org/10.1007/s13277-013-1267-0), [van Leeuwen 2018](https://doi.org/10.1186/s13014-018-1040-z), [Gupta 2020](https://doi.org/10.1016/j.oraloncology.2020.105045). Release workflow references: [renv reproducibility introduction](https://rstudio.github.io/renv/articles/renv.html), [r-lib setup-renv action](https://github.com/r-lib/actions/tree/v2/setup-renv), [r-lib setup-r-dependencies action](https://github.com/r-lib/actions/tree/v2/setup-r-dependencies).

## 9 Ordered implementation plan

### Step 1 Freeze the regression evidence

Start a correction branch from the reviewed commit and retain this report/evidence package. Save current default candidate and agreement CSVs, the R1 invalid-candidate cases, the R5 150-day range case and actual schema-1 files with original results. Record source hashes and versions. This provides an explicit baseline rather than updating expected results to whatever the new code produces.

Done when fixtures are independent of presentation rounding, contain full parameters/anchor/calendar/policies, and can run in a clean R session. Keep scientific model and UI tests distinct so a text change cannot hide a numerical change.

### Step 2 Repair ranking and align scenario summaries

Implement R1 before changing any default or adding a model. Add per-row ranking eligibility and finite-score guarantees, then fix next-dose messages and export semantics. Implement R5's chosen range population. Preserve the PI's invalid-scenario reporting decision and all default numerical outputs. Increment selection version and update new exports accordingly.

Done when the 2.10-Gy valid candidate beats the invalid 2.50-Gy candidate under smaller-change ranking, no model-specific list contains an NA score, 16 criterion/rank edge combinations behave as specified, and the 150-day count/ranges are reconstructible from scenario detail. Re-run the independent 800-case equation and eight-grid oracle; tolerance 1e-8 Gy is adequate for these double-precision checks. Core equation outputs should remain unchanged.

### Step 3 Complete strict replay and explicit migration

Implement R2 scalar/type/range checks at the public replay boundary. Resolve agreement's redundant fields. Implement R3 strict replay versus migration, preserve original provenance and state behavior for old equation versions. Add bounded grid validation and useful field-specific error messages. Do not put JSON coercions inside numerical functions as an accidental substitute for schema validation.

Done when all valid current exports round-trip listed rows, full qualifying rows and scenario detail within 1e-8; malformed flags and inconsistent agreement settings fail before calculation; object-key reordering preserves ID; old files either replay faithfully or undergo explicit recorded reanalysis. R1's selection change must not be silently applied to historical results.

### Step 4 Correct the scientific illustration and active explanations

Implement R4 as a named conditional model derivation with exact numerical fixtures. State whether gamma is assumed on effective EQD2 or transformed, and show central/interval slope effects and the conditional PH conversion. Complete R6-R7 comparator scope, v1 history, default prose and preset provenance. Annotate superseded assertions in the development plan without fabricating a new historical test record.

Done when every probability statement specifies its extra model assumptions, no sentence equates a BED margin with clinical NI, all default narrative signs match the displayed comparator, and every preset field is labeled source/conversion/chosen/inherited. This step can preserve D18 and the selected acute policy exactly.

### Step 5 Improve cross-tab context and clinical inspection

Add a compact context strip to each analysis: evaluated schedule, actual comparator, anchor with T, tumor/acute/late alpha/beta, selected acute/late policy and effective caps, dose grid, gain/loss allowance and saved days. Keep shared settings visibly shared. Add a row-inspection view that carries the selected schedule, its comparator, calendar and anchor together; do not merely overwrite n,d while leaving a different comparison behind. Add visits saved, toxicity slack and source/population details for benchmarks.

Done when moving between candidate, agreement and sensitivity cannot silently change the reference implied by the text. Inspecting a row reproduces its source table exactly. Default computations remain unchanged. Optional later improvements are a standalone sensitivity export and a separately labeled current-custom-scenario row. Robust optimization over all scenarios, TCP/NTCP, chemotherapy corrections and repair-aware delivery are separate scientific features requiring their own specifications.

### Step 6 Restore and test the exact release environment

Complete R9's restore/activation workflow and run the unmodified suite on a clean supported Windows host. Re-run source fixtures, targeted regressions, all selector paths, numeric-invalid states, both exports and the browser procedure below. Investigate any digest precheck problem on that clean host before deciding whether application code needs a portable-read change. Do not ship the diagnostic monkeypatch as an unexplained dependency modification.

Done when the released commit and versions match all logs, source fixtures and exported metadata; every newly repaired behavior has a regression that fails on the reviewed version for the expected reason; all required checks pass without unexplained errors; and the current browser record is complete. Make no clinical validation claim from passing software gates.

## 10 Browser acceptance procedure for the next release

Run against the exact post-fix commit with a fresh session, on a normal desktop width and a narrow laptop/mobile-like width. Record browser/version, app/equation/selection versions and date. A screenshot alone is not an arithmetic test; compare values with the fixtures as well.

1. Open/close each of six sidebar accordions; tab through inputs with the keyboard. Verify labels, units, focus and validation messages stay visible. Enter a blank, out-of-range and fractional fraction count; restore valid values. No stale output should look current after invalidation.
2. Toggle regimen and comparison auto/manual time. Confirm 35 fractions gives 47 days, 20 gives 26 and 15 gives 19; too-short manual T is rejected. Select both HYPNO comparators, confirming 20/2.75/26 and 33/2/38 and appropriate custom-field visibility. Verify slider bounds after each comparison change.
3. Select every tumor preset, then edit a controlled value, an inherited DD value and g. Check preset feedback and explanations. Change Withers r to 5 and 10 and verify 0.84 and 0.72 BED Gy/day. Verify fixed-anchor invariance under comparator changes.
4. Exercise zero DI, zero DD with positive DI, non-positive reference and non-positive candidate cases. Check Results, table, plot, agreement and sensitivity consistently label unavailable values. Run R1's exact 32-fraction reproduction through the controls.
5. In Results, inspect all rows, late sensitivities, Fowler/Fenwick text and benchmark notes. Change acute and late alpha/beta to nonstandard values and check dynamic units, horizontal scrolling and clinical comparator descriptions.
6. In Candidate, select all four criteria and rankings, each dose step, all four acute policies and both late policies. Check conditional excess inputs, custom-Fowler blocking, chemotherapy effect restricted to Fowler, empty-state messages, table sort/pagination and plot legend/axes. Compare default rows to section 7.
7. In Agreement, switch HYPNO/sidebar comparison, allowance and days. Clear candidate-only gain/days/range and confirm agreement remains usable; restore them. Change shared toxicity/step/maximum settings and verify agreement updates. Check all table columns scroll into view and range/count denominators agree.
8. In Sensitivity, confirm the evaluated regimen is the sidebar regimen, not a candidate row; verify 15 scenario rows, inherited DD/g, invalid reasons, summary threshold and late-toxicity sensitivities. If row inspection is implemented, verify its context transfers correctly.
9. Download both ZIPs through the browser for normal, empty and DD-unavailable states where supported. Open each archive, inspect filenames and metadata, and replay in a fresh R process. Compare full-precision listed/qualifying rows and scenario detail, not just rounded screenshots.
10. Open Model and sources. Check MathJax fractions/exponents/symbols, all source links, long provenance tables, version footer and section navigation. Confirm no clipped columns or browser console errors. Repeat essential steps at narrow width and with browser zoom enlarged.

## 11 Verification evidence and how to use it

The accompanying evidence ZIP contains the independent mathematical harness, targeted counterexamples, server/export probes, raw and diagnostic suite logs, source manifests and the clinical derivation notes. It excludes application copies and publication PDFs. Its README gives paths and rerun commands; use an explicit `RepoRoot` pointing at the app being tested. Results must be interpreted against the recorded source hashes, not merely the folder name.

The equation oracle checked 800 cases: 532 feasible calibrations, 220 infeasible and 48 zero-accelerated-repopulation cases, with zero mismatches. Maximum absolute discrepancy was about 2.04e-10 Gy overall and 4.73e-11 Gy for tumor BED; calibration q differed by at most 4.32e-12. Calibration boundary probes correctly reject zero/negative onset and retain arbitrarily small positive feasible onset within numerical precision; for the D18 anchor the critical DD rate is approximately 0.1999442642/day.

Eight independently enumerated grids had zero membership/ranking/value mismatches in their specified settings: D18 reference-acute 0/0 listed/qualifying; HYPNO 3/3; Fowler 11/14; Fenwick 12/17; agreement 3/3; legacy reference 7/9; negative-reference-DI fallback 28/508; nonzero-g agreement 3/3. Maximum numeric difference was about 3.55e-13 Gy. These grids do not contain R1's candidate-invalid ranking condition, which is why the targeted reproductions are required in addition to established regression tests.

The raw full suite ran 101 cases with 1,644 passed assertions, five errors in source hashing, zero failed assertions and zero skips. The documented process-local actual-read diagnostic run had 1,668 passed assertions, zero failures/errors/skips, and two in-suite package-build warnings, plus a separate testthat startup warning. The warnings reflect packages built under a newer R patch release. Preserve both logs and repeat the unmodified suite in the release environment. Twelve candidate/agreement export probes otherwise reproduced their analyses, including empty, unavailable-DD, custom acute, nonzero-g and invalid-scenario detail cases.

The server audit completed 97 control/state cases, all 135 numeric-validator boundary checks, 16 criterion/rank combinations and 28 additional warning/plot-output checks across seven representative states. Valid-change cases produced no unexpected renderer errors. Active blank fields raised local validation messages; hidden inactive fields and independent-tab settings did not block unrelated analyses. The criterion/rank matrix independently reproduced R1. These checks emulate widget update echoes and verify the server layer, not real browser interaction. One-point plot messages and package-build warnings are retained in the log.

All 17 source/configuration/test files in the final software manifest still match the initial reviewed snapshot, and the original HEAD remains the reviewed commit. The manifests and software subreview supplement section 6; the browser acceptance procedure remains explicitly outstanding.

