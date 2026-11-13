# Clinical and source review: frequentist comparison plan and current fitting explanation

Date: 12 November 2026. Reviewed physical files in `C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45`, reported HEAD `79e9c812f91482fa57ffa1dddabe89dc329a6303`, including the working-copy frequentist plan and its audit. No production files were changed and no model was fitted. This report complements the mathematical and software reviews; it does not repeat the unrelated late-toxicity review.

## Assessment

The proposed offline comparison is useful for assessing sensitivity to prior specification under the existing likelihood. It should not be presented as separating biological truth from assumptions, establishing frequentist coverage for the original clinical estimators, validating hypofractionated HNSCC regimens, or selecting a clinically superior model. The existing Bayesian results can be retained. The current explanation has substantially improved its description of calibration, conditional samples, source populations and uncertainty.

A fresh primary-source check found one remaining concrete source-description error in the HYPNO benchmark. The new methods page also displays a tumor hinge different from the one actually fitted. Two plan-level reporting requirements should be made explicit before implementation so that the resulting intervals and prior comparisons cannot be overinterpreted.

## Closure of the previous clinical/source findings

| Prior item | Current evidence | Status |
|---|---|---|
| R1: +5 clinical-risk claim, universal chemotherapy effect, ascertainment count | `early_fit_in_brief.md:10,13–15,30,36,38` and `about.md:79–84` now call +5 a numerical BED margin with no established clinical-risk meaning; identify the chemotherapy coefficient as a coarse pooled adjustment; distinguish 12 families with any post-RT visit from 2 meeting the declared ascertainment rule. | Closed in the inspected text. Preserve these qualifications. |
| R2: lost schedule mixtures and clinical qualifications | `tools/early_tox/export_fit_for_app.py:132–187` now exports all CAIR variants and their weights, schedule notes, assessment gap, population, endpoint time point and structured chemotherapy. `draws/early_fit_trials.csv` includes EORTC's unresolved pattern, DAHANCA's as-treated population, Tandon's assumed equal therapy, HYPNO's four versus five planned courses, and ORO's third post-RT chemotherapy course. | Substantially closed. The new HYPNO source error below is in the underlying record, not a recurrence of the export omission. |
| R7: different endpoint can corroborate an equal count | `tools/early_tox/build_arms.py:120–156` uses explicit endpoint pairings or a source-specific adjudication in `data/early_tox/source/endpoint_equivalence.csv`; free-form presence of the word mucositis no longer suffices. The Bartelink and Poddar exceptions state their basis. | The inspected clinical contract is improved. Mutation-test execution belongs to the software review; I did not independently rerun it. This validates an adjudication contract, not automatic clinical equivalence of every grade/window. |
| R8: unclear active model and editable unused parameters | The implementation record describes hidden Fowler fields in research mode and prior-use records in research-mode exports. The current brief/methods state that fitted early parameters are used and no patient risk is predicted. | Clinical interpretation corrected in inspected text. Runtime and browser closure remain the root/software reviewers' checks. |

Independent CSV checks confirmed 51 primary early arms, 22 families, 6,042 evaluable patients and 2,945 events. The tumor dataset contains 45 rows, 21 fitting units and nominal analyzed N of 8,279, with outcomes at 2, 3 or 5 years. All 45 tumor input rates are time-to-event estimates: 39 actuarial/Kaplan–Meier and 6 cumulative-incidence estimates. These are not 45 observed binomial success counts.

## C12-01 — P2 confirmed: HYPNO ascertainment and analysis population remain incorrectly described

**Locations.** `data/early_tox/source/ascertainment.csv:30` records the population as `randomized and treated`, with blank follow-up fields and the note `benchmark only`. This becomes `follow-up after RT not stated` in `draws/early_fit_trials.csv:66–67` through `tools/early_tox/export_fit_for_app.py:117–129,174–187`.

**Primary evidence.** Bentzen et al., HYPNO, printed p. 391 (PDF p. 5), Outcomes and Data analysis, states weekly assessment during treatment and immediately after its end; the first scheduled post-treatment visit was about two months after RT. Printed p. 394 (PDF p. 8) explicitly calls 753 of 792 randomized patients evaluable for early adverse events: 373 HFX and 380 NFX. The recorded mucositis counts, 190/373 and 208/380, are correct. The analysis retained randomized assignment, but the toxicity denominator is an evaluable subset rather than all randomized or all treated patients. The source also confirms four versus five planned weekly cisplatin courses on p. 391. [Primary record and DOI](https://pubmed.ncbi.nlm.nih.gov/41711624/).

**Effect.** A reader assessing transport and observation-window bias receives false information about a central external benchmark. This does not change the observed counts, chemotherapy scenarios or saved validation scores.

**Fix.** Change the population to `early-adverse-event evaluable patients, analyzed by randomized arm (753/792 overall; HFX 373, NFX 380)`. Use a timing note such as: `Weekly during and immediately after RT; first scheduled post-RT follow-up about two months after RT. Duration/frequency of the immediate post-RT assessments is not fully specified; departmental follow-up policies applied.` Label the selected endpoint `maximum grade 3+ mucositis during and immediately after treatment`.

Do not manufacture an exact eight-week maximum ascertainment window from an approximately two-month routine visit, and do not mark HYPNO adequate under the primary D3 rule merely because follow-up occurred. Its precise early peak-coverage duration is not established by those passages. Preserve its external-benchmark role.

**Acceptance.** A source-to-display test finds the evaluable population, the weekly/immediate-after-RT wording, and the approximate two-month visit; rejects `follow-up after RT not stated`; and confirms unchanged counts 190/373 and 208/380 and unchanged external-only role. The primary 12/22 and 2/22 ascertainment totals remain unchanged because HYPNO is not in that primary set. Regenerate the affected display/provenance artifacts and hashes according to the pipeline's dependency rules. No refit or validation-score recomputation is indicated unless numeric design inputs are deliberately changed.

## C12-02 — P2 confirmed: the methods equation describes an exact tumor hinge that was not fitted

**Locations.** `how_fitted.md:27` writes `max(0, T − Tk)`. The fitted tumor model uses `smooth_max0` at `fit/stan/di.stan:52,80,84`, and the calculator uses it at `R/bed_models.R:92`. The frequentist plan correctly calls it `smax` at lines 37 and 48.

**Effect.** The new reproducibility explanation is not the actual model near onset. The difference is small but deliberate: at T = Tk, the smoothed hinge is ln(2)/3 = 0.231049 days, not zero. With the plan's rounded λ = 0.028 and α = 0.044, the corresponding absolute BED difference is approximately 0.147 Gy. This example is illustrative, not a recomputation of an app result.

**Fix.** In the tumor equation use `smax(T − Tk)` and define `smax(x) = log(1 + exp(3x))/3`, evaluated stably in code. Alternatively label the exact-max equation explicitly as an approximation, but the exact technical formula is preferable on this methods page. Keep the research early-toxicity exact hinge: `how_fitted.md:70` is correct for that model. Do not change either evaluator to make them artificially identical.

**Acceptance.** The rendered tumor methods show the smoothed hinge and its smoothing constant, while the early methods retain their exact hinge. A focused documentation test distinguishes the two models. All fitted results and exported draw files are retained.

## C12-03 — P2 plan gap: intervals need an explicit clinical-estimand limitation, alongside the numerical bounds qualification

**Locations.** Frequentist plan `docs/development/2026_11_12_frequentist_comparison_plan.md:18–21,31,67–94,135`. Its audit calls the comparison a separation of what trials say from priors at line 18. The current `how_fitted.md:29` already discloses rounded tumor counts, but the new plan's interval/reporting specification does not carry this limitation into its mandatory interpretation.

The model comparison intentionally holds the observation likelihood constant. It therefore also holds constant its approximations. Tumor inputs are reconstructed from KM/cumulative-incidence point estimates, with different horizons and failure definitions; their reported N does not supply the censoring-aware sampling variance of those estimators. The early model combines different grading scales, assessment windows, prescription-dose surrogates, systemic-therapy contrasts and evaluable/as-treated/completer populations. A parametric bootstrap from that same model reproduces the assumed observation process; it cannot restore missing censoring histories, unobserved mucositis peaks or patient-level mucosal dose distributions.

**Fix.** Add a required paragraph to the plan and output template: `These intervals measure parameter and contrast uncertainty under the specified likelihood and design. They do not quantify the full uncertainty of the published clinical endpoints or establish transport to a new regimen/population. Holding the same likelihood makes this a sensitivity analysis to estimation and prior assumptions, not a new clinical validation.` Link the source endpoint/clock/ascertainment records and existing targeted sensitivities. Keep the approved mixed-window primary target; do not silently substitute a full-course maximum or an ITT toxicity target.

Call the proposed width labels operational diagnostics, with their thresholds declared analyst choices. Prefer `relatively precise/weakly identified under this likelihood` when explaining the likelihood profile, reserving stronger assertions about prior influence for the full comparison and prior-sensitivity evidence. Mathematical corrections to the label rule, differing support constraints, nuisance integration and the assertion that a point-estimate difference is exactly a prior effect are covered by the mathematical review.

**Acceptance.** Every final interval table names the estimator/likelihood, target contrast, time convention, confidence versus credible interval, and boundary/ridge flags. The report contains the limitations paragraph, and does not infer clinical validation, biological truth or a regimen's patient-level safety from agreement between methods. The source facts 39 KM plus 6 cumulative-incidence rates and the early mixed-window target are reproducible from the current files. No new source model or blanket rerun is required for this correction; any future likelihood change would be a separate analysis with an explicit estimand and dependency list.

## C12-04 — P2 plan gap: the methods-page decision rule can suppress clinically relevant findings

**Locations.** Frequentist plan lines 18, 69–94 ask about every parameter and calculator contrast; D2 at line 136 recommends mentioning the result only if α/β is labeled prior-dominated, and nothing otherwise.

**Effect.** A well-identified α/β does not imply that the onset, time-loss term, overdispersion or a decision-schedule contrast agrees across estimation approaches. Conversely, a weak α/β may leave clinically relevant contrasts stable. Selecting the public update only by one latent parameter is disconnected from the calculator's purpose.

**Fix.** Retain the full offline report whatever the result and add a stable methods link once complete. Predeclare that the brief summary is driven by whether any prespecified schedule contrast or substantive source conclusion changes, with numerical uncertainty and qualifications, not only an α/β label. Do not add an app switch under this plan. State explicitly if the existing point scenario and draw packages are unchanged.

**Acceptance.** A narrative test/example covers (a) weak α/β but stable contrasts and (b) apparently precise α/β with changed time-loss/contrast conclusions. Both produce an accessible completed report; only the substantive conclusions alter the short explanation. There is no automatic promotion of an estimator or selection of favorable results. This changes the reporting plan, not existing fits.

## Recommended order and retained work

1. Revise the plan's statistical definitions using the mathematical review, and add the clinical-estimand and reporting rules above before expensive computation.
2. Correct the HYPNO source description and the tumor methods hinge. Regenerate only the dependent source display/metadata; verify exact count and role preservation. Do not change clinical observations to match a model.
3. Complete the bounded tests and smoke run specified by the revised plan; record approval/execution status accurately. This review did not execute them.
4. If the offline comparison is authorized and its gates pass, run only its declared new fitting units. Existing Bayesian tumor and early fits, historical-label comparisons and external validation scores remain valid for their stated targets and need not be repeated merely because the explanation or metadata changed.
5. Publish the offline comparison report with an explicit retention/change statement, the conditional-uncertainty caveat, and results for all prespecified contrasts. Keep Fowler as the accepted default and the early trial-count fit as the explicitly selected research scenario.

## Verification and coverage

`verification.json` records independent dataset totals, the current CAIR/HYPNO display rows and SHA-256 values of the inspected plan, audit, methods and source artifacts. I read the full frequentist plan/audit, current methods template and model-methods JSON, brief, relevant About sections, export/source mapping logic and the cited source records. This pass freshly extracted and visually checked HYPNO printed pp. 391 and 394 from the local original PDF; those page images and text extracts accompany this report. The primary PubMed HYPNO and Mortensen records were also checked online. The CAIR PubMed request yielded no useful article content, so no new full-text claim is made for that request.

Prior verified source readings for unchanged CAIR, EORTC, DAHANCA, Olmi, Tandon and other trial counts were carried forward transparently. The current source display was checked against those previously established facts. This is not a claim to have reread the entire publication collection this turn. No browser test, expensive fitting, source-data mutation, deployment or general late-toxicity audit was performed by this subreview.
