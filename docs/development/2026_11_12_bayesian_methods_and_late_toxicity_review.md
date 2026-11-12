# Bayesian methods access and late toxicity model review

Igor Shuryak, MD, PhD — is144@cumc.columbia.edu

2026-11-12. Internal review and implementation plan for the head and neck squamous cell carcinoma (HNSCC) biologically effective dose (BED) calculator, version 4.6.0.

## 1. Recommendation

**Implement the unified methods page. Proceed with a late-toxicity feasibility analysis, but revise the proposed statistical plan before fitting a late model.** The existing audit correctly identifies important deployment, documentation, and data limitations. However, several of its proposed corrections are themselves incomplete or incorrect.

The Bayesian explanation is not entirely absent. On the live site, the early-toxicity dialog already explains the model, priors, sampling checks, validation, and limitations. The tumor dialog provides a useful summary and links to a full explainer. The practical problem is fragmented access and an incomplete first explanation of the tumor fitting method. A common destination in the menu and Details should resolve this.

The same **trial collection, treatment schedules, extraction framework, and computation tools** can support a late-toxicity analysis. The existing **acute mucositis counts and tumor-control outcomes cannot serve as late-toxicity outcomes**. New extraction from the existing papers is necessary. Several papers contain useful late data, so feasibility work is justified. However, they do not measure one interchangeable outcome, and they may not identify both late fractionation sensitivity and repair kinetics precisely.

The highest-priority corrections are:

1. Explain the deployed calibrated dose-dependent tumor model, including its dependence on the dose-independent fit. Do not describe the standalone alternative model as the active fit.
2. Replace the proposed combined endpoint, “late mucosal or any-site toxicity at 2–5 years,” with separate endpoint, population, time, and estimator definitions.
3. Match each late-data likelihood to what the publication actually estimates. Do not convert actuarial rates into supposedly exact counts.
4. Evaluate identification of the risk slope, α/β, and repair together. A converged sampler cannot establish that the data identify these quantities.
5. Treat uncertainty as a combination of parameter uncertainty, outcome measurement, schedule reconstruction, and model structure. A wider α/β interval alone is insufficient.

This document supplies the decisions, implementation sequence, and acceptance criteria for those changes. It does not authorize a production refit, source-code modification, commit, or deployment. Those were not performed during this review.

## 2. What was checked

The main source tree is `C:/Users/Igor/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45`. Paths below are relative to that tree unless identified otherwise. The reviewed worktree reference is `52f0d3b1b672524ca245c893dffcf186486a6950`. The review covers the pasted audit and `docs/development/2026_11_12_plan_audit_fitting_page_and_late_refit.md`, the current application routes and model evaluators, active fitting code, model packages, selected primary trial reports, and the live website.

The live site displayed version 4.6 and the default list of 69 dose/fraction combinations. The tumor and early-toxicity fitting dialogs were opened, and the six-item menu was inspected. It has no dedicated fitting-methods destination. This establishes the current visible behavior, but not the live source revision. The local deployment and test markers contain `faf68173fc53d2cad467c6c7133277191c71ea15`, while the README and prior audit name `c80f63b`. Neither marker should be used as proof of the live commit without checking the deployment record.

The dose-list sensitivity results below were recomputed with the actual R application engine and checked against an independent Python implementation. Research early-toxicity parameters were loaded through the application's hash-validated fit loader. This review did not rerun the expensive early or tumor fitting pipelines, extract every late outcome in the library, or certify a complete production release.

## 3. Issue 1: a clear and accessible explanation of Bayesian fitting

### 3.1 Correct the diagnosis

The current tumor dialog says that α/β is mainly determined by its prior, describes trial baselines, and links to the full explainer. Therefore, the audit's statement that it “mentions no prior” is too strong. What is missing is an accessible explanation of what the prior and likelihood mean, how they combine, what was sampled, and which assumptions are inherited by the displayed results.

The early dialog already names Bayesian hierarchical logistic regression, describes the priors and fixed parameters, and reports validation limitations. Preserve that useful content. Consolidation should improve discovery and consistency rather than create another independent numerical summary.

Use the menu and page label **“How the models were fitted (Bayesian)”**. Keep the existing tumor and early contextual links recognizable, but send them to the relevant section of the same page. A contextual link can retain “How was the early-toxicity model fitted?” while opening the common page at its early section. This avoids making a user who asks about mucositis search through the tumor explanation.

### 3.2 Required page content

The page should have a short opening followed by expandable technical sections. A suitable opening is:

> The fitted models combine trial results with stated assumptions about plausible parameter values. The resulting parameter samples describe uncertainty under those assumptions. Tumor control and acute mucositis were fitted separately. Late toxicity currently uses an assumed fractionation-sensitivity range. The calculator's support percentages describe whether model criteria are met, not a patient's probability of tumor control or toxicity.

Use the following content order:

| Section | Required explanation | Source of truth |
|---|---|---|
| What “Bayesian” means | A prior describes assumptions before the analyzed data. The likelihood describes how the observed outcomes could arise under a parameter set. The posterior combines both. Explain why limited fractionation contrasts can leave a parameter largely determined by its prior. | Plain-language narrative, checked against the actual models |
| What was fitted | Distinguish tumor control, acute mucositis, and the currently unfitted late-toxicity model. State outcome definitions, selection rules, fitting units, follow-up differences, and concurrent-treatment scope. | Active data manifests and extraction tables |
| Tumor model | Binomial observation model with the Poisson tumor-control response form, trial-specific baselines, and arm-level extra variation. Explain reconstructed counts and their approximation. Separate dose-independent (DI) and calibrated dose-dependent (DD) models. | `fit/stan/di.stan`, `fit/stan/dd_cal.stan`, active driver/data files and exported draw identities |
| Acute mucositis model | Hierarchical logistic count model using within-trial BED differences, grading-group adjustment, planned chemotherapy indicator, and arm-level extra variation. Explain which parameters were fitted and which were fixed in each scenario. | `tools/early_tox/step5_model.stan`, `step5_fit.py`, active scenario data and `draws/early_fit_manifest.json` |
| Priors and fixed assumptions | Parameter, units, distribution including bounds, source or assumption status, scenario, and sensitivity analyses. Show actual runtime values rather than defaults from unused code. | Serialized fit metadata generated from run inputs |
| How fitting was checked | Define Markov chain Monte Carlo (MCMC), effective sample size (ESS), chain agreement measured by R-hat, divergences, and posterior predictive checks. Explain that computational checks do not establish clinical validity. | Saved diagnostics for retained exported draws |
| What validation supports | Separate in-sample replication, simulated-data calibration, held-out trial prediction, and transport to shorter courses. Include unfavorable findings and outcome-assessment limitations. | Existing validation outputs with model identifiers |
| How samples affect this app | Distinguish point values, tumor-only Refit support, joint toxicity uncertainty, and separate scenario comparisons. Explain shared DI/DD draws and assumed relationships between tumor and toxicity fits. | `R/uncertainty.R`, `R/v4_home.R`, current exports |
| Sources and reproducibility | Link to the fitted trial lists, full explainer, relevant primary references, and a downloadable methods snapshot containing model and package identifiers. | Deployed assets and immutable fit manifests |

Do not expose file-system details in the opening clinical explanation. Put fit identifiers and provenance in the technical section or download.

### 3.3 Correct the tumor-model description

The audit's proposed ban on the word “hierarchical” is not a valid correctness test. The tumor model contains arm effects with an estimated common scale, which are a hierarchical component. The precise distinction is that its trial intercepts have a **fixed common prior specification**, whereas the early model estimates the population distribution of trial baselines.

Recommended tumor wording:

> The tumor fit uses a Bayesian binomial model with a Poisson tumor-control response function. Each fitting unit has a separate baseline with a fixed weak normal prior. Arm-level random effects, with an estimated common spread, allow additional variation between arms. This differs from the acute-mucositis model, which also estimates a common distribution for trial baselines.

The active DD fit is especially important. `dd_cal.stan` takes α, α/β, and background growth from a DI draw, solves its onset constant to agree with DI at the standard regimen, and fits its remaining rate and nuisance parameters conditionally. It also includes an onset penalty. Consequently:

1. The DI and calibrated DD predictions are dependent. Agreement at the anchor is imposed.
2. The standalone `dd.stan` prior on its free onset constant is not a prior on a freely estimated parameter in the deployed calibrated model.
3. The DD prior table must include the actual runtime rate prior and onset-penalty settings from the active driver/data. Reading Stan source alone cannot recover values passed as data.
4. Describe the sequential conditional fitting construction accurately. Do not call it two independent confirmations or an unrestricted joint fit of both tumor models.
5. Use the fitting unit count separately from the publication/trial-family count. A trial split into strata can supply more than one fitting unit.

The current briefs also mix BED and equivalent dose in 2-Gy fractions (EQD2). The conversion at α/β = 10 Gy is `EQD2 = BED / 1.2`. Thus 0.63 Gy of BED/day is approximately 0.53 Gy EQD2/day. These are compatible rounded descriptions, not conflicting estimates. When generating a posterior summary, convert each draw before taking quantiles. Dividing a marginal median by a conversion factor based on a different marginal median need not reproduce the posterior median of the converted quantity.

Two additional distinctions belong on the page. The home list's Refit support uses the **either-tumor-model** rule, while the Uncertainty page has a configurable criterion whose default requires both tumor models. Also, the tumor export's 4,000 paired rows contain 1,000 distinct DI draws, each repeated four times with conditional DD samples. They are not 4,000 independent DI draws. Preserve that grouping when assessing Monte Carlo precision. The early package contains 4,000 joint draws from its own fit. Pairing tumor and toxicity draws independently is an assumption, not an estimated cross-endpoint joint distribution.

### 3.4 Active settings that the page must represent

The following is the implementation reference for the primary scenarios. For a lognormal distribution, the two parameters below are the mean and standard deviation of the natural logarithm. Positive-constrained normal distributions are half-normal where the mean is zero. Standardized random effects have standard normal priors.

| Active component | Verified prior or fixed setting |
|---|---|
| DI tumor α, Gy⁻¹ | Lognormal(log 0.2, 1) |
| DI accelerated regrowth rate λ, day⁻¹ | Lognormal(log 0.1, 1.5) |
| DI onset, days | Normal(28, 10), truncated to 0–60 |
| DI tumor α/β, Gy | Lognormal(log 10, 0.354) |
| DI background growth, day⁻¹ | Uniform over 0.0046–0.021 |
| Tumor trial/stratum baseline | Normal(0, 3), fixed prior scale |
| Tumor arm-effect standard deviation | Half-normal with scale 0.5 |
| Calibrated DD rate λ, day⁻¹ | Lognormal(log 0.5, 1.5), conditional on the DI draw |
| Calibrated DD onset condition | Soft exponential penalty favoring onset at or after day 7 for every fitted arm, with contributions summed and scale 1 day. This is not a hard day-7 lower bound. |
| Early population baseline mean | Normal(0, 2.5) |
| Early between-family baseline standard deviation | Half-normal with scale 1 |
| Early grading-group offset | Normal(0, 1) |
| Early log-odds slope per 10 Gy of BED | Half-normal with scale 1 |
| Early α/β, Gy | Lognormal(log 10, 0.354) |
| Early time-loss rate, Gy of BED/day | Lognormal(log 0.8, 0.5) |
| Early planned-chemotherapy coefficient | Normal(0, 1) |
| Early arm-effect standard deviation | Half-normal with scale 0.5 |
| Primary early fixed settings | Onset 7 elapsed days, repair half-time 3.2 hours, exact time hinge |

These values trace to `di.stan`, `dd_cal.stan`, `step5_model.stan`, their drivers, and the retained fit identities. The retained tumor identities are stored under the sibling `Custom_BED_calculator_v3/fit/output_refit_2026_10_25c_main/` tree. Their Stan hashes agree with the reviewed current model files. Package generation should serialize this information so deployment never depends on that sibling path.

The DI fit used four chains, 1,500 warmup iterations and 2,000 post-warmup draws per chain, target acceptance 0.99, and maximum tree depth 12. The conditional DD process used 1,000 DI draws, initially four chains with 600 warmup and 250 post-warmup draws per chain, with retries and documented repairs. The primary early fit used four chains with 1,000 warmup and 1,000 post-warmup draws per chain and target acceptance 0.95. Report retained-run settings and repair histories rather than implying all fits used one sampling configuration.

The tumor manifest explicitly distinguishes current per-quantity checks from older fits for which only aggregate diagnostics were retained. The early diagnostic code records tree-depth hits but does not include them in its current pass condition. Therefore, the page must not claim that one identical diagnostic gate governed every fit. These are documentation requirements and a separate candidate for a targeted gate-consistency repair, not evidence that the current draws failed convergence.

### 3.5 Implementation sequence for the methods page

**M1. Freeze the source and content contract.** Record the worktree revision, active fit identifiers, draw hashes, data hashes, and the assets that will be deployed. Reconcile the stale deployment markers before making a new live-version claim. Write a content checklist using the sections above. Preserve existing explanatory caveats that remain supported.

**Verify:** every numerical page item is linked to a named data or summary artifact. The page distinguishes the tumor, early, and late model states. No claim depends on a superseded standalone DD fit.

**M2. Create one deployable page and one navigation helper.** Add `how_fitted.md` at the application root. Add a page panel and a single helper for opening it, optionally at the tumor or early section. Route the menu item, Details button, three current tumor links, and the current `early_how` link through that helper. Preserve access to the trial tables and the full explainer.

**Verify:** use a real matching navigation value. Existing panels use full page-title values. An assertion that `section == "fitted"` is only meaningful if that explicit value is added consistently to the panel and navigation calls. Test the emitted navigation target or rendered panel, not merely that clicking caused no error.

**M3. Separate narrative from generated fit metadata.** Reuse the existing briefs as the starting result summaries, but generate shared numerical values and prior tables from active run inputs and exported artifacts. Add a small machine-readable methods metadata file, for example `draws/model_methods.json`, with schema version, model identity, transformations, units, fixed parameters, priors, fit settings, diagnostics, and source hashes. Generate it during packaging, not by fitting models when the page is opened.

**Verify:** compare metadata against both model code and serialized run data. A test must fail when a runtime prior width, truncation bound, fixed repair setting, calibrated-DD penalty, or fit identity is changed without regenerating the metadata. A text search for 0.53/0.63 is not a sufficient anti-drift test.

**M4. Render equations and tables deliberately.** Use the established Shiny rendering approach with math support, `withMathJax(HTML(markdown::mark(text = txt, options = "+latex_math")))`, verified in the installed package version. The `HTML()` wrapper preserves rendered markup. Extend equation/table styling to the new page. Existing equation styling scoped to `#about` will not automatically apply. Avoid nested duplicate headings when embedding the briefs.

**Verify:** desktop and phone views show all equations, subscripts, symbols, table columns, and source links legibly. Keyboard navigation reaches the page, focus moves predictably, and a return action preserves the user's model and margin settings.

**M5. Make deployment fail on absent or empty content.** Update `tools/deploy_preflight.R::check_bundle`, its required-file lists, and test fixtures. `docs/`, `fit/`, and `tools/` are excluded from the deployed bundle, so the page cannot read them at runtime. Include the generated metadata and any required summary assets in the bundle.

**Verify:** remove each required asset in a temporary bundle and confirm refusal. Repeat with a zero-byte or unreadable methods file. Adding a filename to the existence list alone does not exercise the current separate empty-served-asset checks. Start the app from the deployment-shaped bundle and open every methods route.

**M6. Update only necessary ripple items.** The affected locations include `V4_MENU`, `V4_MENU_NOTES`, panel routing, the version fields, `test_v27_colleague_review.R`, `test_v31_stage_e.R`, `test_v45_review.R`, `test_deploy_preflight.R`, the walkthrough, README, and the current user guide. The existing audit omitted the six-destination/menu-name assertions in `test_v45_review.R`. Use the next appropriate patch version after confirming the current release state.

**Verify:** menu membership, version checks, and walkthrough expectations agree. Run the full normal suite after implementation, then inspect the change diff for unrelated edits. If an existing user-guide document is revised, provide its required clean and tracked/diff versions together.

**M7. Release independently of the research refit.** The documentation change does not need to wait for late-model fitting. Complete a local bundle walkthrough, prepare the concrete release diff and test results, and deploy only within the authorized release workflow. After deployment, open the actual live page through the menu, Details, tumor links, and early link, and verify the displayed package identity.

**Verify:** the live page is served from the uploaded bundle, all routes converge on it, and no late model is described as fitted before such a model exists.

## 4. Issue 2: what the existing evidence can support

### 4.1 Reuse the collection, not the outcomes

The raw early-toxicity outcome table has 68 rows of mucositis-related outcomes. That is not the number of independent late observations, nor the size of the active primary early fit. The live primary early summary reports 22 trial families, 51 arms, and 6,042 evaluable patients. Selection, alternative endpoints, and sensitivity records explain why the raw table and fitted analysis have different counts. Keep those levels explicit in the late extraction as well.

The prior audit's search result, 96 of 120 portable document format (PDF) files containing late-effect terms, is a screening yield. It does not establish 96 usable studies or independent trial families. Full reports nevertheless provide concrete evidence that a late extraction is worth doing. The four RTOG 9003 arms below are ordered as standard fractionation, hyperfractionation, accelerated split course, and accelerated concomitant boost:

| Trial/report | Directly observed information | Appropriate use and restriction |
|---|---|---|
| RTOG 9003, Fu 2000, Table 4 | Grade 3/4 mucous-membrane counts are 11/254, 19/253, 19/261, and 29/261 across its four arms. The report's late definition begins more than 90 days after treatment start. | Candidate mucosal dataset after adjudicating the evaluation population and observation window. These are worst-recorded late grades, not automatically fixed-horizon cumulative incidence. [1] |
| RTOG 9003, Beitler 2014, Tables 3 and 4 | Five-year grade 3–5 any-site prevalence is 6/71, 8/84, 12/70, and 13/80. Worst late grade after 180 days uses different denominators and gives 68/234, 65/232, 69/239, and 85/232. | Distinct estimands within the same trial family. Do not count the later report as an independent trial or substitute its prevalence for cumulative risk. [2] |
| DAHANCA 6&7, Mortensen 2012, Table 2 | Peak severe fibrosis counts are 70/632 and 83/617, whereas five-year Kaplan–Meier estimates are 13.0% and 16.7%. | Tissue-specific fibrosis evidence. Counts and actuarial estimates are not interchangeable observation models. [3] |
| RTOG 83-13, Fu 1995 | Dedicated late-effect analysis explicitly accounts for death as a competing event. | Important dose/repair evidence, but extract the estimator and actual schedule information. This contradicts the blanket claim that all late outcomes are counted only among survivors. [4] |
| CHART, Bentzen 1999 | Late laryngeal edema, subcutaneous fibrosis, and telangiectasia are analyzed, with actuarial estimates and uncertainty. Repair estimates depend on assumed/external fractionation sensitivity and dose-response steepness. | Separate tissue outcomes and a repair sensitivity analysis. These are not direct estimates for severe late mucosal injury. Do not reuse the derived repair estimates as independent evidence when fitting the same CHART outcomes. [5] |
| EORTC 22791, analyzed by Stuschke and Thames 1999 | The reported grade ≥2 late α/β estimate is 4.0 Gy with a 3.3–5.0 Gy interval at fixed dose-response steepness γ₅₀ = 5. Smaller assumed steepness gives smaller estimates. | Evidence for a conditional estimate, not a universal empirical late α/β distribution. The source overlaps with a trial that may enter the new likelihood. [6] |

The original EORTC 22791 report also illustrates exposure uncertainty: late effects were evaluated among patients with locoregional control, and reported subcutaneous hot spots differ from prescription doses. Table XI and Figure 2 have a small event-count discrepancy that requires adjudication before extraction. Consequently, prescription BED cannot automatically be interpreted as tissue-specific absorbed-dose BED. [7]

These examples are extraction leads and proof of feasibility, not a finalized pooled dataset. Dedicated reports of consequential late injury also support testing a model beyond direct fractionation effects, but an acute-to-late association does not identify a causal mediation effect from aggregate trial arms. [8]

### 4.2 Choose the scientific target before the likelihood

Two objectives should be kept separate:

**A. Estimate radiobiological fractionation sensitivity for a defined late tissue endpoint.** A coherent fibrosis endpoint may be more informative for α/β than a sparse severe-mucosal endpoint. It still needs compatible grades, exposure definitions, follow-up, and populations. An estimate for fibrosis must remain labeled as fibrosis.

**B. Predict a clinically important late outcome for comparing schedules.** Severe mucosal injury, swallowing dysfunction, persistent feeding-tube dependence, necrosis, and any-site severe toxicity are different outcomes. A model of one cannot be advertised as a general late-complication model.

Recommended decision rule: extract separate tissue-specific datasets first. Select the primary analysis based on prespecified endpoint coherence and usable design information, before viewing fitted α/β values or which model best favors shorter courses. If no endpoint supports a defensible refit, report that result and retain explicit sensitivity scenarios. Do not broaden to “mucosal or any site” merely to increase sample size.

Specify each candidate target as a complete tuple:

> tissue or composite definition; grading system and threshold; time origin and horizon; cumulative incidence, net risk, prevalence, or worst-ever occurrence; eligible population; handling of death, recurrence, salvage treatment, and missing assessments; physical dose definition.

A two-year and a five-year outcome may be modeled jointly only with an explicit time model and appropriate dependence. They must not become independent binomial observations with the same underlying patients. If uncertainty/covariance is unavailable, select one prespecified horizon for the primary analysis and use the other as a separate sensitivity analysis.

## 5. Does changing late α/β matter to the current candidate list?

The actual `v4_regimens()` calculation was repeated at the default tumor point scenario, standard 35 × 2 Gy in 47 inclusive calendar days, default margins of tumor −5/early +5/late +5 Gy of BED, and the current 0.05-Gy dose grid. The acute evaluator converts this to 46 elapsed days. Only late α/β was varied. Both early modes used their actual current application settings. Preserve this clock distinction in the new methods page and late schedule schema.

| Assumed late α/β, Gy | Qualifying doses with Fowler early model | Qualifying doses with research early model |
|---:|---:|---:|
| 1.00 | 55 | 9 |
| 1.50 | 65 | 9 |
| 1.61 | 68 | 9 |
| 1.62 | 69 | 9 |
| 2.00 | 69 | 9 |
| 3.00 | 69 | 9 |
| 4.00 | 69 | 9 |
| 5.00 | 69 | 9 |

These are dose/fraction combinations, not patients or independent schedules from trials. At α/β = 3 Gy, the Fowler list spans 18 fraction counts and the research list spans five.

For complete repair, a candidate with $n$ fractions of $d$ Gy has the following late-BED difference from 70 Gy in 35 fractions, with $a = α/β$:

$$
\Delta B_L = nd - 70\,\mathrm{Gy} + \frac{nd^2 - 140\,\mathrm{Gy}^2}{a}.
$$

For $nd < 75\,\mathrm{Gy}$ and $nd^2 > 140\,\mathrm{Gy}^2$, the current +5-Gy allowance is satisfied when:

$$
a \geq \frac{nd^2-140\,\mathrm{Gy}^2}{75\,\mathrm{Gy}-nd}.
$$

The first existing Fowler candidate is lost below **1.61566265 Gy**, at **29 × 2.30 Gy**. Therefore, “about 2 Gy or lower” is a coarse description, not the actual boundary. The previous examples of toxicity-ceiling ties can also be irrelevant to the final list when the tumor criterion excludes those doses.

Every current research-mode candidate has both total dose below 70 Gy and the complete-repair quadratic dose sum below 140 Gy². Its late BED is therefore below the reference for **every positive α/β**, under this specific late model. Changing α/β alone cannot remove any of those nine candidates.

This does not make late-model improvement unnecessary. It establishes that a narrow project limited to replacing the α/β point estimate is unlikely to change the present default list much. The useful next analyses concern joint repair uncertainty, endpoint-specific clinical interpretation of a BED allowance, uncertainty-based eligibility, other margin settings, and other regimens. A fitted risk model also cannot by itself validate the existing +5-Gy allowance.

## 6. A defensible late-model analysis

### 6.1 Schedule and dose calculation

Use one tested exposure evaluator in fitting, validation, and application prediction. For fraction doses $d_i$, delivery times $t_i$ in hours, and a monoexponential repair half-time $h$ in hours:

$$
D = \sum_i d_i, \qquad
Q(h) = \sum_i d_i^2 + 2\sum_{i<j}d_i d_j\exp\left[-\frac{\ln 2\,(t_j-t_i)}{h}\right], \qquad
B_L(a,h) = D+\frac{Q(h)}{a}.
$$

Here $D$ and $B_L$ have units of Gy, $Q$ has units of Gy², and $a$ has units of Gy. Require nonnegative fraction doses, positive $a$ and $h$, and chronologically ordered delivery times, with $t_j\geq t_i$ for $j>i$. The monoexponential form is a structural model assumption. A more complex repair model should be attempted only if the schedule contrasts can support it.

Do not silently assign a six-hour gap to every twice-daily arm. `tools/early_tox/build_arms.py` currently uses a default when timing is unspecified. A late dataset must distinguish exact reported intervals, minimum protocol intervals, reconstructed schedules, and assumptions. Preserve treatment breaks, boost volumes, mixed fraction doses, and weekday patterns. If an arm contains schedule variants, retain their weights and uncertainty. A nonlinear model generally requires averaging variant-specific probabilities, not applying the risk function to the average BED.

Use repair consistently for the candidate and reference. “Once daily” is not a universal license to drop cross-day repair. As a mathematical sensitivity example, for 70/35 on the Monday–Friday schedule and $a=3$ Gy, a four-hour half-time adds approximately 1.18 Gy to absolute BED compared with complete repair. Whether an approximation is acceptable depends on its effect on the **candidate-minus-reference difference** and the decision boundary. This example is not a proposed clinical repair prior.

### 6.2 Likelihoods must match the reported outcome

For a compatible endpoint with genuine integer event counts, a suitable starting model is:

$$
y_{jk}\sim\operatorname{Binomial}(N_{jk},p_{jk}),\qquad
\operatorname{logit}(p_{jk})=b_j+\kappa\frac{B_{jk}-B_{j0}}{10\,\mathrm{Gy}}+\gamma z_{jk}+u_{jk}.
$$

The indices $j$ and $k$ identify trial family and arm. The reference arm is $0$. The denominator $N_{jk}$ and event count $y_{jk}$ must refer to the declared endpoint and population. The probability $p_{jk}$ refers to that endpoint in that arm. The quantities $B_{jk}$ and $B_{j0}$ are the arm and reference-arm late BED values from the specified evaluator. The baseline $b_j$ and arm effect $u_{jk}$ are on the log-odds scale. The slope $\kappa$ is the log-odds change per 10 Gy of BED. The covariate $z_{jk}$ is a coded planned co-treatment contrast, whose contribution $\gamma z_{jk}$ is on the log-odds scale. The equation is a proposed working model, not a fitted result.

Use a family baseline with a prespecified pooling structure. A fixed-prior baseline sensitivity can assess whether baseline pooling drives results. Add extra arm variation only with explicit prior predictive checks and comparison to the simpler model. Sparse data do not support a large collection of tissue, grading, treatment, and era interactions simultaneously.

Observation rules:

| Published information | Primary treatment |
|---|---|
| Integer events and valid denominator, compatible outcome and observation window | Binomial model, with additional variation if justified. A randomized denominator can estimate cumulative incidence when follow-up to event, death, or horizon is adequate. It is not inherently wrong. |
| Fixed-horizon cumulative incidence with a standard error or confidence interval | An explicitly approximate likelihood for that summary estimator, or an event-history model if sufficient data are available. Preserve the competing-event definition. |
| Kaplan–Meier estimate treating death as censoring | A separate net-risk target or a justified survival model. Do not label it observed cumulative incidence. |
| Landmark toxicity among assessed survivors | A prevalence model conditional on the selected assessment population. Do not interpret the contrast as an unconditional randomized treatment effect. |
| Worst grade during variable follow-up | A separate analysis with its observation-window limitation. Do not silently treat it as a fixed-horizon risk. |
| Percentage without recoverable denominator, uncertainty, or estimator | Narrative evidence or an explicitly labeled approximation sensitivity. Do not create exact events by rounding a guessed effective sample size. |

For compatible published probabilities with uncertainty, a transformed-normal summary likelihood is one possible approximation. Its transformation, standard-error calculation, behavior near zero/one, and treatment of covariance must be specified and tested. Never use a confidence interval for one estimator as uncertainty for a different estimand. Repeated horizons, overlapping publications, multiple endpoints in the same patients, and shared controls require dependence handling or prespecified selection.

### 6.3 Identification is the principal feasibility gate

For fixed repair, the within-trial contrast contains both total-dose difference and quadratic-dose difference. The risk slope and its fractionation component need independent information. Examine the nuisance-adjusted design using these two contrasts, then inspect likelihood profiles and posterior correlations. If repair is estimated, include sensitivity to the repair half-time in the full Jacobian.

Varying dose per fraction within every trial is not mathematically necessary if complementary total-dose and fractionation contrasts across compatible trials identify the common slope. Conversely, several trials with different fraction sizes can still provide weak practical identification. Numerical rank alone is insufficient.

Before a substantive fit, simulate repeated datasets on the **actual proposed design and denominators**, covering several α/β, repair, slope, baseline-risk, and heterogeneity scenarios. Evaluate bias, interval coverage, interval width, convergence, and decision accuracy. Also generate misspecified data with outcome selection, unmodeled treatment-time effects, and schedule error. Report Monte Carlo error for simulation summaries.

The audit's claim that omitting repair necessarily biases α/β upward is too strong. That direction is plausible in some hyperfractionation designs, but depends on contrasts, slope estimation, and nuisance effects. A single simulated dataset cannot prove bias or nominal interval coverage. The simulation acceptance test must not require a predetermined bias sign.

Simulation-based calibration checks whether the fitting implementation recovers its own prior-generated model. It does not establish that the data distinguish α/β from repair or that the model applies clinically. Perform both implementation calibration and design-specific recovery.

### 6.4 Priors and realistic uncertainty

Keep four sources of uncertainty distinct in both outputs and explanations:

1. **Parameter uncertainty:** joint draws of α/β, repair, slope, baseline, and heterogeneity within a named model.
2. **Measurement and extraction uncertainty:** grading crosswalks, denominators, follow-up, uncertain event classification, and source discrepancies.
3. **Exposure uncertainty:** intraday intervals, actual treatment duration, dose hot spots, tissue volumes, and schedule mixtures.
4. **Structural and transport uncertainty:** direct versus consequential injury, tissue specificity, chemotherapy context, competing-event interpretation, and applying historical techniques to contemporary regimens.

Do not hide the last three inside an artificially broad α/β prior. Use named alternative datasets/models or explicit latent exposure models when justified. Preserve joint draws rather than independently sampling marginal intervals. Unless justified model probabilities exist, report scenarios separately instead of pooling them into a single apparently calibrated probability.

The current late range is an analyst-assumed range, not an empirical 95% interval from all late tissues. Stuschke and Thames' conditional result cannot be copied into a prior while its EORTC 22791 data also enter the likelihood. The same issue applies to CHART-derived repair estimates. Track overlap at the underlying trial-family level, including derived analyses and later publications.

Recommended initial strategy: retain the current α/β = 3 Gy calculation as the benchmark, use transparent alternative α/β and repair scenarios for feasibility, and choose proper prior distributions only after the endpoint and identifiable parameters are established. Do not carry the acute repair half-time of 3.2 hours into the late model simply because the evaluator can accept it.

If a broad illustrative lognormal prior is used, label it as a sensitivity assumption and report its actual quantiles. A lognormal with median 3 Gy and central 95% limits 1 and 10 Gy is not mathematically possible. Limits of 1 and 10 Gy imply a median of approximately 3.162 Gy. A median of 3 Gy with symmetric logarithmic limits 1 and 9 Gy is a different distribution. Neither is established by the existing late literature as a universal prior.

Choose slope and baseline priors through predictions of plausible risks across the actual dose contrasts. The early model's slope prior is not automatically appropriate for a rarer late endpoint. If α/β or repair is weakly identified, retain fixed named scenarios or a prior-sensitive research posterior. Do not report a narrow interval conditional on a fixed repair assumption as total late-model uncertainty.

### 6.5 Treatment time, consequential injury, and other context

Identical fraction number and dose do not guarantee identical late biological effect if calendars change repair. A useful model check compares arms matched on total dose, fraction size, relevant repair, co-treatment, and endpoint definition. It tests the direct late-damage model under those conditions.

A planned early-BED contrast can be added as a structural sensitivity for consequential injury. Check its collinearity with treatment duration and late exposure, and propagate early-model uncertainty. Observed acute toxicity is a post-treatment variable. Adding it changes the interpretation and can introduce mediation or selection bias. An aggregate arm association does not identify an individual causal pathway.

Keep planned chemotherapy contrasts separate from delivered chemotherapy, interruptions, salvage surgery, or other post-randomization variables. Adjustment for delivered treatment can remove part of a schedule effect or induce bias. Older radiotherapy-only cohorts, modern chemoradiation, postoperative treatment, nasopharyngeal disease, and frail-patient split courses should have explicit inclusion or validation roles rather than being absorbed into an unexplained pooled “late” parameter.

## 7. Step-by-step late-analysis implementation plan

The file names below are proposed additions. Existing early-fit data and scripts remain the reference implementation until a specific shared component is deliberately extracted and regression-tested.

| Step | Work and concrete output | Verification and stopping rule |
|---|---|---|
| L0. Lock scope | Create `docs/development/late_toxicity_protocol.md` with the two scientific objectives, endpoint tuples, primary-selection rule, population exclusions, causal interpretation, and proposed decision outputs. | Every proposed outcome has a tissue, grade, time origin/horizon, estimator, population, and dose definition. No combined “mucosal or any-site” target. |
| L1. Build a source-family inventory | Create `data/late_tox/sources.csv` and `trial_families.csv`. Start with RTOG 9003, RTOG 83-13, EORTC 22791, CHART, and DAHANCA. Add other already-held reports based on usable endpoint content. | Map every report, update, subgroup, and derived parameter paper to its parent cohort. Record PDF path/hash and exact table/page. A term-search hit is not an eligible observation. |
| L2. Extract independently | Create separate `outcomes.csv`, `schedules.csv`, and `extraction_discrepancies.csv`. Extract counts and summary estimates into different fields. Adjudicate discrepancies against the page image. | Reject noninteger exact counts instead of truncating `float` to `int`. Validate events ≤ denominator, grades, time units, and estimator fields. Source conflicts stay flagged until resolved. |
| L3. Freeze compatible analysis sets | Produce `analysis_sets.csv` and a human-readable inventory: endpoint families, compatible contrasts, exclusions, denominator completeness, follow-up, and dose quality. | Select the primary endpoint without using fitted results. Use one record per overlapping cohort/endpoint/horizon unless a dependence model is specified. If no coherent set exists, stop at evidence synthesis. |
| L4. Implement the exposure evaluator | Add `tools/late_tox/evaluator.py`, or a tested shared evaluator with an explicit late interface. Produce arm-level D, Q(h), reference contrasts, and timing-assumption flags. | Check complete-repair and no-repair limits, two-fraction analytic examples, mixed doses, breaks, overnight/weekend repair, and invariance to shifting all delivery times. Compare R and Python. |
| L5. Assess identification | Add `tools/late_tox/identify.py`, a design report, profiles, and recovery simulations for the final proposed observation model. | Demonstrate which parameter combinations the data identify. If α/β and repair cannot be separated, fix named repair scenarios or stop the joint refit. Do not solve nonidentification by reporting only sampler convergence. |
| L6. Specify likelihoods and priors | Create `fit/late/model_spec.yaml`, prior predictive outputs, overlap exclusions, and the exact model equations. Register all assumed and literature-derived parameters. | Prior predictions are compatible with the declared clinical endpoint. No prior uses the same trial outcomes as independent likelihood information. All lognormal bounds/quantiles and units agree. |
| L7. Implement and fit | Add `fit/late/` model code and a `tools/late_tox/run_fit.py` orchestrator. Preserve raw chains, diagnostics, inputs, code hashes, and resumable run manifests. | Syntax check, independent equation/evaluator checks, small end-to-end run, simulation calibration, then full runs. Check every retained chain/unit for diagnostics. No silent removal of unfavorable draws or failed conditional units. |
| L8. Validate by trial family | Add `tools/late_tox/validate.py`. Refit leaving out whole trial families. Compare with fixed-α/β, complete-repair, and no-dose-effect baselines using the same outcome likelihood and validation protocol. | Group all related reports and source-derived priors with the held-out family. Report predictive calibration, log scores, uncertainty, and failures by endpoint and schedule class. No reliance on in-sample fit alone. |
| L9. Evaluate decision relevance | Produce `late_decision_impact.csv` with both early modes, margins, late scenarios, point and posterior criteria, affected candidates, and extrapolation flags. | Reproduce the current 69/9 fixtures first. Compare candidate-minus-reference differences using paired parameter draws. Show which uncertainty source changes eligibility. Do not convert the default +5-Gy margin into an unvalidated clinical-risk limit. |
| L10. Decide integration level | Write `late_fit_decision.md`: evidence synthesis only; sensitivity-scenario update; research fit; or a justified candidate for default adoption. | Default replacement requires defensible endpoint transport, stable prediction, and meaningful decision behavior. An inconclusive or prior-dominated fit is a legitimate result. Research-only is the appropriate initial integration level if fitting is feasible. |
| L11. Package an accepted research model | Only after the preceding gates, add a late-fit draws file and manifest, loader, evaluator identifier, endpoint/horizon labels, methods-page section, export fields, and replay checks. | Hash mismatch or missing endpoint metadata fails closed. All app surfaces use the same selected late model. Existing defaults remain reproducible, and full uncertainty draws preserve dependence within each model. |

The minimum outcome schema should include: `trial_family_id`, `report_id`, `arm_id`, `endpoint_id`, tissue/composite definition, grade scale/version/threshold, event definition, time origin, horizon, estimator type, analysis population, numerator, denominator, estimate, standard error or interval, confidence level, death/recurrence/salvage handling, assessment completeness, overlap group, source page/table, PDF hash, extraction status, and eligibility reason. Missing values remain missing. Use a unique key that includes endpoint, horizon, estimator, and population, not only arm and endpoint.

The schedule schema should separately record prescription and tissue dose, boost/field changes, individual or repeated fraction patterns, delivery-time units, interval status (exact/minimum/assumed), breaks, mixture strata/weights, and reconstruction uncertainty. When tissue dose is unavailable, label the result a prescription-dose surrogate model and assess its limitations explicitly.

Long-running fitting and simulation scripts need fixed seeds, bounded workers, capped numerical-library threads, memory logging, checkpoints, a smoke-test option, and separate aggregate/shareable outputs. Before full execution, review the code, parse every script, and run the end-to-end smoke test. These are acceptance conditions for implementation, not work already completed by this review.

## 8. Disposition of the existing audit

| Existing item | Disposition | Required amendment |
|---|---|---|
| F1, existing links and unified page | Retain with clarification | The tumor brief mentions prior dependence and links to full details. The early explanation is already substantial. Retarget `early_how` as well as the tumor links, and test real navigation behavior. |
| F2, “only the early fit is hierarchical” | Replace | Describe exactly which levels are pooled. The tumor arm effects have an estimated common spread. Do not ban the word with a regex. Use `dd_cal.stan` and active drivers for the deployed DD model. |
| F3, deployment exclusion | Retain and extend | Modify `tools/deploy_preflight.R`, include metadata, test empty as well as missing content, and test a deployment-shaped bundle. |
| F4, avoid duplicate numbers | Retain and strengthen | Generate summaries from active artifacts. Runtime prior values, bounds, transformations, and conditional-DD settings are not all recoverable from Stan text. |
| F5, jargon and reader-first check | Retain, change verification | Use a clinical explanation and actual page tests. The generic manuscript checker currently demands an Abstract/Results/Discussion structure and has empty inventories. Its unadapted exit status is not an appropriate page acceptance test. |
| F6, menu/version ripple | Retain and extend | Include `test_v45_review.R` menu assertions, current deployment identity, early-context routing, and phone math/table checks. |
| L1, data not extracted | Retain with distinction | Reuse trial collection and schedules. Extract new late outcomes. Do not equate 68 raw acute rows, 51 primary acute arms, 120 PDFs, or 96 search hits with a late sample size. |
| L2, late α/β around 2 Gy matters | Replace with verified conditional result | Use Section 5. Both early modes were computed. `early_cal` is calendar information, not a missing statistical calibration step. |
| L3, repair bias and simulation | Retain concern, replace certainty | Bias direction is design-dependent. Fit and predict with the same repair evaluator. Use repeated design-specific simulation and report uncertainty, not one successful recovery. |
| L4, survivors and effective sample size | Major revision | Separate cumulative incidence, net risk, prevalence, and worst-ever toxicity. Choose estimator-specific likelihoods. Randomized denominators are not automatically invalid. |
| L5, same-dose/time-only check | Retain conditionally | Equal n and d imply equal direct late BED only if repair and other exposure components agree. Consequential injury is a structural sensitivity, not an identified mediation analysis. |
| L6, Stuschke prior | Retain overlap warning, correct source | Conditional 4.0 Gy, interval 3.3–5.0 at fixed γ₅₀ = 5, from EORTC 22791. It is not an independent universal 3.3–4.0 Gy prior. Check overlap for repair sources too. |
| L7, identification | Upgrade to a major gate | Examine risk slope, α/β, and repair jointly. Require practical recovery and prior-sensitivity evidence, not only varying fraction size. |
| Decision 4, mucosal or any-site at 2–5 years | Reject as the primary endpoint definition | Replace with the complete endpoint tuple and selection rule in Section 4.2. |

The two workstreams can proceed independently after their scope is accepted. There is no scientific need to delay read-only late-data inventory until the methods-page patch is deployed. A separate branch or isolated work area is useful for preventing an exploratory model from entering a documentation release.

## 9. Audit results and remaining readiness work

### 9.1 Methods verdict

| Dimension | Verdict | Reason |
|---|---|---|
| Conceptual validity | Partial | The reuse/feasibility approach is sound. The proposed pooled late endpoint is not a single estimand. |
| Implementation correctness | Partial | Current engine sensitivity results are reproducible. Proposed navigation and prior-table tests need correction, and late extraction cannot reuse all early assumptions. |
| Parameterization and reporting | Needs revision | Active DD conditioning, late prior provenance, joint identification, and uncertainty types must be explicit. |
| Reproducibility | Partial | Fit manifests and evaluators are useful foundations. Generic evidence/routing records and deployment markers are not current release certification. |
| Fit for purpose | Partial | A unified methods page meets the user need. A late fit can support research only to the extent its endpoint and exposure information support transport. |
| Strategic alignment | Pass with narrower scope | Prioritize methods access and endpoint-specific late feasibility. Do not promise that a new α/β estimate will materially expand or shrink the default list. |

**Final methodological recommendation: needs revision, not fundamental redesign.** The main weaknesses are fixable through a better endpoint contract, correct active-model documentation, and identification gates before fitting.

### 9.2 Recorded auditor triage

| Auditor | Routing and result for this request |
|---|---|
| Methods audit | RUN: endpoint, likelihood, identification, and implementation review above. |
| Paper/code audit | RUN: current explanations and proposed plan compared with active model/evaluator/navigation code. |
| Pipeline audit | RUN, scoped: runtime/export file and row-count checks, timestamp checks, actual-engine sensitivity, and isolated deployment-preflight tests. Three timestamp advisories remain. Full release suite and fitting pipeline not rerun. |
| Stale-constant audit | RUN: strict registry check on a scratch snapshot passed for 31 registered constants. Coverage/exclusions limit this result. |
| Manuscript audit | Not applicable: no manuscript or journal submission is being edited. The new internal review received citation, numeric, and equation checks appropriate to its content. |
| Tough grant review | Not applicable: no grant or proposal. |
| Biomedical data audit | Not applicable: the reviewed fitting inputs are aggregate trial-arm data, not patient-level extracts. Aggregate schema and denominator checks are addressed separately. |
| Radiation oncology review | RUN: tissue specificity, dose exposure, treatment context, consequential injury, and clinical interpretation reviewed against primary reports. |
| Surprising-result skeptic | RUN: “late α/β will change the list” and repair-bias assertions tested and qualified. |
| Evidence audit | RUN: existing strict ledger is incomplete for current `about.md`; see below. This review's source and calculation evidence is separately documented. |
| Simulation audit | RUN: schema validates, secondary literal scan flags unresolved registrations, and formal run-provenance certification is unavailable without the required trace files. |
| Subagent scan | RUN: independent code, statistical/methodological, and clinical-source reviews consolidated, including conflicting findings. |
| Reader-first audit | RUN: live page and explanations inspected. Generic manuscript-oriented executable is unsuitable without adaptation, and its inventories are empty. |
| Statistician review | RUN: estimator/estimand alignment, summary-data likelihoods, repeated outcomes, prior sensitivity, and identifiability reviewed. |
| Radiation epidemiology review | RUN: competing events, survivor selection, prescription versus tissue dose, and treatment-era transport reviewed. |
| Machine-learning practitioner review | RUN for relevant fitted-model validation issues: family-level separation, source leakage through priors, extrapolation, and distinction between computation checks and predictive validation. |

### 9.3 Executable results and limitations

The isolated deployment-preflight suite completed **13 test blocks and 47 expectations with no failures or errors**. This is a scoped result. The attempted broader application-suite run did not yield a completed result, so this review does not claim a full-suite pass.

The scoped pipeline runner found the declared runtime/export files and expected row counts, with no errors and three freshness advisories: the early exporter has a slightly later modification time than its three exported assets, all within the same displayed minute. This is not proof of numerical drift, particularly in a worktree, but it must be reconciled using source identities or a controlled rebuild before release. The runner had no custom code-quality patterns, so it does not replace the separate code review.

The existing strict evidence-ledger run reported **767 findings**. These include uncovered numeric tokens, obsolete anchors, references absent from the ledger, and missing external-proof mappings. This does not mean 767 scientific values are wrong. It means the ledger cannot currently certify `about.md`. Supplying the existing local PDF paths would address the proof-environment failures but not the coverage and obsolete-anchor problems.

The simulation registry passed schema validation. Its secondary scan reported 11 unregistered literals in the selected calibration script, some of which are seeds, indices, or tuning constants rather than scientific parameters. The formal run check could not find `simulation_input_trace.csv` and `simulation_provenance.json` in the inspected calibration output directory. Existing calibration outputs may still be useful, but this engine cannot certify their parameter provenance without those records.

The active application tree lacks a local `audit_manifest.yaml`. The parent project's manifest remains a template with generic declarations. The gate therefore cannot certify the active tree. These are specific readiness tasks for implementation: bind the intended deliverables and code to a current manifest, update evidence coverage and proof mappings, register simulation parameters and traces, and adapt reader checks to the application page. Do not manufacture passing records or replace source-level checks with an old `.suite_passed` marker.

No existing project scaffolds were overwritten during this review. For a plan-only request, listing the exact repairs is preferable to silently replacing project configuration. Journal-compliance and Word-layout requirements are not relevant to the requested internal Markdown review or the live application page. They become applicable if those materials are later converted into a journal manuscript or Word deliverable.

## 10. Acceptance checklist

The methods-page work is complete when all routes reach one deployed page, the active DI/calibrated-DD/early models are described correctly, priors and diagnostics trace to the retained fits, mathematical content renders on desktop and phone, and missing/empty deployment assets fail before upload.

The late feasibility work is complete when every included record has a source, endpoint tuple, valid denominator or estimator uncertainty, family overlap mapping, and exposure definition, and the design report says which parameters are identifiable. An explicit finding that the library cannot support a reliable α/β refit is a complete and useful feasibility result.

A late research fit is complete when computation checks, design recovery, family-held-out prediction, sensitivity analyses, source-overlap controls, and decision-impact calculations all have reproducible outputs. It is not ready for default adoption merely because its posterior is narrower than the existing assumption range.

Before integration, record a single decision using the prespecified gates in Section 7. This avoids repeated open-ended reviews: unresolved source discrepancies, nonidentification, leakage, failed diagnostics, or unsupported transport each have a defined outcome rather than triggering ad hoc model changes.

## 11. Reproducibility and source notes

The accompanying `2026_11_12_bayesian_methods_and_late_toxicity_evidence.zip` contains the sensitivity scripts and tables, source/output hashes, audit-engine logs, and domain-review notes. It contains aggregate calculations and provenance only, not patient-level data or licensed reference PDFs. Copies of the seven cited PDFs already held locally are co-located under `References/late_toxicity_review_2026_11_12/`, with a source manifest. Beitler 2014 was checked in primary full text, but its PDF could not be retrieved: the PMC endpoint returned an HTML access response and the Europe PMC fallback returned HTTP 403. This missing PDF is recorded explicitly, rather than represented by an HTML file with a PDF extension.

The research early package used for the sensitivity calculation has manifest SHA-256 `e80ecc45015651cf975d6ae7da43d26080ca7d48df539ecccdca93b4fa1a4355` and draws SHA-256 `a6fde4626c88a1872a84408d66e78ddd1460cdabb98734d06427ce186f07b4fd`. The R check used version 4.5.1 with `jsonlite` and `digest`, called `load_early_fit()` followed by `fitted_early_params()`, and then the same `v4_regimens()` function used by the home page. An independent implementation reproduced the baseline counts and sensitivity results.

Equation walk: the EQD2 conversion uses a dimensionless divisor; the complete-repair BED contrast and eligibility boundary have units Gy and reproduce the 29 × 2.30 boundary; the repair kernel is dimensionless, approaches complete repair as h approaches zero for distinct delivery times, and gives Q approaching D² as h grows without bound; the proposed binomial/logistic model has a dimensionless linear predictor. Every symbol is defined with its equation. These equations are source Markdown mathematics, not a claim of Word/PDF rendering verification.

## Verified References

1. Fu KK, Pajak TF, Trotti A, et al. A Radiation Therapy Oncology Group (RTOG) phase III randomized study to compare hyperfractionation and two variants of accelerated fractionation to standard fractionation radiotherapy for head and neck squamous cell carcinomas: first report of RTOG 9003. *International Journal of Radiation Oncology, Biology, Physics*. 2000;48(1):7–16. [PubMed](https://pubmed.ncbi.nlm.nih.gov/10924966/). Local proof: `References_refit_trials/eric_2026_10_22/Part 3/Fu_KK_2000_10924966_A_Radiation_Therapy_Oncology_Group_RTOG_phase_III_randomized_study_to.pdf`, relative to the parent HNSCC project. Relevant pages: PDF 8–9.

2. Beitler JJ, Zhang Q, Fu KK, et al. Final results of local-regional control and late toxicity of RTOG 9003: a randomized trial of altered fractionation radiation for locally advanced head and neck cancer. *International Journal of Radiation Oncology, Biology, Physics*. 2014;89(1):13–20. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4664465/). Relevant sources: Tables 3–4. Local PDF missing: primary full text was checked, but PDF retrieval was blocked as described in Section 11.

3. Mortensen HR, Overgaard J, Specht L, et al. Prevalence and peak incidence of acute and late normal tissue morbidity in the DAHANCA 6&7 randomised trial with accelerated radiotherapy for head and neck cancer. *Radiotherapy and Oncology*. 2012;103(1):69–75. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22398313/). Local proof: `References_refit_trials/eric_2026_10_22/Part 3/Mortensen_HR_2012_22398313_Prevalence_and_peak_incidence_of_acute_and_late_normal_tissue_morbidit.pdf`. Relevant source: Table 2.

4. Fu KK, Pajak TF, Marcial VA, et al. Late effects of hyperfractionated radiotherapy for advanced head and neck cancer: long-term follow-up results of RTOG 83-13. *International Journal of Radiation Oncology, Biology, Physics*. 1995;32(3):577–588. [PubMed](https://pubmed.ncbi.nlm.nih.gov/7790242/). Local proof: `References_refit_trials/eric_2026_10_22/Part 3/Fu_KK_1995_7790242_Late_effects_of_hyperfractionated_radiotherapy_for_advanced_head_and_n.pdf`. Relevant method: PDF 4.

5. Bentzen SM, Saunders MI, Dische S. Repair halftimes estimated from observations of treatment-related morbidity after CHART or conventional radiotherapy in head and neck cancer. *Radiotherapy and Oncology*. 1999;53(3):219–226. [PubMed](https://pubmed.ncbi.nlm.nih.gov/10660202/). Local proof: `References_refit_trials/eric_2026_10_22/Part 3/Bentzen_SM_1999_10660202_Repair_halftimes_estimated_from_observations_of_treatment_related_morb.pdf`. Relevant sources: methods and Table 1.

6. Stuschke M, Thames HD. Fractionation sensitivities and dose-control relations of head and neck carcinomas: analysis of the randomized hyperfractionation trials. *Radiotherapy and Oncology*. 1999;51(2):113–121. [PubMed](https://pubmed.ncbi.nlm.nih.gov/10435801/). Local proof: `References_refit_trials/pdfs/Stuschke_1999_review_10435801.pdf`. Relevant sources: PDF 5–6 and abstract.

7. Horiot JC, Le Fur R, N'Guyen T, et al. Hyperfractionation versus conventional fractionation in oropharyngeal carcinoma: final analysis of a randomized trial of the EORTC cooperative group of radiotherapy. *Radiotherapy and Oncology*. 1992;25(4):231–241. [PubMed](https://pubmed.ncbi.nlm.nih.gov/1480768/). Local proof is identified in the clinical-source inventory in the evidence archive. Relevant sources: late-effects population, Table XI, Figure 2, and tissue hot-spot description.

8. Denham JW, Peters LJ, Johansen J, et al. Do acute mucosal reactions lead to consequential late reactions in patients with head and neck cancer? *Radiotherapy and Oncology*. 1999;52(2):157–164. [PubMed](https://pubmed.ncbi.nlm.nih.gov/10577701/). Local proof: `References_refit_trials/eric_2026_10_22/Part 3/Denham_JW_1999_10577701_Do_acute_mucosal_reactions_lead_to_consequential_late_reactions_in_pat.pdf`.
