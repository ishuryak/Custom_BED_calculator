# HNSCC BED calculator: frequentist plan and completed Bayesian refit review

**12 November 2026 · application 4.6.1 · reviewed HEAD 79e9c812f91482fa57ffa1dddabe89dc329a6303**

## Assessment and recommended disposition

**Retain the completed Bayesian tumor and early-toxicity fits. Revise the frequentist plan before implementing its engine or running its large bootstrap.** The offline comparison is worthwhile, but the present plan and its internal audit overstate its readiness. Several proposed tests would reject correct code, one recovery requirement is statistically impossible to guarantee, and the proposed ordinary fixed-node quadrature can materially distort the likelihood on the project's actual trial sizes.

The Bayesian implementation has improved substantially since the 11 November review. The earlier custom-calendar error and changed-posterior replay problem are fixed in the checked cases. The underlying early posterior, likelihood and primary numerical design are unchanged. There are four remaining groups of low-cost corrections: incomplete numeric/schema validation at export/runtime, a methods-record check that does not check actual draw-file bytes, incorrect HYPNO observation/population metadata, and the tumor equation on the methods page.

Keep Fowler as the early-toxicity default and the trial-count fit as an explicitly selected research scenario. Keep the frequentist comparison offline. A likelihood-only fit would not by itself supply a clinically calibrated toxicity limit, establish transport to hypofractionation, or justify changing the app's default.

This document consolidates the findings into an execution plan. Each issue states the evidence, effect, exact repair, dependencies and acceptance criteria. The later checklist is the completion gate; it does not require another general review after those criteria pass.

## Scope, versions and evidence boundaries

All application-relative paths below refer to **Custom_BED_calculator_v45**, inside `C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator`. Despite the folder name, this is version 4.6.1. All **1,542 tracked files** matched the canonical Git revision and the physical checkout. There are 120 changed tracked files since the previously reviewed 5384f60 revision. The recorded local origin branch agrees with HEAD; this is not a claim of a fresh remote fetch.

The two current frequentist documents are physical working-copy files **absent from that Git tree**. They were captured separately, so this review includes the latest plan rather than only the committed release:

| File under docs/development | Reviewed identity |
|---|---|
| 2026_11_12_frequentist_comparison_plan.md | Revision 2; SHA-256 `44d2324f5e7820d9b9155cef43c17f8d5d32bc09e5d25623f3c369c95a0cf292` |
| 2026_11_12_frequentist_comparison_plan_audit.md | SHA-256 `c1b445873cdb7c7026d727b960c247216326c5e1321dbce567f8002defd95deb` |

“Plan” in the findings means the first file. “Audit” means the second. Their line references apply to these captured versions.

The live site was freshly opened. It identifies **4.6.1**, with deployed revision **7ec2fafadde0747f8175163569959abcb828cef4** on its methods page. Main showed 69 qualifying dose/fraction points. The early-fitting link opened the correct section of the Bayesian methods page, with rendered mathematics, the current summaries and the source link. This was a targeted live check, not a repetition of the project's recorded 73-item walkthrough.

This review made no production, source-data, Git, deployment or app-default changes. Bounded scripts, deliberately invalid fixtures, a source snapshot and test outputs were kept in the review workspace. No new MCMC, full frequentist refit campaign or bootstrap was run. One inexpensive reoptimization of the 21 nuisance intercepts at the already archived tumor parameters independently checked the historical likelihood value.

## What is correct and what can be retained

### Model and data contracts

The planned tumor likelihood has the correct link orientation: for tumor control probability `TCP = exp(−exp(η))`, the complementary log-log link is on **failure**, not control. Its dose term uses the published fraction components, not an invented uniform-fraction replacement. The tumor time convention is inclusive duration. The research early model uses elapsed time between first and last fraction, an exact onset hinge and explicit fraction times for repair.

The current early exposure is:

\[
E=D+Q/b-r\max(0,T_{\mathrm{elapsed}}-T_k),\qquad
Q=\sum_i d_i^2+2\sum_{i<j}d_i d_j\exp[-\log(2)(t_j-t_i)/h],
\]

where b is α/β, fraction times and h are hours, and r is BED units/day. The primary scenario fixes onset at 7 elapsed days and h at 3.2 hours. Trial-centered exposure enters a logistic count model with a trial-baseline distribution, grading-group shift, coarse planned-therapy contrast and arm variation. Integrating the random effects is the correct frequentist counterpart for this model; optimizing one unrestricted latent effect for each arm is not that marginal likelihood.

The intended count designs remain **45 tumor arms/strata in 21 fitting units** and **51 early arms in 22 families**. Early denominators sum to **6,042**, with **2,945 events**. Tumor nominal analyzed N sums to **8,279**. All 45 tumor input rates are time-to-event estimates—39 actuarial/Kaplan–Meier and six cumulative-incidence estimates—rounded into pseudo-counts. They are not observed binomial success totals. A frequentist calculation retaining this likelihood must retain this limitation too.

### Closure of the preceding review

| Earlier concern | Current evidence and disposition |
|---|---|
| Unsupported +5 risk statement and universal chemotherapy multiplier | Corrected in inspected app/brief/About text. The margin has no established clinical-risk meaning, and γ is a coarse pooled strategy adjustment. Retain. |
| Lost schedule mixtures and source qualifications | CAIR variants/weights, EORTC's unresolved pattern, DAHANCA's gap/as-treated status, structured therapy and Tandon's assumed equality are now preserved. Retain; B03 corrects a different underlying HYPNO record. |
| Research repair ignored gaps or used the candidate calendar for its comparator | Fixed in independent R checks. The skipped-day example gives 31.11404545 BED; the six-day comparator gives 29.73698842 BED. A manual duration inconsistent with its dates is unavailable. |
| Exporter accepted sensitivity fits under the primary label | U04/U15/U17/U19 are now refused; valid U02 is accepted. Retain that restriction. B01 closes remaining invalid numeric metadata. |
| Reordered early draws could replay an old exact scenario | Fixed: early draw/package digests are recorded, and the original scenario refuses the changed file. A valid custom-calendar scenario replayed successfully. |
| Per-chain energy and cache integrity gaps | NaN/Inf in each chain position and wrong chain count now refuse. Changed NPZ bytes now refuse. All 52 current cache hashes and identities match their manifest. |
| A grade-1 mucositis row could corroborate a severe endpoint | Now refused by explicit endpoint equivalence rules. Preserve source-specific adjudications and positive controls. |
| Active early model and unused Fowler controls unclear | Current code and live methods/banner identify the active model; the implementation hides the inactive fields and records prior use. Related maintained R checks are included in this review's test run. |

The retained chain-recheck table has **114 fits checked, 111 passing**. The three failures are the already disclosed binomial U01/U05/U07, not the primary overdispersed fit. All 114 row identities match the corresponding current 62 step-5 units and 52 step-6 caches; each records four finite energy diagnostics. All seven step-6 report-file hashes match. This is verification of the retained record and artifact identities; this review did not reopen the external Linux raw-chain archive.

The current primary early draw CSV is byte-identical to the previously verified export. The primary exposure evaluator, Stan model, fraction and outcome tables are also unchanged. The complete rebuilt primary Stan input dictionary was directly compared with the previously verified generation and is identical, despite changes to validation code. The earlier independent held-out score and conditional-interval reproduction can therefore be retained for that same scientific representation; another wholesale score replay or sampling campaign is unnecessary merely because validation code or text changed.

## Findings and priorities

**P1 below blocks the proposed frequentist production run, not continued use of the unchanged Bayesian fit.** P2 identifies a required repair to the affected feature or scientific description. Hypothetical corrupt fixtures are clearly distinguished from the actual shipped artifacts, whose inspected numeric values are valid.

| ID | Priority | Required correction |
|---|---|---|
| F01 | P1 | Replace the unsupported fixed-quadrature accuracy claim and repair the independent R comparison contract |
| F02 | P1 | Correct likelihood sign, design/hash identity and rounded-output reference tests |
| F03 | P1 | Replace the fixed-cluster, inflated-N recovery requirement |
| F04 | P2 | Define the comparison as estimation/prior sensitivity with matched support; repair interpretation labels |
| F05 | P1 | Specify valid profile/bootstrap procedures, boundaries and failed-replicate treatment |
| F06 | P2 | Specify reproducible restart, complete test gates and a measured compute budget |
| B01 | P2 | Validate finite typed diagnostics and exact runtime model definitions |
| B02 | P2 | Make methods-record freshness check actual draw files and safely handle invalid records |
| B03 | P2 | Correct HYPNO's ascertainment and evaluable-population descriptions |
| B04 | P2 | Display the actual smoothed tumor hinge on the methods page |

## F01. Use demonstrably accurate marginal integration and genuinely equivalent cross-checks

**Locations:** Plan §§3.1–3.2, especially lines 40–58; T3/T6 at lines 113 and 116; proposed fit/frequentist/likelihood.py and glmer_check.R.

**Confirmed problem.** Ordinary Gauss–Hermite quadrature with 30 nodes is called “essentially exact”; early likelihood integration is assigned 20 × 20 nodes without a convergence rule. Accuracy depends on the posterior concentration and location of the latent effects, not just the dimension of the integral. The real arm sizes already defeat that claim:

| Bounded numerical fixture | Ordinary fixed-node log-likelihood error versus accurately scaled adaptive integration |
|---|---:|
| Tumor arm N=750, y=525, TCP at zero arm effect=0.7, σ=0.25, 30 nodes | −0.03395 |
| Same arm, σ=1, 30 nodes | −6.87156 |
| Early RTOG 0129 family, N=361/360, y=141/119, τ=1, arm spread=0, 20 outer nodes | −7.97402 |

These errors can alter an optimum and swamp a 1.92 profile cutoff. They are counterexamples to the planned numerical method, not newly fitted clinical results. A positive control using mode/curvature-adapted quadrature matched the independent integration in the tested single-arm and outer-integral fixtures within approximately 2e−14. The primary methods literature and the independent calculation both support an adaptive, accuracy-checked implementation. [Bilodeau, Stringer and Tang](https://arxiv.org/html/2202.07864).

**Fix.**

1. Implement normalized log marginal likelihoods using log-sum-exp and stable binomial/cloglog terms. Include the normal density and change-of-variable factor correctly. Provide analytic zero-variance branches.
2. Use adaptive quadrature or another independently validated integration method. For F-E, the outer log integrand must contain the **sum of inner log integrals**, not a product of independently averaged family effects or a mean of log likelihoods.
3. Register a convergence ladder, for example 20→40→80 adaptive nodes, extending when needed. Validate the entire family and total likelihood, not only isolated arms. Test concentrated, displaced, large-variance and near-boundary cases, including the three counterexamples above.
4. Accept a reported optimum/profile endpoint only after a higher-accuracy evaluation changes total log likelihood and the profile log-likelihood drop by less than 1e−6, changes finite reported contrasts by less than 1e−4 in their natural units, and preserves crossing/boundary classifications. These are proposed engineering tolerances, well below the 1.92 cutoff; record them and revise only with a documented numerical reason. Cross-check representative integrals against an independently scaled adaptive integrator to 1e−9 absolute log-likelihood error.
5. Check finite-difference or independent derivatives at interior points and projected gradients/KKT conditions at bounds. A positive-definite full Hessian is not required for a legitimate boundary optimum; a flat/rank-deficient profile is a scientific diagnostic, not permission to force a curvature estimate.

The F-T R formula also needs correction. It says α/β is fixed, then fits independent dD and dQ coefficients. At fixed b and Tk, use a single predictor **B = dD + dQ/b**, plus the specified hinge and centered time. Its coefficient is −α. Alternatively describe a free-ratio cross-check and apply the same sign/support constraints. An unconstrained glmer/glm optimum cannot be required to equal a constrained Python optimum when their feasible sets differ. For F-E's zero-random-effect GLM check, choose known interior fixtures satisfying the coefficient signs; test boundary/non-identification behavior separately.

Use R's adaptive quadrature for the eligible single-random-effect comparison, but evaluate the project's own full likelihood at R's returned coefficients before comparing objectives: lme4 documents convention differences in reported deviance across approximations. The two-random-effect Laplace result remains descriptive, including possible differences in fixed effects as well as variances. [Official glmer documentation](https://lme4.github.io/lme4/reference/glmer.html).

For coefficient comparisons use a registered combined absolute/relative tolerance, such as 1e−8 absolute plus 1e−4 relative, so a valid zero does not make a relative-error test undefined. Keep objective and integration tolerances separate from coefficient tolerances.

**Acceptance:** both independent positive controls and deliberately low-accuracy controls behave as specified; total/profile likelihood convergence is logged; R and Python compare the same model, support and likelihood convention. No full optimization/bootstrap starts before this gate passes. Existing Bayesian fits are unaffected.

## F02. Repair the plan's reference values and distinguish source, design and run identities

**Locations:** Plan T1/T2/T4, lines 111–114; Audit findings 1–2; docs/refit_explainer/aicc_check.py:29–33 and main_free.json; fit/fit_refit.py:93–139; early U02 draws.csv.gz.

**T2 has the wrong sign.** The stored 5342.406617400148 is **negative** log likelihood without the binomial coefficient. The plan tells the implementer to add the coefficient to this positive value. Correct values, independently calculated from all 45 rows, are:

| Quantity | Value |
|---|---:|
| Σ log C(N,y) | 5211.115489844190 |
| Archived negative log likelihood, coefficient omitted | 5342.406617400148 |
| Full negative log likelihood | **131.291127555958** |
| Full log likelihood | **−131.291127555958** |
| Incorrect “add to positive NLL” target | 10553.522107244338 |

Reoptimizing only the intercepts at the archived shared parameters reproduced the omitted-constant NLL to about 4.3e−8. Thus this is a plan/audit arithmetic error, not an error in the historical optimizer. Replace T2 with `full_nll = stored_nll − sum_log_choose`, or compare `full_ll = −stored_nll + sum_log_choose`. Name the sign convention in code and output.

**T1 compares unlike hashes.** The CSV file hash, serialized Stan-data hash and early U02 run identity hash different objects. The latter also includes settings, model, sampling configuration and software versions. A likelihood-design hash cannot equal that full identity. Check these layers separately: exact input bytes; canonical data arrays and ordered unit/arm IDs; reconstructed scientific/run identity. Test changing one count, one row order, one time convention and one setting independently. Preserve both the file digest and the canonical design digest in the new manifest.

**T4 exceeds the saved early precision.** The retained U02 columns are written at six significant digits. Recomputing log odds ratios from their rounded parameter rows differs from the saved rounded generated quantities by up to **1.51e−5**; **10,641 of 16,000** values exceed the proposed 1e−6 tolerance. A correct evaluator will fail. Use original full-precision chain values when available, or use a documented bound on decimal input/output rounding. For this unchanged six-significant-digit fixture, a 2e−5 absolute tolerance covers the measured discrepancy; the independent unrounded formula tests should retain tighter tolerances. Do not globally loosen all numerical checks.

The tumor half of T4 currently passes the proposed 0.1-Gy median check, with maximum observed difference **0.09225 Gy**. However, the app's 4,000 paired rows contain **1,000 distinct DI draws repeated four times**, not 4,000 independent draws thinned from 8,000. Treat the ten summary comparisons as a versioned distribution sanity check, not an exact cross-language identity or a guaranteed tolerance for future exports. Compare R/Python values on the **same rows** for numerical correctness; compute and label posterior contrast intervals from the chosen draw artifact. The existing tumor contrast summary has q10/q90, not the requested 95% endpoints.

**Acceptance:** T2 reproduces the full NLL above within 1e−3 and refuses a sign/constant mutation; T1 distinguishes file/design/run changes; early rounded and unrounded tests use separate justified tolerances; tumor same-row calculations agree tightly and the unique-draw count is stated correctly. All existing posterior files are retained.

## F03. Replace the recovery test that confuses more patients with more independent random effects

**Locations:** Plan T5, line 115; Audit finding 4's sigma-bias experiment; proposed test_frequentist.py.

The plan multiplies N by 1,000 on the same trial/arm design and requires every identified parameter to be recovered within 2%. That reduces binomial observation noise but leaves the same finite number of realized family and arm effects. It does not consistently recover their population variances from one simulated dataset. Weak dose/time contrasts and the tumor's many nuisance intercepts remain too.

Even in the easier case of 22 exactly observed independent normal effects, the simple ML spread has an approximate central 95% ratio to the true spread of **0.684–1.270**; its chance of falling within 2% is only about **10.1%**. Requiring the more complicated mixed model to meet a deterministic 2% target is not a valid code-consistency test.

**Fix the testing hierarchy.**

1. Use deterministic likelihood, derivative and contrast identities plus zero-random-effect cases as unit tests. Choose parameters away from a ridge, with full-rank predictor design and interior constrained solutions.
2. Use repeated simulated datasets on the real design to estimate bias, optimizer success and coverage. Record actual Monte Carlo uncertainty. Do not require one random realization to equal the generating value within a tiny relative tolerance.
3. If studying asymptotic recovery, explicitly increase the relevant independent information—families, arms/contrasts and within-unit information as appropriate—and explain the nuisance-intercept regime. Increasing family count alone does not automatically remove the tumor variance's incidental-parameter problem.
4. Separate the proposed sigma-bias study from a required green test. Report its measured direction and uncertainty; do not require a predetermined downward result in every seed/configuration. A posterior-fixed σ sensitivity is a declared hybrid sensitivity, not a proof that likelihood-only inference has been corrected.

**Acceptance:** a wrong sign or missing dose term is caught by deterministic tests; correctly generated finite-cluster variability is not mislabeled an implementation failure; repeated-simulation summaries include seeds, truth, replicate counts, bias/coverage and their Monte Carlo errors. T5 no longer claims that inflated N alone identifies all mixed-model parameters.

## F04. Define a fair comparison and avoid treating arbitrary search bounds as evidence about priors

**Locations:** Plan Q1/Q2 and §§3.1, 3.3, 4, 8; particularly lines 18–21, 42, 62, 85–94 and 135–138. Audit findings 3, 5 and 7.

The statement that a maximum-likelihood estimate and posterior median differ “by exactly the prior effect” is too strong. The procedures also differ in nuisance profiling versus integration and in mode versus median summaries. Both retain the normal random-effect distributions and the observation-model assumptions. Changing parameter support between fits adds another difference. This is useful estimation/prior sensitivity analysis, not a pure decomposition of biology into “data” and “assumptions.”

The primary proposal restricts Tk to the observed 12–54 days and calls that a physical constraint, while the secondary uses the Bayesian 0–60 interval. The observed minimum/maximum duration are design features, not physical onset limits. Excluding onset below 12 days also hides the very low-onset ridge the plan intends to study. Conversely, an onset later than all observations can be scientifically unidentified rather than disproved.

**Recommended replacement contract:**

- Primary comparison: the same parameter support as the Bayesian target, while removing the population-parameter priors and explicitly retaining the normal mixing distributions. State which nuisance parameters are profiled and which latent effects are integrated.
- Secondary sensitivity: expanded scientifically defensible support, with computational bounds distinguished from physical limits. Expand bounds when the likelihood has not demonstrably closed; report unbounded/non-identifiable regions rather than manufacture finite endpoints.
- Describe the result as “likelihood-based and Bayesian estimates under the same observation model,” with matched-support and expanded-support results separate.
- Replace automatic “data-dominated/prior-dominated” conclusions based solely on interval-width ratios or a search limit. Report the ratio as an **analyst-defined descriptive measure**, alongside profile shape, prior-sensitivity results, bound dependence and the actual change in schedule contrasts. A numerical search failure cannot establish prior dominance.
- Specify transformations separately. A log-width for a positive parameter can be a chosen descriptive scale, but κ has a half-normal prior, not a lognormal one; zero/boundary values make naive logs undefined. Signed contrasts cannot be put on a log scale. Derived contrasts have no independent assigned prior: their prior distribution must be generated from the full declared joint prior and the exact contrast function, or be marked unavailable. Do not invent a prior interval by combining marginal endpoint limits.
- Report `(λ+g)/α` in BED/day, and its draw/estimate-specific EQD2 conversion when useful clinically, alongside λ+g in inverse days. Keep the units distinct.

It is reasonable to leave calibrated DD out of this scoped comparison. Replace “no likelihood-only analogue” with “not included in this plan.” A conditional plug-in calibration and repeated two-stage bootstrap could be defined; that would be additional work and is not the same as the already available free-DD ML fit. No DD refit is needed for this review's proposed scope.

**Clinical qualification required in the output:**

> These intervals describe parameter and contrast uncertainty under the specified likelihood and design. They do not quantify the full uncertainty of the published clinical endpoints or establish transport to a new regimen or population. Holding the likelihood constant makes this a sensitivity analysis to estimation and prior assumptions, not a new clinical validation.

The binomial bootstrap cannot restore censoring-aware uncertainty in the 45 reconstructed tumor pseudo-counts, missing early mucositis peaks, delivered chemotherapy histories or patient mucosal dose distributions.

**Acceptance:** matched-support and expanded-support analyses have separate identities and labels; numerical limits are never used as a clinical conclusion; every contrast has defined units and an uncertainty source; D2 reporting covers **all prespecified schedule contrasts**, not only whether α/β receives one label. Publish the offline report regardless of direction, and link it from the methods page when complete. Keep the app's numerical packages unchanged under this plan.

## F05. Specify profiles, boundary inference and failed-replicate handling before fitting

**Locations:** Plan §4, lines 69–83; T8, line 117; implementation steps 4–6. These are currently incomplete specifications, not observed defects in a frequentist implementation, which does not yet exist.

**Profiles.** Saying “a 95% profile interval for every quantity” is insufficient for functions such as ΔBED and κΔE. For each target q, define the constrained optimization `max ℓ(θ) subject to q(θ)=q0`, nuisance reoptimization, multiple starts, feasible branches and endpoint bracketing. Box-constrained L-BFGS-B alone does not enforce that nonlinear equality. Use an analytical reparameterization where safe or a constrained optimizer, then verify the target residual and objective independently. Allow disconnected confidence sets, one-sided/unbounded regions and undefined ratios at α or κ approaching zero. Do not collapse them into a reassuring finite min/max interval.

**Boundaries and ridges.** The 1.92 cutoff is the regular one-parameter approximation. Adding an “approximate” flag does not establish coverage at a boundary, separation, or the λ–g ridge. A flat profile or nonidentified onset is a valid completed result, provided numerical accuracy is demonstrated. The available boundary theory gives specific counterexamples to routine bootstrap validity; it is not a universal guarantee that every percentile bootstrap here fails or succeeds. [Andrews, publisher abstract](https://doi.org/10.1111/1468-0262.00114).

For the proposed λ=0 likelihood-ratio test, generate from the **fitted null**; refit both null and alternative for each simulated dataset, including the alternative's full Tk search; use the same accurate integration, constraints and optimizer rules. The ordinary bootstrap generated at the alternative MLE is not this null distribution. State whether this additional test is necessary for the offline comparison and budget it separately. If omitted, do not report a regular χ² p-value as a substitute.

**Bootstrap algorithm.** State what is regenerated: F-T draws fresh independent arm effects with fitted fixed unit intercepts and then binomial counts; F-E draws fresh family and arm effects and then counts, keeping the declared design fixed. Record this target explicitly. A different conditional/resampling scheme answers a different question.

Replace “intervals with and without failed replicates.” A failed optimization has no defensible estimate to include in a percentile interval. Save every attempt, retry deterministic failed replicates under a predeclared stronger strategy, distinguish numerical failures from legitimate boundary solutions, and retain any unresolved failures with their reason. Never silently drop them or replace them with zero/the starting value. Report successful-replicate summaries only as conditional/provisional when failures remain, with sensitivity bounds where possible. Exceeding the declared 5% unresolved-failure limit blocks a final interval claim and triggers the specified numerical remediation, not an open-ended new general review. A rate below 5% is not by itself proof that selection is harmless.

Specify any variance bias correction mathematically and report its uncertainty and boundary behavior; do not imply that correcting one σ point estimate automatically repairs the contrast intervals. Keep the requested σ=posterior-median and upper-quantile scenarios clearly labeled as externally fixed sensitivities.

**Acceptance:** synthetic interior profiles recover known likelihood intervals; genuine boundary/ridge fixtures report their nonregular state; derived-contrast constraints are met within 1e−7 in their stated units and remain stable under tighter integration; restart/retry tests preserve replicate identities; the null LR test, if retained, simulates under the null and refits both models; every reported interval names its method and qualifications. Golden tests should pin objectives and identifiable predictions, not demand identical arbitrary coordinates along a ridge.

## F06. Make the planned run restartable and its completion criteria complete

**Locations:** Plan §5, line 101; tests/order at lines 111–120; §7 and stopping condition line 142. The planned append-and-skip-if-present bootstrap would reintroduce a failure mode already repaired in the Bayesian pipeline.

**Required implementation contract:**

1. Before writing a result, make a run identity covering input-file and canonical-design hashes, ordered IDs, likelihood/evaluator code, model/support, quadrature accuracy policy, optimizer/start policy, target/profile definitions, RNG, software environment and seeds. Preserve the reviewed plan revision and subsequent approved changes.
2. Store each replicate as an atomic completed record, with seed, generating target, estimates/objective, convergence status, active bounds, numerical checks, attempts and output digest. An incomplete or malformed record is not a completed replicate. Use per-replicate files or another transactional scheme; do not treat a partially appended CSV row as successful work.
3. Resume only when identity and completion checks match. A numerical failure may be retried according to policy; a valid nonidentified/boundary outcome should be retained as such. Never rerun completed healthy replicates merely because the process restarted.
4. Publish assembled tables/figures and one generation manifest only after all required pieces verify. Reject mixed old/new generations. A changed dataset, model, integration tolerance or constraint set starts a new generation.
5. The final gate must include **T8 as well as T1–T7**, plus restart, malformed-output, wrong-generation and publication-interruption tests. Existing T1–T6 wording does not close the newly required inference controls. An lme4 skip is permissible in the ordinary app suite, but it does not count as a passed required independent cross-check for this analysis.
6. Measure the pilot's CPU/wall time, memory, node escalation, optimizer starts and retry rate. Include profiles, sigma-bias simulations, null-LR simulations if retained, and all bootstrap refits in the budget. The current “1–3 worker-CPU hours” is unverified. Keep B=1000 tumor/500 early as provisional targets; also report Monte Carlo uncertainty of interval endpoints and increase B only for a prespecified precision need after authorization.
7. Put pure numerical/identity tests in a non-sampling CI job. Tests of run identity should use a controlled software-version/installation-path fixture rather than require a working CmdStan installation when no sampling is requested. The current recorded release had two such dependency failures; keep them distinguishable from model failures.

**Acceptance:** uninterrupted and interrupted/resumed pilots produce the same accepted replicate records and summaries; truncated rows, edited bytes and stale identities refuse; a deliberately failed optimizer is counted and retried without losing its seed; zero replacement sampling occurs on a valid no-change resume. The measured run manifest accounts for every planned expensive operation. These are code/pilot requirements before the long run, not reasons to rerun the Bayesian fits.

## B01. Complete typed, finite validation at the source-export and runtime boundaries

**Locations:** tools/early_tox/export_fit_for_app.py:57–78; R/uncertainty.R:195–224, especially :201 and :207–208.

The per-chain diagnostic gate is improved, but summaries crossing later boundaries are still accepted under invalid types/values. A private source-unit probe with rehashed diagnostics accepted R-hat=−Inf and ESS/E-BFMI=Inf. Independently, the R loader accepted each of these copied-package mutations:

- an unknown `early-bed-fit-10` evaluator because `startsWith(..., "early-bed-fit-1")` matched;
- ESS values set to the text `unavailable`;
- R-hat set to −1;
- a missing model ID.

The valid package was accepted, and a false gate was refused, as intended. The shipped diagnostic values are finite and valid; these tests do **not** establish that the current posterior is defective.

**Fix.** Validate each required scalar's JSON type, length, finiteness and scientific support before comparing thresholds. Reject booleans masquerading as numbers, strings, nulls, arrays, NaN and infinities. Require a finite positive R-hat, finite positive ESS/E-BFMI and the exact declared integer divergence/gate contract. Reuse one documented threshold policy without accidentally tightening/changing the actual sampling target. Use an exact evaluator ID field separate from its description, or parse the current ID token exactly; do not prefix-match version 10 to version 1. Require a nonempty well-formed model/source identity and supported row-order/point rules. Return a clear unavailable-package result at runtime.

**Acceptance:** valid U02 and correctly rounded point values still load/export; each malformed value at each boundary refuses before publishing or calculating; unknown versions and absent identities refuse; the current numeric draw CSV is unchanged. No new fit or raw-chain rerun is required when the retained valid diagnostics meet the strengthened boundary validator.

## B02. Validate the actual posterior bytes when declaring the methods record current

**Locations:** tools/build_model_methods.py:143–166, 194–214, 230–259; R/methods_page.R:16–31; tests/python/test_model_methods.py; methods-page deployment checks.

The generator does valuable checks: it reconstructs scientific input/identity hashes and refuses a changed prior/model. It then takes exported-draw hashes from manifests without confirming the actual CSV bytes. In private copies, changing a tumor DI α row or an early r row to 100, while leaving manifests untouched, still made `--check` print that model_methods.json was current. The runtime methods loader similarly compares record hashes to manifest hashes rather than a validated loaded draw object. The ordinary calculation loaders have additional byte checks, so this finding is a false methods/provenance attestation, not proof that the app would calculate from those damaged files.

**Fix.** Before building or checking the methods record, verify actual DI, DD and early draw-file SHA-256 values, required columns/rows, paired-row IDs and source/package identity consistency. Validate the record schema and types before dereferencing nested fields; invalid metadata should render a clear unavailability message rather than an uncaught error or misleading technical table. Prefer using the same validated artifacts/objects as the calculator. Add the builder freshness check to release preflight/CI, not just a nonempty-file check. Publish the record from the same validated generation as its draws.

**Acceptance:** unchanged files reproduce the current record; independently changing any draw file without its manifest causes a named refusal; changing the source scientific identity also refuses; malformed/missing fields fail safely; regenerating a valid record does not alter posterior values. Tests must mutate real temporary files, not only the JSON description.

## B03. Correct HYPNO's source population and follow-up description

**Locations:** data/early_tox/source/ascertainment.csv:30; draws/early_fit_trials.csv:66–67; export_fit_for_app.py:117–129,174–187.

The record says `randomized and treated` and exports `follow-up after RT not stated`. The primary report explicitly states weekly early-reaction assessments during and immediately after RT, with the first scheduled post-RT visit about two months later (printed p. 391). It identifies 753/792 patients as evaluable for early adverse events—373 HFX and 380 NFX (p. 394). The current mucositis counts **190/373 and 208/380 are correct**. [HYPNO primary report](https://pubmed.ncbi.nlm.nih.gov/41711624/).

**Replace the population with:** `Early-adverse-event evaluable patients, analyzed by randomized arm (753/792 overall; HFX 373, NFX 380).`

**Replace the timing note with:** `Weekly during and immediately after RT; first scheduled post-RT follow-up about two months after RT. Duration/frequency of the immediate post-RT assessments is not fully specified; departmental follow-up policies applied.`

Label the selected outcome as maximum grade ≥3 mucositis during and immediately after treatment. Do not infer complete peak ascertainment or manufacture an exact eight-week observation window from the approximate routine visit. Preserve HYPNO's external-benchmark role and the existing therapy sensitivity scenarios.

**Acceptance:** source and exported display contain the corrected population/timing, unchanged event counts and external-only role. Primary totals and the 12/22-any-visit versus 2/22-adequate distinction stay unchanged. Rebuild dependent display/provenance artifacts and their hashes together; do not resample or rescore unchanged numerical inputs merely because this descriptive field changed.

## B04. Make the public tumor equation match the fitted equation

**Location:** how_fitted.md:27, compared with fit/stan/di.stan:52,80,84; fit/stan/common_functions.stanfunctions; R/bed_models.R:92. The mismatch is visible on the live methods page.

The page writes an exact `max(0, T−Tk)`, whereas the tumor fit and app use `smax(x) = log(1+exp(3x))/3`, evaluated stably. At x=0 this is ln(2)/3 = 0.231049 days rather than zero. With the illustrative rounded λ=.028 and α=.044, that is about 0.147 BED units in the time-loss term. This is a documentation discrepancy near onset, not evidence that the already fitted numerical evaluator changed.

**Fix and acceptance:** show smax in the tumor equation and define its constant and stable evaluation. Keep the **exact** hinge in the research early equation. A focused text/render test must distinguish the two; the live page must display the corrected formula after any separately authorized documentation release. Retain every fitted value and draw.

## Decisions: recommended answers and reasons

| Plan decision | Recommended answer |
|---|---|
| D0: offline check or app switch | Offline check. Preserve current app choices; interval/support semantics differ and clinical validation is unchanged. |
| D1: primary constraints | Matched Bayesian support first; expanded-support sensitivity second. Do not call observed-duration limits physical constraints. |
| D2: public methods update | Link the complete report whatever its result. Summarize meaningful changes/stability across all prespecified schedule contrasts, not only an α/β label. |
| D3: bootstrap sizes | Retain 1000/500 as provisional budgets, with measured pilot cost, explicit extra simulations and endpoint Monte Carlo precision. Do not start them under the present uncorrected plan. |
| D4: calibrated DD | Exclude from this scope, with accurate “additional two-stage analysis would be needed” wording. Existing calibrated DD remains unchanged. |

The review itself does not authorize the expensive run or a deployment. Once the user has authorized implementation/execution and the concrete gates below pass, proceed within that authorization; no additional general-review approval is prescribed.

## Ordered execution roadmap

### Phase 1 — Freeze and correct the specification

1. Preserve the reviewed source revision, captured plan/audit, posterior packages and retained recheck/cache/report manifests. Commit the revised plan through the project's normal process so the actual run is attributable to a specific specification.
2. Apply F01–F06 to the plan and its output template. Correct the old audit's sign, inflated-N, quadrature and prior-attribution conclusions explicitly rather than leaving it as a conflicting governing document.
3. Resolve D0–D4 using the recommendations above. Define matched/expanded supports, target contrasts, interval methods, quadrature accuracy, retry policy, boundary outcomes and computational budget before seeing new fit results.

**Exit:** the plan states an implementable mathematical and numerical contract, and no test compares different hashes, incompatible estimands or unattainable precision. No expensive fitting is needed for this phase.

### Phase 2 — Close current Bayesian packaging and source-description issues

Implement B01/B02 guards first, then B03/B04 source/text corrections. Extend the existing export, model-methods and v4.6/v4.6.1 R tests with the deliberate negative cases above. Regenerate only affected displays/manifests/metadata. Keep the early draw CSV and tumor draw CSVs byte-identical.

**Exit:** every malformed package or changed actual draw file refuses; actual valid packages still work; corrected source/formula displays match the model; default Main remains 69 points and research Main nine points. This is independent of whether the frequentist experiment proceeds.

### Phase 3 — Build and validate the likelihood engine

Implement F-T/F-E pure numerical functions, stable integration, zero-variance limits and exact contrast functions. Use one shared canonical design import; preserve both source and scientific identity checks. Add corrected T1–T6 and T8 controls before recording a golden optimum. Perform the independent R check in the designated analysis environment where lme4 is present.

**Exit:** accurate full likelihoods and derived contrasts agree with independent references; the fixed-GH counterexamples are detected; valid constrained interior and boundary fixtures have appropriate outcomes. No failure is “fixed” by changing the clinical counts to match an optimizer.

### Phase 4 — Implement inference and restart; run the bounded pilot

Implement constrained profiles, bootstrap generation, deterministic replicate seeds, saved attempt states, boundary classifications and atomic generations. Run the declared B=20 pilot after unit gates pass. Interrupt and resume it deliberately in a temporary output location. Measure time and memory, numerical escalation and failures. Test a forced failure and a genuine nonidentified outcome separately.

**Exit:** resumed/uninterrupted accepted output matches, scientific identities reject stale work, and cost/failure estimates cover all planned operations. A pilot may legitimately identify a flat profile; that is a reportable result, not a demand for more arbitrary starts until a finite interval appears.

### Phase 5 — Execute only the authorized, measured new work

Run the specified frequentist fits, profiles and bootstrap. Add null-LR or repeated-simulation jobs only if retained in the approved contract and budget. Reuse validated completed replicate records. Escalate a numerical method only for a named failed accuracy/convergence criterion, creating the appropriate new generation or declared attempt record.

**Exit:** every result has source/run identity, numerical checks, convergence/boundary status and an interval or explicit reason no regular interval is justified. Both successful and unsuccessful attempts remain auditable.

### Phase 6 — Report, verify integration and close

Produce one report containing every prespecified parameter and schedule contrast, matched/expanded-support distinctions, confidence/credible interval labels, profile plots, the clinical-estimand caveat, simulation/failure summaries and retained-versus-changed app behavior. Include meaningful negative or inconclusive results. Run the maintained tests once on the final revision, plus the targeted new cross-checks. Update the methods link/text only within separate authorization; verify the published revision and the specific corrected source/formula views if deployed.

The work is complete when the checklist below passes or records a scientifically justified nonidentified outcome under the declared rules. Reopen only for a failed criterion, changed target or new evidence—not for another routine general review.

### Maintained verification commands

Run these in the final app worktree using the project's analysis Python environment, not an arbitrary interpreter missing pytest/CmdStanPy. These commands do not launch new sampling:

```powershell
Set-Location -LiteralPath 'C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator\Custom_BED_calculator_v45'
python tools/build_model_methods.py --check
python -m pytest tests/python -m 'not sampling' -q
$env:BED_LIBRARY = 'C:\Users\ishur\AppData\Local\R\cache\R\HNSCC_BED_calculator\library-2026-10-15-R4.5.1'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e "Sys.setlocale('LC_ALL','English_United States.utf8'); source('tools/run_tests.R')"
```

Inspect each exit status before continuing. The standard R runner must write its own release marker on the actual final committed revision; do not create it manually or suppress errors. For the new frequentist runner, document separate validate-only, bounded-pilot, resume and full-run invocations as part of Phase 3/4, and record the exact invocation in each run manifest. Its commands are not yet available, and this review does not pretend to have executed them. Run the required lme4 cross-check in its dedicated analysis environment even if the ordinary app test suite legitimately skips it.

## Retain, recalculate or refit

| Existing or planned artifact | Required treatment |
|---|---|
| Current Bayesian tumor DI and calibrated DD | Retain. No new defect in their fitted numerical inputs established here. |
| Primary early U02 and its 4,000 exported joint rows | Retain byte-for-byte. Fix boundary validation around them. |
| Earlier held-out early scores/conditional intervals and structural sensitivities | Retain for unchanged numeric inputs/model. Their prior independent reproduction remains applicable. |
| 114-fit recheck, 52 cache identities and seven report hashes | Retain verified records. Do not repeat all raw-chain diagnostics merely to reject impossible JSON types downstream. |
| HYPNO counts and chemotherapy scenarios | Retain. Correct descriptive population/ascertainment and dependent hashes. |
| Methods page/JSON, trial table and manifests affected by B01–B04 | Rebuild/recheck cheaply; preserve scientific source lineage and exact-replay policy. |
| Frequentist results | None exist under this plan yet. Generate only after the revised specification and numerical gates pass. |
| Any genuinely changed selected count, exposure, likelihood, prior or scientific support | Enumerate dependent units and start a distinct analysis; do not silently mix identities. This is not required by the current metadata fixes. |
| Tumor/early app defaults, late and organ models | No numerical/default change follows from this review. Late-toxicity refitting is outside this review's requested scope. |

## Completion checklist

- [ ] Plan and audit revisions are consistent, versioned and identifiable; decisions D0–D4 are recorded.
- [ ] F01 full likelihood and profile drops pass the adaptive-integration convergence contract and independent controls.
- [ ] F02 full tumor NLL is 131.29112756 for the archived reference configuration; file/design/run identities are distinct; rounded early fixtures use justified precision.
- [ ] F03 recovery distinguishes deterministic code checks from finite-cluster bias/coverage; there is no inflated-N 2% guarantee.
- [ ] F04 matched/expanded support, prior-sensitivity interpretation, units and derived-prior handling are explicit.
- [ ] F05 all constrained profiles, boundary outcomes, null simulation if used, failed-attempt handling and confidence-interval qualifications are implemented.
- [ ] F06 restart, corruption, stale-identity and interrupted-publication tests pass; T8 is in the final gate; the full budget is measured.
- [ ] B01 source/runtime invalid-number, wrong-version and missing-identity tests refuse; valid primary exports retain their numbers.
- [ ] B02 methods checks verify actual posterior bytes and safely reject malformed records.
- [ ] B03 HYPNO displays the evaluable denominator and actual assessment schedule without an invented peak window.
- [ ] B04 tumor smax and early exact max are shown distinctly and correctly.
- [ ] Maintained Python non-sampling tests and R tests pass in their declared environments; a skipped required cross-check does not count as passed.
- [ ] Every prespecified contrast is reported, including instability/non-identification; bootstrap endpoints have Monte Carlo precision information.
- [ ] Posterior files and default app results are retained unless a separately justified scientific change explicitly supersedes them.
- [ ] Any authorized publication/deployment is verified on its actual revision, with no additional broad review required.

## Verification performed and remaining limits

**Evidence archive:** `HNSCC_BED_frequentist_Bayesian_review_evidence_2026-11-12.zip`. It contains the captured plan/audit, domain reports, scripts, numerical outputs and SHA-256 inventory. It excludes copied application source, installed libraries and the external raw-chain archive. The evidence README explains reproduction dependencies.

**Independently executed:** canonical/physical tracked-file comparison; data totals and source/identity checks; current Bayesian guard/cache/extraction positive and negative controls; current methods-builder positive check and draw/model mutation checks; historical NLL/binomial-constant verification; all ten current tumor contrast median comparisons; rounded early contrast checks; ordinary/adaptive quadrature counterexamples and positive controls; finite-cluster recovery illustration; R calendar, package and exact-replay probes; and the full maintained R test directory on an isolated source snapshot.

**R suite:** the full run recorded **6,781 expectations, zero assertion failures, one error, zero skips and zero warnings**, exit status 1. All **285 expectations** in test_v46_app_switch.R, test_v46_early_research.R and test_v461_methods_page.R passed. The sole error was test_v31_stage_e.R:17: digest's file-permission precheck reported the source CSV unreadable despite direct reads and an independently verified correct SHA-256. This is the previously established environment limitation, not an observed numerical regression; subsequent assertions in that one test were not executed by testthat. A separate probe verified the dataset digest and all four recorded gate Booleans that follow that point. A full clean-environment suite pass remains required for release certification. The wrapper used Windows UTF-8 and the pinned release library, with a byte-identical DT copy in the workspace to avoid the earlier dependency-directory access problem; a fresh comparison verified all 368 DT files. It did not create a release-pass marker. The startup's unsupported inherited C.UTF-8 warnings precede the explicit locale setup and are retained in the log; the stated zero warnings is the testthat result count.

**Python:** the native review environment lacks pytest/CmdStanPy, so no fresh full pytest pass is claimed. Direct probes use NumPy/pandas/SciPy and explicitly prevent sampling where imports need a stub. The unmodified methods builder's `--check` succeeded. The project's 4.6.1 record reports 324 Python tests passed, one skipped and two CmdStan-installation-dependent identity tests failed; that is historical project evidence, not this review's own complete test run.

**Live check:** observed 4.6.1 at deployed revision 7ec2faf, Main's 69 points and the early link's correct methods-section destination, rendered formula and current summaries. The inaccurate exact tumor hinge was visible there. No full fresh browser walkthrough of every app workflow is claimed.

**Source reading:** freshly checked HYPNO printed pp. 391/394 against the local original PDF and the primary report record. Inspected current source mappings, ascertainment/therapy records, equations, generated metadata and explanations. Previously verified unchanged CAIR/EORTC/DAHANCA/Olmi/Tandon and other counts were carried forward. This is not a claim to have reread every page of the publication library. Quadrature and R-method documentation were checked against primary sources; only the accessible publisher abstract was used for the Andrews boundary citation.

No raw-chain bytes were independently reopened, and no scientific estimator can be guaranteed defect-free by this finite review. The checks and explicit completion criteria above cover the identified issues without requiring another broad review to fill predictable implementation gaps.
