# Mathematical and statistical review: frequentist comparison plan, 12 November 2026

## Scope and conclusion

Reviewed the **physical revision 2** of `docs/development/2026_11_12_frequentist_comparison_plan.md` and its accompanying audit in `C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45`. These files are uncommitted additions relative to reported HEAD `79e9c812f91482fa57ffa1dddabe89dc329a6303`; the plan, not merely HEAD, is the review target. Its SHA-256 is `44d2324f5e7820d9b9155cef43c17f8d5d32bc09e5d25623f3c369c95a0cf292`; the audit hash is `c1b445873cdb7c7026d727b960c247216326c5e1321dbce567f8002defd95deb`. Copies of both reviewed files are preserved beside this report.

**The offline comparison is worthwhile, but revision 2 needs material changes before implementation and the full bootstrap.** Its principal model equations and random-effect integration structure are correct. The audit's verdict “accurate after minor fixes” is too favorable: it introduced or retained a likelihood-constant sign error, an invalid high-N recovery criterion, an inadequately specified quadrature calculation, mismatched cross-check models, and unjustified attribution of all MLE/posterior differences to priors. Boundary and failure handling remain insufficiently specified for execution.

These are mostly defects in a proposed analysis, not evidence against the completed Bayesian fit. The current early evaluator, Stan model, fraction table, outcome table and 4,000 exported posterior rows are unchanged from the independently verified 11 November generation. Current `draws/early_fit_draws.csv` still has SHA-256 `a6fde4626c88a1872a84408d66e78ddd1460cdabb98734d06427ce186f07b4fd`. Retain its fitted parameters, held-out results and source draws. No new model fits or production edits were performed for this review.

Priorities below mean **P1: correct before executing the planned analysis** and **P2: complete before relying on or reporting its intervals/labels**. They are not statements that the deployed calculator has a corresponding current defect.

## What is mathematically sound

The tumor model in plan lines 35–41 agrees with `fit/stan/di.stan`: with control probability `TCP=exp(-exp(eta))`, the complementary log-log link applies to failure probability. Mean BED and time centering can be absorbed by one intercept per fitting unit in the unpenalized likelihood. The Gaussian arm-effect density must remain part of the random-effects likelihood; “without priors” should mean removing priors on structural/fixed parameters, not removing the standard-normal distribution defining the random effects.

For tumor F-T the marginal likelihood factorizes into independent one-dimensional arm integrals. For early F-E it factorizes by family as

\[
 L_j=\int\phi(u)\prod_{i\in j}\left[\int\phi(z_i)\,
 \operatorname{Binom}\{y_i;N_i,\operatorname{logit}^{-1}(x_i\beta+\tau u+\sigma z_i)\}\,dz_i\right]du.
\]

Integrating the product over the shared family effect, rather than multiplying separate family-effect integrals per arm, preserves within-family dependence. The plan describes this ordering correctly. The exact-variance-zero GLM limit is an appropriate oracle once the parameter constraints and predictor columns are aligned.

Keeping the work offline and keeping posterior-draw support distinct from frequentist resampling frequencies are good decisions. The existing wide-α/β Bayesian sensitivity and early U02 prior/repair sensitivities remain useful evidence and should accompany, rather than be replaced by, a new likelihood analysis.

## M1 — P1: fixed-node quadrature is not “essentially exact”; specify adaptation and convergence

**Locations:** plan lines 42, 56 and 113; audit's conceptual/implementation pass does not address this issue.

The plan specifies 30-node Gauss–Hermite for F-T and nested 20×20 for F-E without saying whether nodes adapt to the conditional integrand. Ordinary fixed nodes around the Gaussian prior can miss a sharply concentrated likelihood. Testing isolated toy arms at unspecified N and predictors is not enough. Increasing N by 1,000 in T5 makes this numerical problem worse.

**Independent evidence, with actual arm denominators/counts and permitted variance values.** The probe uses stable log probabilities, mode-scaled adaptive `scipy.integrate.quad` as an independent integral, and ordinary Hermite nodes. It does not fit a model.

| Fixture | Correct log marginal likelihood | Fixed-node result minus correct result |
|---|---:|---:|
| Tumor arm N=750, controls=525, eta=log(-log(0.7)), sigma=0.25, 30 nodes | −4.799961516 | −0.033952656 |
| Same arm, sigma=1, 30 nodes | −6.153702024 | −6.871555069 |
| Early DAHANCA7_6F, N=398, events=269, eta=logit(269/398), sigma=0.25, 20 nodes | −4.085205090 | −0.003356541 |
| Same early arm, sigma=1, 20 nodes | −5.392732536 | −3.621221976 |
| RTOG0129 family, N=(361,360), y=(141,119), fixed predictors=observed logits, tau=1, sigma_arm=0, outer 20 nodes | −8.812489172 | −7.974022741 |

These are valid likelihood-evaluation fixtures, not claims about the eventual fitted parameter values. The family fixture is a required boundary submodel of the nested integral. Errors of several log-likelihood units exceed the proposed 1.92 profile cutoff and can change maxima, apparent ridges and confidence sets.

**Positive control:** shifting/scaling the nodes to each integrand's mode and curvature makes adaptive GH30 for the six single-arm fixtures and adaptive GH20 for the three family-boundary fixtures agree with the independent integral within **2×10⁻¹⁴**. The problem is fixable without changing the model. This behavior accords with the primary analysis of fixed versus adaptive numerical integration by [Bilodeau, Stringer and Tang](https://arxiv.org/html/2202.07864).

**Fix:** explicitly implement adaptive quadrature or an independently validated adaptive integrator at both nesting levels, with log-sum-exp accumulation, exact zero-variance branches, and error control. Check integrals, objective changes and gradients, not just scalar likelihood values at the starting point. Record the accepted node/integration policy in every result identity.

**Acceptance:** all fixtures above pass the specified log-likelihood tolerance (1e-9 is attainable here). Add nonzero-inner/nonzero-outer family fixtures, extreme/all/zero events and large N. At the fitted point, all profile endpoints and selected bootstrap fits, increasing quadrature order/tightening integration changes total log likelihood by <1e-6, changes reported contrasts by <1e-4 in their natural units, and leaves profile crossing classifications unchanged. Treat these numerical thresholds as proposed registered criteria, not tests already passed on a future implementation. A deliberately nonadaptive 20/30-node implementation must fail the fixtures above. Re-estimate the computational budget only after this algorithm is fixed.

## M2 — P1: T5's high-N 2% recovery rule is invalid for these random-effects models

**Locations:** plan line 115; implementation gate lines 120, 125–128.

Multiplying the patients in the same 45 arms/21 units or 51 arms/22 families does not create additional independent arm or family effects. As N grows, the observed arm probabilities become more precise, but the particular realized Gaussian effects do not average away. Variance-component estimation and structural coefficients remain subject to that finite-cluster variation. A correct program will frequently fail “every identified parameter recovered within 2%.”

**Analytical numerical control:** even if the 22 early-model family effects were observed perfectly and independently as Normal(mu,tau²), their MLE standard deviation satisfies `22*(tau_hat/tau)^2 ~ chi-square(21)`. Its central 95% ratio range is **0.683670–1.269912**, and the chance of recovering tau within 2% is only **10.12%**. Exact observation of those effects is already more informative than binomial counts with huge denominators. This is a counterexample to T5's criterion, not a calculation of the eventual model's actual variance uncertainty.

**Fix:** separate deterministic correctness from statistical recovery. Use zero-random-effect/high-N fixtures or exact expected-score checks for coefficient signs and algebra. Test positive random-effect integration using independent numerical oracles. For statistical recovery, simulate many independent replicated trial/arm designs (or evaluate bias and coverage across repeated finite-design datasets) and use predeclared Monte Carlo error bands, not a 2% pass rule for one realization. Do not gate code correctness on an asserted direction of sigma bias in 200 simulations; report that as an empirical result, with uncertainty.

**Acceptance:** a correct finite-cluster simulation may produce a non-2% estimate without being called an implementation failure. Deliberately reversing a time coefficient fails deterministic or expected-score checks. A replicated-cluster recovery experiment reports its number of independent units, bias, Monte Carlo standard errors, boundary frequency and coverage. Tests distinguish optimizer failures from a valid boundary maximum or a nonidentifiable parameter.

## M3 — P1: the proposed “exact” R cross-checks can fit different parameter spaces

**Locations:** plan lines 48, 59–60 and 116.

For F-T at a **fixed** b=α/β, the mean contains `-alpha*(dD+dQ/b)`. Fitting free coefficients for `dD` and `dQ` does not hold b fixed. It fits a larger model, whose implied ratio is `b_dD/b_dQ`. Furthermore, unconstrained glmer can return signs forbidden by the Python model for alpha, lambda or g. At an interior joint global optimum an unconstrained linear fit at Tk may happen to agree, but that does not make it an exact fixed-(b,Tk) oracle or cover the boundary cases this plan expects.

For F-E with tau=sigma_arm=0, plain glm is equivalent only if its unrestricted solution satisfies `b_D>=0`, `b_Q>0` and `b_T<0` with meaningful nonzero denominators. The nonlinear constrained parameterization does not match arbitrary unrestricted GLM coefficients. The Laplace sensitivity can also differ in fixed-effect estimates, not only variance components; excluding tau and sigma_arm from numerical comparison does not isolate the approximation error.

**Fix:** at fixed (b,Tk), use a single `dBED=dD+dQ/b` column for F-T. Compare the same constrained region or use an explicitly interior synthetic fixture and verify signs before treating an unrestricted R estimate as an oracle. Alternatively describe the free-dD/free-dQ check correctly as a fixed-Tk free-b check, with matching admissibility conditions. For F-E, require an interior matching GLM fixture or compare objective values at identical admissible coefficient vectors instead of insisting on equality of potentially different optima. Evaluate both programs' likelihood at the same parameters to isolate numerical integration from model restrictions.

The [official glmer documentation](https://lme4.github.io/lme4/reference/glmer.html) supports adaptive GH for a single scalar random-effects term. It also warns that internal deviance conventions can differ by approximation; do not confuse a constant convention with a coefficient discrepancy. The dedicated offline gate must actually run its required R check. A release-suite `skip_if_not_installed` is acceptable for packaging but cannot count as a passed scientific validation gate.

**Acceptance:** the one-column fixed-b formulation matches the Python likelihood and interior maximum; an intentionally free ratio or wrong link fails the designated oracle; boundary solutions are compared under matching constraints, without requiring invalid parameter ratios. Use mixed absolute/relative tolerances so a valid zero coefficient does not make “1e-4 relative” undefined.

## M4 — P1: T2 has the binomial-constant sign reversed; T4 asks rounded draws for unavailable precision

**Locations:** plan lines 45 and 112–114; audit findings 1 and 2.

`aicc_check.py` returns a positive **negative** log likelihood without the combinatorial constant. If C=sum(log choose(N,y)), then

`full negative log likelihood = old negative log likelihood − C`,

`full log likelihood = −old negative log likelihood + C`.

The audit's instruction to add C to positive 5342.4066 is wrong. The parent review independently calculated C=**5211.11548984419**, so the correct full NLL is **131.291127555958**, or LL −131.291127555958. Addition instead gives 10553.522107244338. Fix the test's target and explicitly name LL versus NLL before any implementation is made to pass it.

For early T4, U02's retained selected-draw CSV stores six significant digits (`step5_fit.py` selected-draw output). Recomputing its four `lor_dec` columns from its rounded kappa, b and r cannot meet a universal 1e-6 absolute or relative tolerance. Direct current-data reproduction gives:

| Decision schedule | Maximum absolute disagreement | Rows exceeding 1e-6 absolute |
|---|---:|---:|
| 15×3.5 Gy | 1.50794e-5 | 3,397/4,000 |
| 20×2.75 Gy | 1.32642e-5 | 3,091/4,000 |
| 25×2.4 Gy | 6.39965e-6 | 2,648/4,000 |
| 30×2.2 Gy | 3.43016e-6 | 1,505/4,000 |

These differences are rounding, not a new model error. Use full-precision retained raw-chain fixture values for a tight formula test, or derive a forward-error allowance from six-significant-digit rounding of every input and output. An arbitrary parameter change to satisfy 1e-6 would be a regression. A ~2e-5 absolute allowance covers these current rows, but a propagated rounding bound is more transferable. Keep a separate independently generated full-precision oracle at a tight tolerance.

For tumor T4, the parent review found all current differences from the original 8,000-row contrast medians below 0.1 Gy (maximum **0.092252 Gy**), so do not report a present failure. However, the 4,000 app rows are **1,000 distinct DI draws repeated four times**, not 4,000 independent/thinned DI draws. A fixed 0.1 Gy threshold is a property of this saved selection, not a general precision guarantee. Compare deterministic formulas on the same rows; preserve the full-draw reference and selected-index provenance separately. Tumor `summary.json` stores q10/q90 for contrasts, not a 95% interval: compute the requested 2.5%/97.5% values from an identified posterior set and label which set was used.

**Acceptance:** corrected T2 passes with explicit sign convention; dropping the constant fails. Full-precision early formulas pass tightly; stored-rounded fields pass only within their registered rounding bounds; dropping Q/b fails both. Tumor selected-draw comparison reports distinct-draw count and does not confuse 80% with 95% intervals. No model refit is needed for any of these corrections.

## M5 — P1: the scientific question and parameter domains do not isolate “the prior effect”

**Locations:** plan lines 18–21, 25–27, 44, 88–94 and 135; audit finding 3 still retains the same attribution error.

An MLE and a posterior median do not differ by “exactly the prior effect.” They use different summaries and different treatments of nuisance parameters. Even for a likelihood `L(lambda)=lambda*exp(-lambda)` with a constant prior density on lambda>0, the MLE and posterior mode are 1 while the posterior median is **1.678347**. No nonconstant prior weighting produced that median/mode difference. Marginalizing nuisance parameters versus profiling them adds another difference.

The proposed primary tumor domain also changes support relative to Stan: g≥0 instead of [0.0046,0.021], and Tk∈[12,54] instead of [0,60]. Restricting an onset to observed course durations is not a physical law: onset may precede the shortest or follow the longest course. This primary is not uniformly less restrictive than the Bayesian support. Calling it “physical constraints only” and declaring that “only priors differ” is inaccurate. The existing primary Bayesian DI comparison in line 25 has b estimated (median 10.6195), not fixed at 10; its approximate numerical proximity to the fixed-b no-overdispersion MLE is encouraging but is not a matched-model equality demonstration.

**Fix:** make the **matched-support marginal-likelihood fit** the primary comparison for Q1, retaining Stan's bounds and the same structural/random-effect distribution. Label a wider scientific support fit a separate sensitivity and predeclare its support; if a finite maximum does not exist, report that fact. Rewrite Q1 as “likelihood information and sensitivity to prior/constraint choices,” not exact causal attribution to priors. Use the existing genuine prior-width sensitivities as the direct evidence of prior sensitivity. Do not promise the new fits' results in advance based on the older 35-arm expected-information calculation; the present early design has 51 arms/22 families and additional heterogeneity.

Width ratios can be reported as descriptive measures on declared scales, but “data-dominated/prior-dominated” should not be treated as an established statistical diagnosis. Reaching an arbitrary numerical search edge means **unresolved numerical extent** until the grid is expanded or a limiting model is evaluated. A finite profile narrower than a prior need not establish negligible prior influence. Kappa's prior is half-normal, not lognormal, and kappa=0 has no finite log value. Derived contrasts require their induced joint prior distribution; the plan does not specify how to derive that distribution or classify negative/zero contrasts.

**Acceptance:** primary support matches the Bayesian fit; sensitivities identify every changed bound. Demonstrate support expansion or a limiting model before declaring unboundedness. Record profile extent, uncertainty and prior-width ratio descriptively; reserve “prior-sensitive” for an observed prespecified prior perturbation. Prior ranges honor truncation/support and contrast ranges are calculated from joint priors. A flat-prior skew example must not be labeled as a measured nonconstant-prior effect.

## M6 — P2: boundary, profile and bootstrap procedures need executable definitions

**Locations:** plan lines 46, 69–94, 101, 117 and 142; audit findings 4–7.

The added flags identify real concerns but do not specify a reliable computation or cure coverage problems:

- A single 1.92-drop contour has its usual nominal interpretation only under regular scalar likelihood-ratio conditions. A boundary, unidentified Tk under lambda=0, ratio singularity, or multiple profile components cannot be fixed by adding the word “approximate.” [Andrews' primary boundary-bootstrap result](https://doi.org/10.1111/1468-0262.00114) supports caution about unrestricted bootstrap inference; it does not mean all null-imposed parametric bootstrap tests are invalid.
- L-BFGS-B handles boxes, not arbitrary equality constraints such as a fixed BED contrast. The plan needs an actual profiling algorithm for every derived quantity: a validated reparameterization/elimination where denominators are safe, or equality-constrained optimization with checked residuals.
- b→infinity, alpha→0 or kappa→0 can create unbounded or undefined ratios and contrasts. A positive-definite Hessian and a finite estimate are not guaranteed outcomes. Approximate lambda/g ridges are expected under an onset well before all courses. The proposed Tk=5 ridge fixture is outside the proposed [12,54] primary domain; put it explicitly in the broader-support test.
- “Intervals with failed replicates” is undefined if a replicate has no valid estimate. Including arbitrary failed optimizer outputs is not a sensitivity interval. With 5% unresolved failures, even the entire 2.5% tail may be missing; a 5% threshold alone does not protect percentile endpoints.
- The specified lambda=0 LR bootstrap must simulate from the **fitted null**, refit null and alternative for each replicate, search Tk under the alternative, and use an explicit Monte Carlo p-value/critical-value convention. The ordinary bootstrap under the unconstrained fitted alternative does not automatically provide these samples. This is additional computation not included by merely naming B=1,000/500.
- The finite-design sigma bias direction/magnitude and its effect on every contrast are not universal consequences of subtracting parameter counts. A plug-in bootstrap bias correction can remain zero at a boundary and does not supply calibrated coverage. Fixing sigma at posterior values is a Bayesian-informed conditional sensitivity, not a likelihood-only correction.

**Fix:** define estimator status codes (interior maximum, admissible boundary maximum, nonunique ridge, supremum at infinity, numerical failure), profile constraints, adaptive bracketing and endpoint status. Use projected/KKT gradients at box boundaries. Require every profile objective to be no better than the claimed global maximum within tolerance; an improvement triggers global optimization repair. Return disconnected confidence sets as such. Parameterize inverse-b when appropriate to evaluate the b=infinity limit explicitly, and define denominators/undefined estimands without forcing an arbitrary finite answer.

For bootstrap generation, hold the design and N fixed; regenerate the stipulated arm and family Gaussian effects and then binomial counts, conditional on fixed unit intercept estimates for F-T. Do not resample fitted latent effects as if they were the population. Archive the generated data/seed. Retry a numerically failed dataset with a stronger deterministic strategy on the **same dataset**, not a replacement draw. Boundary solutions that satisfy the optimizer contract are valid outcomes, not failures. Suppress primary percentile intervals if unresolved failures remain; describe success-only intervals as conditional descriptive summaries, or supply worst-case tail-rank bounds showing what the failures could do. Keep the >5% stop as a resource-warning rule, not an inferential guarantee.

B=500 gives an expected 12.5 observations below a true 2.5th percentile, with a central binomial count range approximately **6–20**; B=1,000 gives 25 with range **16–35**. These are reasonable initial budgets for exploratory intervals, but report endpoint Monte Carlo uncertainty and extend only when it prevents the declared conclusion. Do not require an arbitrary larger B without a precision need.

**Acceptance:** explicit positive/zero-variance and lambda=0 fixtures return the intended status; deliberate wrong sign/failed optimizer is caught; every profile constraint residual is within a registered tolerance; likelihood increases from a profile trigger a better global fit; no failed numerical result enters an ordinary percentile as if valid; a forced restart reproduces the same datasets/results; null LR simulations are distinguishable by identity from alternative-estimator bootstrap simulations. Include T8 in the completion gate—line 142 currently stops at T1–T7.

## M7 — P2: identity, restart, numerical budget and reproducibility must be specified before long runs

**Locations:** plan lines 101, 105, 111 and 124–142.

“Skip if present” plus a per-replicate seed is insufficient. A changed likelihood, quadrature policy, constraints, retry strategy, N or design can otherwise silently mix generations. T1 also equates different hashes: an input CSV byte hash, a canonical design-array hash and a whole fit identity are different objects and should not be asserted equal. Early U02's identity includes model/settings/sampler/environment, not merely its design arrays.

**Fix:** freeze a manifest covering source revision plus actual uncommitted file hashes, data bytes, canonical design/order and rounding, likelihood version, parameter support, quadrature and optimizer settings, starts, profiles, seed mapping, bootstrap DGP and software. Compare like with like in T1: source bytes against source-data hash; canonical arrays against canonical arrays/Stan-data hash; the full recorded identity via its own serializer. Validate all expected outputs and columns before resuming. Write one complete replicate record atomically with a checksum or use a transactional storage format, then derive CSV reports. Refuse altered identity, duplicate/missing/truncated records and nonfinite entries masquerading as completed work. Preserve partial failures rather than silently skipping them.

The 1–3 CPU-hour estimate is not established. Adaptive nested integration, 40 starts, profile scans, retries and a separate null LR bootstrap may dominate cost. Measure these individually in the smoke run; state which starts are repeated for every bootstrap, how warm starts are used, and whether the best solution is independently checked. A B=20 smoke is an operational check, not calibration of 95% tails or boundary coverage. Golden tests should pin likelihood and identifiable contrasts/status within tolerance, not insist that a nonunique ridge select identical coordinates across library versions.

**Acceptance:** uninterrupted and interrupted/resumed smoke runs have identical per-replicate seeds, generated datasets, identities and results. Mutating a bound, node policy, design arm or code hash refuses reuse. A crash during append does not create a “completed” row. The dedicated scientific R gate actually runs and passes on its pinned library. The work estimate explicitly includes profile and null-test costs before a long run begins.

## M8 — P2 documentation: leaving calibrated DD out is a scope decision, not a mathematical impossibility

**Locations:** plan lines 11 and 64; audit “DD omitted” dismissal.

A frequentist analogue of the calibrated procedure need not be a free-DD fit. One could estimate DI, hold its structural values, solve the same anchor constraint, estimate conditional DD parameters, and repeat both stages on bootstrap datasets to propagate estimation uncertainty. This would be a defined two-stage/constrained procedure whose inferential properties need separate work; it would not automatically equal a joint unrestricted maximum-likelihood model or posterior draw mixture.

Leaving it out of this project phase is nevertheless sensible: it expands the estimand, uncertainty and validation scope and is unnecessary for the specified DI/early likelihood-information check. Replace “has no likelihood-only analogue” with “a calibrated frequentist two-stage analogue is outside this comparison and would require a separate uncertainty/validation design.” No DD fits are recommended in the present plan.

## Recommended decisions and execution order

1. **D0: offline comparison, no app switch.** Retain this scope. State what the inference does and does not represent; no calculator draw replacement or promotion follows automatically.
2. **D1: matched support primary; broader support sensitivity.** Remove the [12,54] “physical” assertion. Define the broader domain or limiting cases before results are seen.
3. Revise the question/labels, correct T2 and T4, specify constraint-aware R checks, and replace T5. Preserve the existing data/model code and fitted outputs. Include all required tests, including T8, in one final gate.
4. Implement/test the likelihood **first** against exact variance-zero cases, independent scaled integration and adaptive quadrature fixtures. Positive/negative controls above must run before any optimization result is treated as scientific.
5. Implement global fitting, profiles, admissible-boundary statuses and manifests. Test known interior and boundary/ridge fixtures. Validate array equality, source provenance and actual scientific R checks. A release test skipped for lack of lme4 is not a pass for this stage.
6. Run the bounded operational smoke with manifests and restart testing. Measure integration, 40-start fits, profiles, bootstrap retries and any null LR test separately. Decide the executable budget from this measurement.
7. **D3: retain 1,000/500 as initial exploratory estimator-bootstrap budgets**, with endpoint Monte Carlo diagnostics and a separate explicitly scoped null-LR budget if that test is retained. Run only after the preceding numerical gates pass. Resolve numerical failures on the same simulated datasets. Report boundaries/unbounded sets rather than manufacturing finite intervals.
8. Generate the full offline report whether results suggest strong or weak likelihood information. **D2:** a neutral methods summary should not be conditional on a desired “prior-dominated” label; profile figures can be chosen for scientific relevance. **D4:** omit DD and late-toxicity fitting for the stated scope reasons.
9. Completion checklist: model/domain contract frozen; known LL/NLL conventions correct; actual-design hashes verified; independent integration and gradient checks pass; matched R check run; profiles and boundary statuses verified; bootstrap outputs complete or explicitly noninferential; restart corruption controls pass; nominal/Monte Carlo intervals distinguished; all reported Bayesian and frequentist quantities use identified source sets and interval levels; no app/deployment change under this plan. No new general review is required after these concrete criteria pass unless new evidence changes the task.

## Retain, recalculate and fit only when justified

- **Retain:** current tumor posterior/draw pairing; completed Bayesian early U02, held-out scores, HYPNO scenarios, corrected observed-scale intervals, repair sensitivities and SBC summaries. The probes find no reason to refit them. Current primary early numerical inputs are unchanged from the prior verified review.
- **Cheap recalculation/document changes:** plan/audit text, T2 constant conversion, T4 rounding-aware fixtures, posterior/prior 95% contrast intervals from specified draws, source/design identities, descriptive likelihood/prior-width labels, numerical and restart validation. A prior-range simulation does not require fitting a new outcome model.
- **New computation justified by this proposal:** the new frequentist fits, profiles and registered bootstrap after their implementation is validated. Additional Bayesian fits are justified only by an actual changed input/model/prior question or a named unresolved diagnostic issue, not by ordinary MLE/posterior disagreement.
- **Do not infer:** clinical superiority, an alternative clinical threshold, a posterior probability from bootstrap frequency, or complete identification from a full-rank design. The clinical/source review separately addresses tumor pseudo-count outcome uncertainty and early endpoint applicability.

## Evidence and reproduction

`plan_probes.py` is a portable, read-only-source script using Python 3.12, NumPy, SciPy and pandas. It writes files only beside itself; CmdStan compilation/sampling imports are replaced by an explicit raising stub. It evaluates numerical integrals and closed-form statistical controls, **not model fits**.

```powershell
& 'C:/Users/ishur/AppData/Local/Programs/Python/Python312/python.exe' -B `
  'C:/Users/ishur/Documents/Codex/2026-10-18/rev-2/work/review_2026_11_12/math/plan_probes.py' `
  'C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45'
```

Outputs: `plan_probes.json`/`.log`, exact reviewed plan/audit copies, nine fixed-versus-adaptive integration controls, early T4 comparisons, the finite-cluster recovery counterexample, flat-prior summary counterexample, bootstrap percentile order uncertainty, and current/previous verified-artifact comparisons. The previous-artifact comparison expects the prior review snapshot in this workspace; the numerical plan checks themselves use the project supplied as the argument.

Parent evidence `work/review_2026_11_12/plan_numeric_checks.json` supplies the independently reproduced tumor constant and selected-draw median checks quoted here. Primary sources consulted were the official glmer documentation, Bilodeau–Stringer–Tang's full accessible mathematical text, and Andrews' publisher record/abstract. Andrews' complete PDF could not be retrieved through the browsing tool, so no full-paper reading claim is made. No bootstrap, MCMC, whole-model maximum-likelihood run, production write or deployment was performed in this review.
