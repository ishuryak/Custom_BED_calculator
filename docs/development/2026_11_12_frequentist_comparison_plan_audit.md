# Audit: plan for likelihood-only refits of the tumor DI and acute-mucositis models

Date: 2026-11-12. Audited file: `2026_11_12_frequentist_comparison_plan.md` (first draft, same day).
Framework: `~/.claude/methods-audit/METHODS_AUDIT_INSTRUCTIONS.md` (six factors), in the actionable format (finding, fix, test, verification, order). Single-context review, no subagents.

**Scope.** I checked the following against the repository:

- the plan's model equations against `fit/stan/di.stan` and `tools/early_tox/step5_model.stan`,
- its "already exists" numbers against `docs/refit_explainer/main_free.json` and `fit/output_refit_2026_10_25c_main/{di_primary,di_wide_ab}/summary.json`,
- its test inputs against `draws/tumor_di_primary.csv`, `draws/early_fit_draws.csv` and the U02 `draws.csv.gz`,
- the likelihood convention of `docs/refit_explainer/aicc_check.py`,
- the R package availability in the Windows R 4.5.1 library.

**Not checked:** no code was written or run, so the bootstrap cost, the behavior of the optimizer and the actual intervals are unknown. Magnitudes below are labeled as derived or as judgment.

## Six factors

**Factor 1, conceptual validity: pass.** Comparing a likelihood-only fit with the Bayesian fit under the same likelihood is the right way to separate what the trials say from what the priors say. Keeping it out of the app (D0) follows from the identification facts the plan cites, and those facts check out:

- `di_wide_ab` moves the α/β upper limit from 19.1 to 31.5 Gy while α, λ, Tk and *g* are unchanged,
- the existing σ = 0 maximum-likelihood fit at α/β = 10 Gy (α 0.041, λ 0.033, Tk 38.4, *g* 0.0094) is close to the posterior medians (0.044, 0.028, 37.2, 0.0117).

One conceptual slip is in Q2 (finding 3).

**Factor 2, implementation correctness: partial.** The models match the Stan code term by term:

- the cloglog link on failures is correct, because TCP = exp(−e^η),
- centering within units is a reparameterization absorbed by c_j,
- the F-E linear form gives κ = b_dD, α/β = b_dD / b_dQ and r = −b_dTe / b_dD.

The problems are with the tests and the inference: findings 1, 2, 4, 5 and 6.

**Factor 3, parameterization and reporting: partial.** The interpretation rule uses "the log scale" for every parameter. That is wrong for Tk (normal prior) and *g* (uniform prior) (finding 7). Likelihood-ratio intervals are used near bounds without qualification (finding 5).

**Factor 4, reproducibility: pass, with minor fixes.** The plan includes per-replicate seeds, checkpoints, registry entries and golden values. The R invocation path is not specified (finding 9).

**Factor 5, fit for purpose: pass.** The plan answers how much the priors drive α/β and the contrasts, which is what decides D2. It does not produce anything a user could select, and it says so.

**Factor 6, strategic alignment: pass.** The work is about a day plus background compute. It reuses the existing builders, so it does not create a second design pipeline to maintain.

## Findings

### 1. Test T4 cannot run as written for the tumor (blocker for T4, verified)

**Evidence.** `draws/tumor_di_primary.csv` has only `alpha, ab, lambda, Tk, g, di_draw` (4000 rows, thinned from 8000). `di.stan`'s `delta_eval` is not stored per draw anywhere. `fit/output_refit_2026_10_25c_main/di_primary/` holds only `identity.json` and `summary.json`, which report the median, q10 and q90 of each contrast. The early U02 `draws.csv.gz` does have the four `lor_dec` columns, so T4 works for F-E.

**Fix.** Change T4 (tumor) to compare the Python contrast with the app's own R contrast (`R/bed_models.R`) on the same 20 parameter rows, to 1e-9. Add a distribution check: applied to all 4000 rows, the medians match `summary.json` `contrasts_vs_anchor` to Monte Carlo precision (0.1 Gy, judgment: thinning 8000 to 4000 moves a median by much less than this).

**Test.** `tests/python/test_frequentist.py::test_tumor_contrast_matches_app`. Red: drop the Q/(α/β) term.

**Verify.** `python3 -m pytest tests/python/test_frequentist.py -k contrast` passes.

### 2. T2 will fail for a reason unrelated to the code (major, verified)

**Evidence.** `aicc_check.py` computes `-(y log P + (N − y) log(1 − P))` without the binomial coefficient, while Stan's `binomial_lpmf` includes log C(N, y). The target 5342.4066 is on the first convention. A correct F-T that includes the constant will miss it by Σ log C(N_i, y_i), and someone could "fix" that mismatch the wrong way.

**Fix.** State that F-T and F-E return the log-likelihood with the constant, and that T2 adds Σ log C(N, y) (computed once) to the `aicc_check` value before comparing. Also state that T2 uses `aicc_check`'s *g* bounds (0.0046–0.021, logistic transform) and σ = 0, and that both scripts round y half-to-even.

**Test.** T2 as amended. Red: drop the constant from one side.

**Verify.** T2 green, with the constant printed in the log.

### 3. Q2 cannot be answered by comparing with the posterior (major, inferred from the design)

**Evidence.** The plan says the posterior comparison checks whether "an independent implementation reproduces the fits". A maximum-likelihood estimate and a posterior median differ by the effect of the priors, so a gap between them cannot be read as a code error. That is exactly the quantity Q1 measures.

**Fix.** Rewrite Q2 so that code correctness is answered by T1–T6 only: known values, the same design, the R cross-checks and recovery. The posterior comparison then answers only Q1.

**Test / verify.** The text change is checked by reading. No code test is needed.

### 4. ML σ is biased low with fixed unit intercepts, which narrows the contrast intervals (major, derived in direction, size unknown)

**Evidence.** F-T has 21 unit intercepts plus 6 shared parameters for 45 arms. That leaves about 18 residual degrees of freedom for the arm-level variance (derived: 45 − 21 − 6). Maximum likelihood with many nuisance intercepts underestimates a variance component (the Neyman–Scott effect), and a smaller σ narrows every profile and bootstrap interval. The Stan fit has the same fixed-prior intercepts but averages over σ (95% interval 0.004–0.17).

**Fix.**

- Report σ both as estimated and with a parametric-bootstrap bias correction.
- Add a sensitivity: re-profile the contrasts with σ fixed at the posterior median (0.065) and at the posterior 97.5% point (0.17).
- Report the contrast-interval widths for all three.

**Test.** `test_sigma_bias_recovery`: on 200 simulations from the real design at σ = 0.1, the mean ML σ is below 0.1. This documents the bias rather than gating on it. Red: simulate with σ = 0.

**Verify.** The write-up table has three columns of contrast intervals.

### 5. Likelihood-ratio and bootstrap intervals are not valid at bounds or along ridges (major, standard theory)

**Evidence.** σ, *g* (≥ 0), κ (≥ 0), λ (≥ 0) and Tk can sit at a bound. If λ is at 0, Tk is not identified at all, because the hinge term vanishes. When Tk is below every arm's time, smax(T − Tk) is linear in T within each unit and λ and *g* lie on a ridge. In these cases the χ²₁ cut-off of 1.92 and percentile bootstrap intervals have no guaranteed coverage. The bootstrap is inconsistent when the true value is on the boundary (a standard result, Andrews 2000 in Econometrica, cited from memory and not re-checked for this internal audit).

**Fix.**

- Label such intervals "approximate (at a bound)" in the table.
- Use the profile shape, not the interval, for the data-dominated/prior-dominated label.
- Report λ + *g* (identified on the ridge) beside each term.
- Test λ = 0 by parametric bootstrap of the likelihood ratio, not by χ².

**Test.** `test_ridge_detected`: simulate with Tk = 5 days and confirm that the fit flags the λ–*g* ridge, using the Hessian condition number above a registered threshold. Red: Tk = 35.

**Verify.** Every bounded row in the output table has its flag.

### 6. Bootstrap refits that fail must be counted, not dropped (minor, judgment)

**Evidence.** With 27 parameters, boundaries and a possible ridge, some bootstrap refits will fail to converge. Dropping them silently biases the percentile interval toward well-behaved samples.

**Fix.** Record convergence per replicate in `bootstrap.csv`. Report the failure rate and the intervals with and without failed replicates. Stop and report if more than 5% fail.

**Test.** `test_bootstrap_records_failures`: force a failure with `maxiter = 1` and confirm that the row is written with `converged = 0`. Red: filter before writing.

**Verify.** The failure rate appears in the write-up.

### 7. The "log scale" in the interpretation rule fits only lognormal-prior parameters (minor, verified)

**Evidence.** Priors: Tk normal(28, 10), *g* uniform(0.0046–0.021), and the early τ and σ_arm half-normal. A log-scale width is meaningless for a uniform prior that starts near 0, and for Tk.

**Fix.** Compare widths on the prior's natural scale: log for α, α/β, λ, r and κ, and identity for Tk, *g*, γ, τ and σ. Write the rule as width(likelihood 95%) / width(prior 95%) on that scale. Keep the labels and thresholds (below 0.5 data-dominated, above 1 or reaching a search limit prior-dominated).

**Test.** `test_label_rule`: a hand-computed example for each scale. Red: apply log to *g*.

**Verify.** The label table states the scale used.

### 8. The F-E Laplace cross-check needs a stated purpose (minor, judgment)

**Evidence.** The plan already gates on the exact `glm()` check and calls the `glmer` Laplace comparison descriptive. That is correct. The risk is that a reader treats a Laplace difference as a failure.

**Fix.** Add one sentence: Laplace is known to bias binomial variance components with small counts, so differences in τ and σ_arm are expected, and only κ, α/β and r are compared numerically.

**Test / verify.** Reading only.

### 9. The R invocation path is not specified (minor, verified)

**Evidence.** `Rscript` is not on the WSL path. Windows `Rscript.exe` exists at `/mnt/c/Program Files/R/R-4.5.1/bin/`, and `lme4` is in `C:/Users/ishur/AppData/Local/R/win-library/4.5`, not in the `renv` release library. Windows R cannot resolve `/mnt/c/...` arguments.

**Fix.** Run `glmer_check.R` from its directory with bare filenames, or pass `wslpath -w` paths. Record the `lme4` version in the output. Add `tests/testthat/test_frequentist_glmer.R` to the suite only behind a `skip_if_not_installed("lme4")`, so the release suite is unchanged.

**Test.** The skip fires when `lme4` is absent. Red: remove the skip in a clean library.

**Verify.** `Rscript.exe tools/run_tests.R` gives the same count of failures (0) as before.

## Investigated and dismissed

- **cloglog orientation.** Correct: TCP = exp(−e^η) makes P(failure) = 1 − exp(−e^η).
- **Centered versus uncentered predictors** (di.stan versus aicc_check). This is the same likelihood, because the centering constants are absorbed by c_j.
- **"21 units" versus "19 trials" elsewhere in the docs.** Fitting units include strata, and `main_free.json` gives n_units = 21.
- **Patient-level data.** Every input is published trial-arm aggregates, so no git-safety issue arises.
- **DD omitted.** Justified: the calibrated DD has no likelihood-only analogue, and the free DD is already covered by `aicc_check.py`.

## Implementation order

1. Claude: apply fixes 1–9 to the plan text (done in the same pass, see the revision note in the plan).
2. Igor: rule on D0–D4. This blocks everything below.
3. Claude: tests for findings 1, 2 and 4–7 with T1–T6, each shown red.
4. Then smoke, tests again, and the full run, as in plan Section 7.

## Bottom line

The plan is sound in concept and scoped correctly. Two of its tests would have failed or misled as written (findings 1 and 2), its second question was mis-assigned (finding 3), and its intervals needed bias and boundary qualifications (findings 4 and 5). After those fixes it is **accurate after minor fixes**.

**Stopping condition:** once fixes 1–9 are in the plan and Igor has ruled on D0–D4, implementation proceeds. Reopen only if a test cannot be made red, or the deployed fits change.
