# Plan: likelihood-only (frequentist) refits of the tumor DI and acute-mucositis models

Date: 2026-11-12. Status: DRAFT for approval, revision 2 (audit findings 1–9 applied). No code has been written or run for this plan.
Author: Igor Shuryak (plan drafted with Claude). Audit: `2026_11_12_frequentist_comparison_plan_audit.md`.

## 0. Decision already taken, to be confirmed (D0)

Igor asked whether to add frequentist mixed-effects refits "so the user can choose either". The recommendation, which this plan implements, is to run them as an **offline check recorded in `docs/development/`, not as a user switch in the app**. Reasons:

- α/β (tumor and mucosa) and background growth *g* are weakly informed by the trials, so a likelihood-only fit either fixes them (then it is a scenario, not a competing fit) or reports very wide intervals.
- The calibrated DD model has no likelihood-only analogue (Section 3.3), so a switch could cover DI only.
- "Refit support" and the Uncertainty-page shares are counts over posterior draws. A frequentist version needs bootstrap or simulated draws, which mean something different under the same label.

**D0 (Igor): confirm "check, not switch".** If Igor wants the switch anyway, this plan is still step 1 of that work, and a separate plan covers the app.

## 1. Questions

- **Q1. How much do the priors drive the results?** For each parameter and for each calculator contrast, compare the likelihood-only interval with the prior range and the posterior interval.
- **Q2. Is the independent implementation correct?** This is answered only by tests T1–T6 (same design, known values, R cross-checks, recovery). It is not answered by agreement with the posterior, because a maximum-likelihood estimate and a posterior median differ by exactly the prior effect that Q1 measures (audit finding 3).

Not a question here: which model is "better". Both fits use the same likelihood, and only the priors differ.

## 2. What already exists (not repeated)

- `fit/mle_check.py` and `docs/refit_explainer/aicc_check.py`: maximum likelihood of DI and DD with one intercept per fitting unit, **no overdispersion term**, α/β fixed at 10 Gy, *g* fixed (0.0154 per day) or free within 0.0046–0.021. On the current main set with *g* free (`docs/refit_explainer/main_free.json`, 45 arms, 21 units): DI α 0.041 per Gy, λ 0.033 per day, Tk 38.4 days, *g* 0.0094 per day, Pearson χ²/df 1.66. The DI posterior medians (`fit/output_refit_2026_10_25c_main/di_primary/summary.json`) are α 0.044, λ 0.028, Tk 37.2, *g* 0.0117. **Point estimates at α/β = 10 Gy therefore already agree.** The new information this plan adds is α/β itself, the arm-level overdispersion, and intervals.
- `fit/output_refit_2026_10_25c_main/di_wide_ab`: the Bayesian DI fit with a wider α/β prior. The α/β posterior moves from 10.6 (6.6–19.1) to 11.7 (6.1–31.5) Gy while α, λ, Tk and *g* barely move. This already shows that the upper tail of α/β is set by the prior.
- `tools/early_tox/step4_identification.py`: design-only identification for the early model (Fisher information and simulated recovery on the 35-arm primary set). Expected posterior SD of log α/β was 0.54 to 0.87 times the prior SD, and of log *r* 0.20 to 0.40 times. So the mucosal α/β is expected to be only partly identified, and *r* well identified.

## 3. Models

All fits use the **same data files, the same arm-level terms and the same likelihood** as the deployed Bayesian fits. The design is imported from the existing builders, not rewritten: `fit/fit_refit.py` (tumor `stan_data`) and `tools/early_tox/step5_fit.py` (`load_rows`, `stan_data`). A test asserts that the design is identical (Section 6, T1).

### 3.1 Tumor DI (model F-T)

Arm *i* of fitting unit *j*, with failures *f* = N − y:

cloglog P(failure) = η_i = c_j − α(BED_i − BED̄_j) + λ·smax(T_i − Tk) + g(T_i − T̄_j) + σ z_i,  z_i ~ N(0, 1),  BED_i = D_i + Q_i / (α/β).

This is `fit/stan/di.stan` without priors. TCP = exp(−e^η), so the complementary log-log link applies to **failures**, not controls.

- **Fixed effects:** c_j (one per unit, as in Stan, where c_j has a fixed N(0, 3) prior and is not pooled). **Random effect:** independent arm-level z_i.
- **Marginal likelihood:** each arm's integral over z_i is one-dimensional and independent of the others, so it is computed by Gauss–Hermite quadrature (30 nodes). This is essentially exact and avoids the Laplace approximation.
- **Parameters estimated jointly:** c, α, α/β, λ, Tk, g, σ.
- **Constraints (decision D1):** primary uses the physical constraints only (α > 0, α/β > 0, λ ≥ 0, σ ≥ 0, g ≥ 0, Tk within the observed range of overall times, 12–54 days). The secondary uses the Stan bounds (g in 0.0046–0.021, Tk 0–60).
- **Likelihood convention:** both F-T and F-E return the full log-likelihood, including log C(N, y), as Stan does. `aicc_check.py` omits that constant, so T2 adds it before comparing (audit finding 2).
- **Optimizer:** L-BFGS-B from 40 starting points (the multistart scheme of `aicc_check.py`), then convergence checks: gradient norm, a positive-definite Hessian away from bounds, and best-of-starts agreement.

Linear form for the R cross-check: given (α/β, Tk), η is linear in (D − D̄_j), (Q − Q̄_j), smax(T − Tk) and (T − T̄_j). So the call `lme4::glmer(cbind(N − y, y) ~ 0 + unit + dD + dQ + hinge + dT + (1 | arm), family = binomial("cloglog"), nAGQ = 25)` at the maximum-likelihood (α/β, Tk) must reproduce the Python coefficients. Then α = −b_dD and α/β = b_dD / b_dQ. glmer allows `nAGQ > 1` because there is a single scalar random effect.

### 3.2 Acute mucositis (model F-E), matching unit U02

logit p_i = m0 + δ·grp_j + τ u_j + κ(dD_i + dQ_i / (α/β) − r·dTe_i) + γ·chemo_i + σ_arm z_i,  u_j, z_i ~ N(0, 1).

This is `tools/early_tox/step5_model.stan` as run for U02: cumulative chemotherapy coding, CTCAE scale group, arm overdispersion, Tk 7 days elapsed, repair half-time 3.2 h, exact hinge.

- **Marginal likelihood:** nested Gauss–Hermite. The outer integral is over the family effect u_j and the inner one over each arm's z_i given u_j (20 × 20 nodes). Families are independent, so the cost is small.
- **Constraints:** κ ≥ 0 (as the half-normal prior implies), α/β > 0, r > 0, τ ≥ 0, σ_arm ≥ 0.
- **Cross-checks:**
  - **Exact:** with τ = σ_arm = 0 the model is a plain binomial GLM. `glm()` in R (linear form: coefficients of dD, dQ, dTe, giving κ = b_dD, α/β = b_dD / b_dQ, r = −b_dTe / b_dD) must match the Python fit with the random effects switched off.
  - **Descriptive:** `glmer(cbind(y, N − y) ~ dD + dQ + dTe + chemo + grp + (1 | family) + (1 | arm), binomial)`, which uses the Laplace approximation because it has two random-effect terms. Differences from the quadrature fit are reported, not gated. Laplace is known to bias binomial variance components with small counts, so only κ, α/β and r are compared numerically, and differences in τ and σ_arm are expected (audit finding 8).

### 3.3 Not refitted, and why

- **Calibrated DD:** its onset constant is solved per DI draw, and α, α/β and *g* are copied from that draw. A likelihood-only analogue would be a free DD fit, which did not sample reliably in Stan and which `aicc_check.py` already covers by maximum likelihood. The calculator does not use it.
- **Late toxicity:** not fitted (protocol L0 pending).

## 4. Outputs and comparison

For F-T: α, α/β, λ, Tk, g, λ + g (the late-time regrowth rate, which is better identified than either term), σ, and the effective-BED contrast against 35 × 2 Gy in 47 days for every regimen in `fit/fit_refit.py::REGIMENS` (the quantity `delta_eval` in `di.stan`).

For F-E: κ (and L* = 10κ), α/β, r, γ, τ, σ_arm, and the log odds ratio of each decision schedule against 70 Gy in 35 fractions (`lor_dec`).

For each quantity:

1. the maximum-likelihood estimate,
2. a 95% profile-likelihood interval (a drop of 1.92 in log-likelihood), marked "reaches the search limit" when it does,
3. a 95% parametric-bootstrap percentile interval (decision D3: B = 1000 for the tumor and 500 for the early model),
4. the prior 95% range and the posterior median with its 95% interval, from the deployed fits.

**Qualifications (audit findings 4–6):**

- **Flags.** An interval for a parameter at a bound (σ, g, κ, λ, Tk) or on the λ–g ridge is flagged "approximate (at a bound)". The test of λ = 0 uses a parametric bootstrap of the likelihood ratio, not χ², because Tk is not identified when λ = 0.
- **σ bias.** ML σ is biased low with 21 fixed unit intercepts (about 18 residual degrees of freedom). σ is reported as estimated and with a bootstrap bias correction. The contrast intervals are also re-profiled with σ fixed at the posterior median (0.065) and at the posterior 97.5% point (0.17).
- **Bootstrap failures.** Every replicate's convergence is recorded. The failure rate is reported, and intervals are given with and without failed replicates. The run stops for review if more than 5% fail.

**Interpretation rule, declared before any fit is run:**

Widths are compared on each prior's natural scale: log for α, α/β, λ, r and κ, and identity for Tk, g, γ, τ and σ (audit finding 7). The label uses the shape of the profile, not an interval flagged at a bound.

- **Data-dominated:** the likelihood interval is finite and narrower than half the prior's central 95% range.
- **Prior-dominated:** the likelihood interval reaches a search limit, or is wider than the prior range.
- **Partly informed:** anything between the two.

Each parameter and contrast gets one of these labels in the write-up.

**Figures:** profile log-likelihood curves for tumor α/β, Tk and g, and for mucosal α/β and r, each with the prior density and the posterior interval overlaid.

## 5. Files

- `fit/frequentist/likelihood.py`: F-T and F-E negative log marginal likelihoods, quadrature, and contrast functions (no I/O).
- `fit/frequentist/fit_ml.py`: multistart fits, profiles, and the bootstrap. Bootstrap replicates are written one row at a time to `fit/output_frequentist_2026_11_xx/<model>/bootstrap.csv`, with skip-if-present on resume, a per-replicate seed (base seed + replicate index, PCG64), and a progress line per replicate (Rule 12).
- `fit/frequentist/glmer_check.R`: the R cross-checks, run with Windows `Rscript.exe` from the script's directory with bare filenames (Windows R cannot resolve `/mnt/c` paths). The `lme4` version is recorded in the output, and the testthat file uses `skip_if_not_installed("lme4")` so the release suite is unchanged (audit finding 9). `lme4` is in the user library, not the release `renv` library, so this script never runs in the app.
- `tests/python/test_frequentist.py` and `tests/testthat/test_frequentist_glmer.R`.
- `docs/development/2026_11_xx_frequentist_comparison.md`: results, labels and figures. Linked from `docs/FIXES.md` only if anything outside `fit/` and `docs/` changes.
- `constant_registry.yaml`: quadrature nodes, bootstrap B, grids, seeds and the Tk search range. `sim_provenance.yaml`: the parameter values used in the recovery test (T5).

No change to `app.R`, `R/`, `draws/` or the deployed site under this plan.

## 6. Tests (each to be made to fail on purpose)

- **T1. Identical design.** The F-T and F-E designs equal the arrays the Stan fits used, and their sha256 matches the data hash in `draws/model_methods.json` (tumor) and the U02 manifest identity (early). Red: drop one arm.
- **T2. Known value.** With σ = 0, α/β = 10 and *g* within 0.0046–0.021 (the `aicc_check.py` logistic transform), F-T reproduces `main_free.json` DI after Σ log C(N, y) is added to its 5342.4066 (to 1e-3), together with α, λ, Tk and *g* (to 1e-3 relative). Red: swap controls and failures, or drop the constant from one side.
- **T3. Quadrature.** For single arms, at σ ∈ {0.05, 0.25, 1}, Gauss–Hermite matches `scipy.integrate.quad` to 1e-9. At σ = 0 it equals the binomial likelihood exactly. Red: use 3 nodes at σ = 1.
- **T4. Contrast function.** `delta_eval` is not stored per draw (`draws/tumor_di_primary.csv` has only the parameters), so the tumor test is done in two parts. First, the F-T contrast at 20 rows of that file equals the app's R contrast (`R/bed_models.R`) to 1e-9. Second, applied to all 4000 rows, its medians match `summary.json` `contrasts_vs_anchor` to 0.1 Gy. For F-E, `lor_dec` reproduces the four stored columns of the U02 `draws.csv.gz` to 1e-6. Red: drop the Q / (α/β) term (audit finding 1).
- **T5. Recovery.** Simulate from known parameters on the real designs with N multiplied by 1000. Every identified parameter is recovered to within 2% (a consistency test of the code). Red: wrong sign on the time term.
- **T6. R cross-checks.** glmer F-T at fixed (α/β, Tk) and glm F-E with τ = σ_arm = 0 match Python to 1e-4 relative. Red: use the logit link in F-T.
- **T8. Audit tests.** `test_sigma_bias_recovery` (finding 4), `test_ridge_detected` (finding 5, Tk = 5 against Tk = 35), `test_bootstrap_records_failures` (finding 6) and `test_label_rule` (finding 7), each made red as described in the audit.
- **T7. Golden value.** The F-T and F-E maximum-likelihood estimates are pinned, with scipy, numpy and lme4 versions and a 1e-4 relative tolerance.

Order: T1–T6 written first and each shown red, then a smoke run (40 starts, B = 20), then T1–T6 again, then the full bootstrap.

## 7. Implementation order

1. Igor rules on D0–D4 (Section 8).
2. Claude: `likelihood.py` with T1–T5 (about 4 hours).
3. Claude: `fit_ml.py` fits and profiles, then `glmer_check.R` and T6 (about 3 hours).
4. Gate: T1–T6 green and each confirmed red when broken, and the combined smoke run (B = 20) green.
5. Claude: full bootstrap, run in the background with checkpoints (estimate 1 to 3 worker-CPU hours, to be measured on the smoke run).
6. Claude: write-up, labels and figures, with T7 pinned.
7. Igor: reads the write-up and decides D2 (methods page).

## 8. Decisions for Igor

- **D0.** Check, not switch. Recommended: yes.
- **D1.** Constraints in the primary tumor fit. Recommended: physical constraints only, with the Stan bounds as the secondary. The primary fit shows what the data say, and the secondary isolates the effect of the priors from the effect of the bounds.
- **D2.** Methods page, after the results. Recommended: one sentence and the α/β profile figure if α/β is labeled prior-dominated, and nothing otherwise.
- **D3.** Bootstrap size. Recommended: B = 1000 (tumor), 500 (early).
- **D4.** Leave calibrated DD out. Recommended: yes (Section 3.3).

## 9. Stopping condition

Done when T1–T7 pass on the revision that produced the results, every parameter and contrast has a label, and the write-up is committed. Reopen only if a test fails or the deployed fits change.
