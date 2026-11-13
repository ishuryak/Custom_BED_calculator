# Plan: likelihood-based and Bayesian estimates under the same observation model

Date: 2026-11-12. Status: **revision 3, the governing specification.** Revision 2 (SHA-256 `44d2324f…`, committed unchanged in `aa6c438`) was reviewed externally the same day: `docs/development/2026_11_12_frequentist_review/` (review findings F01–F06). This revision applies all six findings. It was written and committed before any fit of real data was run. The code that implements it is `fit/frequentist/`.

Author: Igor Shuryak (plan drafted with Claude). The audit of revision 2 (`2026_11_12_frequentist_comparison_plan_audit.md`) is superseded: see the correction note at its top.

## 0. Decisions (recorded 2026-11-12)

Igor Shuryak asked for the review to be implemented without further questions. The review's recommended answer was taken for each decision.

| ID | Decision | Answer |
|---|---|---|
| D0 | Offline check or app switch | **Offline check.** App choices and defaults stay as they are. |
| D1 | Supports | **Matched Bayesian support first, expanded support second** (Section 3). Observed-duration limits are not called physical constraints. |
| D2 | Public methods update | **Link the complete report whatever its result.** It covers all prespecified schedule contrasts, not only α/β. Any link on the site waits for a separately authorized release. |
| D3 | Bootstrap sizes | **B = 1000 (tumor) and 500 (early) as provisional budgets.** They run only after the measured pilot (B = 20), with Monte Carlo precision reported for every interval endpoint. |
| D4 | Calibrated DD | **Not included in this plan.** A conditional plug-in calibration with a repeated two-stage bootstrap could be defined, but that would be additional work. |

## 1. Question

How do likelihood-based estimates and Bayesian estimates compare under the same observation model? This is an estimation and prior-sensitivity analysis. It is **not** a decomposition of the results into "data" and "assumptions". The two procedures also differ in other ways:

- nuisance parameters are profiled in one and integrated in the other,
- the summaries differ (mode versus median),
- where the supports differ, that is a further difference.

Both procedures keep the normal random-effect distributions and the binomial observation model.

Code correctness is established only by the tests in Section 6, never by agreement with the posterior.

**Clinical qualification, printed in the report:**

> These intervals describe parameter and contrast uncertainty under the specified likelihood and design. They do not quantify the full uncertainty of the published clinical endpoints or establish transport to a new regimen or population. Holding the likelihood constant makes this a sensitivity analysis to estimation and prior assumptions, not a new clinical validation.

The binomial bootstrap cannot restore censoring-aware uncertainty in the 45 tumor pseudo-counts, which are rounded from 39 actuarial and 6 cumulative-incidence estimates. Nor can it restore missing mucositis peaks, delivered chemotherapy histories or patient mucosal dose distributions.

## 2. Models

Designs are imported from the builders of the deployed fits (`fit/fit_refit.py`, `tools/early_tox/step5_fit.py`). A run refuses unless the retained Bayesian identities reproduce: the DI identity `2c828af6…` and its Stan-data hash, and the U02 identity `4747b297…`.

- **Tumor, 45 arms in 21 fitting units, N = 8,279.** η_i = c_j − α(B_i − B̄_j) + λ smax(T_i − Tk) + g(T_i − T̄_j) + σ z_i. Here B = D + Σd²/(α/β), TCP = exp(−e^η), and the cloglog link is on failures. The unit intercepts c_j are profiled. The arm effect z_i ~ N(0, 1) is integrated.
- **Early (U02), 51 arms in 22 families, 6,042 patients, 2,945 events.** logit p = m0 + δ·grp + τu_j + κ(dD + dQ/(α/β) − r·dTe) + γ·chemo + σ_arm z_i. The terms dD, dQ and dTe are taken at onset 7 elapsed days and repair half-time 3.2 h, with the exact hinge. m0, δ and γ are profiled. u_j and z_i ~ N(0, 1) are integrated, z nested in u.

All log-likelihoods are full, including log C(N, y).

## 3. Supports (D1)

Each bound is either a **support** edge (the Bayesian model's own, or a stated physical assumption) or a **computational** edge. An estimate or interval end reaching a computational edge is reported as **open**, never as a finite limit.

| Parameter | Matched (primary) | Expanded (tumor, secondary) |
|---|---|---|
| α | (0, 2], with 2 computational | same |
| α/β | [0.05, 1000], both computational (Stan: > 0) | same |
| λ | [0, 2], with 2 computational | same |
| Tk | [0, 60] days, both support | [0, 120], with 120 computational. Onset after every observed time is not identified |
| g | [0.0046, 0.021] per day, both support | [0, 0.2]: 0 is a physical assumption (no net clonogen loss without treatment), 0.2 is computational |
| σ | [0, 5], with 5 computational | same |
| c_j | [−30, 30] computational | same |
| Early κ, τ, σ_arm, r | ≥ 0 (support), upper ends 1, 10, 10, 20 computational | none: the matched support is already the scientifically defensible one |
| Early α/β | [0.05, 1000] computational | none |
| m0, δ, γ | [−20, 20] computational | none |

The two tumor supports are separate runs with separate identities and are reported separately.

## 4. Targets, intervals and comparison

**Tumor targets:**

- α, α/β, λ, Tk, g and σ,
- λ + g (per day),
- (λ + g)/α (Gy of BED per day), and its EQD2 conversion at each estimate's α/β (Gy EQD2 per day),
- the effective tumor-BED contrast of each of the ten regimens in `fit/fit_refit.py::REGIMENS` against 35 × 2 Gy in 47 days (Gy of BED at that estimate's α/β).

**Early targets:**

- κ, α/β, r, γ, τ and σ_arm,
- the log odds ratio of each of the four decision schedules (15 × 3.5, 20 × 2.75, 25 × 2.4 and 30 × 2.2 Gy) against 70 Gy / 35.

**Profile likelihood.** For a parameter, it is fixed by its bounds and every other parameter is reoptimized. For a derived target q, the profile solves max ℓ(θ) subject to q(θ) = q0, using SLSQP with the equality constraint from the warm start plus two perturbed restarts. The constraint residual must be below 1e−7 × max(1, |q0|) and is checked independently.

The profile is walked outward from the estimate in both directions. A walk stops when the log-likelihood has fallen 0.5 below the cutoff and stayed below it for four more doubling steps, or when it reaches a limit. Crossings are refined by root finding.

- The cutoff 1.92 (χ²₁ 95% / 2) is the regular one-parameter approximation.
- The confidence set is every region above the cutoff, so disconnected sets are reported as such.
- An end that reaches a support or computational limit is reported **open**, with the kind of limit.
- A profile that never drops by 1.92 is reported **flat**: that is a valid result, not a failure.
- A profile point above the estimate's log-likelihood means the estimate was not global. The estimate is refitted from that point and the profiles are rerun.

Where the estimate sits on a support edge, or on the λ–g ridge, intervals are labeled **nonregular**. No coverage is claimed for them.

**Parametric bootstrap (matched support only).**

- **Tumor:** fresh independent arm effects with the fitted unit intercepts, then binomial counts on the fixed denominators.
- **Early:** fresh family and arm effects, then counts.

The design is fixed. Each replicate is refitted with the same integration, constraints and optimizer rules. Intervals are the 2.5% and 97.5% percentiles. Each endpoint gets a distribution-free 95% Monte Carlo interval from order statistics.

Bias is reported as mean(θ*) − θ̂ with its Monte Carlo standard error, as a descriptive finite-cluster quantity. No bias-corrected interval is constructed.

**Failure policy:**

- Every attempt is saved.
- A failed replicate is retried, at most twice, with 8 random starts and doubled nodes. Its seed and simulated counts never change.
- A legitimate boundary optimum counts as a success and is labeled as one.
- Unresolved failures are kept with their reason.
- Any unresolved failure makes the bootstrap intervals **provisional**. More than 5% blocks a final interval claim and triggers numerical remediation.

**Not run:** the null-simulation likelihood-ratio test of λ = 0. It is not needed for this comparison, and no regular χ² p-value is reported for it. A separate repeated-simulation coverage study is also not run: the bootstrap already gives finite-cluster bias, and coverage is not claimed.

**Hybrid sensitivity (externally fixed, labeled as such).** The ten tumor contrasts are re-profiled with σ fixed at the posterior median 0.065 and at the posterior 97.5% point 0.17.

**Comparison table.** For each target the report shows, side by side:

- the maximum-likelihood estimate and its profile confidence set,
- the bootstrap percentile interval,
- the posterior median and 95% credible interval, computed from the deployed draws with the same contrast function. For the tumor these are 4,000 rows holding 1,000 distinct DI draws, so the quantiles have the precision of 1,000 draws.
- the prior 95% interval.

For a derived target, the prior interval is simulated from the full declared joint prior through the exact contrast function (seed 20261112, 200,000 draws). It is never formed from marginal endpoints.

The ratio of the likelihood-interval width to the prior-interval width is reported as an **analyst-defined descriptive measure**. It uses a log scale for α, α/β, λ and r (lognormal priors) and the natural scale for all other targets: Tk, g, σ, κ (half-normal prior), γ, τ, σ_arm and every contrast.

The report shows this ratio beside the profile shape, the bound dependence and the change in each schedule contrast. It assigns no "data-dominated" or "prior-dominated" label, and a numerical search limit is never read as prior dominance.

## 5. Numerical contract

- **Integration.** Adaptive Gauss–Hermite quadrature centered at each integrand's mode (bracketed Newton, which is safe because the integrands are log-concave) and scaled by its curvature. Zero variances use the closed form. The early outer integrand is the sum of the inner log integrals of a family's arms.
- **Accuracy.** Every reported optimum is re-evaluated at doubled nodes. It is accepted only if the total log-likelihood changes by less than 1e−6, and otherwise nodes escalate and the fit is repeated. Contrasts must be stable to 1e−4 in natural units. Independent references (scipy QUADPACK, and a dense double grid for the nested integral) must agree to 1e−9.
- **Gradients.** Analytic, by the Fisher identity at the quadrature weights. They are checked against central differences at interior points. Convergence is the projected-gradient (KKT) criterion, max |projected ∂ℓ/∂x × scale| < 1e−3. A boundary optimum needs no positive-definite Hessian.
- **Starts.**
  - Maximum-likelihood fit: the posterior medians plus 15 random starts (seeds `[20261112, 1]`).
  - Bootstrap replicate: the fitted estimate, the posterior medians and 2 random starts.
  - Retry: 8 random starts.
  - Starts within 1e−3 of the best log-likelihood count as ties. The spread of every target across ties is recorded as a ridge indicator.
- **Run identity** (`fit_ml.run_identity`). It covers:
  - the input-file and canonical-design hashes, and the ordered IDs,
  - model, support, cutoff and quadrature policy,
  - optimizer and start policy, and the seed,
  - the code hashes and software versions.

  A changed dataset, model, tolerance or support starts a new generation folder, `fit/output_frequentist/<model>_<support>_<identity16>/`.
- **Records.** Each bootstrap replicate is one JSON record, written atomically by temporary file and rename, with:
  - its identity, seed, attempts and status (and why),
  - the hash of its simulated counts,
  - its estimates, log-likelihood, nodes, gradient check, active bounds and timing,
  - a sha256 of its own content.

  A resume accepts only verifying records of the same identity, and never reruns a completed one. `summary.json` and then `MANIFEST.json` are written only after every record verifies. A generation without its manifest is incomplete.

## 6. Tests (`tests/python/test_frequentist.py`, each made red on purpose before it was accepted)

- **T1. Identities.** The Bayesian identities reproduce. The design totals are correct: 45/21/8,279 for the tumor, and 51/22/6,042/2,945 for the early model. Changing one count, the row order, one time convention or one setting moves the right layer: file and design for data changes, the run identity only for a setting.
- **T2. Known value.** The full tumor negative log-likelihood of the archived reference configuration (σ = 0, α/β = 10, g within the 2018 bounds, intercepts reoptimized) is **131.29112756**. That is the archived 5342.406617400148 minus Σ log C(N, y) = 5211.115489844190. Adding the constant, as revision 2 said, is refused.
- **T3. Integration.**
  - the review's counterexamples, where ordinary fixed-node quadrature is off by −6.87 and −7.97 log units, are detected,
  - single-arm, nested-family and whole-likelihood values match independent references to 1e−9,
  - the zero-variance branches are exact,
  - the ladder converges at large variances,
  - the mode search does not cycle.
- **T4. Contrasts.**
  - Tumor: the medians over the 4,000 deployed rows (1,000 distinct draws) are within 0.1 Gy of the stored posterior summary. This is a versioned sanity check, not an identity.
  - Early: `lor_dec` from the six-significant-digit stored rows agrees to 2e−5, the measured rounding bound.
  - Unrounded formula identities are checked to 1e−12.
- **T5.** Deterministic identities only: likelihood, gradient against finite differences, and contrast formulas. The inflated-N 2% recovery requirement of revision 2 is withdrawn, because one dataset cannot recover population variances from 22 realized effects.
- **T6. R cross-checks** (`fit/frequentist/glmer_check.R`, run by Windows `Rscript.exe` with `lme4`). A missing `lme4` fails this required check rather than passing it.
  - **Tumor:** `glmer(..., binomial("cloglog"), nAGQ = 25)` at the fitted (α/β, Tk), with the single predictor B (coefficient −α), the hinge, centered T and unit fixed effects. It is compared with the Python fit on the same fixed (α/β, Tk) and the same unconstrained support. The Python likelihood is evaluated at R's coefficients and must agree within 1e−4. Coefficients must agree within 1e−6 + 1e−4 relative.
  - **Early:** `glm()` with τ = σ_arm = 0, on synthetic counts simulated from known interior parameters on the real design. The profile interval of κ must equal R's profile interval of the dD coefficient.
  - The two-random-effect `glmer` Laplace fit is descriptive only.
- **T7. Golden values.** The maximum-likelihood objectives and identifiable predictions are pinned with the numpy and scipy versions. Arbitrary coordinates along a ridge are not pinned.
- **T8. Inference and restart controls.**
  - A quadratic fixture recovers known profile intervals for a parameter and a derived target.
  - A boundary fixture reports an open end at its support limit.
  - A flat fixture reports flat, and a bimodal fixture reports a disconnected set.
  - Resume reruns nothing. An edited or truncated record is rerun with the same seed and counts. A stale identity is refused.
  - A failed record is retried with its seed unchanged.
  - A generation without `MANIFEST.json`, or with a record of another identity, is refused by the report assembler.

Gate order: T1–T8 green → combined pilot (both models, B = 20, interrupted and resumed once) → T1–T8 again → measured budget → full run.

## 7. Outputs

- `fit/output_frequentist/…`: estimates, profiles, bootstrap records, summaries and manifests. These are aggregate, trial-level quantities with no patient data.
- `docs/development/2026_11_12_frequentist_comparison_results.md`: every prespecified parameter and contrast, including inconclusive and non-identified ones, with:
  - the comparison table,
  - profile figures for tumor α/β, Tk, g and λ + g, and for early α/β and r, each with its prior density and the posterior interval,
  - the bootstrap failure summary and the measured cost.

  It is published whatever its direction.
- No change to `app.R`, `R/`, `draws/` or any default.

## 8. Stopping condition

Done when T1–T8 pass on the revision that produced the results and the report is committed. It must contain every prespecified target with its interval method, qualifications and Monte Carlo precision. Reopen only for a failed criterion, a changed target, or new evidence.
