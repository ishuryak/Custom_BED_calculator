# Likelihood-based and Bayesian estimates under the same observation model: results

Date: 2026-11-12. Specification: `2026_11_12_frequentist_comparison_plan.md`, revision 3 (commit `c373941`), written before any fit of real data. Run at revision `4a639ce` by `fit/frequentist/run_full.sh` (log `fit/output_frequentist/run_full.log`). Generations:

- `tumor_matched_f0e067bed39d276d`
- `tumor_expanded_5affd91c876d9b96`
- `early_early_matched_48a6a4f0039e68f0`

Every number below is generated from the verified records by `fit/frequentist/compare.py` and `report_tables.py`, not typed by hand.

> These intervals describe parameter and contrast uncertainty under the specified likelihood and design. They do not quantify the full uncertainty of the published clinical endpoints or establish transport to a new regimen or population. Holding the likelihood constant makes this a sensitivity analysis to estimation and prior assumptions, not a new clinical validation.

The binomial bootstrap cannot restore censoring-aware uncertainty in the 45 tumor pseudo-counts (39 actuarial and 6 cumulative-incidence estimates), missing mucositis peaks, delivered chemotherapy histories or patient mucosal dose distributions.

## Summary

**The schedule comparisons the calculator uses are informed by the trials, not mainly by the priors.**

- **Early model.** The four decision log odds ratios have likelihood intervals close to the posterior intervals. For 15 × 3.5 Gy against 70 Gy in 35 fractions, the estimate is 1.39 (profile 0.84 to 1.98) against a posterior of 1.39 (0.84 to 1.96).
- **Tumor model.** The likelihood intervals of the ten effective-BED contrasts are about as wide as the posterior intervals. They are shifted, though, in a pattern that follows the higher repopulation rate under maximum likelihood (1.04 against 0.90 Gy of BED per day at the posterior median). The 19- to 26-day courses lose 1.1 to 1.7 Gy more tumor BED. Accelerated courses (6 fractions a week, twice daily, or 65 Gy in 30 fractions over 40 days) gain 0.9 to 1.6 Gy more. 30 × 2 Gy in 40 days and 66 Gy in 33 fractions over 45 days lose 1.6 and 0.5 Gy less. For HYPNO 55 Gy in 20 fractions the estimate is −3.0 Gy (profile −6.7 to 1.8) against a posterior of −1.9 (−5.6 to 3.2).

**Several individual parameters are not bounded by the trials, and there the published posterior range comes from the prior.**

- **Tumor α/β.** The estimate is 13.4 Gy and the profile interval runs from 6.6 Gy to the 1,000 Gy computational limit. The log-likelihood is still within 1.92 of its maximum at 1,000 Gy. The posterior upper limit of 18.9 Gy is therefore the prior's (5 to 20 Gy). The lower limits agree (6.6 against 6.9 Gy).
- **Mucosal α/β.** The likelihood rises steadily toward large α/β, so the estimate is the 1,000 Gy computational limit. The profile set runs from 11.3 Gy upward. The posterior (12.9 Gy, 6.9 to 24.0) reaches below 11.3 Gy only because the prior is centered at 10 Gy. A Laplace `glmer` fit, run as an independent check, also estimates a negative quadratic-dose coefficient.
- **Background growth g.** The profile is flat across the whole prior range (0.0046 to 0.021 per day), so g is not identified. The sum λ + g is well identified: 0.041 per day (0.029 to 0.056) against a posterior of 0.040 (0.025 to 0.059).
- **Mucosal time-loss rate r.** The estimate is 1.75 Gy of BED per day (profile 1.08 to 4.87). The posterior upper limit of 2.71 is held down by the prior (95% range 0.30 to 2.13).
- **Tumor arm spread σ.** The estimate is exactly 0, on the support edge (posterior median 0.065). That makes σ, and every interval resampled at σ = 0, nonregular.

**The tumor contrasts are sensitive to the arm-level spread.** With σ held at the posterior 97.5% point (0.17), a hybrid sensitivity that is not a likelihood-only result, the HYPNO profile interval widens from −6.7 to 1.8 Gy to −7.9 to 9.6 Gy. The estimate moves by 1.3 Gy, to −1.8 Gy.

**Practical relevance.**

- **No default changes.** Nothing here changes a default, a draw or a number in the app.
- **The research early-toxicity scenario is robust to its α/β prior.** Its schedule comparisons do not depend on the prior on mucosal α/β.
- **The tumor model's losses for the shortest courses may be understated.** Under maximum likelihood, the 19- to 26-day courses would lose 1.1 to 1.7 Gy more than the posterior-median values the calculator uses. That shift is smaller than the calculator's default ±5 Gy margins. Whether that matters is a decision about the margins, not something this analysis settles.
- **Quote the α/β values as assumptions.** Neither α/β upper limit should be quoted as an estimate from the trials.

## How the results were checked

- **Bayesian identities.** Both reproduced from the current data and code before the run (`fit_ml.py validate`).
- **Maximum-likelihood fits.**
  - Every fit used 16 starts. The 12 (tumor) and 13 (early) starts that reached the best log-likelihood agree on every target to about 1e−5, so the optimum is unique and not a ridge.
  - The projected gradient is below 1e−5, and the quadrature change at doubled nodes is below 1e−12.
  - No profile point exceeded the estimate's log-likelihood, so every estimate is global.
- **Bootstraps.** The tumor bootstrap has 1,000 of 1,000 replicates converged, the early bootstrap 500 of 500, with no unresolved failures. Both generations load only through the manifest check.
- **R cross-checks (T6).** Both passed again on the final revision (`fit/output_frequentist/r_check/r_check.json`).
- **Tests.** All 37 tests (`tests/python/test_frequentist.py`) pass on the run revision. Each was made red on purpose: 13 mutations, `2026_11_12_frequentist_review/mutation_checks.json`.
- **Cost.** 3.4 CPU-seconds per tumor replicate and 32.5 per early replicate. The whole run took 1 h 33 min of wall time on 20 workers, against a pilot-based estimate of about 30 minutes. The difference came from the derived-target profiles, which took up to 20 minutes each against about 2.5 minutes in the pilot timing.

## Qualifications and nonregular results

- **σ = 0 (tumor).** The estimate is on the support edge, and 998 of the 1,000 bootstrap replicates also put σ at 0, because they are generated at σ = 0. The bootstrap interval for σ (0 to 0) carries no information. The profile set, 0 to 0.047, is one-sided and nonregular.
- **g.** It hits a prior bound in 456 of the 1,000 bootstrap replicates. Its profile is flat, so neither interval is a confidence statement about g.
- **α/β (tumor).** The bootstrap upper percentile (115 Gy) is unstable: its 95% Monte Carlo interval is 80 to 752 Gy, and 15 replicates reached the 1,000 Gy computational limit. Only the lower end (6.6 Gy, Monte Carlo 6.4 to 6.9) is a usable bootstrap limit.
- **α/β (mucosa).** 231 of the 500 replicates reached the 1,000 Gy limit, so the bootstrap upper limit is that limit by construction.
- **68 × 1.2 Gy twice daily.** This contrast has the same overall time as the reference, so it depends only on α/β. Its upper profile end is open: it approaches its value as α/β → ∞, about 11.6 Gy.
- **One endpoint not resolved.** In the hybrid sensitivity with σ = 0.17, the lower end for 33 × 2 Gy at 6 fractions a week is listed as "open". Root finding could not bracket it, because the refinement restarted from the overall estimate and found a worse local optimum. The walked profile places the crossing between −0.43 and +0.43 Gy: the log-likelihood drops by 2.22 at −0.43 and by 1.74 at +0.43. Fixing the refinement would change the code and therefore the run identity of every generation, so it is recorded as an open item for the next run and not patched here.
- **Width ratios** are descriptive only and are not used to label anything. For the tumor contrasts they are small (0.03 to 0.06), because the declared joint prior implies extremely wide contrast ranges: its upper ends of 28 to 292 Gy come from small α under the lognormal(log 0.2, 1) prior.
- **Expanded tumor support.** It gives the same estimates. With g allowed down to 0 and Tk up to 120 days, most contrast intervals widen by up to about 1.4 Gy at one end (for example HYPNO, lower end −7.2 against −6.7 Gy). The upper end for g closes at 0.0216 per day.

## Open items

1. Warm-start the profile endpoint refinement from the neighbouring walked points. That fixes the one unresolved endpoint above, and the run identity of the next generation will record the change.
2. `about.md` line 153 ("Every max(0, x) is computed as …") predates the exact-hinge research early model. It needs a separate text pass.
3. A methods-page link to this report waits for a separately authorized release (decision D2).

## Tables

### Tumor DI, matched support

Log-likelihood at the estimate −131.097885; estimate on a limit: sigma at its lower support limit.

| Target | Units | ML estimate | Profile 95% set | Bootstrap 2.5–97.5% | Posterior median (95%) | Prior 95% | Width ratio |
|---|---|---|---|---|---|---|---|
| α | per Gy | 0.03977 | 0.0279 to 0.0551 | 0.03 to 0.0551 | 0.0443 (0.0301 to 0.0602) | 0.0281 to 1.42 | 0.17 |
| α/β | Gy | 13.43 | 6.57 to 1000 (open) | 6.59 to 115 | 10.6 (6.87 to 18.9) | 5 to 20 | – |
| λ | per day | 0.03455 | 0.0171 to 0.0494 | 0.0206 to 0.0559 | 0.0277 (0.0102 to 0.0485) | 0.00528 to 1.88 | 0.18 |
| Tk | days | 38.29 | 32.7 to 42.2 | 29.2 to 44.8 | 37.3 (23.7 to 42.8) | 8.88 to 47.4 | 0.25 |
| g | per day | 0.006735 | 0.0046 (open) to 0.021 (open) (flat) | 0.0046 to 0.021 | 0.0121 (0.00495 to 0.0204) | 0.00501 to 0.0206 | – |
| σ (arm spread) | cloglog scale | 0 | 0 (open) to 0.0474 | 0 to 0 | 0.065 (0.00387 to 0.174) | 0.0156 to 1.12 | – |
| λ + g | per day | 0.04128 | 0.0294 to 0.0556 | 0.0313 to 0.0679 | 0.0395 (0.0254 to 0.0592) | 0.0168 to 1.9 | 0.014 |
| (λ + g)/α | Gy of BED per day | 1.038 | 0.702 to 1.55 | 0.693 to 1.7 | 0.903 (0.547 to 1.46) | 0.0344 to 17 | 0.05 |
| (λ + g)/α in EQD2 | Gy EQD2 per day | 0.9035 | 0.599 to 1.36 | 0.605 to 1.47 | 0.756 (0.464 to 1.24) | 0.0283 to 14.1 | 0.054 |
| ΔBED, HYPNO 55 Gy/20 (26 d) | Gy of tumor BED | −3.041 | −6.65 to 1.77 | −6.35 to 1.29 | −1.95 (−5.62 to 3.19) | −13.4 to 266 | 0.03 |
| ΔBED, HYHOPE-type 46.5 Gy/15 (19 d) | Gy of tumor BED | −10.88 | −14.6 to −6.1 | −14.2 to −5.91 | −9.18 (−13.2 to −3.46) | −22.5 to 292 | 0.027 |
| ΔBED, 17 × 3.05 Gy (23 d) | Gy of tumor BED | −5.169 | −9.33 to 0.0202 | −8.78 to −0.129 | −3.64 (−7.7 to 2.02) | −16 to 281 | 0.031 |
| ΔBED, 60 Gy/25 (33 d) | Gy of tumor BED | 0.2342 | −2.67 to 4.93 | −2.88 to 3.54 | 0.0934 (−2.99 to 4.56) | −9.23 to 200 | 0.036 |
| ΔBED, 65 Gy/30 (40 d) | Gy of tumor BED | 2.326 | −0.0248 to 5.73 | −0.414 to 4.86 | 1.22 (−1.13 to 4.47) | −4.72 to 107 | 0.051 |
| ΔBED, 66 Gy/33 (45 d) | Gy of tumor BED | −2.52 | −3.43 to −1.4 | −3.4 to −1.28 | −2.97 (−3.7 to −1.83) | −5.12 to 28.1 | 0.061 |
| ΔBED, 30 × 2 Gy (40 d) | Gy of tumor BED | −4.225 | −7.05 to −0.805 | −7.41 to −1.9 | −5.8 (−8.19 to −2.53) | −12.6 to 100 | 0.055 |
| ΔBED, 33 × 2 Gy, 6 per week (38 d) | Gy of tumor BED | 4.393 | 1.2 to 7.85 | 0.73 to 6.67 | 2.75 (0.0511 to 6.46) | −4.71 to 137 | 0.047 |
| ΔBED, 70 Gy/35, 6 per week (40 d) | Gy of tumor BED | 7.264 | 4.91 to 10.7 | 4.52 to 9.78 | 6.14 (3.79 to 9.39) | 0.195 to 112 | 0.051 |
| ΔBED, 68 × 1.2 Gy twice daily (47 d) | Gy of tumor BED | 8.467 | 5.2 to 11.3 (open) | 5.21 to 11.2 | 7.62 (5.47 to 9.38) | 3.19 to 9.49 | – |

Bootstrap: 1000 replicates, 1000 converged, 0 unresolved failures, 998 with at least one parameter on a limit; 3.41 CPU-seconds per replicate.

Hybrid sensitivity (σ held at an externally chosen value; not a likelihood-only result):

| Regimen | σ = 0.065 (posterior median): estimate, profile 95% set | σ = 0.17 (posterior 97.5%): estimate, profile 95% set |
|---|---|---|
| HYPNO 55 Gy/20 (26 d) | −2.75: −7.15 to 2.97 | −1.77: −7.93 to 9.58 |
| HYHOPE-type 46.5 Gy/15 (19 d) | −10.9: −15.1 to −4.97 | −10.1: −15.9 to 1.32 |
| 17 × 3.05 Gy (23 d) | −5.05: −9.86 to 1.2 | −3.98: −10.7 to 8.18 |
| 60 Gy/25 (33 d) | 0.914: −3.08 to 6.04 | 1.83: −3.78 to 12.1 |
| 65 Gy/30 (40 d) | 2.3: −0.335 to 6.55 | 2.91: −0.689 to 11.7 |
| 66 Gy/33 (45 d) | −2.51: −3.53 to −1.09 | −2.47: −3.79 to 0.459 |
| 30 × 2 Gy (40 d) | −4.21: −7.37 to 0.0508 | −3.93: −8.01 to 4.78 |
| 33 × 2 Gy, 6 per week (38 d) | 4.7: 1.08 to 8.73 | 5.34: open to 14 |
| 70 Gy/35, 6 per week (40 d) | 7.24: 4.61 to 11.5 | 7.84: 4.25 to 16.6 |
| 68 × 1.2 Gy twice daily (47 d) | 8.55: 4.71 to open | 7.88: 1.69 to open |

### Tumor DI, expanded support

Log-likelihood at the estimate −131.097885; estimate on a limit: sigma at its lower support limit.

| Target | Units | ML estimate | Profile 95% set | Posterior median (95%) | Prior 95% | Width ratio |
|---|---|---|---|---|---|---|
| α | per Gy | 0.03977 | 0.0255 to 0.0551 | 0.0443 (0.0301 to 0.0602) | 0.0282 to 1.43 | 0.2 |
| α/β | Gy | 13.43 | 6.57 to 1000 (open) | 10.6 (6.87 to 18.9) | 4.97 to 20.1 | – |
| λ | per day | 0.03455 | 0.0171 to 0.0512 | 0.0277 (0.0102 to 0.0485) | 0.00538 to 1.93 | 0.19 |
| Tk | days | 38.29 | 32.5 to 42.2 | 37.3 (23.7 to 42.8) | 8.78 to 47.5 | 0.25 |
| g | per day | 0.006734 | 0 (open) to 0.0216 | 0.0121 (0.00495 to 0.0204) | 0.00501 to 0.0206 | – |
| σ (arm spread) | cloglog scale | 0 | 0 (open) to 0.0474 | 0.065 (0.00387 to 0.174) | 0.0154 to 1.12 | – |
| λ + g | per day | 0.04128 | 0.0293 to 0.0556 | 0.0395 (0.0254 to 0.0592) | 0.0169 to 1.94 | 0.014 |
| (λ + g)/α | Gy of BED per day | 1.038 | 0.702 to 1.64 | 0.903 (0.547 to 1.46) | 0.0336 to 17.8 | 0.053 |
| (λ + g)/α in EQD2 | Gy EQD2 per day | 0.9035 | 0.599 to 1.51 | 0.756 (0.464 to 1.24) | 0.0277 to 14.7 | 0.062 |
| ΔBED, HYPNO 55 Gy/20 (26 d) | Gy of tumor BED | −3.041 | −7.15 to 1.78 | −1.95 (−5.62 to 3.19) | −13.4 to 275 | 0.031 |
| ΔBED, HYHOPE-type 46.5 Gy/15 (19 d) | Gy of tumor BED | −10.88 | −15.7 to −6.1 | −9.18 (−13.2 to −3.46) | −22.5 to 306 | 0.029 |
| ΔBED, 17 × 3.05 Gy (23 d) | Gy of tumor BED | −5.169 | −10.3 to 0.0202 | −3.64 (−7.7 to 2.02) | −15.9 to 295 | 0.033 |
| ΔBED, 60 Gy/25 (33 d) | Gy of tumor BED | 0.2342 | −2.67 to 5.9 | 0.0934 (−2.99 to 4.56) | −9.23 to 210 | 0.039 |
| ΔBED, 65 Gy/30 (40 d) | Gy of tumor BED | 2.326 | −0.0248 to 6.43 | 1.22 (−1.13 to 4.47) | −4.72 to 113 | 0.055 |
| ΔBED, 66 Gy/33 (45 d) | Gy of tumor BED | −2.52 | −3.43 to −1.03 | −2.97 (−3.7 to −1.83) | −5.13 to 29.8 | 0.069 |
| ΔBED, 30 × 2 Gy (40 d) | Gy of tumor BED | −4.225 | −7.05 to 0.459 | −5.8 (−8.19 to −2.53) | −12.7 to 106 | 0.063 |
| ΔBED, 33 × 2 Gy, 6 per week (38 d) | Gy of tumor BED | 4.393 | 1.2 to 9.23 | 2.75 (0.0511 to 6.46) | −4.71 to 144 | 0.054 |
| ΔBED, 70 Gy/35, 6 per week (40 d) | Gy of tumor BED | 7.264 | 4.91 to 11.4 | 6.14 (3.79 to 9.39) | 0.194 to 118 | 0.055 |
| ΔBED, 68 × 1.2 Gy twice daily (47 d) | Gy of tumor BED | 8.467 | 5.2 to 11.3 (open) | 7.62 (5.47 to 9.38) | 3.13 to 9.5 | – |

### Acute mucositis (U02)

Log-likelihood at the estimate −179.775687; estimate on a limit: ab at its upper computational limit.

| Target | Units | ML estimate | Profile 95% set | Bootstrap 2.5–97.5% | Posterior median (95%) | Prior 95% | Width ratio |
|---|---|---|---|---|---|---|---|
| κ (slope) | log odds per Gy of early BED | 0.0442 | 0.0147 to 0.0755 | 0.0166 to 0.075 | 0.0467 (0.0245 to 0.0768) | 0.00312 to 0.224 | 0.28 |
| α/β | Gy | 1000 | 11.3 to 1000 (open) | 4.15 to 1000 | 12.9 (6.94 to 24) | 4.98 to 20 | – |
| r | Gy of BED per day | 1.749 | 1.08 to 4.87 | 1.07 to 4.38 | 1.55 (1 to 2.71) | 0.301 to 2.13 | 0.77 |
| γ (chemotherapy) | log odds | 0.5838 | 0.263 to 0.901 | 0.266 to 0.895 | 0.553 (0.2 to 0.904) | −1.97 to 1.97 | 0.16 |
| τ (trial spread) | log odds | 0.7284 | 0.526 to 1.05 | 0.446 to 0.956 | 0.785 (0.557 to 1.14) | 0.0312 to 2.24 | 0.24 |
| σ_arm (arm spread) | log odds | 0.1822 | 0.0888 to 0.31 | 0 to 0.244 | 0.24 (0.134 to 0.392) | 0.016 to 1.12 | 0.2 |
| log OR, 15 × 3.5 Gy | log odds ratio | 1.393 | 0.839 to 1.98 | 0.902 to 2.04 | 1.39 (0.84 to 1.96) | −0.398 to 5.29 | 0.2 |
| log OR, 20 × 2.75 Gy | log odds ratio | 0.9609 | 0.511 to 1.43 | 0.523 to 1.43 | 0.877 (0.438 to 1.3) | −0.826 to 3.28 | 0.22 |
| log OR, 25 × 2.4 Gy | log odds ratio | 0.6405 | 0.341 to 0.955 | 0.347 to 0.944 | 0.571 (0.278 to 0.859) | −0.604 to 2.13 | 0.22 |
| log OR, 30 × 2.2 Gy | log odds ratio | 0.3646 | 0.232 to 0.506 | 0.24 to 0.508 | 0.344 (0.208 to 0.482) | −0.12 to 1.28 | 0.2 |

Bootstrap: 500 replicates, 500 converged, 0 unresolved failures, 274 with at least one parameter on a limit; 32.5 CPU-seconds per replicate.


## Figure

![Profile log-likelihoods](2026_11_12_frequentist_comparison_results/frequentist_profiles.png)

**Profile log-likelihoods (blue) of six parameters, with the 95% cutoff (dashed), the posterior 95% interval and median (orange) and the prior density (green, rescaled so its peak is at 0).** Top row: tumor α/β (log axis) and onset Tk. Middle row: background growth g, which is flat across its prior range, and λ + g. Bottom row: mucosal α/β (log axis) and time-loss rate r. Matched support. Where a profile stays above the cutoff to the edge of the axis, the trials do not bound that parameter on that side.
