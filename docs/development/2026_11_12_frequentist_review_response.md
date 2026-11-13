# 2026-11-12: implementation of the external review of the frequentist plan and the Bayesian packaging

Igor Shuryak asked for the review (`2026_11_12_frequentist_review/`) to be checked, analyzed and implemented autonomously. This record lists what was done for each finding, the test that pins it, the evidence that the test fails without the change, and what each change moves. Branch `frequentist_review_2026_11_12`, from `79e9c81` (4.6.1 as deployed). Nothing was pushed or deployed.

**Changes published numbers: no.**

- Every draw file is byte-identical, and the tumor and early defaults are unchanged.
- The default main list still has 69 points and the research list 9.
- The frequentist results are a new offline analysis and change nothing in the app.

## Verification of the review's claims before acting

Each claim the fixes depend on was reproduced first:

- **Plan file hashes.** The reviewed plan and audit hash to `44d2324f…` and `c1b44587…`, so the review's line references apply. Both files were committed unchanged (`aa6c438`) so the corrections are attributable.
- **F02 (sign error).** Σ log C(N, y) = 5211.115489844190, and the full negative log-likelihood is 131.2911275 by an independent intercept refit (review: 131.29112756). Revision 2's instruction to add the constant was wrong.
- **F01 (quadrature).** Ordinary 30-node Gauss–Hermite misses the review's tumor fixture by −6.87156 and its RTOG 0129 outer fixture by −7.97402. Both are now pinned in T3.
- **B01.** The R loader accepted a prefix-matched `early-bed-fit-10` evaluator, a text ESS and an R-hat of −1. The exporter accepted −Inf R-hat and +Inf ESS. Each was observed before the fix.
- **B02.** Changing a draw value without touching any manifest left `build_model_methods.py --check` reporting the record current. Observed before the fix.
- **B03.** HYPNO printed pp. 391 and 394, from the review's page extracts of the user's PDF, state 753/792 evaluable for early adverse events (HFX 373, NFX 380). They also state weekly scoring during and immediately after RT and a first follow-up about two months after RT. The counts 190/373 and 208/380 match the data.
- **B04.** `how_fitted.md` showed max(0, T − Tk) for the tumor fit, which uses smax.

## Phase 2: version 4.6.2 (commit `c156fbc`)

| ID | Change | Test (observed red before the change) | What it moves |
|---|---|---|---|
| B01 | Exporter: every diagnostic must be a finite, non-boolean number with the right type and support before the gate comparison (`diagnostic_problem`). R loader: the evaluator ID is the token before the colon, matched exactly. Model and source identities must be well formed. Gate fields are type-checked (`early_gate_problem`). | `tests/python/test_early_tox_export.py::test_non_finite_or_mistyped_diagnostics_refuse` (13 cases, all red before). `tests/testthat/test_v462_review.R` B01 blocks (6 failures and 3 errors before). | Nothing for valid packages: U02 still loads, with the same point values. |
| B02 | `build_model_methods.py` hashes the actual DI, DD and early draw files and checks columns, rows, the DI/DD row pairing and the distinct-draw count. `load_model_methods()` hashes the files and validates every field the page reads. It fails safely on malformed records. `check_methods_record()` is added to `tools/deploy.R`. | `tests/python/test_model_methods.py::test_a_changed_draw_file_is_refused_even_with_its_manifest_untouched` (3 files) and `::test_unpaired_or_short_tumor_draws_are_refused` (red before). `test_v462_review.R` B02 blocks (4 failures and 1 error before). | `draws/model_methods.json` source hashes only. |
| B03 | HYPNO population and assessment schedule corrected in `data/early_tox/source/ascertainment.csv`, using a new optional column `assessment_schedule` (empty for every other trial). `endpoint_mapping.csv` and `early_fit_trials.csv` were regenerated: the diff touches only the two HYPNO rows. | `tests/python/test_early_tox_export.py::test_hypno_population_and_assessment_match_the_report` (red before). | The trial-table text. The package manifest changes only in the trial-table hash. 4.6.0/4.6.1 research scenarios still replay exactly (`EARLY_FIT_EQUIVALENT_PACKAGES`, bound to the unchanged draws hash). `test_v462_review.R` B03 block was red before. |
| B04 | The tumor equation shows smax(x) = log(1 + e^{3x})/3, its constant and its value at onset (0.23 days). The early equation shows its exact max. | `test_v462_review.R` "B04: the page shows smax…" (the old text contains `\lambda\,\max(0, T_i`, which the test refuses). | Text only. |

Version 4.6.2 (`DESCRIPTION`, `APP_VERSION` and both version tests).

**Verification on the Phase 2 tree:**

- `Rscript tools/run_tests.R` (via `testthat::test_dir`): 6,860 expectations, 0 failed, 0 errors, 0 skipped. The review's run had 6,781 with one environment error, which did not recur here.
- `pytest -m "not sampling" tests/python`: 346 passed.
- `stale_constant_audit.py --strict`: 0 errors.
- `build_model_methods.py --check`: current.

**Found and not changed (outside the review's scope):** `about.md` line 153 says "Every max(0, x) is computed as ln(1 + e^{3x})/3". The research-mode early model (4.6.0) uses an exact hinge, so "every" no longer holds there. It is left for a separate text pass.

## Phases 1, 3 and 4: the frequentist specification and engine

**Specification.** Plan revision 3 (`c373941`) was committed before any fit of real data. Decisions D0–D4 take the review's recommended answers. The revision-2 audit is marked superseded, with its four errors stated.

**Engine** (`fit/frequentist/`, commit `81fb39e`):

- full marginal log-likelihoods with adaptive Gauss–Hermite quadrature, nested for the early model,
- analytic gradients by the Fisher identity,
- designs imported from the deployed builders, with file, design and run identities,
- multistart fits with projected-gradient (KKT) and node-doubling checks,
- profile likelihoods: box-fixed for parameters, SLSQP equality-constrained for derived targets, with a full-range scan for disconnected sets,
- a parametric bootstrap with atomic, digest-checked, identity-bound records and a manifest written last,
- the required R cross-checks.

**Defects found in the engine itself while testing it (both fixed before any result was used):**

1. **The mode search cycled.** The first draft used a step-clipped Newton for the quadrature mode. On a real arm (398 patients, 269 events, η0 = −2, spread 0.96) it oscillated between z = 0 and z = 5 and returned an integral 114 log units too low. The ladder test caught it: 10 and 20 nodes agreed while 80 differed by 1.3e−4. The fix is a bracketed Newton with bisection, which converges because every integrand is log-concave (Prékopa's theorem for the outer one). It is pinned by `test_T3_inner_mode_search_does_not_cycle` and `test_T3_full_early_likelihood_matches_a_grid_reference_at_large_spreads`.
2. **Disconnected sets could be missed.** The profile walk doubled its step once below the cutoff, so it could step over a second region above the cutoff. A 25-point scan over each target's whole range is merged with the walk. It is pinned by `test_T8_a_bimodal_profile_gives_a_disconnected_set` (red with the scan removed).

Two reference defects in the tests were found and fixed in the same way. Nested `scipy.quad` overflowed at extreme outer values, so it was replaced by a dense double-grid reference. One test also added log C(N, y) twice.

**Tests: `tests/python/test_frequentist.py`, 37 tests, all green.** Each was made red on purpose. The mutation record (`2026_11_12_frequentist_review/mutation_checks.json`) breaks the engine 13 ways, and every one turns a test red:

| Mutation | Caught by |
|---|---|
| cloglog orientation swapped | T2/T3 |
| binomial constant dropped | T2 |
| mode search removed | T3 |
| tumor contrast without Q/(α/β) | T4 |
| early contrast without dQ/(α/β) | T4 |
| sign of the Tk gradient | T5 |
| τ gradient dropped | T5 |
| wrong profile cutoff | T8 |
| full-range scan removed | T8 |
| record digest not checked | T8 restart |
| run identity without the design hash | T1 |
| smax sharpness 2.9 | T7/T4 |
| R check on the logit link | T6 |

The first M09 run selected no test because its filter was wrong, and pytest's "no tests" exit looked red. It was rerun with the correct filter and is genuinely red.

**R cross-check (T6, Windows R 4.5.1, lme4 1.1.37).**

- **Tumor.** α, λ and g agree with `glmer(binomial("cloglog"), nAGQ = 25)` within 1e−4 relative. The Python log-likelihood at R's coefficients is within 5e−10 of the Python optimum. Both put σ at 0.
- **Early GLM.** On synthetic counts, the coefficients agree to about 1e−10. The profile interval of κ, (0.0396084, 0.0703220), matches R's profile interval of the dD coefficient, (0.0396084, 0.0703219).
- **Laplace fit (descriptive only).** The two-random-effect `glmer` estimates a negative dQ coefficient (−0.0094). That independently corroborates the open upper end of the mucosal α/β.

**Pilot (B = 20 per model).** Each pilot was run uninterrupted, and again killed by SIGKILL and resumed:

- tumor: killed after 6 records, then 14 run on resume,
- early: killed after 7, then 13 run.

The accepted records of the interrupted and uninterrupted runs are identical: same simulated-count hashes, same status, targets within 1e−9. No manifest was written at the kill.

Measured cost:

- tumor: 3.7 CPU-seconds per replicate, 106 MB,
- early: 31 CPU-seconds per replicate, 144 MB,
- one tumor profile: about 35–140 s,
- one early profile: about 200 s.

The full run was estimated at about 6.5 CPU-hours, within the D3 budget.

## Results

See `2026_11_12_frequentist_comparison_results.md`.
