# Stages B to E: corrected trial data, fit runner, refit and app (version 3.1.0)

Plan: `docs/development/PLAN_codex_review_2026_10_25.md` (approved by Igor Shuryak on 2026-10-25, all 12 recommendations). Stage A is recorded in `2026_10_25_stage_A_app_fixes.md`. Base revision 677ce21. All work below is uncommitted on that revision until Igor asks for a commit.

## Outcome

- The app (version 3.1.0) uses the verified-trial refit: DI posterior draws and DD calibrated to DI, with DD row i fitted on DI row i. D18 stays the default.
- Every fit passes the release gate of decision 9. DI: 0 divergences, R-hat ≤ 1.006, bulk and tail ESS ≥ 907, E-BFMI ≥ 0.3, in all 12 fits. Calibrated DD: 3600 of 3600 conditional fits with 0 divergences, 0 tree-depth hits, R-hat ≤ 1.01, E-BFMI ≥ 0.3, after 35 were re-fitted with longer chains (below).
- Python tests: 38 passed. R suite (`Rscript tools/run_tests.R`, R 4.5.1, release library): 5265 expectations, 0 failed, 0 errors, 0 warnings, on 2026-10-26.

## Stage B: data corrections (decisions 4, 6, 7)

Canonical builder: `fit/build_refit_dataset.py` (moved into the repository. The old builder in `References_refit_trials/search/extraction/` is unchanged as the record of the 10-22 to 10-25 files). Main file `data/refit_trials_2026_10_25c_main.csv`: 45 arms, 19 trials, 21 fitting units, 8,279 patients, d 1.15–2.55 Gy, sha256 `51ed5d8a…8409`.

| Item | Change | Source |
|---|---|---|
| KROG 0201 HYPO arm | Mixture (65 × 28 + 9 × 30)/74 = 28.243 fractions. 9 T2 patients, not 8 | Moon 2014 Table 1 (T2a 7 + T2b 2) |
| RTOG 9003 | Rates and N from the same report: N 268/263/274/268, 2-year LRC 46.0/54.4/47.5/54.5% | Fu 2000 (PMID 10924966), Tables 1 and 2. Beitler 2014 rates are the sensitivity set `--rtog9003-beitler` (its per-arm N is not reported, so it borrows the 2000 N) |
| Mixed fraction sizes | Kept as components. BED_LQ = Σd + Σd²/(α/β) | ARTSCAN AF 23 × 2.0 + 20 × 1.1. Ghoshal 25 × 1.8 + 15 × 1.5. RTOG 88-09 30 × 1.8 + 11 × 1.5. RTOG 9003 CB 30 × 1.8 + 12 × 1.5. Vienna V-CHART 1 × 2.5 + 32 × 1.65 (each checked in the full text) |
| Carried fields | Estimator (39 KM, 6 cumulative incidence), endpoint definition, exact source percentage, source and nominal overall time | Extraction |
| Refusals | Builder exits on an N that differs from the extraction without a stated reason, on components that do not sum, on a missing nominal schedule | Tests below |

The diff of the new main file against the 10-25 main file is only the KROG HYPO row and the RTOG 9003 rows. ARTSCAN stays out of the main set (opt-in `--artscan`).

## Stage C: fit runner

- `fit/stan/{di,dd,dd_cal}.stan`: data `dsum`/`d2sum` in place of `d`/`nfx`. Intercept prior sd as data. Generated `delta_eval[R]` (each reference regimen's effective BED minus the anchor's) and, in `dd_cal`, `onset_eval[R]`.
- `fit/fit_refit.py`: run identity (sha256 of data, Stan sources, runner, priors, sampler settings, seed, software versions). A finished fit is reused only when its identity matches, else the run refuses. Full diagnostics over every sampled quantity including `c[]`, `z[]` and `delta_eval[]` (82 for main DI). Leave-one-study-out grouped by trial.
- `fit/fit_dd_calibrated.py`: 4 chains × (600 warmup, 250 samples) per DI draw, retried up to 3 times with doubled length (decision 10 asked for 1000 outer draws. The inner length is a deviation, below). Folder identity binds the parent DI identity, its draws hash and the picks.
- `fit/export_draws.py`: writes the app's draws only when every fit passes the gate. Manifest `draws_version` 3.1.

## Stage D: fits

Orchestrator `fit/run_refit_2026_10_25c.sh`, log `fit/run_refit_2026_10_25c.log`. Ran 2026-10-25 12:57 to 2026-10-26 02:10 (13.2 h) on 4 cores, serial, while another job used about 16 of the 22 cores. Main set: 7.3 h (four DI/DD variants, leave-one-study-out, calibrated DD 1000 draws for two variants). Each 200-draw sensitivity set: 34–50 min.

### Main results (primary variant)

| Quantity | Median (95% interval) |
|---|---|
| α (1/Gy) | 0.044 (0.030–0.060) |
| α/β (Gy) | 10.6 (6.9–18.9). Close to the prior |
| DI λ/α (Gy per day at each draw's α/β) | 0.63 (0.21–1.19) |
| DI onset (inclusive day) | 37 (24–43) |
| g/α (Gy per day) | 0.28 (0.11–0.48) |
| Overdispersion σ | 0.065 (0.004–0.17) |
| Calibrated DD λ (1/day) | 0.40 (0.12–1.03). Prior after the onset preference 0.65 (0.11–10.4) |
| Arms inside their 90% predictive interval | 45 of 45 |
| DI vs DD order agreement, 55 regimen pairs | 0.921 (outer-draw SE 0.001) |

HYPNO 55 Gy/20 (26 d) against 35 × 2 Gy (47 d): DI −2.0 (80% −4.3 to +1.0), DD −3.3 (−5.9 to −0.6), same sign in 87.5% of paired draws.

Free DD fits (record only, not used by the app) fail the gate: 74, 64, 35 and 28 divergences (primary, binomial, wide α/β, paper restrictions). Leave-one-study-out over 24 predicted arms, conditional plug-in score: DI −99.12 (0 divergences), free DD −102.68 (415 divergences across folds). Not a firm model comparison.

The DD onset preference is soft: 0.5% of pooled DD draws (13,506 of 2,629,000), in 983 of 1000 outer fits, have some fitted arm with onset below day 7 (minimum 5.1 days). No evaluated reference regimen has a DD onset below day 7.

### Sensitivity sets (primary variant, calibrated DD 200 draws each)

| Set | Arms | HYPNO DI median (80%) | HYPNO DD median (80%) | Order agreement | α/β median |
|---|---|---|---|---|---|
| Main | 45 | −2.0 (−4.3, +1.0) | −3.3 (−5.9, −0.6) | 0.921 | 10.6 |
| Effective N halved (decision 5) | 45 | −1.9 (−4.5, +1.5) | −3.1 (−5.8, −0.3) | 0.916 | 10.5 |
| Nominal clock (decision 6) | 45 | +0.3 (−2.4, +3.7) | −2.1 (−5.1, +0.9) | 0.883 | 10.7 |
| KM-only units | 39 | −1.9 (−4.3, +1.2) | −3.4 (−6.1, −0.5) | 0.915 | 10.1 |
| RTOG 9003 Beitler 2014 rates | 45 | −1.9 (−4.3, +1.1) | −3.3 (−5.8, −0.4) | 0.917 | 10.6 |
| + ARTSCAN | 47 | −2.7 (−5.0, +0.4) | −3.7 (−6.0, −0.8) | 0.926 | 10.6 |
| + three trials added 2026-10-24 | 51 | −1.7 (−3.8, +1.0) | −2.9 (−5.3, −0.5) | 0.918 | 10.8 |
| + same-chemotherapy trials | 51 | −1.0 (−3.5, +2.2) | −2.4 (−5.2, +0.2) | 0.913 | 10.2 |
| Intercept prior sd 10 | 45 | −1.9 (−4.3, +1.1) | −3.3 (−6.1, −0.7) | 0.915 | 10.6 |

All DI fits in these sets pass the gate (0 divergences, R-hat ≤ 1.006, bulk ESS ≥ 907). The nominal clock moves the DI result most: with the planned rather than the delivered time, HYPNO's DI change is near zero.

A KM-only set that also drops cumulative-incidence units leaves 39 arms. A CIF-only set is not feasible (3 units), so it was not run (deviation from the plan's list).

### Calibrated-DD fits that needed longer chains

35 of 3600 conditional fits ended attempt 3 (4 × 1000 samples) with R-hat 1.010–1.017, worst on `sigma_raw`, `z[]` or `c[]`, with no divergence. Main primary 9, main paper restrictions 12, effective N halved 4, ARTSCAN 4, chemotherapy 2, RTOG 9003 Beitler 2, KM-only 1, prior sd 10 1, nominal clock 0, later trials 0.

The gate was not changed. `fit/repair_dd_calibrated.py` continues the runner's attempt sequence for those fits only (attempt 4 at 8× length, attempt 5 at 16×, same seed rule, data, model and gate). It refuses a folder whose identity does not match the current data, Stan model, runner or parent DI fit, keeps each replaced file in `repaired_originals/`, lists every repair in `repair.json` and recomputes `summary.json`. The runner itself is unchanged, so the identity of the finished folders stays true. All 35 passed at attempt 4 (R-hat 1.0013–1.0059, 0 divergences), in 16 minutes (`fit/run_repair_2026_10_25c.sh`, log `fit/run_repair_2026_10_25c.log`).

The run's pooled summaries weight a retried fit by its larger number of rows. Recomputed with 250 rows per outer draw: primary λ_DD median 0.404 (0.122–1.033) vs pooled 0.401, order agreement 0.9185 vs 0.9205, HYPNO DD −3.28 vs −3.25. The exported app draws use 4 rows per outer draw and are not affected.

## Stage E: app

- Draws: `draws/tumor_{di,dd}_{primary,paper_restrictions}.csv`, 4000 rows each (1000 distinct DI draws × 4 paired DD draws), from `fit/output_refit_2026_10_25c_main`. The manifest records the dataset and its sha256, per-file diagnostics with `gate_pass`, the pairing, the binomial-approximation note and `outer_fits_repaired` (9 and 12). `load_draws()` refuses a file without `gate_pass` TRUE.
- Refit preset (`draws/refit_point_values.json`): DI medians of identified quantities with α/β 10 Gy (α 0.0443, λ 0.0280, Tk 37.3, g 0.0123). DD λ 0.408 calibrated at the anchor (C/α 68.4 Gy, onset 38.3 days).
- Units: "Gy₁₀" only when every draw has α/β 10 Gy, otherwise "Gy, α/β varies by draw". MCSE is described as approximate over distinct DI draws.
- Labels: "DD (calibrated to DI)". The DD mode where the user types DD values reads "DD values entered directly (not calibrated)" in the banner.
- About page: refit section rewritten with the numbers above. Organ wording matches stage A (maximum doses converted, mean doses physical with a flag when the fractionation differs, uploads saved as a temporary file and deleted at once). Data range 1.15–2.55 Gy.
- APP_VERSION and DESCRIPTION 3.1.0.

## Tests and red/green evidence

| Test | Red without the change | Green |
|---|---|---|
| `tests/python/test_refit_data_and_identity.py` (KROG, RTOG 9003, components, carried fields, calendar, equal-split invariance, builder refusals, `load` components, identity refused/reused/sensitive) | Dataset tests fail on the 10-25 file (`BED_TEST_DATA`). Runner tests fail on the HEAD `fit_refit` ("DID NOT RAISE", KeyError d2sum). Builder checks mutated: red | Yes |
| `tests/python/test_export_draws.py` (pairing, gate refusal, DD `gate_pass` present) | Export gate mutated: red. `gate_pass` lines removed: KeyError | Yes |
| `tests/python/test_repair_dd_calibrated.py` (only failing fits re-fitted, recorded, originals kept, idempotent. Identity refusal for data, runner, Stan, parent) | Selection mutated to re-fit every draw: red. Identity check disabled: 4 of 4 red | Yes |
| `tests/testthat/test_v31_stage_e.R` (manifest 3.1 and dataset sha, gate, DI/DD row pairing with anchor equality < 1e-6, paper restrictions α/β 10, gate refusal, data range, units, point values, superseded text) | Red on the 10-21 draws (manifest 3.0, no pairing). Red on the first 3.1 export because DD diagnostics lacked `gate_pass` | Yes |
| `test_v3_app.R::fitted DD` banner | Old banner text restored: red | Yes |
| `test_v3_uncertainty.R` golden values | Repinned to −2.28476945801911, 1.44023467546532 (was 1.749, 5.260 on the 10-21 draws). R 4.5.1 ucrt, tolerance 1e-10 (exact arithmetic on fixed draws) | Yes |
| `test_v3_uncertainty.R` About numbers, ESS attribute. `test_v27`, `test_v25`, `test_v26` versions | Red on the new draws and 3.1.0 until the About page and DESCRIPTION were updated | Yes |

## Deviations from the plan

1. Inner DD fits use 4 × 250 samples (600 warmup) with up to 3 doublings, plus the repair step's attempts 4–5, instead of a single fixed length. Inner fits are gated on divergences, tree depth, R-hat and E-BFMI. ESS is not part of the inner gate (each inner fit contributes only 4 rows to the app).
2. The builder moved into the repository (old builder kept unchanged).
3. No CIF-only sensitivity set (3 units).
4. The repair step was added after the run. The gate is unchanged.

## Open items

- Igor to confirm stage A decision 3 (mean-dose limits on physical dose when the fractionation differs from the source's. Flagged per row).
- Stage F, done 2026-10-26 (Igor approved all three decisions as recommended, including decision 12 on physical mean-dose limits):
  - S6: `check_bundle()` refuses any bundled file not tracked by git (`deployed_revision.txt` excepted, written by the deploy). On its first real use it found `.pytest_cache/` in the upload list, now in `.rscignore`.
  - Full suite on `e4ff22f`: 5269 expectations, 0 failed, 0 errors; pytest 38 passed.
  - Local walkthrough 34 of 34 (first run 31 of 34: three checks still held 3.0 strings, kept as `walkthrough/20261026_071429_local_stale_3_0_checks_fail`).
  - Bundle: 40 files (app, R, data, draws, DESCRIPTION, LICENSE, renv.lock).
  - Deployed `e4ff22f` to `custom_bed_calculator_v3` only. Live walkthrough 34 of 34, including an exact replay of the uncertainty download. v1 and v2 respond unchanged. Rollback: redeploy `f951cd3`.
  - Not pushed to GitHub (no decision taken on push).
