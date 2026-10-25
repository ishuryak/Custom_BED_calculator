# Calibrated DD fit and DI–DD regimen agreement (2026-10-24)

Nothing in the app changed. This record extends `2026_10_22_refit_verified_trials.md`.

## Why

The free DD fit is not identified: its repopulation rate and onset constant reproduce their priors, so DD and DI cannot be compared as competing fits. Igor Shuryak (2026-10-24) set a different purpose. The question is not which model fits better. It is where the two models agree or disagree about which regimens to recommend. If they mostly agree, the recommendations do not depend on the details of repopulation modeling.

## Decisions (Igor Shuryak, 2026-10-24)

- **Calibrate DD to DI at an anchor regimen, as the app does.** The onset constant is solved so that DD repopulation equals DI repopulation at 35 × 2 Gy in 47 days, which removes the free onset constant C. Anchoring at HYPNO was rejected, because it would make DD equal DI at HYPNO by construction.
- **Keep repopulation onset at or after day 7**, the normal mucosal repopulation time. This is applied as a smooth exponential penalty, −Σ exp((7 − onsetᵢ)/1 day) over the fitted arms, rather than a hard bound.

## Method

`fit/stan/dd_cal.stan` and `fit/fit_dd_calibrated.py`.

For each of 200 evenly spaced draws of the DI primary refit (`fit/output_refit_2026_10_22/di_primary`), the fit takes α, α/β and g from the draw and solves the onset at the anchor. It then fits λ_DD, the trial intercepts and the arm-level overdispersion: 4 × 250 draws per DI draw, 200,000 pooled draws in all.

- **Prior on λ_DD:** the same as in the free DD fit, lognormal(log 0.5, 1.5).
- **Same draws for both models:** every regimen's tumor-dose change is computed under DI and calibrated DD on identical pooled draws, so the two are directly comparable.
- **Data reuse:** DD borrows α, α/β and g from DI, which was fitted to the same trials. This is the app's calibration, not an independent DD fit.

## Results

All results are from `fit/output_refit_2026_10_22/dd_calibrated_onset7/summary.json`.

**Fit diagnostics:**
- 0 divergent transitions, max R-hat for λ_DD 1.013.
- λ_DD has a 95% interval of 0.11–1.11, against a prior 95% interval of 0.03–9.5.
- The onset at the anchor has a median of 38.8 days (80% interval 29.5–43.6).
- The earliest fitted arm onset has a median of 13.4 days (80% interval 10.2–15.1).

Tumor-dose change against 35 × 2 Gy in 47 days, in Gy₁₀, median (80% interval):

| Regimen | DI | Calibrated DD | Same sign |
|---|---|---|---|
| HYHOPE-type 46.5 Gy/15 (19 d) | −9.8 (−12.3 to −6.2) | −10.1 (−12.8 to −6.4) | 0.99 |
| 30 × 2 Gy (40 d), last week dropped | −5.8 (−7.6 to −3.7) | −6.0 (−8.0 to −3.5) | 0.98 |
| 17 × 3.05 Gy (23 d) | −4.2 (−7.0 to −0.5) | −5.0 (−7.8 to −1.4) | 0.97 |
| 66 Gy/33 (45 d) | −2.9 (−3.5 to −2.0) | −2.4 (−3.4 to −0.9) | 0.96 |
| HYPNO 55 Gy/20 (26 d) | −2.6 (−5.1 to +0.5) | −3.6 (−6.3 to −0.7) | 0.91 |
| 60 Gy/25 (33 d) | −0.3 (−2.3 to +2.3) | −2.2 (−4.3 to +0.2) | 0.67 |
| 65 Gy/30 (40 d) | +1.3 (−0.5 to +3.3) | −1.0 (−2.2 to +0.5) | 0.33 |
| 33 × 2 Gy, 6 per week (38 d) | +2.7 (+0.7 to +4.8) | −0.4 (−1.6 to +1.2) | 0.40 |
| 70 Gy/35, 6 per week (40 d) | +6.2 (+4.4 to +8.3) | +1.9 (+1.1 to +2.9) | 1.00 |
| 68 × 1.2 Gy twice daily (47 d) | +7.5 (+5.9 to +9.0) | +4.5 (+1.9 to +6.7) | 0.97 |

The two models give the same ordering for 92.6% of regimen pairs (the anchor is included as a regimen), averaged over draws.

**Where they agree.** Hypofractionation to 15–20 fractions, dropping the last week, acceleration and twice-daily treatment all point the same way under both models.

**Where they disagree.** The models differ on moderately accelerated 2 Gy courses of 38–40 days. DI's onset (Tk 38.6 d) falls after most of those courses end, so DI credits them with the full gain from shortening. DD's onset falls earlier as BED per day rises, to about day 25 for these schedules, so its gain is smaller. The differences are 2–4 Gy₁₀.

## Checks

**Penalty effect.** A run without the penalty (`dd_calibrated_onset-10000`, 40 DI draws) lets the lower tail of the earliest fitted onset reach 4.4 days (10th percentile). It moves no regimen median by more than 0.5 Gy₁₀ relative to the penalized run on the same 40 draws.

**Number of DI draws.** A first run on 40 DI draws gave a DI HYPNO median of −3.3, against −2.48 over all 8000 draws. That gap is about 2 Monte Carlo standard errors, so the main run uses 200 draws, which give −2.6.

**Tests.** `tests/python/test_dd_calibrated.py` has 8 tests:
- DD repopulation equals DI at the anchor.
- DD differs from DI at HYPNO.
- The onset and C/α match the app's R `dd_calibration` to 1e-9.
- `inv_smax` inverts `smax`.

When I removed the kill-probability factor from the calibration on purpose, 5 tests failed. The runner also asserts in every draw that DD repopulation at the anchor equals R0.

## Limits

- The trial data cover 1.15–2.55 Gy per fraction, so HYPNO, 17 × 3.05 and HYHOPE are extrapolations under both models.
- Arms are not adjusted for stage mix.

## Files

- **Committed:** `summary.json` and the run logs.
- **Ignored as regenerable:** the per-DI-draw CSVs (`di_draw_*.csv`) and `draws_params.csv`, following the existing `fit/output` convention.
- **Unused draw files:** the onset-7 folder also holds 38 per-draw files from the 40-draw run that the 200-draw summary does not use.
