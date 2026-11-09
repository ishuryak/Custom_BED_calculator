# Jackson 1997 accelerated arm: dose-encoding check at the fitted parameters (2026-11-09)

Review of 2026-11-09 (`HNSCC_BED_v44_independent_review_and_action_plan_2026-11-09.md`), item 2.1, step 3, first half only.

## What was found

`fit/build_refit_dataset.py` gives all 41 twice-daily patients of Jackson et al. 1997 (PMID 9165135) 66 Gy in 33 × 2 Gy. The trial paper says otherwise (p.41, Methods, and p.43, Discussion): "The last ten patients in the accelerated arm received a total dose of 62.8 Gy because the second of the daily treatments was reduced to 1.8 Gy." That is 17 × 2 + 16 × 1.8 Gy, Σd² = 119.84 Gy². The paper reports one pooled 3-year actuarial locoregional recurrence-free survival per arm (49.1% twice daily, 44.3% once daily, p.39 and Fig. 1) and mean delivered overall times of 24.25 and 47.6 days (Table 3, p.41). Checked in `../References_refit_trials/fulltext_text/9165135.txt` and the PDF on 2026-11-09.

The fitted-trial list note now states this (commit "Fitted-trial list: correct CAIR, Jackson and PMH schedule notes"). The fitting data and the draws are unchanged.

## What was computed

`fit/jackson_exposure_check.py` (no Stan run) evaluates the accelerated arm under three encodings at the DI parameters, with the once-daily arm unchanged:

- (a) uniform: 66 Gy, Σd² 132 Gy² (as fitted).
- (b) weighted exposure: 31/41 and 10/41 prescription weights, 65.219512 Gy and 129.034146 Gy² (a mean exposure inside a nonlinear model, so an approximation).
- (c) prescription mixture: 31/41 × TCP(66 Gy) + 10/41 × TCP(62.8 Gy) on the probability scale.

It uses the fit's own link, log(−log TCP) = c − α(B − B̄) + λ s(T − Tₖ) + g(T − T̄) (`fit/stan/di.stan`, `fit/fit_refit.py` lines 98–112). The unit intercepts are not exported with the draws, so the shared Jackson baseline is set by holding the once-daily arm at its fitted count rate, 18/41 = 43.9%, with arm deviations at zero. The accelerated arm then follows from the within-unit difference, in which the centering cancels.

At α/β = 10 Gy the weighted encoding lowers the arm's LQ BED from 79.200 to 78.122927 Gy (the review's numbers reproduce).

## Results

| Encoding | Point values: accelerated TCP | Accelerated minus once daily | 1,000 DI draws: accelerated TCP, median (95% interval) |
|---|---|---|---|
| (a) uniform 66 Gy | 62.97% | +19.07 points | 62.85% (57.14–69.36) |
| (b) weighted exposure | 61.56% | +17.66 points | 61.46% (56.13–67.74) |
| (c) prescription mixture | 61.51% | +17.60 points | 61.42% (56.09–67.68) |

Per draw, (b) minus (a) is −1.38 points (95% interval −1.82 to −0.94) and (c) minus (a) is −1.43 points (−1.93 to −0.96). On the log(−log) scale the point-value contrast moves from −0.5765 to −0.5287 (about 8% smaller). The weighted-exposure approximation and the mixture agree to within 0.1 point.

Reading: at fixed parameters the amendment lowers the predicted accelerated-arm control by about 1.4 percentage points and shrinks the modeled Jackson contrast from about +19 to about +18 points. The observed contrast is +4.8 points (49.1 − 44.3), so the corrected encoding moves the model slightly toward the data. The change is small next to the arm's binomial uncertainty at N = 41 (one patient is 2.4 points).

## What was not done, and the decision left open

- **No targeted DI refit was run.** Whether the amended encoding shifts the posterior (α, λ, Tₖ, g), the HYPNO contrast, threshold or support crossings, or model agreement is not known from this check. Running it (review item 2.1, step 3, second half) is a separate decision for Igor, with the dataset outside Jackson, the priors and the gate held fixed.
- Delivered individual doses and outcomes by prescription period are not reported, so the 31/41 and 10/41 weights are illustrative.
- Promoting a revised input (review item 2.1, step 4) would version the dataset and require refitting primary DI and conditional DD and regenerating every downstream artifact. Not done.

## Command

    python3 fit/jackson_exposure_check.py
