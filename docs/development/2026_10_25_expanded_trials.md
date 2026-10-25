# Trial-set review: main fit, sensitivity sets and source checks (2026-10-25)

Nothing in the app changed. This record extends `2026_10_22_refit_verified_trials.md` and `2026_10_24_dd_calibrated_agreement.md`.

## Summary

- **The new main fit** is the 10-22 trial set with ARTSCAN removed: 45 arms in 21 units.
- **Sensitivity sets:**
  - three newly found trials (Vienna V-CHART, ELAN-RT, Gupta 2015).
  - ARTSCAN.
  - trials comparing fractionation under the same chemotherapy.
- **Main-fit HYPNO result:** the tumor-dose change of HYPNO against 35 × 2 Gy in 47 days is −2.0 Gy₁₀ under DI and −3.3 under calibrated DD.
- **Agreement:** the two models give the same pairwise regimen ordering 91.6% of the time.
- **Source checks:** four values flagged in the reference check were each confirmed against their primary PDF, so no value changed.

## Decisions (Igor Shuryak, 2026-10-25)

1. **The three newly found trials go to sensitivity, not the main fit.** Their endpoints are weaker:
   - Vienna and Gupta report crude control.
   - ELAN-RT reports a competing-risk cumulative incidence in frail elderly patients.
2. **ARTSCAN goes to sensitivity.** About 37% of its patients had preoperative radiotherapy followed by surgery.
3. **Flagged values are corrected only from the actual PDFs.** This check found nothing to correct (see below).

## Source checks of flagged values

Each was read from the primary full text (`References_refit_trials/fulltext_text/` and the PDFs).

| Trial | Flag | What the PDF shows | Action |
|---|---|---|---|
| EORTC 22791 | 38/56 against 40/59 in the fit | The 1992 final analysis (PMID 1480768) gives 5-year control of 40% (CF) vs 59% (HF), 159/166 analyzed, median overall time 51 days in both arms (Table IX). The 38/56 is from the 1990 preliminary analysis (PMID 2145894). The flag wrongly cited it as Horiot 1996. | None. The fit uses the final analysis. |
| TROG 91.01 | Conventional arm 47 or 49 days against 50 in the fit | The primary paper (PMID 11439206) table gives a delivered median of 50 days (range 9–77) for CRT and 24 for ART. 49 days is the planned time. The 47 is from a toxicity companion (Denham 1999). | None. The rule uses the delivered median. |
| CHART | Randomized 552/366 against 550/365 in the fit | In Dische 2010 (PMID 20394851), 550/365 are the numbers at risk at time 0 on the locoregional relapse-free curve, the endpoint used. 552/366 are on the overall-survival curve. | None. Every trial uses the count analyzed for its endpoint. |
| RTOG 83-13 | 5–14% nasopharyngeal | Cox 1990 (PMID 2180866), Table 2 shows nasopharynx in 8%, 10% and 10% of the 67.2, 72.0 and 76.8 Gy arms. Analyzed counts are 63/58/116, and 2-year control is 25/37/42% (Table 3). Control is reported only for each arm as a whole, not by site. | None possible. Kept in the main fit. |

## Newly found trials (reference check, 2026-10-24)

The downloaded references were searched again (three agents plus my own quote check against each PDF). A trial is used only when its arm-level values are stated in a downloaded full text.

| Trial | Arms | Value source | Caveat |
|---|---|---|---|
| Vienna V-CHART (PMID 11054514) | CF 2 Gy × 35, 49 d, n 81, control 31%. V-CHART 55.3 Gy/33, 17 d, n 78, control 32% | Crude locoregional control | Crude |
| ELAN-RT | HSC-RT 2.75 Gy × 20, 43 d, n 101, failure 27.3%. SF-RT 2 Gy × 35, 52 d, n 100, failure 24.7% | 3-year cumulative incidence, ITT, Table 2 p. 8 | Frail elderly patients, split course, competing risk |
| Gupta 2015 | Arm A 2 Gy × 33, 6 per week, 38.5 d, n 55, 52.7%. Arm B 2 Gy × 33, 45.5 d, n 54, 44.4% | Crude control | Crude, median follow-up 43 months |

**Same-chemotherapy trials, for the chemo sensitivity set:** HYPNO, ART DECO and Shivananjappa.

**Excluded, with reasons:**
- JCOG0701: no overall time in the paper held.
- Kachhwaha: DFS = OS, no local-control endpoint.
- Lyhne 2015: same patients as DAHANCA 6 and 7.
- ORO 93-01: DFS only.
- GORTEC 99-02: only one radiotherapy-alone arm.
- CHA-CHA: endpoint is not tumor control.
- Bjordal: quality of life only.
- Pinto: crude control not stated.
- Poddar: under 2 years of follow-up and per-protocol.

## Data and runs

- **Builder:** `References_refit_trials/search/extraction/build_refit_dataset.py`, outside the repository.
  - Every tumor-control value is asserted against a quote-checked extraction row.
  - The 10-22 and 10-24 files rebuild byte-identical.
- **Datasets:**

  | File | Builder flags | Arms / units |
  |---|---|---|
  | `data/refit_trials_2026_10_25_main.csv` | `--drop-artscan` | 45 / 21 (the 10-22 file minus the two ARTSCAN rows, checked by merge) |
  | `data/refit_trials_2026_10_25_chemo.csv` | `--drop-artscan --chemo` | 51 / 24 |

- **Earlier fits reused as sensitivity sets:**
  - main + ARTSCAN = the 10-22 fit.
  - main + new trials = `output_refit_2026_10_24_drop_artscan`.
  - main + both = `output_refit_2026_10_24`.
  - 10-24 + ARTSCAN without ELAN = `output_refit_2026_10_24_drop_elan`.
- **Runs:** `fit/fit_refit.py` (all variants, plus leave-one-trial-out on the main set), then `fit/fit_dd_calibrated.py` on 200 DI draws with the day-7 onset penalty.
  - Code revision: 6618fec.
  - Outputs: `fit/output_refit_2026_10_25_{main,chemo}/`.
- **Diagnostics:**
  - every DI primary fit and every calibrated-DD fit had 0 divergences.
  - max R-hat was 1.003 (DI) and 1.020 (λ_DD).
  - in every set, every arm falls inside its 90% posterior predictive interval.

## DI parameters

Median [95% interval]:

| Data | α (Gy⁻¹) | α/β (Gy) | λ (day⁻¹) | Tk (days) |
|---|---|---|---|---|
| **Main (10-25)** | 0.043 [0.029, 0.060] | 10.7 [6.6, 19.2] | 0.029 [0.011, 0.052] | 37.4 [23.9, 43.0] |
| Main + ARTSCAN (10-22) | 0.042 [0.028, 0.059] | 10.3 [6.2, 19.1] | 0.028 [0.008, 0.054] | 38.6 [24.2, 44.3] |
| Main + new trials | 0.041 [0.028, 0.056] | 10.7 [6.7, 19.1] | 0.029 [0.012, 0.049] | 37.7 [26.1, 43.0] |
| Main + both | 0.040 [0.027, 0.054] | 10.4 [6.4, 18.9] | 0.029 [0.011, 0.051] | 38.8 [25.0, 44.1] |
| Main + same-chemotherapy trials | 0.040 [0.026, 0.057] | 10.2 [6.3, 18.5] | 0.027 [0.010, 0.049] | 36.7 [22.6, 43.5] |

**Leave one trial out, main set:**
- Summed elpd over 21 units: −98.4 for DI against −102.2 for the free DD fit.
- The free DD leave-out fits had 228 divergent transitions, consistent with DD not being identified. DI had 0.

**External prediction of the new trials.** The 10-22 fit, which never saw them, was given each trial's first arm and asked for its second. The intervals include parameter uncertainty and overdispersion, not binomial sampling.

| Trial | Given arm | Predicted other arm (80% interval) | Observed |
|---|---|---|---|
| Vienna | CF 30.9% | V-CHART 25.6% (18.4–33.5) | 32.1% |
| ELAN-RT | HSC-RT 72.3% | SF-RT 77.2% (72.5–81.4) | 75.0% |
| Gupta 2015 | Arm A 52.7% | Arm B 44.3% (38.2–50.1) | 44.4% |

## Regimen results, main fit

Tumor-dose change against 35 × 2 Gy in 47 days, Gy₁₀, median (80% interval), from `fit/output_refit_2026_10_25_main/dd_calibrated_onset7/summary.json`:

| Regimen | DI | Calibrated DD | Same sign |
|---|---|---|---|
| HYHOPE-type 46.5 Gy/15 (19 d) | −9.4 (−11.9 to −6.2) | −9.7 (−12.3 to −6.5) | 1.00 |
| 30 × 2 Gy (40 d), last week dropped | −5.5 (−7.3 to −3.5) | −5.5 (−7.6 to −2.9) | 0.99 |
| 17 × 3.05 Gy (23 d) | −3.8 (−6.3 to −0.1) | −4.6 (−7.5 to −1.6) | 0.94 |
| 66 Gy/33 (45 d) | −2.8 (−3.5 to −2.2) | −2.3 (−3.2 to −0.7) | 0.95 |
| HYPNO 55 Gy/20 (26 d) | −2.0 (−4.4 to +1.0) | −3.3 (−6.0 to −0.5) | 0.88 |
| 60 Gy/25 (33 d) | +0.3 (−2.0 to +2.7) | −2.0 (−4.1 to +0.5) | 0.54 |
| 65 Gy/30 (40 d) | +1.4 (−0.2 to +3.6) | −1.0 (−2.1 to +0.6) | 0.30 |
| 33 × 2 Gy, 6 per week (38 d) | +3.0 (+1.0 to +4.9) | −0.3 (−1.6 to +1.3) | 0.39 |
| 70 Gy/35, 6 per week (40 d) | +6.3 (+4.7 to +8.5) | +1.8 (+1.1 to +2.8) | 1.00 |
| 68 × 1.2 Gy twice daily (47 d) | +7.7 (+6.3 to +8.8) | +4.4 (+1.7 to +6.4) | 0.97 |

- **λ_DD:** 95% interval 0.13–1.09.
- **Onset:** the median earliest fitted onset is 13.2 days. No regimen has a DD onset before day 7 in any draw.
- **Ordering agreement:** 91.6%.

**HYPNO across all sets:**

| Data | DI | Calibrated DD | Same sign | Ordering agreement |
|---|---|---|---|---|
| **Main** | −2.0 (−4.4 to +1.0) | −3.3 (−6.0 to −0.5) | 0.88 | 91.6% |
| Main + ARTSCAN | −2.6 (−5.1 to +0.5) | −3.6 (−6.3 to −0.7) | 0.91 | 92.6% |
| Main + new trials | −1.1 (−3.3 to +2.0) | −2.4 (−5.0 to +0.3) | 0.82 | 90.9% |
| Main + both | −1.9 (−4.2 to +1.2) | −2.9 (−5.7 to −0.3) | 0.89 | 92.2% |
| Main + same-chemotherapy trials | −0.9 (−3.1 to +2.4) | −2.3 (−5.0 to +0.8) | 0.80 | 90.2% |

## Interpretation

- **The parameters hardly depend on which trials are included.** α/β stays at 10.2–10.7 Gy and Tk at 36.7–38.8 days.
- **The new trials agree with the existing fit.** All six new arms are inside their predictive intervals, both before they are added and after.
- **HYPNO's deficit is consistent in sign but not well determined.**
  - Under DI, the median ranges from −0.9 to −2.6 Gy₁₀ across sets, and every 80% interval includes zero.
  - ARTSCAN pulls HYPNO down (−2.0 to −2.6). The new trials pull it up (−2.0 to −1.1).
  - Calibrated DD puts HYPNO 1–1.4 Gy₁₀ lower than DI in every set.
- **Recommendations agree on most regimens.** The two models agree on hypofractionation to 15–20 fractions, dropping the last week, strong acceleration and twice-daily treatment.
- **The disagreement is limited and unchanged.** It is on 2 Gy courses of 38–40 days, by 2–3 Gy₁₀. DI's onset near day 37 credits these courses with the full gain from shortening, and DD's earlier onset does not.
- **The same-chemotherapy trials lower α slightly** (0.043 to 0.040), the direction expected if chemotherapy adds an effect the model does not represent.

## Limits

- **Extrapolation:** the main-set data cover 1.15–2.55 Gy per fraction, so HYPNO, 17 × 3.05 and HYHOPE are extrapolations under both models.
- **Stage mix:** arms are not adjusted for it.
- **RTOG 83-13** includes 8–10% nasopharyngeal patients, which cannot be removed from the reported arm-level control.
