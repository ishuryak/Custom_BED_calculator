# Refit to the verified trial data (2026-10-22)

This record covers the refit of the dose-independent (DI) and dose-dependent (DD) repopulation models to randomized trials whose full text is in hand. It follows the version 3 refit (`2026_10_21_v3/fit_report.md`), which used the values printed in Table 1 of Shuryak et al. 2018. Nothing in the app changed: its defaults and draws still come from the version 3 fit.

## Bottom line

1. **Hypofractionated schedules lose their apparent advantage.** Relative to 35 × 2 Gy in 47 days, the DI refit gives HYPNO 55 Gy in 20 fractions −2.5 Gy (80% range −4.8 to +0.5), against +1.4 in version 3. The 17–19 fraction schedules come out 4–5 Gy lower. All three sensitivity fits agree in direction (−1.0 to −2.2 for HYPNO), and the ranges still include zero in most of them.
2. **Dropping the last week of 35 × 2 Gy loses tumor dose in every fit** (DI −5.6 to −5.8 Gy), as in version 3.
3. **The DI model is better determined than before.** Repopulation onset Tk is 38.6 days (95% interval 24.2–44.3) and its rate λ is 0.029 per day, against 31.8 days (13.0–45.4) and 0.015 in version 3. α/β stays at 10.3 (6.2–19.1).
4. **The DD model is still not identified and does not predict better.** Leave one trial out over 22 units gives DD minus DI −3.96 (standard error 4.08). The DD fits sample poorly (76 divergent transitions, bulk ESS 250 in the primary fit), and their λ and onset constant reproduce their priors.
5. **No single trial drives the shift.** Removing any one of the 22 units moves the HYPNO estimate by at most 0.8 Gy. The shift comes mainly from the corrected values of the 2018 trials, not from the new trials (see Robustness).

## Data

**Sources.** 135 papers were in hand: 99 PDFs that Eric Wang downloaded through Columbia on 2026-10-22, 21 open-access or local PDFs, and 15 PubMed Central full texts. Every file passed an identity check (DOI on the first two pages or the title on page 1, and one French original checked by eye). The 49 listed papers not in hand are not available at Columbia and were not used (`References_refit_trials/references_availability.md`).

**Extraction.** Arm-level data were extracted from 68 trials (292 rows). Every value carries a verbatim quote, and a script confirmed each quote is in the paper's text and contains the value. The 33 trials that could enter a fit were then extracted a second time, blinded to the first. The first and second reads agree on 79 of the 80 primary rows. The 80th is a tie between two Leclerc dose levels with the same value. All files are in `References_refit_trials/search/extraction/`.

**Decisions (Igor Shuryak, 2026-10-22).** The main fit uses definitive radiotherapy of non-nasopharyngeal HNSCC, early glottic included. The 2018 Table 1 trials are re-extracted from their papers. Leclerc 2013 and DAHANCA split-course (Hansen 1997) are dropped as not randomized, and GORTEC 94-01 and Wang 1985 as not available. Failure-only trials enter as 100 − failure, and figure readings are used. Postoperative and nasopharyngeal (NPC) trials enter only as sensitivity fits.

**Conventions**, as in the 2018 paper:

- dose per fraction = total dose ÷ fractions;
- a range of fraction counts uses its midpoint;
- overall time is the delivered median or mean when reported, otherwise the planned time.

Two trials whose arms mix schedules by stratum are split into their randomized strata, each fitted with its own intercept:

- **CAIR:** the 2.0 Gy and 1.8 Gy periods.
- **Yamazaki:** tumor length up to or over two thirds of the cord.

Data file: `data/refit_trials_2026_10_22.csv` (47 arms, 22 fitting units, 20 trials, sha256 `d980ddc3…`). It is built by `References_refit_trials/search/extraction/build_refit_dataset.py`, which refuses any value that is not in a checked extraction row.

| Unit | Arms | Endpoint, time point | Source |
|---|---|---|---|
| ARTSCAN | 2 | locoregional control, 5 y | Zackrisson 2015 |
| CAIR 2.0 Gy, CAIR 1.8 Gy | 2 + 2 | local control, 5 y | Skladowski 2006, Table 2 |
| CAIR vs concomitant boost | 2 | locoregional control, 5 y | Skladowski 2013 |
| CHART | 2 | locoregional relapse-free survival, 5 y, read from figure | Saunders 2010 |
| DAHANCA 6&7 | 2 | locoregional control, 5 y | Overgaard 2003 |
| EORTC 22791 | 2 | local control, 5 y | Horiot 1992 |
| EORTC 22811 (misonidazole arm excluded) | 2 | locoregional control, 5 y, read from figure | Van den Bogaert 1995 |
| EORTC 22851 | 2 | locoregional control, 5 y | Horiot 1997 |
| Ghoshal (PGIMER) | 2 | locoregional control, 2 y | Ghoshal 2008 |
| Hliniak | 2 | locoregional control, 3 y | Hliniak 2002 |
| IAEA-ACC | 2 | locoregional control, 5 y | Overgaard 2010 |
| Jackson (Vancouver) | 2 | locoregional recurrence-free survival, 3 y | Jackson 1997 |
| KROG-0201 (glottic) | 2 | local progression-free survival, 5 y | Moon 2014 |
| Princess Margaret | 2 | locoregional control, 5 y (competing risk) | O'Sullivan 2020 |
| RTOG 83-13, scheme A | 3 | locoregional control, 2 y | Cox 1990 |
| RTOG 88-09 | 2 | 100 − locoregional failure, 2 y | Fu 1995 |
| RTOG 9003 | 4 | locoregional control, 2 y | Beitler 2014 (n from Fu 2000) |
| RTOG 9512 (T2 glottic) | 2 | locoregional control, 5 y | Trotti 2014 |
| TROG 91.01 | 2 | locoregional control, 5 y | Poulsen 2001 |
| Yamazaki (T1 glottic), two strata | 2 + 2 | local control, 5 y | Yamazaki 2006 |

Rows built on an assumption are flagged in the file's `time_rule` and `notes` columns:

- **EORTC 22811:** the control arm allowed 70 Gy in 7 weeks or 75 Gy in 9 weeks without reporting the split, so 70 Gy is used. The three-times-daily arm's 6–7 weeks is entered as 45.5 days.
- **RTOG 83-13:** overall time is derived from 2 fractions a day, 5 days a week. This is consistent with the paper's "5½ weeks for 67.2 Gy" and "less than 6½ weeks for 76.8 Gy".
- **KROG-0201:** fractions are weighted by the T1/T2 mix.

## How the 2018 Table 1 data changed

Dropped:

- **Not randomized:** Leclerc 2013 (three sequential dose levels) and Hansen 1997, whose split and continuous arms came from two consecutive trials that randomized only a sensitizer.
- **Not available at Columbia:** GORTEC 94-01 (Bourhis 2006) and Wang 1985.
- **No control at 2–5 years by arm:** Pinto 1991.

Values that differ between the 2018 table and the papers (2018 value → value used here):

| Trial | Arms | 2018 Table 1 | Read from the paper |
|---|---|---|---|
| EORTC 22791 | conventional / hyperfractionated | 62.9 / 63.9 | 40 / 59 (5-y local control) |
| CHART | conventional / CHART | 53.0 / 51.8 | 50 / 45 (5-y, figure) |
| DAHANCA 6&7 | 5 / 6 fractions per week | 64.3 / 71.3 | 60 / 70 |
| TROG 91.01 | conventional / accelerated | 49.1 / 57.6 | 47 / 52 |
| EORTC 22851 | conventional / accelerated | 48.6 / 60.3 | 46 / 59 |
| RTOG 9003 | standard / HFX / AFX-S / AFX-C | 42.5 / 52.5 / 42.3 / 52.2 | 45.7 / 53.3 / 47.8 / 53.8 |
| CAIR | 5 / 7 days per week | 32.7 / 78.4 (1.9 Gy) | 30 / 73 (2.0 Gy) and 35 / 78 (1.8 Gy) |
| Yamazaki | 2.0 / 2.25 Gy | 83.0 / 93.5 | 76 / 91 and 78 / 92 by stratum |

The 2018 paper does not say which endpoint or time point each Table 1 value represents, so these differences cannot be traced further.

## Fits

Sampler and model: `fit/fit_refit.py` and `fit/stan/` exactly as in version 3 (4 chains × 2000 draws, seed 20261021). The script now takes the data file and output folder from `BED_FIT_DATA` and `BED_FIT_OUT`, with the version 3 paths as defaults. Outputs are in `fit/output_refit_2026_10_22*/`.

- **DI:** all four variants have 0 divergent transitions, R-hat ≤ 1.004, bulk ESS ≥ 985, and 47 of 47 arms inside their 90% predictive interval.
- **DD primary:** 76 divergent transitions, R-hat 1.016, bulk ESS 250.
- **DD binomial:** does not mix (R-hat 1.118, bulk ESS 27) and is not interpreted.

Tumor-dose change relative to 35 × 2 Gy in 47 days, DI primary fit, posterior median (80% range). Effective tumor BED is BED − (repopulation + gT)/α, computed for every draw:

| Regimen | Version 3 | Main refit | Without figure readings | Plus postoperative | Plus postoperative and NPC |
|---|---|---|---|---|---|
| HYPNO 55 Gy/20 (26 d) | +1.4 (−3.6 to +10.2) | −2.5 (−4.8 to +0.5) | −2.2 (−5.3 to +2.6) | −1.0 (−4.0 to +3.8) | −1.0 (−3.9 to +4.2) |
| 17 × 3.05 Gy (23 d) | +0.4 | −4.1 (−6.7 to −0.8) | −3.8 | −2.4 | −2.3 |
| 30 × 2 Gy (40 d), last week dropped | −6.0 (−8.2 to −1.8) | −5.8 (−7.6 to −3.6) | −5.6 | −5.6 | −5.8 |
| 33 × 2 Gy, 6 per week (38 d) | +2.9 | +2.6 (+0.7 to +5.0) | +2.7 | +2.8 | +2.7 |
| 68 × 1.2 Gy twice daily (47 d) | +7.4 | +7.5 (+5.9 to +8.8) | +7.6 | +7.8 | +7.8 |

The DD columns, parameter tables and leave-one-trial-out scores for every run are in `References_refit_trials/search/refit/comparison_all_runs.md`.

## Robustness

**Drop-one influence** (`influence_di.json`): DI refit 22 times, each time without one unit (4 × 500 draws). The HYPNO median ranged from −3.0 (without RTOG 9003) to −1.6 (without EORTC 22791), against −2.4 with all units at this sampler setting. Tk stayed between 37.4 and 39.3 days.

**Where the shift comes from** (DI primary, HYPNO median):

| Data | HYPNO median |
|---|---|
| 2018 Table 1 as published (version 3) | +1.4 |
| 2018 Table 1 without Leclerc | +1.4 |
| 2018 Table 1 without all five dropped trials | +0.9 |
| The same trials with their values re-read from the papers | −1.8 |
| All verified trials (main refit) | −2.5 |

Most of the change therefore comes from correcting the 2018 values, not from the new trials or the removed entries.

## Limitations

- Endpoints differ between trials: locoregional or local control, Kaplan–Meier or competing-risk estimates, 2 to 5 years. The per-trial intercept absorbs a constant difference between trials but not a difference in how an endpoint responds to dose or time.
- The arms are not adjusted for stage mix. This was item 4 of the version 3 open items and is still open.
- Three trials contribute values read off figures (CHART, EORTC 22811, and Ang in the postoperative set). Two independent readings agree within about 1 point.
- Overall time is derived for RTOG 83-13 and approximated for EORTC 22811.
- Several of the arms used are small: CAIR strata of 21–28 patients, Jackson 41 per arm, Yamazaki strata of 31–61.
- HYPNO itself cannot be used. It gave cisplatin to 75% of patients and does not report control by arm for those without it.

## Open items

1. Igor: whether the refit should replace the version 3 fit in the app (decision 2a: only after validation). Nothing has been changed in the app.
2. JCOG0701 was dropped despite the "include" decision: no paper in hand gives its overall time or weekly schedule, which are in the unavailable primary report (Kodaira 2018, PMID 29401241).
3. Peters 1993 (postoperative) was left out of the sensitivity fits: two of its dose cells show identical values (89%, 51 patients), which needs a check against the PDF table.
4. OCAT and Teo 2000 were left out of the sensitivity fits: one arm of each has no stated overall time.
5. The DD model needs either more informative data or a reparameterization before its posteriors mean anything.
6. The `fit_refit.py` change and its test (`tests/python/test_fit_refit_paths.py`, observed failing without the change and passing with it) are on the local branch `refit_trials_2026_10_22`, not committed.
