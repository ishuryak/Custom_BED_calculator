---
title: "HNSCC BED calculator, version 4.3: how the tumor model was refitted to randomized trials"
author: "Igor Shuryak, MD, PhD – is144@cumc.columbia.edu"
date: "November 4, 2026"
---

# Summary

This document explains the tumor-model refit behind the head and neck squamous cell carcinoma (HNSCC) biologically effective dose (BED) calculator, version 4.3.0 (https://ishuryak.shinyapps.io/custom_bed_calculator_v4/). It answers five questions. How was the refit done? Which trial arms differ from the 2018 paper? How does each trial enter the fit? How do the two repopulation models compare by the Akaike information criterion with small-sample correction (AICc)? How and why is the fit Bayesian?

- **Data.** The refit uses 45 arms from 19 randomized trials (8,279 patients). Every value was read from the trial's own full text and checked against a verbatim quote. The 2018 paper (Shuryak, Hall and Brenner) used 35 arms from 16 publications. Eleven of its trials are kept and re-read, five are dropped, and eight trials are new.
- **Design.** Each trial has its own intercept, so only the randomized comparisons within a trial inform the radiobiological parameters. Differences between trials in stage, site, endpoint and era do not.
- **Models.** Two repopulation models are fitted. In the dose-independent (DI) model, accelerated regrowth starts on a fixed day. In the dose-dependent (DD) model, it starts once a set amount of cell kill is reached.
- **Model comparison.** The 2018 paper found an AICc 40 units lower for DD. With one intercept per trial, that preference disappears. On the 2018 data DD is lower by only 3.9 units. On the current data DI is lower by 9.0 units, or 4.4 after allowing for overdispersion (more variation between arms than binomial sampling allows). Leave-one-trial-out prediction points the same way but within its noise. Neither model is decisively preferred, and the DD repopulation parameters are not identified by these trials.
- **Why Bayesian.** About 24 within-trial contrasts must inform 5–6 structural parameters, one of which (tumor α/β) the trials cannot inform at all. The trial data are also overdispersed. A Bayesian fit makes each assumption an explicit prior, models the overdispersion, and gives joint parameter draws. The calculator turns those draws into an uncertainty interval for any schedule.
- **Main result.** The reference is 70 Gy in 35 fractions of 2 Gy over 47 days. Against it, HYPNO's 55 Gy in 20 fractions loses 2.0 Gy₁₀ of effective tumor BED under DI (80% interval −4.3 to +1.0) and 3.3 Gy₁₀ under calibrated DD (−5.9 to −0.5). Effective tumor BED is the linear-quadratic BED minus the dose lost to tumor regrowth during treatment, and Gy₁₀ means BED computed with α/β = 10 Gy. Calibrated DD is the DD model with its onset tied to DI at the standard schedule, because DD cannot be fitted on its own. The 2018-based fit had put it at +1.4 Gy₁₀. Most of that shift comes from correcting the 2018 values, not from the new trials.

**Where the refit is used on the site.** The main page lists schedules under the default tumor parameters (scenario D18), which the developer chose and did not fit to outcome data. The refit enters in three places. It is a preset on the Model settings page. Under Details on the main page, its parameter draws give the "Refit support" column. It also supplies the draws for the Uncertainty page. The draws shipped with version 4.3.0 (draws version 3.2) come from the main fit described here.

# 1. How the refitting analysis was done

## 1.1 Sources

135 papers were in hand as full text: 99 PDFs downloaded through Columbia on 2026-10-22, 21 open-access or local PDFs, and 15 PubMed Central full texts. Each file passed an identity check (the digital object identifier (DOI) on the first two pages or the title on page 1). Forty-nine listed papers were not available at Columbia and were not used. No value was taken from an abstract or from memory.

## 1.2 Extraction and its checks

Arm-level data were extracted from 68 trials (292 rows). Every value has a verbatim quote, and a script confirmed that each quote occurs in the paper's text and contains the value. The 33 trials that could enter a fit were then extracted a second time, blinded to the first. The two reads agree on 79 of the 80 primary rows. The 80th is a tie between two dose levels of Leclerc 2013, which was dropped as not randomized. A builder script assembles the fitting file and refuses any value that is not in a checked extraction row.

## 1.3 Which trials enter the main fit

The main fit includes randomized trials of definitive radiotherapy for non-nasopharyngeal HNSCC, early glottic cancer included, that compare fractionation schedules with tumor control reported by arm at 2 to 5 years. Five further rules apply:

- A trial reporting only failure enters as 100 minus the failure percentage.
- Values read off a published figure are used (two independent readings agreed within about 1 point).
- Postoperative and nasopharyngeal trials enter only sensitivity fits.
- ARTSCAN enters only a sensitivity fit, because about 37% of its patients had preoperative radiotherapy followed by surgery.
- Three trials found later (Vienna V-CHART, ELAN-RT, Gupta 2015) and three trials comparing fractionation under the same chemotherapy (HYPNO, ART DECO, Shivananjappa) also enter only sensitivity fits, because their endpoints are crude or competing-risk estimates in frail patients, or because chemotherapy is not modeled.

## 1.4 Conventions

- **Overall time** is the delivered median or mean when the paper reports it, otherwise the planned time.
- **Mixed fraction sizes** within an arm (for example a concomitant boost) are kept as components, so the linear-quadratic (LQ) BED is BED~LQ~ = Σd + Σd²/(α/β) over all fractions.
- **Number controlled** is the reported control rate times the number analyzed, rounded. Most rates are Kaplan–Meier estimates (39 arms) and six are cumulative-incidence estimates, so this count is an approximation.
- **Effective tumor BED** of a schedule is BED~LQ~ minus the BED lost to repopulation and background growth, (repopulation + gT)/α. Results are reported as the change against 35 × 2 Gy in 47 days, in Gy₁₀ when α/β is 10 Gy.

## 1.5 The two models

Tumor control probability (TCP) follows the Poisson form TCP = exp(−K S), where K is the initial number of clonogens and S their surviving fraction. Taking log(−log TCP) turns this into a linear predictor. For arm *i* of fitting unit *j*:

$$\log(-\log \mathrm{TCP}_i) = c_j - \alpha\,(B_i - \bar{B}_j) + R_i + g\,(T_i - \bar{T}_j) + \sigma z_i$$

Here *B*ᵢ is the arm's BED~LQ~ and *T*ᵢ its overall time in days, while B̄ⱼ and T̄ⱼ are the means within unit *j*. The unit's intercept is *c*ⱼ, α is the linear radiosensitivity and *g* the background growth rate per day. The term σ*z*ᵢ, with *z*ᵢ ~ N(0, 1), is the arm-level overdispersion. The repopulation term *R*ᵢ differs between the models:

- **DI:** $R_i = \lambda\, s(T_i - T_k)$, with a fixed onset day *T*ₖ and rate λ per day.
- **DD:** $R_i = \lambda_{DD}\, P_{kill,i}\, s\!\left(T_i - \frac{C}{\alpha}\,\frac{T_i}{B_i}\right)$, with $P_{kill,i} = 1 - e^{-\alpha B_i / T_i}$, the mean fraction of clonogens killed per day. Onset comes once a set BED, C/α, has been delivered, and the rate grows with the daily cell kill.

The function s(x) = log(1 + e³ˣ)/3 is a smooth version of max(0, x), the same hinge the calculator uses. Centering BED and time within each unit removes the correlation between the intercepts and α, which otherwise made the sampler struggle. The repopulation equations are those of the 2018 paper. The main structural change is that its three stage-specific clonogen numbers are replaced by one intercept per unit (Section 3). The refit also adds the overdispersion term, estimates α/β and *g* instead of fixing them, and uses the smooth hinge where 2018 truncated at zero.

## 1.6 Fitting and quality checks

The DI and free DD models are fitted in Stan (CmdStan, 4 chains, 1500 warm-up and 2000 sampling iterations, fixed seed). The calibrated DD fits (Section 5.3) are smaller conditional fits, 4 chains of 600 warm-up and 250 sampling iterations, lengthened up to 16-fold until they pass. A fit is accepted only if every sampled parameter and every regimen contrast passes a fixed quality gate. The gate requires no divergent transitions, no hits of the maximum tree depth, R-hat at most 1.01, and bulk and tail effective sample size (ESS) at least 400. It also requires an estimated Bayesian fraction of missing information (E-BFMI, an energy diagnostic) of at least 0.3. The primary DI fit passes with maximum R-hat 1.002 and minimum bulk ESS 1,489, and all 45 arms fall inside their 90% posterior predictive interval. Four variants were fitted: primary, binomial (no overdispersion), the 2018 restrictions (α/β fixed at 10 Gy and *g* at 0.0154 per day), and a doubled-width α/β prior.

## 1.7 Sensitivity fits

HYPNO 55 Gy/20 against 35 × 2 Gy/47 d, tumor BED change in Gy₁₀, median (80% interval):

| Data set | Arms | DI | Calibrated DD | α/β median (Gy) |
|---|---|---|---|---|
| Main | 45 | −2.0 (−4.3 to +1.0) | −3.3 (−5.9 to −0.5) | 10.6 |
| Effective N halved | 45 | −1.9 (−4.5 to +1.5) | −3.1 (−5.8 to −0.4) | 10.5 |
| Planned instead of delivered overall time | 45 | +0.3 (−2.4 to +3.7) | −2.1 (−5.2 to +1.0) | 10.7 |
| Kaplan–Meier units only | 39 | −1.9 (−4.3 to +1.2) | −3.3 (−6.1 to −0.5) | 10.1 |
| RTOG 9003 rates from the 2014 update | 45 | −1.9 (−4.3 to +1.1) | −3.3 (−5.9 to −0.4) | 10.6 |
| Plus ARTSCAN | 47 | −2.7 (−5.0 to +0.4) | −3.7 (−6.1 to −0.8) | 10.6 |
| Plus three later trials | 51 | −1.7 (−3.8 to +1.0) | −2.9 (−5.3 to −0.5) | 10.8 |
| Plus same-chemotherapy trials | 51 | −1.0 (−3.5 to +2.2) | −2.4 (−5.2 to +0.2) | 10.2 |
| Intercept prior sd 10 instead of 3 | 45 | −1.9 (−4.3 to +1.1) | −3.4 (−6.1 to −0.8) | 10.6 |

The calibrated DD fits for the sensitivity sets used 200 DI draws each, against 1,000 for the main set. The overall-time rule matters most: with the planned time in place of the delivered time, HYPNO's DI change is close to zero.

# 2. Trial arms compared with the 2018 paper

The 2018 paper fitted 35 arms from 16 publications (its Table 1, 7,283 patients). The main refit fits 45 arms from 19 trials, with 21 fitting units and 8,279 patients.

## 2.1 Trials kept, with values re-read from the papers

| Trial | 2018 arms | Refit arms | Refit endpoint, time point |
|---|---|---|---|
| CHART | 2 | 2 | locoregional relapse-free survival, 5 y (figure) |
| DAHANCA 6&7 | 2 | 2 | locoregional control, 5 y |
| EORTC 22791 | 2 | 2 | local control, 5 y |
| EORTC 22851 | 2 | 2 | locoregional control, 5 y |
| Hliniak | 2 | 2 | locoregional control, 3 y |
| Jackson (Vancouver) | 2 | 2 | locoregional recurrence-free survival, 3 y |
| RTOG 9003 | 4 | 4 | locoregional control, 2 y |
| TROG 91.01 | 2 | 2 | locoregional control, 5 y |
| CAIR (Skladowski 2006) | 2 | 4 (2.0 Gy and 1.8 Gy periods) | local control, 5 y |
| CAIR vs concomitant boost (Skladowski 2013) | 2 | 2 | locoregional control, 5 y |
| Yamazaki (T1 glottic) | 2 | 4 (two tumor-length strata) | local control, 5 y |

CAIR and Yamazaki each report results for strata that used different schedules, so each stratum is a separate fitting unit with its own intercept.

Values that differ between the 2018 table and the papers:

| Trial | Arms | 2018 Table 1 (%) | Read from the paper (%) |
|---|---|---|---|
| EORTC 22791 | conventional / hyperfractionated | 62.9 / 63.9 | 40 / 59 |
| CHART | conventional / CHART | 53.0 / 51.8 | 50 / 45 |
| DAHANCA 6&7 | 5 / 6 fractions per week | 64.3 / 71.3 | 60 / 70 |
| TROG 91.01 | conventional / accelerated | 49.1 / 57.6 | 47 / 52 |
| EORTC 22851 | conventional / accelerated | 48.6 / 60.3 | 46 / 59 |
| RTOG 9003 | standard / hyperfractionated / accelerated split / concomitant boost | 42.5 / 52.5 / 42.3 / 52.2 | 46.0 / 54.4 / 47.5 / 54.5 |
| CAIR | 5 / 7 days per week | 32.7 / 78.4 (one pair at 1.9 Gy, apparently both periods pooled) | 30 / 73 (2.0 Gy) and 35 / 78 (1.8 Gy) |
| CAIR vs concomitant boost | CAIR / concomitant boost | 65.9 / 66.3 | 63 / 65 |
| Jackson | once / twice daily | 43.9 / 48.8 | 44.3 / 49.1 |
| Yamazaki | 2.0 / 2.25 Gy | 83.0 / 93.5 | 76 / 91 and 78 / 92 by stratum |

The 2018 paper does not say which endpoint or time point each Table 1 value represents, so these differences cannot be traced further.

## 2.2 Trials dropped

| 2018 entry | Reason |
|---|---|
| Leclerc 2013 (3 arms) | not randomized (three sequential dose levels) |
| Hansen 1997, split course (2 arms) | not randomized: the arms came from two consecutive trials that randomized only a sensitizer |
| Bourhis 2006, GORTEC 94-01 (2 arms) | full text not available at Columbia |
| Wang 1985 (2 arms) | full text not available at Columbia |
| Pinto 1991 (2 arms) | control at 2–5 years not reported by arm |

## 2.3 Trials added

| Trial | Arms | Contrast | Endpoint, time point |
|---|---|---|---|
| EORTC 22811 | 2 | three fractions a day, split course, vs conventional (misonidazole arm excluded) | locoregional control, 5 y (figure) |
| Ghoshal (Postgraduate Institute of Medical Education and Research) | 2 | concomitant boost vs conventional | locoregional control, 2 y |
| IAEA-ACC | 2 | 6 vs 5 fractions per week | locoregional control, 5 y |
| KROG 0201 (glottic) | 2 | 2.25 Gy vs 2.0 Gy | local progression-free survival, 5 y |
| Princess Margaret (O'Sullivan 2020) | 2 | hyperfractionated vs 2.55 Gy daily | locoregional control, 5 y (cumulative incidence) |
| RTOG 83-13, scheme A | 3 | 67.2, 72.0 and 76.8 Gy at 1.2 Gy twice daily | locoregional control, 2 y |
| RTOG 88-09 | 2 | split-course vs concomitant boost acceleration | 100 − locoregional failure, 2 y |
| RTOG 9512 (T2 glottic) | 2 | hyperfractionated vs once daily | locoregional control, 5 y (cumulative incidence) |

These eight trials add 17 arms. The 28 other arms come from the eleven kept trials (24 arms in 2018, plus 4 from splitting CAIR and Yamazaki into strata).

## 2.4 Where the change in result comes from

The DI result for HYPNO moved from +1.4 Gy₁₀ (fit to the 2018 table) to about −2 Gy₁₀. To separate the causes, the DI model was refitted step by step. These runs used an earlier version of the trial set, which still included ARTSCAN, took RTOG 9003 rates from the 2014 update and treated mixed fraction sizes as one average size. The re-read values behind the last two rows therefore differ slightly from those in Section 2.1:

| Data | HYPNO median (Gy₁₀) |
|---|---|
| 2018 Table 1 as published | +1.4 |
| 2018 Table 1 without Leclerc | +1.4 |
| 2018 Table 1 without all five dropped trials | +0.9 |
| The same trials with values re-read from the papers | −1.8 |
| All verified trials | −2.5 |

Most of the change therefore comes from correcting the 2018 values, not from the new trials or the removed entries.

# 3. How each trial gets its own result

**Each arm's value comes from its own paper.** The endpoint and time point are the ones that trial reports (Section 2), and the number analyzed is the count for that endpoint.

**Each trial has its own intercept.** The 2018 paper explained differences in tumor control between trials with three shared clonogen numbers, one per stage group (T1–2, T3, T4). Refitting that structure needs each trial's stage mix, which Table 1 does not give and which many papers report only in part. The refit instead gives each fitting unit its own intercept *c*ⱼ. The intercept absorbs anything that shifts a whole trial's control up or down: stage mix, tumor site, the endpoint and its time point, patient selection and treatment era. α, *g*, the repopulation parameters and σ are shared by all trials. As a result, a trial informs them only through how its randomized arms differ from one another, which is the comparison randomization protects. The cost is that a trial with one usable contrast contributes one piece of information, however many patients it has.

**Each trial's fit is checked against its own data.** The table below compares each arm's observed control with the primary DI fit's posterior predictive interval. The prediction uses that trial's fitted intercept, so this is an in-sample check with limited power to detect misfit. All 45 arms fall inside their 90% interval.

| Trial | Arm | N | Observed (%) | Predicted median (90% interval) |
|---|---|---|---|---|
| CAIR, 2.0 Gy period | CAIR | 23 | 73.9 | 60.9 (39.1–78.3) |
| | control | 21 | 28.6 | 42.9 (19.0–66.7) |
| CAIR, 1.8 Gy period | CAIR | 28 | 78.6 | 64.3 (46.4–82.1) |
| | control | 28 | 35.7 | 46.4 (28.6–67.9) |
| CAIR vs concomitant boost | CAIR | 172 | 62.8 | 63.4 (54.7–72.1) |
| | concomitant boost | 171 | 64.9 | 64.3 (55.0–72.5) |
| CHART | CHART | 550 | 45.1 | 45.5 (37.1–53.3) |
| | conventional | 365 | 49.9 | 49.6 (40.8–57.5) |
| DAHANCA 6&7 | 6 per week | 750 | 70.0 | 68.8 (63.1–74.1) |
| | 5 per week | 726 | 60.1 | 61.4 (55.2–67.5) |
| EORTC 22791 | conventional | 159 | 40.3 | 45.3 (35.2–55.3) |
| | hyperfractionated | 166 | 59.0 | 54.2 (44.0–63.3) |
| EORTC 22811 | conventional | 168 | 26.8 | 25.0 (16.7–33.9) |
| | three per day | 163 | 22.7 | 24.5 (16.6–33.7) |
| EORTC 22851 | conventional | 255 | 45.9 | 41.2 (31.4–49.8) |
| | accelerated | 257 | 59.1 | 64.2 (56.0–72.0) |
| Ghoshal | concomitant boost | 143 | 73.4 | 68.5 (59.4–76.9) |
| | conventional | 142 | 54.2 | 59.2 (48.6–69.0) |
| Hliniak | conventional | 193 | 76.2 | 76.2 (68.4–82.9) |
| | accelerated | 188 | 80.9 | 80.9 (74.5–86.7) |
| IAEA-ACC | 6 per week | 452 | 42.0 | 40.9 (33.6–48.5) |
| | 5 per week | 448 | 29.9 | 31.0 (23.9–38.8) |
| Jackson | once daily | 41 | 43.9 | 36.6 (22.0–53.7) |
| | twice daily | 41 | 48.8 | 56.1 (39.0–73.2) |
| KROG 0201 | 2.0 Gy | 82 | 78.0 | 81.7 (70.7–90.2) |
| | 2.25 Gy | 74 | 87.8 | 85.1 (75.7–91.9) |
| Princess Margaret | hyperfractionated | 169 | 49.1 | 46.7 (36.7–56.2) |
| | 2.55 Gy daily | 162 | 40.1 | 42.6 (32.7–52.5) |
| RTOG 83-13 | 67.2 Gy | 63 | 25.4 | 30.2 (19.0–42.9) |
| | 72.0 Gy | 58 | 36.2 | 34.5 (22.4–48.3) |
| | 76.8 Gy | 116 | 42.2 | 39.7 (29.3–50.9) |
| RTOG 88-09 | split course | 38 | 55.3 | 47.4 (28.9–63.2) |
| | concomitant boost | 32 | 43.8 | 56.2 (37.5–71.9) |
| RTOG 9003 | standard | 268 | 45.9 | 45.1 (36.9–53.4) |
| | hyperfractionated | 263 | 54.4 | 56.7 (48.7–64.3) |
| | accelerated, split | 274 | 47.4 | 45.6 (37.2–53.6) |
| | concomitant boost | 268 | 54.5 | 55.2 (47.8–62.7) |
| RTOG 9512 | once daily | 119 | 67.2 | 66.4 (55.5–75.6) |
| | twice daily | 120 | 73.3 | 75.0 (65.8–82.5) |
| TROG 91.01 | conventional | 171 | 46.8 | 48.5 (38.0–57.9) |
| | accelerated | 172 | 51.7 | 50.6 (40.7–59.9) |
| Yamazaki, smaller tumors | 2.0 Gy | 31 | 77.4 | 83.9 (67.7–96.8) |
| | 2.25 Gy | 31 | 90.3 | 83.9 (71.0–96.8) |
| Yamazaki, larger tumors | 2.0 Gy | 57 | 77.2 | 84.2 (71.9–93.0) |
| | 2.25 Gy | 61 | 91.8 | 86.9 (75.4–93.4) |

The observed percentages are the rounded counts divided by N, so they can differ from the published rate by up to 1.4 points in the smallest arms.

**Each trial's out-of-sample result.** In leave-one-trial-out cross-validation, each of the 19 trials was removed in turn and the DI and free DD models were refitted to the other 18. In each fitting unit of the removed trial, the observed control of the first arm fixed the intercept, and the other arms (24 in total) were predicted. The score is this conditional log predictive density, where higher is better. It does not include the uncertainty of the first arm's observed rate:

| Trial | DI | DD | DD − DI |
|---|---|---|---|
| CAIR (both periods) | −12.58 | −15.61 | −3.03 |
| CAIR vs concomitant boost | −3.07 | −3.31 | −0.24 |
| CHART | −4.35 | −4.80 | −0.46 |
| DAHANCA 6&7 | −4.75 | −5.70 | −0.95 |
| EORTC 22791 | −5.02 | −4.63 | +0.40 |
| EORTC 22811 | −3.45 | −3.42 | +0.03 |
| EORTC 22851 | −7.41 | −3.73 | +3.68 |
| Ghoshal | −5.29 | −6.00 | −0.71 |
| Hliniak | −2.82 | −3.12 | −0.30 |
| IAEA-ACC | −4.14 | −5.10 | −0.96 |
| Jackson | −3.72 | −2.87 | +0.85 |
| KROG 0201 | −2.83 | −3.46 | −0.63 |
| Princess Margaret | −3.65 | −3.69 | −0.03 |
| RTOG 83-13 | −6.61 | −6.50 | +0.12 |
| RTOG 88-09 | −4.01 | −3.53 | +0.48 |
| RTOG 9003 | −11.62 | −12.25 | −0.63 |
| RTOG 9512 | −2.86 | −2.88 | −0.02 |
| TROG 91.01 | −3.50 | −4.06 | −0.56 |
| Yamazaki (both strata) | −7.44 | −8.06 | −0.61 |
| **Total** | **−99.12** | **−102.68** | **−3.56 (standard error 5.37)** |

DI predicts 13 of the 19 trials better, although three of these differences (EORTC 22811, Princess Margaret, RTOG 9512) are below 0.05 and amount to ties. CAIR favors DI most and EORTC 22851 favors DD most. The free DD fits in this exercise had 415 divergent transitions across folds (175 of them in the fold without RTOG 83-13) against none for DI, so the DD scores are less reliable.

**No single trial drives the HYPNO result.** This check used the earlier trial set described in Section 2.4 (47 arms in 22 units, ARTSCAN included), where the HYPNO median with all units is −2.4 Gy₁₀. When the DI model was refitted 22 times, each time without one fitting unit, the HYPNO median ranged from −3.0 Gy₁₀ (without RTOG 9003) to −1.6 (without EORTC 22791). The onset day *T*ₖ stayed between 37.4 and 39.3 days.

# 4. Model-versus-model AICc

The 2018 paper compared DI and DD by AICc after maximum-likelihood fits. The Bayesian refit compares them by leave-one-trial-out prediction (Section 3), which is the better test of how a model extrapolates to a new trial. To answer the AICc question directly, the 2018-style comparison was also repeated by maximum likelihood with one intercept per fitting unit. That calculation was done for this document. Its script and outputs are in the calculator's public code repository (https://github.com/ishuryak/Custom_BED_calculator/tree/v4/docs/refit_explainer). It uses a binomial likelihood, α/β fixed at 10 Gy and *g* at 0.0154 per day as in 2018, and no overdispersion term.

| Data and model structure | Arms / units | Parameters per model | −log L, DI | −log L, DD | AICc(DI) − AICc(DD) | Favors |
|---|---|---|---|---|---|---|
| 2018 Table 1, stage-based clonogen numbers (published) | 35 / – | – | – | – | +40 (reported) | DD |
| 2018 Table 1, one intercept per trial | 35 / 16 | 19 | 4713.10 | 4711.14 | +3.9 | DD, weakly |
| Current main set, one intercept per unit | 45 / 21 | 24 | 5342.82 | 5347.31 | −9.0 | DI |
| Current main set, *g* also estimated | 45 / 21 | 25 | 5342.41 | 5346.59 | −8.4 | DI |

The quasi-likelihood AICc (QAICc) corrects for overdispersion by dividing each log-likelihood by the larger of the two models' Pearson χ²/df. With it, the difference is +1.1 for the 2018 data (χ²/df 3.65) and −4.4 for the current data (χ²/df 2.05). With *g* estimated it is −4.0.

How to read these numbers:

- **The two models have the same number of parameters in every row.** The AICc difference is therefore exactly twice the log-likelihood difference. It does not depend on the small-sample correction or on whether the sample size is counted in arms or patients. The absolute AICc values are omitted because they include constants that cancel.
- **The published 40-unit preference for DD does not survive per-trial intercepts.** The per-arm log-likelihood differences printed in the 2018 Table 1 sum to 20.2, consistent with the reported 40. The 2018 models explained between-trial differences in control through three shared stage-based clonogen numbers plus their time and dose terms. Two explanations are plausible, and neither was tested separately. First, much of DD's advantage may have come from fitting differences between trials, which randomization does not protect. Second, the 2018 DI fit pushed its clonogen numbers to their lower limit of 10, as that paper notes, so a bound may have constrained DI. With per-trial intercepts only within-trial contrasts remain, no such bound applies, and the difference falls to 3.9 units.
- **On the corrected and expanded data, the sign reverses.** DI is preferred by 9.0 units, or about 4 once overdispersion is allowed for. That is below the threshold of 10 that the 2018 paper used for significant support.
- **The other checks agree.** Leave-one-trial-out favors DI by 3.6 log units, with a standard error of 5.4, so the difference is within noise. The free Bayesian DD fit fails the quality gate (28–74 divergent transitions per variant), and its rate and onset constant mostly reproduce their priors (onset constant *C* median 10.5 against a prior median of 10).
- **The calculator therefore keeps both models.** It treats DI as the fitted model and DD as a variant calibrated to DI (Section 5.3). It reports where the two disagree instead of choosing one.

In the DD fit with *g* estimated, *g* sits at its upper bound of 0.021 per day, so that row's DD likelihood may be slightly limited by the bound.

# 5. The Bayesian model: how and why

## 5.1 Why a Bayesian fit

1. **Few contrasts for the parameters.** 45 arms minus 21 intercepts leave about 24 within-unit contrasts to inform α, α/β, λ, *T*ₖ, *g* and σ. Maximum likelihood with this little information gives point estimates at parameter bounds and confidence intervals that depend on approximations. Priors keep each parameter in its plausible range and state openly what is assumed.
2. **One parameter is not informed by the data at all.** The trials span only 1.15–2.55 Gy per fraction, too narrow a range to estimate tumor α/β. Its posterior (median 10.6 Gy, 95% interval 6.6–19.1) essentially reproduces its prior (95% 5–20 Gy). A Bayesian fit passes that prior uncertainty into every schedule's result. A fixed α/β of 10 Gy, as in 2018, would hide it.
3. **Overdispersion.** Arms vary more than binomial sampling allows. On the 2018 data, Pearson χ²/df was 3.65 with per-trial intercepts. Ignoring this makes intervals too narrow and model comparisons overconfident. The arm-level term σ*z*ᵢ absorbs it.
4. **Joint uncertainty, delivered as draws.** The calculator needs, for any schedule a user types, an interval for the change in tumor BED. That change depends on α, α/β, λ, *T*ₖ and *g* together, and those parameters are correlated. Posterior draws give the answer for any schedule by evaluating it once per draw. The calculator's "Refit support" column, the share of paired DI and calibrated-DD draws in which a schedule meets the tumor condition (within 5 Gy₁₀ of the standard schedule under DI or DD), is computed this way. Profile-likelihood intervals for single parameters, as in 2018, cannot provide it.
5. **Built-in checks.** Posterior predictive intervals for every arm (Section 3) and leave-one-trial-out scores come directly from the draws.

Stan never runs inside the app. The fitted draws are exported once, checked against the quality gate and their content hashes, and shipped with the site.

## 5.2 Priors

In this table, sd is the standard deviation and log-sd the standard deviation on the log scale.

| Parameter | Prior | Reason |
|---|---|---|
| α/β (Gy) | lognormal, median 10, 95% 5–20 | the range the 2018 paper calls plausible |
| α (per Gy) | lognormal, median 0.2, log-sd 1 | weakly informative |
| λ, DI (per day) | lognormal, median 0.1, log-sd 1.5 | weakly informative, spans about 0.005–1.9 |
| *T*ₖ, DI (days) | normal, mean 28, sd 10, limited to 0–60 | centered on the classic onset near 4 weeks |
| λ~DD~ (per day) | lognormal, median 0.5, log-sd 1.5 | weakly informative |
| *C*, DD onset constant | lognormal, median 10, log-sd 1 | weakly informative |
| *g* (per day) | uniform 0.0046–0.021 | the 2018 paper's plausible range |
| *c*ⱼ | normal(0, 3) on the centered scale | weak (sd 10 tested as a sensitivity fit) |
| σ | half-normal(0, 0.5) | weak |

Priors from earlier pooled analyses of the same trials (for example Withers or Fowler estimates) were deliberately not used, because several trials would then be counted twice.

## 5.3 Calibrated DD

The free DD fit is not identified: its rate and onset constant reproduce their priors, and it fails the quality gate. The calculator therefore uses a calibrated DD, fitted conditionally on the DI posterior. For each of 1,000 DI draws:

1. α, α/β and *g* are taken from the DI draw.
2. The DD onset constant is solved so that DD repopulation at the anchor schedule (35 × 2 Gy in 47 days) equals DI repopulation there. Both models therefore give the same effective tumor BED for the standard schedule.
3. λ~DD~, the intercepts and σ are fitted to all 45 arms, with a soft penalty against any arm's onset falling before day 7.
4. Each conditional fit must pass the same gate (4 chains, lengthened and retried when needed). Four DD draws are kept per DI draw, giving 4,000 paired rows.

The two models agree by construction at the standard schedule and differ only in how they extrapolate to shorter or longer courses, which is the question the calculator asks. Their agreement at the anchor is therefore not independent confirmation. All 1,000 conditional fits pass the gate (0 divergences, R-hat at most 1.010, bulk ESS at least 400). λ~DD~ has median 0.41 per day (95% interval 0.12–1.04), narrower than its prior.

## 5.4 Posterior estimates (primary DI fit)

| Quantity | Median (95% interval) |
|---|---|
| α (per Gy) | 0.044 (0.030–0.060) |
| α/β (Gy) | 10.6 (6.6–19.1), close to the prior |
| λ (per day) | 0.028 (0.011–0.049) |
| Onset *T*ₖ (days) | 37.2 (23.7–42.8) |
| *g* (per day) | 0.012 (0.005–0.020) |
| λ/α, BED lost per day after onset (Gy) | 0.63 (0.21–1.19) |
| *g*/α, BED lost per day throughout (Gy) | 0.28 (0.11–0.48) |
| σ (arm-level overdispersion) | 0.065 (0.004–0.174) |

The small α (about 0.04 per Gy) is a population-level slope. It reflects patient-to-patient variation within each arm, not cell-line radiosensitivity. It does not change BED differences, but it makes any absolute tumor-control prediction shallow.

## 5.5 What the refit says about common schedules

Tumor BED change against 35 × 2 Gy in 47 days, Gy₁₀, median (80% interval), paired DI and calibrated-DD draws:

| Schedule | DI | Calibrated DD | Same sign in |
|---|---|---|---|
| HYHOPE-type 46.5 Gy/15 (19 d) | −9.2 (−11.9 to −5.8) | −9.6 (−12.3 to −6.3) | 100% |
| 30 × 2 Gy (40 d), last week dropped | −5.8 (−7.3 to −3.8) | −5.8 (−7.7 to −3.2) | 99% |
| 17 × 3.05 Gy (23 d) | −3.6 (−6.3 to −0.2) | −4.5 (−7.4 to −1.5) | 94% |
| 66 Gy/33 (45 d) | −3.0 (−3.5 to −2.3) | −2.4 (−3.3 to −0.9) | 97% |
| HYPNO 55 Gy/20 (26 d) | −1.9 (−4.4 to +1.2) | −3.3 (−5.9 to −0.5) | 87% |
| 60 Gy/25 (33 d) | +0.1 (−2.0 to +2.8) | −2.0 (−4.1 to +0.3) | 58% |
| 65 Gy/30 (40 d) | +1.2 (−0.4 to +3.2) | −1.0 (−2.1 to +0.5) | 31% |
| 33 × 2 Gy, 6 per week (38 d) | +2.8 (+0.9 to +5.0) | −0.3 (−1.5 to +1.2) | 39% |
| 70 Gy/35, 6 per week (40 d) | +6.1 (+4.5 to +8.1) | +2.0 (+1.1 to +2.8) | 100% |
| 68 × 1.2 Gy twice daily (47 d) | +7.6 (+6.1 to +8.9) | +4.6 (+2.0 to +6.5) | 98% |

The DI column here uses the 1,000 DI draws paired with DD. The full DI posterior gives HYPNO −2.0 (−4.3 to +1.0). Across all 55 pairs of these schedules, the two models rank the pair the same way 91.9% of the time. They agree on hypofractionation to 15–20 fractions and on dropping the last week. On strong acceleration and twice-daily treatment they agree in sign, but DI predicts a gain 3–4 Gy₁₀ larger. The sign disagreements are mainly for 2 Gy courses of 38–40 days and for 60 Gy in 25 fractions, by 2–3 Gy₁₀. DI's onset near day 37 credits such courses with the full gain from shortening, and DD's earlier onset does not. Schedules above 2.55 Gy per fraction (HYPNO, 17 × 3.05 Gy, HYHOPE) are extrapolations under both models.

# 6. Limitations

- **Endpoints differ between trials.** They include local or locoregional control, Kaplan–Meier or competing-risk estimates, and 2 to 5 years of follow-up. The per-trial intercept absorbs a constant difference between trials, but not a difference in how an endpoint responds to dose or time.
- **Approximate likelihood.** The binomial likelihood is applied to rounded counts built from actuarial rates. The "effective N halved" sensitivity fit changes HYPNO by 0.1 Gy₁₀.
- **No stage adjustment.** Arms within a trial are assumed comparable by randomization, and arms are not adjusted for stage mix.
- **RTOG 83-13** includes 8–10% nasopharyngeal patients, who cannot be removed from the arm-level control.
- **Values read off figures** come from CHART and EORTC 22811. Overall time is derived for RTOG 83-13 and approximated for EORTC 22811.
- **Small arms** come from CAIR (strata of 21–28 patients), Jackson (41 per arm) and Yamazaki (strata of 31–61).
- **Overall-time rule.** Using the planned rather than the delivered overall time moves HYPNO's DI change to about zero.
- **Chemotherapy is not modeled.** The parameters come from radiotherapy-alone trials, and adding same-chemotherapy trials moves HYPNO toward zero (−1.0 Gy₁₀ under DI).
- **Tumor α/β is set by its prior,** and the DD repopulation parameters are not identified by these trials.
- **The AICc comparison** in Section 4 omits overdispersion except through the QAICc correction, and it fixes α/β and *g*.

# Sources

- Shuryak I, Hall EJ, Brenner DJ. Dose dependence of accelerated repopulation in head and neck cancer: supporting evidence and clinical implications. Radiother Oncol 2018;127:20–26. https://pubmed.ncbi.nlm.nih.gov/29534828/
- Trial sources used in the main fit, with PubMed identifiers (PMID):
    - CAIR, Skladowski 2006 (https://pubmed.ncbi.nlm.nih.gov/17011446/)
    - CAIR vs concomitant boost, Skladowski 2013 (https://pubmed.ncbi.nlm.nih.gov/22836063/)
    - CHART, Saunders 2010 (https://pubmed.ncbi.nlm.nih.gov/20394851/)
    - DAHANCA 6&7, Overgaard 2003 (https://pubmed.ncbi.nlm.nih.gov/14511925/)
    - EORTC 22791, Horiot 1992 (https://pubmed.ncbi.nlm.nih.gov/1480768/)
    - EORTC 22811, Van den Bogaert 1995 (https://pubmed.ncbi.nlm.nih.gov/7569030/)
    - EORTC 22851, Horiot 1997 (https://pubmed.ncbi.nlm.nih.gov/9288839/)
    - Ghoshal 2008 (https://pubmed.ncbi.nlm.nih.gov/18343310/)
    - Hliniak 2002 (https://pubmed.ncbi.nlm.nih.gov/11830307/)
    - IAEA-ACC, Overgaard 2010 (https://pubmed.ncbi.nlm.nih.gov/20382075/)
    - Jackson 1997 (https://pubmed.ncbi.nlm.nih.gov/9165135/)
    - KROG 0201, Moon 2014 (https://pubmed.ncbi.nlm.nih.gov/24161568/)
    - Princess Margaret, O'Sullivan 2020 (https://pubmed.ncbi.nlm.nih.gov/32055717/)
    - RTOG 83-13, Cox 1990 (https://pubmed.ncbi.nlm.nih.gov/2180866/)
    - RTOG 88-09, Fu 1995 (https://pubmed.ncbi.nlm.nih.gov/7790243/)
    - RTOG 9003, Fu 2000 (https://pubmed.ncbi.nlm.nih.gov/10924966/)
    - RTOG 9512, Trotti 2014 (https://pubmed.ncbi.nlm.nih.gov/25035199/)
    - TROG 91.01, Poulsen 2001 (https://pubmed.ncbi.nlm.nih.gov/11439206/)
    - Yamazaki 2006 (https://pubmed.ncbi.nlm.nih.gov/16169681/)
- Fitting data: `data/refit_trials_2026_10_25c_main.csv` in the calculator's public code repository, https://github.com/ishuryak/Custom_BED_calculator (branch v4, SHA-256 51ed5d8a…8409), with a per-value provenance table beside it. Fit summaries: `fit/output_refit_2026_10_25c_main/`.
