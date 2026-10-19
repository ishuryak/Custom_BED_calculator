---
title: "Tumor BED under the DI and DD repopulation models for standard and clinically tested head and neck regimens"
subtitle: "DD maximum repopulation rate matched to the DI model (0.666 per day) and at the current default (0.5 per day)"
---

Igor Shuryak, MD, PhD

is144@cumc.columbia.edu

October 18, 2026

# Purpose

This note shows how the two tumor repopulation models in the head and neck squamous cell carcinoma (HNSCC) biologically effective dose (BED) calculator behave on the standard regimen and on regimens that have been tested clinically. It compares two settings of the dose-dependent (DD) model. In the first, its maximum repopulation rate is chosen so that the DD model reproduces the dose-independent (DI) model at the standard 70 Gy in 35 fractions (0.666 per day). The second is the calculator's current default (0.5 per day). The question is whether the repopulation onset times, doubling times and daily BED losses that the models imply are clinically plausible.

# Models and settings

Both models use the linear-quadratic (LQ) BED, $K = nd\,(1 + d/(\alpha/\beta))$, for $n$ fractions of dose $d$ given over $T$ days.

- **DI model.** Repopulation starts at a fixed day $T_k$ and proceeds at a fixed rate $\lambda_{DI}$: $BED_{DI} = K - \lambda_{DI}\max(0, T - T_k)/\alpha$.
- **DD model** (Shuryak, Hall and Brenner 2018 [1]). Repopulation starts when cumulative cell killing reaches a fixed level $C$, at day $C\,T/(\alpha K)$, and proceeds at rate $\lambda_{DD}P_{kill}$, where $P_{kill} = 1 - e^{-\alpha K/T}$ is the mean fraction of tumor cells killed per day: $BED_{DD} = K - \lambda_{DD}P_{kill}\max(0, T - C\,T/(\alpha K))/\alpha$.
- **Calibration.** $C$ is solved so that $BED_{DD} = BED_{DI}$ for the reference regimen (70 Gy in 35 fractions over 47 days).
- **Matching.** Here, in addition, $\lambda_{DD} = \lambda_{DI}/P_{kill,ref}$. With this choice the two models agree at the reference not only in BED but also in onset (day 28) and rate (0.2 per day). They differ only where daily or cumulative cell killing differs from the reference.

Each regimen has one repopulation doubling time, $\ln 2$ divided by its repopulation rate, and one daily BED loss, the rate divided by $\alpha$, in Gy~10~ per day after onset.

| Parameter | Value |
|---|---|
| Tumor α | 0.2 Gy^−1^ (calculator default) |
| Tumor α/β | 10 Gy |
| DI rate λ~DI~ | 0.2 per day (doubling time 3.5 days) |
| DI onset T~k~ | 28 days |
| DD maximum rate λ~DD~ | 0.666 per day (matched to DI, Tables 1 and 2) or 0.5 per day (current default, Tables 3 and 4) |
| DD onset level C/α | 50.0 Gy~10~ (matched) or 38.8 Gy~10~ (default) |
| Reference regimen | 70 Gy in 35 fractions, 47 days |

**Clinical benchmarks used for the plausibility assessment.** Repopulation onset about 4 ± 1 weeks (Withers et al. 1988 [2]) or 21–32 days (Fowler 2010 [3]). Doubling time during repopulation about 4 days (Withers [2]) or 3 days (Fowler [3]). Daily BED loss about 0.65–0.77 Gy~10~: Withers estimated 0.6 Gy per day in 2-Gy fractions, equal to 0.72 Gy~10~ [2], Fowler gave 0.66 Gy~10~ and "no more than about 0.7 Gy/day" [3], and the HYPNO design assumed 0.65 Gy per day [4]. These estimates come mainly from 5–7-week schedules.

# Table 1. Repopulation behavior of each regimen, DD rate matched to DI (0.666 per day)

Overall times $T$ count the first to the last treatment day inclusive. For rows F–M, dose per fraction, number of fractions and $T$ are the average values tabulated in Shuryak et al. 2018, Table 1 [1]. Rows A–E assume a Monday start, with 6 fractions per week for the HYPNO control arm. For the DI model, onset, doubling time and daily loss are the same for every regimen (day 28, 3.5 days, 1.00 Gy~10~ per day), and regimens shorter than 28 days incur no DI repopulation. Abbreviation: fx, fractions.

| Regimen | d (Gy) | n | T (d) | DD onset (day) | DD doubling time (d) | DD loss (Gy~10~/day) | BED~DI~ / BED~DD~ (Gy~10~) | Plausibility of the DD values |
|---|---|---|---|---|---|---|---|---|
| A. Standard, 70 Gy/35 fx (reference) | 2.0 | 35 | 47 | 28.0 | 3.5 | 1.00 | 65.0 / 65.0 | Identical to DI by construction. Onset and doubling time plausible. Daily loss above the 0.65–0.77 Gy~10~ estimates, inherited from the calculator's DI default. |
| B. UK, 65 Gy/30 fx | 2.17 | 30 | 40 | 25.3 | 3.2 | 1.09 | 67.1 / 63.1 | Plausible onset and doubling time. |
| C. HYPNO control, 66 Gy/33 fx, 6 fx/week | 2.0 | 33 | 38 | 24.0 | 3.1 | 1.13 | 69.2 / 63.3 | Plausible onset and doubling time. |
| D. HYPNO test, 55 Gy/20 fx [4] | 2.75 | 20 | 26 | 18.6 | 2.5 | 1.39 | 70.1 / 59.8 | Onset at the early edge of clinical estimates. Doubling time shorter and daily loss about twice the clinical estimates. |
| E. HYHOPE, 46.5 Gy/15 fx [5] | 3.1 | 15 | 19 | 15.6 | 2.2 | 1.57 | 60.9 / 55.6 | Onset earlier than any clinical estimate. Fast doubling. Not testable here (de-escalation trial in favorable human papillomavirus (HPV)-positive disease). |
| F. DAHANCA 6&7, 6 fx/week | 2.0 | 33.5 | 39 | 24.3 | 3.1 | 1.12 | 69.4 / 63.8 | Plausible onset and doubling time. |
| G. CAIR, 7 fx/week | 1.9 | 40 | 40 | 22.1 | 2.9 | 1.21 | 78.4 / 68.8 | Plausible onset and doubling time. |
| H. GORTEC very accelerated, 2 fx/day | 2.0 | 31.5 | 22 | 14.6 | 2.1 | 1.65 | 75.6 / 63.3 | Onset in week 3 and doubling time near 2 days, both outside clinical estimates. |
| I. TROG accelerated, 2 fx/day | 1.8 | 33 | 24 | 17.1 | 2.4 | 1.47 | 70.1 / 60.0 | Onset slightly early, fast doubling. |
| J. CHART, 3 fx/day | 1.5 | 36 | 12 | 9.7 | 1.6 | 2.15 | 62.1 / 57.1 | Implausible: onset in the second week, doubling every 1.6 days, more than 2 Gy~10~ lost per day. |
| K. RTOG 9003 concomitant boost | 1.7 | 42 | 43 | 25.8 | 3.2 | 1.07 | 68.5 / 65.1 | Plausible. |
| L. RTOG 9003 hyperfractionation, 2 fx/day | 1.2 | 68 | 50 | 27.4 | 3.4 | 1.02 | 69.4 / 68.3 | Plausible. |
| M. EORTC 22791 hyperfractionation, 2 fx/day | 1.2 | 70 | 51 | 27.1 | 3.4 | 1.03 | 71.1 / 69.6 | Plausible. |

# Table 2. Predicted and observed differences from each trial's own control arm, DD rate 0.666 per day

Each experimental arm is compared with the conventional arm of the same trial, using that arm's average schedule from Shuryak et al. 2018, Table 1 [1]. Observed tumor control is the value tabulated there (2–5 years, as reported by each trial). For HYPNO, the published 3-year locoregional control differed by at most 1.4 percentage points between arms [4].

| Trial arm (patients, experimental/control) | Control arm schedule | ΔBED~DI~ (Gy~10~) | ΔBED~DD~ (Gy~10~) | Observed Δ tumor control (points) | Which model matches the direction |
|---|---|---|---|---|---|
| D. HYPNO test (395/397) | 66 Gy/33 fx, 38 d | +0.9 | −3.5 | within ±1.4 | DI closer |
| F. DAHANCA 6&7 (750/726) | 2.0 Gy × 33.5, 46 d | +7.0 | +0.5 | +7.0 | DI |
| G. CAIR (51/49) | 1.9 Gy × 40, 54 d | +14.0 | +1.2 | +45.7 | Both positive, DI closer |
| H. GORTEC very accelerated (137/129) | 2.0 Gy × 35, 48 d | +11.6 | −1.6 | +18.2 | DI (DD has the wrong sign) |
| I. TROG accelerated (172/171) | 2.0 Gy × 35, 50 d | +8.1 | −4.8 | +8.5 | DI (DD has the wrong sign) |
| J. CHART (552/366) | 2.0 Gy × 33, 45 d | −0.1 | −5.7 | −1.2 | DI closer |
| K. RTOG 9003 concomitant boost (268/268) | 2.0 Gy × 35, 50 d | +6.5 | +0.3 | +9.7 | DI |
| L. RTOG 9003 hyperfractionation (263/268) | 2.0 Gy × 35, 50 d | +7.4 | +3.5 | +10.0 | Both positive, DI closer |
| M. EORTC 22791 hyperfractionation (166/159) | 2.0 Gy × 35, 51 d | +10.1 | +4.8 | +1.0 | DD closer |

# Assessment of the matched setting

- **Schedules of 5–7 weeks.** With the matched rate, the DD model's onsets (22–28 days) and doubling times (2.9–3.5 days) fall within the clinical estimates (rows A–C, F, G, K–M). The DD and DI models differ by up to about 10 Gy~10~ even here, mainly because DD starts repopulation earlier when the daily cell kill is higher.
- **Schedules shorter than about 4 weeks.** The DD model predicts earlier onset (days 10–19) and faster repopulation (doubling every 1.6–2.5 days and 1.4–2.2 Gy~10~ lost per day, rows D, E and H–J). This is the model's central hypothesis rather than an error. However, no clinical estimate supports repopulation this early or this fast. The published estimates come mainly from longer schedules, so they neither confirm nor exclude it.
- **Direction of the trial results.** With these parameter values, the DD model predicts a tumor BED loss for two accelerated arms that improved tumor control (GORTEC +18.2 points, TROG +8.5 points), and a 3.5 Gy~10~ deficit for HYPNO, which showed no difference. The DI model matches the direction of eight of the nine comparisons. The exception is EORTC 22791, where it overpredicts the gain.
- **What this does and does not show.** This is a check of direction, not a fit. Tumor control was measured at different times and stages across trials, several trials are small (CAIR has 100 patients), and the average schedules are approximate. The 2018 analysis [1] found that the DD model described the same trials better than the DI model, but it fitted α, the rates, the onset parameter and stage-specific clonogen numbers to the data. The comparison here uses the calculator's default parameters, which were not fitted. It therefore tests these parameter values, not the DD model in general.
- **Implication for the calculator.** The matched rate gives clinically plausible DD behavior for conventional and moderately accelerated schedules. For short courses, the DD values should be presented as the hypothesis of a dose-dependent onset, not as an estimate supported by clinical data. The DI daily loss of 1.00 Gy~10~ remains above the published estimates, and both models inherit it at the reference.

# Table 3. Repopulation behavior of each regimen, current default DD rate (0.5 per day)

Same regimens, conventions and DI values as Table 1. With λ~DD~ = 0.5 per day, the calibration places DD repopulation in the reference regimen at day 21.7 rather than day 28, at a slower rate (0.150 rather than 0.200 per day), so that the total BED lost over the course still equals the DI value.

| Regimen | d (Gy) | n | T (d) | DD onset (day) | DD doubling time (d) | DD loss (Gy~10~/day) | BED~DI~ / BED~DD~ (Gy~10~) | Plausibility of the DD values |
|---|---|---|---|---|---|---|---|---|
| A. Standard, 70 Gy/35 fx (reference) | 2.0 | 35 | 47 | 21.7 | 4.6 | 0.75 | 65.0 / 65.0 | Onset, doubling time and daily loss all within clinical estimates. DD starts 6 days earlier than DI at a slower rate, with the same total loss. |
| B. UK, 65 Gy/30 fx | 2.17 | 30 | 40 | 19.6 | 4.2 | 0.82 | 67.1 / 62.5 | Onset slightly before day 21, otherwise plausible. |
| C. HYPNO control, 66 Gy/33 fx, 6 fx/week | 2.0 | 33 | 38 | 18.6 | 4.1 | 0.85 | 69.2 / 62.7 | Onset at the early edge of clinical estimates, otherwise plausible. |
| D. HYPNO test, 55 Gy/20 fx [4] | 2.75 | 20 | 26 | 14.4 | 3.3 | 1.04 | 70.1 / 58.0 | Onset in week 3, earlier than clinical estimates. Doubling time plausible, daily loss above estimates. |
| E. HYHOPE, 46.5 Gy/15 fx [5] | 3.1 | 15 | 19 | 12.1 | 2.9 | 1.18 | 60.9 / 52.8 | Onset in week 2, earlier than any clinical estimate. Doubling time plausible. |
| F. DAHANCA 6&7, 6 fx/week | 2.0 | 33.5 | 39 | 18.8 | 4.1 | 0.84 | 69.4 / 63.4 | Onset at the early edge, otherwise plausible. |
| G. CAIR, 7 fx/week | 1.9 | 40 | 40 | 17.2 | 3.8 | 0.91 | 78.4 / 69.7 | Onset slightly early, otherwise plausible. |
| H. GORTEC very accelerated, 2 fx/day | 2.0 | 31.5 | 22 | 11.3 | 2.8 | 1.24 | 75.6 / 62.3 | Onset in week 2, outside clinical estimates. Doubling time plausible. |
| I. TROG accelerated, 2 fx/day | 1.8 | 33 | 24 | 13.3 | 3.1 | 1.11 | 70.1 / 58.2 | Onset early, doubling time plausible. |
| J. CHART, 3 fx/day | 1.5 | 36 | 12 | 7.5 | 2.2 | 1.61 | 62.1 / 54.8 | Implausible onset (day 7.5) and a daily loss of 1.6 Gy~10~. |
| K. RTOG 9003 concomitant boost | 1.7 | 42 | 43 | 20.0 | 4.3 | 0.80 | 68.5 / 65.0 | Plausible, onset at the early edge. |
| L. RTOG 9003 hyperfractionation, 2 fx/day | 1.2 | 68 | 50 | 21.2 | 4.5 | 0.77 | 69.4 / 69.4 | Plausible. |
| M. EORTC 22791 hyperfractionation, 2 fx/day | 1.2 | 70 | 51 | 21.0 | 4.5 | 0.77 | 71.1 / 71.0 | Plausible. |

# Table 4. Predicted and observed differences from each trial's own control arm, DD rate 0.5 per day

Same comparisons as Table 2. The DI column does not depend on λ~DD~ and is repeated for convenience.

| Trial arm (patients, experimental/control) | ΔBED~DI~ (Gy~10~) | ΔBED~DD~ (Gy~10~) | Observed Δ tumor control (points) | Which model matches the direction |
|---|---|---|---|---|
| D. HYPNO test (395/397) | +0.9 | −4.7 | within ±1.4 | DI closer |
| F. DAHANCA 6&7 (750/726) | +7.0 | +0.5 | +7.0 | DI |
| G. CAIR (51/49) | +14.0 | +1.2 | +45.7 | Both positive, DI closer |
| H. GORTEC very accelerated (137/129) | +11.6 | −2.6 | +18.2 | DI (DD has the wrong sign) |
| I. TROG accelerated (172/171) | +8.1 | −6.6 | +8.5 | DI (DD has the wrong sign) |
| J. CHART (552/366) | −0.1 | −7.3 | −1.2 | DI closer |
| K. RTOG 9003 concomitant boost (268/268) | +6.5 | +0.2 | +9.7 | DI |
| L. RTOG 9003 hyperfractionation (263/268) | +7.4 | +4.6 | +10.0 | Both positive, DI closer |
| M. EORTC 22791 hyperfractionation (166/159) | +10.1 | +6.2 | +1.0 | DD closer |

# Comparison of the two DD settings

- **Where the implausibility sits.** The matched rate (Table 1) keeps onsets within clinical estimates for 5–7-week schedules but makes repopulation fast in short courses (doubling every 1.6–2.5 days). The default rate (Table 3) gives slower, more plausible repopulation (doubling every 2.2–4.6 days, 0.77–1.61 Gy~10~ lost per day, and 0.75 at the reference, close to the published 0.65–0.77). In exchange, onset moves earlier in every regimen: days 17–22 for 5–7-week schedules and days 7.5–14 for courses under 4 weeks. The lower rate trades an implausible rate for an implausible onset.
- **Effect on tumor BED.** The default rate gives a lower DD BED than the matched rate for every regimen except the reference (equal by calibration), CAIR and the two hyperfractionation arms (G, L and M), where the slower rate outweighs the earlier onset. The deficit grows with acceleration, for example HYPNO test 58.0 versus 59.8 Gy~10~ and CHART 54.8 versus 57.1 Gy~10~.
- **Direction of the trial results.** The pattern is unchanged and somewhat stronger: the DD model still predicts losses for the GORTEC and TROG accelerated arms, which improved tumor control, and a larger deficit for HYPNO (−4.7 versus −3.5 Gy~10~). Neither setting removes the mismatch for accelerated schedules. It comes from the structure of the DD model combined with the unfitted calculator parameters (α and the calibrated onset level), not from the choice of λ~DD~ alone.

# Limitations

- Schedules with two or three fractions per day are evaluated with the same formulas as once-daily ones. Incomplete repair between fractions given on the same day is not included, which overstates their BED slightly.
- The average fraction numbers and times from the 2018 compilation (for example 33.5 and 31.5 fractions) stand for trial arms whose individual schedules varied.
- Baseline tumor growth before onset, included in the 2018 models as a background growth term, is not included in the calculator.
- Concurrent chemotherapy (given in about 75% of HYPNO patients and in HYHOPE) is not part of any BED here.

# References

1. Shuryak I, Hall EJ, Brenner DJ. Dose dependence of accelerated repopulation in head and neck cancer: supporting evidence and clinical implications. Radiother Oncol. 2018;127(1):20-26. https://doi.org/10.1016/j.radonc.2018.02.015 (PMID 29534828)
2. Withers HR, Taylor JMG, Maciejewski B. The hazard of accelerated tumor clonogen repopulation during radiotherapy. Acta Oncol. 1988;27(2):131-146. https://doi.org/10.3109/02841868809090333 (PMID 3390344)
3. Fowler JF. 21 years of biologically effective dose. Br J Radiol. 2010;83(991):554-568. https://doi.org/10.1259/bjr/31372149 (PMID 20603408)
4. Bentzen SM, Gupta T, Jacinto AA, et al. Hypofractionated versus normofractionated accelerated radiation therapy with or without cisplatin for locally advanced head and neck squamous cell carcinoma (HYPNO): a randomized, open-label, phase 3, noninferiority trial. Int J Radiat Oncol Biol Phys. 2026;125(2):387-398. https://doi.org/10.1016/j.ijrobp.2025.12.056 (PMID 41711624)
5. Moon DH, Avkshtol V, Hughes RS, et al. HYHOPE: a phase 1 study of deintensified hypofractionated radiation therapy for human papillomavirus-associated oropharynx cancer. Int J Radiat Oncol Biol Phys. 2025;123(3):694-702. https://doi.org/10.1016/j.ijrobp.2025.05.009 (PMID 40383163)
