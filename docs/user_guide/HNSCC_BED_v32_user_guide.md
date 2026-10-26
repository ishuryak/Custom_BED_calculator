---
title: "HNSCC BED calculator, version 3.2: user guide"
---

Igor Shuryak, MD, PhD

is144@cumc.columbia.edu

26 October 2026. Live app: https://ishuryak.shinyapps.io/custom_bed_calculator_v3/

## 1. What the calculator does

The calculator compares radiotherapy schedules for head and neck squamous cell carcinoma (HNSCC). For an evaluated regimen and one comparison regimen it computes the biologically effective dose (BED) of the tumor, the acute oral and pharyngeal mucosa and late-responding tissue, and checks organ-at-risk (OAR) doses against published constraints. Tumor BED is computed under two repopulation models from Shuryak, Hall and Brenner [1]. In the dose-independent (DI) model, accelerated repopulation starts on a fixed day at a fixed rate. In the dose-dependent (DD) model, onset and rate depend on how fast cells are killed. Equivalent dose in 2-Gy fractions (EQD2) is BED/(1 + 2/(α/β)).

It is a tool for generating trial hypotheses. Its outputs are model quantities, not predictions of tumor control or toxicity for a patient, and not treatment recommendations.

Every tab uses the same evaluated regimen, comparison regimen and tumor threshold, shown in a banner at the top. The banner also flags a dose per fraction outside the range of the tumor data (1.15–2.55 Gy).

## 2. Sidebar menus

| Menu | Item | What it sets | Default |
|---|---|---|---|
| Regimen to evaluate | Number of fractions, dose per fraction | The evaluated regimen | 18 × 2.9 Gy (52.2 Gy), an example schedule, not a published regimen |
| | Fractions per week, first treatment day | Calendar (5 or 6 per week) | 5, Monday |
| | Days without treatment | Holidays or gaps as calendar day numbers, which lengthen overall time | none |
| | Compute overall time from this calendar | Overall time from the calendar, or typed by hand | on (24 days) |
| Comparison and tumor threshold | Compare with | The comparison used by every tab: a custom regimen, HYPNO's test arm (55 Gy in 20 fractions) or HYPNO's control arm (66 Gy in 33 fractions, 6 per week) [2] | custom 35 × 2 Gy in 47 days |
| | Tumor threshold | Minimum tumor BED change from the comparison that counts as acceptable (negative values allow a loss) | +1 Gy |
| Tumor (DI and DD models) | Parameter preset | Sets α, α/β, repopulation rates and onset from a source: version 1 values, the Withers rule [3], Hartley et al. [4], the HYHOPE protocol [5], Scenario D18, or the refit medians (section 4) | Scenario D18 (an analyst scenario, not fitted) |
| | α, α/β, DI rate, DD maximum rate, DI onset, background growth g | The individual tumor parameters | D18: α 0.35 per Gy, α/β 10 Gy, DI rate 0.23 per day, DD rate 0.40 per day, onset day 28, g 0 |
| | DD model | Calibrated to DI (default) or values entered directly | calibrated |
| DD calibration anchor | Anchor fractions and dose | The regimen at which DD is calibrated to give the same BED as DI | 35 × 2 Gy |
| Early toxicity (acute mucosa) | α/β, α, doubling time, onset (elapsed days) | Acute mucosal BED parameters of Fowler et al. [6] | 10 Gy, 0.35 per Gy, 2.5 days, day 7 |
| | Concurrent chemotherapy | Lowers only the Fowler tolerance zone used as a limit, by 3–5 Gy₁₀, and changes no BED | off |
| Late toxicity | α/β, late-responding tissue | Late BED at the prescription dose (a coarse screen, not an organ dose) | 3 Gy |

## 3. Tabs

**Report card.** The evaluated regimen against the comparison on one page: tumor BED under DI and DD, acute and late BED, their differences, and the refit's median and 80% range of the tumor BED difference. Organ rows appear when doses are entered on the Organs at risk tab, each with its status, its result over the cord α/β range and, for mean-dose limits, the note that the source limit comes from conventional fractionation. Notes below flag extrapolation and state that organ checks cover selected metrics only.

**Results.** Full numbers for both regimens: tumor BED (DI and DD), acute BED, late BED at α/β 2, 3 and 4 Gy, and EQD2. It also shows the evaluated regimen against the Fowler acute tolerance zones (59–61 Gy₁₀ for radiotherapy alone, 54–58 Gy₁₀ with chemotherapy) [6] and the Fenwick et al. boundary on cell-kill BED [7], and a table of published regimens recomputed with the current parameters: HYPNO's two arms [2] and HYHOPE's 46.5 Gy in 15 fractions [8].

**Compare regimens.** The evaluated regimen, three further regimens (default 17 × 3.05, 15 × 3.1 and 30 × 2 Gy) and the published benchmarks side by side, with fraction sizes above the benchmark range flagged.

**Search.** Lists shorter once-daily regimens that meet the tumor threshold and the chosen toxicity limits. Controls: which model the threshold applies to (both models by default), the ranking (the smaller of the two changes by default), the range of fraction numbers, minimum days saved, the dose grid, and the early and late toxicity limits. The default early limit is the acute BED of HYPNO's 55 Gy in 20 fractions, a regimen with randomized toxicity data. An empty list is a valid result. The results and the full scenario can be downloaded as a ZIP file that replays exactly.

**Both-model screen.** Regimens that meet the threshold under DI and DD in each of 15 parameter settings, using the global comparison and threshold unless overridden on the tab (overrides are highlighted). Agreement between the models is not two independent confirmations, because they share α and are calibrated to agree at the anchor.

**Sensitivity.** The tumor BED difference of the evaluated regimen under each preset and α/β, and the acute and late differences. The settings are not a probability sample.

**Uncertainty.** Applies 4000 parameter draws from the refit (section 5) to the evaluated regimen and the comparison together. It reports medians and 80% and 95% ranges of the tumor BED difference under each model, the share of draws meeting the threshold and every toxicity limit, and a model-conditional tumor control illustration from an entered control rate for the comparison. Controls: refit variant (primary, or α/β and g fixed as in the 2018 paper), tumor criterion, toxicity parameters (prior distributions or the sidebar values), the comparison's control rate and the selection rule for an optional search over regimens. The analysis downloads as a replayable ZIP file.

**Organs at risk.** Organ doses entered as percentages of the prescription, as Gy, or from an uploaded cumulative dose-volume histogram (DVH; files with patient identifiers are refused). Maximum-dose constraints are converted to the evaluated fractionation with the linear-quadratic model. Mean-dose constraints are compared on physical dose, as stated in their sources, with a flag that the effect of changed fractionation is not validated. Two libraries: QUANTEC [9, 10, 11, 12, 13] and the HYHOPE protocol tables for 12, 15 and 20 fractions [5]. Statuses are within, within 5%, exceeds, indeterminate at the DVH resolution, or not assessed, and the cord and brainstem are also checked over α/β 0.87–2 Gy. The limits can optionally restrict the searches. The table downloads as a CSV file.

**Model and sources.** Equations, parameter sources, the refit description, limitations, version and the deployed revision.

## 4. Model fitting

**Data.** 45 arms of 19 randomized trials of definitive fractionation in HNSCC, mainly without cytotoxic chemotherapy (8,279 patients), with local or locoregional control re-extracted from each publication (table below) and checked against the full texts. Arms with mixed fraction sizes enter as the sum of their fractions. Overall time is the delivered time where reported. Within some arms the population is mixed (for example nasopharyngeal and sinus cases in RTOG 83-13, nimorazole in DAHANCA 6/7), and twice-daily arms are fitted without a correction for incomplete repair.

| Trial | Publication | Arms | Patients |
|---|---|---|---|
| CAIR | Skladowski 2006 [14] | 4 | 100 |
| CAIR vs concomitant boost | Skladowski 2013 [15] | 2 | 343 |
| CHART | Saunders 2010 [16] | 2 | 915 |
| DAHANCA 6/7 | Overgaard 2003 [17] | 2 | 1,476 |
| EORTC 22791 | Horiot 1992 [18] | 2 | 325 |
| EORTC 22811 | Van 1995 [19] | 2 | 331 |
| EORTC 22851 | Horiot 1997 [20] | 2 | 512 |
| Ghoshal (India) | Ghoshal 2008 [21] | 2 | 285 |
| Hliniak (Poland) | Hliniak 2002 [22] | 2 | 381 |
| IAEA-ACC | Overgaard 2010 [23] | 2 | 900 |
| Jackson (Vancouver) | Jackson 1997 [24] | 2 | 82 |
| KROG 02-01 | Moon 2014 [25] | 2 | 156 |
| Princess Margaret Hospital | O'Sullivan 2020 [26] | 2 | 331 |
| RTOG 83-13 | Cox 1990 [27] | 3 | 237 |
| RTOG 88-09 | Fu 1995 [28] | 2 | 70 |
| RTOG 90-03 | Fu 2000 [29] | 4 | 1,073 |
| RTOG 95-12 | Trotti 2014 [30] | 2 | 239 |
| TROG 91.01 | Poulsen 2001 [31] | 2 | 343 |
| Yamazaki (Osaka) | Yamazaki 2006 [32] | 4 | 180 |

**Model.** The DI model of Shuryak, Hall and Brenner [1] was refitted in Stan, with the linear-quadratic cell kill K = Σd + Σd²/(α/β):

log(−log TCP) = c(trial) − α (K − K̄(trial)) + λ s(T − T_k) + g (T − T̄(trial)) + σz

Here TCP is the arm's tumor control probability and c(trial) is an intercept for each of the 21 within-trial comparisons, so the randomized contrasts inform the parameters. s is a smooth form of max(0, x), T is overall time in inclusive days and T_k the repopulation onset. g is background growth and σz an arm-level random effect. The likelihood is binomial on the reported control rate times the number of patients, an approximation because most rates are Kaplan–Meier estimates. Priors: α/β lognormal with 95% between 5 and 20 Gy, α and λ broad lognormal, onset normal (28 ± 10 days), g uniform 0.0046–0.021 per day, intercepts normal (sd 3). Four chains of 2000 draws after 1500 warm-up gave 8000 DI draws.

**DD model.** A DD model fitted on its own did not sample reliably. The app therefore uses DD calibrated to DI: for each of 1000 evenly spaced DI draws, the DD onset constant is set so that both models give the same repopulation for 35 × 2 Gy in 47 days, and only the DD rate, intercepts and overdispersion are then fitted given that draw (a modular analysis in which DD does not reweight DI). A soft penalty discourages DD onsets before day 7.

**Convergence.** Every fit passes one gate: no divergent transitions, R-hat at most 1.01, bulk and tail effective sample sizes of at least 400 and an energy Bayesian fraction of missing information (E-BFMI) of at least 0.3. Version 3.2 refitted with longer chains the 748 conditional DD fits that had not met the effective-sample-size criterion.

**Results.** In Gy at each draw's α/β, the DI repopulation loss λ/α is 0.63 per day (95% interval 0.21–1.19) after an onset near day 37 (24–43). α is about 0.04 per Gy, an effective population-level slope, not a cellular radiosensitivity. Tumor α/β (median 10.6 Gy) and background growth are not informed by these data. The DD rate is 0.41 per day (0.12–1.04). For HYPNO's 55 Gy in 20 fractions against 35 × 2 Gy, the median tumor BED change is −2.0 Gy under DI (80% range −4.3 to +1.0) and −3.3 Gy under DD (−5.9 to −0.5). Eight sensitivity analyses (added trials, Kaplan–Meier rates only, an updated RTOG 9003 report, halved patient numbers, planned instead of delivered overall time, a wider intercept prior) moved the DI median between −2.7 and +0.3 Gy and the DD median between −3.7 and −2.1 Gy. The planned overall time moved the DI median to about zero. Every arm lies within its 90% predictive interval, an in-sample check with little power.

## 5. Uncertainty estimation

**Tumor parameters.** The Uncertainty tab uses 4000 draws: 1000 DI posterior draws, each paired with 4 DD draws fitted on that DI draw. Both regimens are evaluated with the same draw, so their difference carries the correlations between parameters. A second variant fixes α/β at 10 Gy and g at 0.0154 per day, as in the 2018 paper.

**Toxicity parameters.** By default, acute, late and organ parameters are drawn from stated prior distributions, independently of the tumor draws: acute α/β lognormal around 10 Gy (95% 5–20), mucosal α uniform 0.30–0.40 per Gy, doubling time uniform 2–4 days, onset uniform 7–19 elapsed days (from Fowler [6] and the Pardo-Montero refit [33]), late α/β lognormal around 3 Gy (95% 2–4.5) [34], cord and brainstem α/β 0.87–2 Gy [9, 2] and other organs 2–4.5 Gy. The sidebar values can be used instead.

**Summaries.** For each model the tab reports the median and the 80% and 95% ranges of the tumor BED difference, and the share of draws meeting the threshold and every toxicity limit, counted draw by draw. Resampling the 1000 DI draws gives a Monte Carlo standard error of about 0.1 Gy for medians, up to 0.25 Gy for 95% limits and about 1.5 percentage points for shares, so BED summaries are rounded to 0.5 Gy.

**Tumor control illustration.** With TCP = exp(−K S), a BED difference Δ changes tumor control from a stated comparison value TCP₀ to TCP₀^exp(−αΔ), using each draw's α. This depends on the entered TCP₀ and on the model, and does not establish efficacy or non-inferiority.

**Search under uncertainty.** For each number of fractions, the dose with the best value of the chosen rule (the highest 10th percentile of the smaller of ΔDI and ΔDD, or the highest share meeting all criteria) is screened on 1000 draws among doses whose share within all toxicity limits reaches the required level (90% by default), then re-checked on all 4000 draws and shown with its grid neighbors.

**What is not included.** The intervals cover parameter uncertainty within this model and data. They do not cover model choice (DI versus DD beyond the pair shown), endpoint and population differences between the trials and a new trial, chemotherapy, HPV status, target volumes or the binomial approximation.

## References

1. Shuryak I, Hall EJ, Brenner DJ. Dose dependence of accelerated repopulation in head and neck cancer: supporting evidence and clinical implications. Radiother Oncol. 2018;127(1):20-26. https://doi.org/10.1016/j.radonc.2018.02.015
2. Bentzen SM, Gupta T, Jacinto AA, et al. Hypofractionated versus normofractionated accelerated radiation therapy with or without cisplatin for locally advanced head and neck squamous cell carcinoma (HYPNO): a randomized, open-label, phase 3, noninferiority trial. Int J Radiat Oncol Biol Phys. 2026;125(2):387-398. https://doi.org/10.1016/j.ijrobp.2025.12.056
3. Withers HR, Taylor JM, Maciejewski B. The hazard of accelerated tumor clonogen repopulation during radiotherapy. Acta Oncol. 1988;27(2):131-146. https://doi.org/10.3109/02841868809090333
4. Hartley A, Sanghera P, Kazi W, et al. Correlation of currently used radiobiological parameters with local control and acute and late mucosal toxicity in randomised studies of altered fractionation for locally advanced head and neck cancer. Clin Oncol (R Coll Radiol). 2011;23(1):29-33. https://doi.org/10.1016/j.clon.2010.08.007
5. Moon D (PI). HYHOPE: Phase I study of de-intensified hypofractionated radiation therapy for human papilloma virus-associated oropharynx cancer. Protocol STU 2020-1079, version 4.2 (amendment of 2024). UT Southwestern Simmons Comprehensive Cancer Center. ClinicalTrials.gov NCT04580446. https://clinicaltrials.gov/study/NCT04580446
6. Fowler JF, Harari PM, Leborgne F, Leborgne JH. Acute radiation reactions in oral and pharyngeal mucosa: tolerable levels in altered fractionation schedules. Radiother Oncol. 2003;69(2):161-168. https://doi.org/10.1016/s0167-8140(03)00231-7
7. Fenwick JD, Lawrence GP, Malik Z, Nahum AE, Mayles WP. Early mucosal reactions during and after head-and-neck radiotherapy: dependence of treatment tolerance on radiation dose and schedule duration. Int J Radiat Oncol Biol Phys. 2008;71(2):625-634. https://doi.org/10.1016/j.ijrobp.2008.02.008
8. Moon DH, Avkshtol V, Hughes RS, et al. HYHOPE: a phase 1 study of deintensified hypofractionated radiation therapy for human papillomavirus-associated oropharynx cancer. Int J Radiat Oncol Biol Phys. 2025;123(3):694-702. https://doi.org/10.1016/j.ijrobp.2025.05.009
9. Kirkpatrick JP, van der Kogel AJ, Schultheiss TE. Radiation dose-volume effects in the spinal cord. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S42-S49. https://doi.org/10.1016/j.ijrobp.2009.04.095
10. Mayo C, Yorke E, Merchant TE. Radiation associated brainstem injury. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S36-S41. https://doi.org/10.1016/j.ijrobp.2009.08.078
11. Deasy JO, Moiseenko V, Marks L, Chao KS, Nam J, Eisbruch A. Radiotherapy dose-volume effects on salivary gland function. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S58-S63. https://doi.org/10.1016/j.ijrobp.2009.06.090
12. Rancati T, Schwarz M, Allen AM, Feng F, Popovtzer A, Mittal B, Eisbruch A. Radiation dose-volume effects in the larynx and pharynx. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S64-S69. https://doi.org/10.1016/j.ijrobp.2009.03.079
13. Bhandare N, Jackson A, Eisbruch A, Pan CC, Flickinger JC, Antonelli P, Mendenhall WM. Radiation therapy and hearing loss. Int J Radiat Oncol Biol Phys. 2010;76(3 Suppl):S50-S57. https://doi.org/10.1016/j.ijrobp.2009.04.096
14. Skladowski K, et al. Continuous accelerated 7-days-a-week radiotherapy for head-and-neck cancer: long-term results of phase III clinical trial. Int J Radiat Oncol Biol Phys. 2006. https://doi.org/10.1016/j.ijrobp.2006.05.026
15. Skladowski K, et al. Radiation-free weekend rescued! Continuous accelerated irradiation of 7-days per week is equal to accelerated fractionation with concomitant boost of 7 fractions in 5-days per week: report on phase 3 clinical trial in head-and-neck cancer patients. Int J Radiat Oncol Biol Phys. 2013. https://doi.org/10.1016/j.ijrobp.2012.06.037
16. Saunders MI, et al. Mature results of a randomized trial of accelerated hyperfractionated versus conventional radiotherapy in head-and-neck cancer. Int J Radiat Oncol Biol Phys. 2010. https://doi.org/10.1016/j.ijrobp.2009.04.082
17. Overgaard J, et al. Five compared with six fractions per week of conventional radiotherapy of squamous-cell carcinoma of head and neck: DAHANCA 6 and 7 randomised controlled trial. Lancet. 2003. https://doi.org/10.1016/s0140-6736(03)14361-9
18. Horiot JC, et al. Hyperfractionation versus conventional fractionation in oropharyngeal carcinoma: final analysis of a randomized trial of the EORTC cooperative group of radiotherapy. Radiother Oncol. 1992. https://doi.org/10.1016/0167-8140(92)90242-m
19. Van den Bogaert W, et al. The EORTC randomized trial on three fractions per day and misonidazole (trial no. 22811) in advanced head and neck cancer: long-term results and side effects. Radiother Oncol. 1995. https://doi.org/10.1016/0167-8140(95)01538-r
20. Horiot JC, et al. Accelerated fractionation (AF) compared to conventional fractionation (CF) improves loco-regional control in the radiotherapy of advanced head and neck cancers: results of the EORTC 22851 randomized trial. Radiother Oncol. 1997. https://doi.org/10.1016/s0167-8140(97)00079-0
21. Ghoshal S, et al. Concomitant boost radiotherapy compared with conventional radiotherapy in squamous cell carcinoma of the head and neck--a phase III trial from a single institution in India. Clin Oncol (R Coll Radiol). 2008. https://doi.org/10.1016/j.clon.2008.01.011
22. Hliniak A, et al. A multicentre randomized/controlled trial of a conventional versus modestly accelerated radiotherapy in the laryngeal cancer: influence of a 1 week shortening overall time. Radiother Oncol. 2002. https://doi.org/10.1016/s0167-8140(01)00494-7
23. Overgaard J, et al. Five versus six fractions of radiotherapy per week for squamous-cell carcinoma of the head and neck (IAEA-ACC study): a randomised, multicentre trial. Lancet Oncol. 2010. https://doi.org/10.1016/s1470-2045(10)70072-3
24. Jackson SM, et al. A randomised trial of accelerated versus conventional radiotherapy in head and neck cancer. Radiother Oncol. 1997. https://doi.org/10.1016/s0167-8140(97)01944-0
25. Moon SH, et al. A prospective randomized trial comparing hypofractionation with conventional fractionation radiotherapy for T1-2 glottic squamous cell carcinomas: results of a Korean Radiation Oncology Group (KROG-0201) study. Radiother Oncol. 2014. https://doi.org/10.1016/j.radonc.2013.09.016
26. O'Sullivan B, et al. Durable therapeutic gain despite competing mortality in long-term follow-up of a randomized hyperfractionated radiotherapy trial for locally advanced head and neck cancer. Clin Transl Radiat Oncol. 2020. https://doi.org/10.1016/j.ctro.2020.01.003
27. Cox JD, et al. Dose-response for local control with hyperfractionated radiation therapy in advanced carcinomas of the upper aerodigestive tracts: preliminary report of radiation therapy oncology group protocol 83-13. Int J Radiat Oncol Biol Phys. 1990. https://doi.org/10.1016/0360-3016(90)90054-n
28. Fu KK, et al. Randomized phase I/II trial of two variants of accelerated fractionated radiotherapy regimens for advanced head and neck cancer: results of RTOG 88-09. Int J Radiat Oncol Biol Phys. 1995. https://doi.org/10.1016/0360-3016(95)00078-d
29. Fu KK, et al. A Radiation Therapy Oncology Group (RTOG) phase III randomized study to compare hyperfractionation and two variants of accelerated fractionation to standard fractionation radiotherapy for head and neck squamous cell carcinomas: first report of RTOG 9003. Int J Radiat Oncol Biol Phys. 2000. https://doi.org/10.1016/s0360-3016(00)00663-5
30. Trotti A, et al. Randomized trial of hyperfractionation versus conventional fractionation in T2 squamous cell carcinoma of the vocal cord (RTOG 9512). Int J Radiat Oncol Biol Phys. 2014. https://doi.org/10.1016/j.ijrobp.2014.04.041
31. Poulsen MG, et al. A randomised trial of accelerated and conventional radiotherapy for stage III and IV squamous carcinoma of the head and neck: a Trans-Tasman Radiation Oncology Group Study. Radiother Oncol. 2001. https://doi.org/10.1016/s0167-8140(01)00347-4
32. Yamazaki H, et al. Radiotherapy for early glottic carcinoma (T1N0M0): results of prospective randomized study of radiation fraction size and overall treatment time. Int J Radiat Oncol Biol Phys. 2006. https://doi.org/10.1016/j.ijrobp.2005.06.014
33. Pardo-Montero J, Parga-Pazos M, Fenwick JD. Classification of tolerable/intolerable mucosal toxicity of head-and-neck radiotherapy schedules with a biomathematical model of cell dynamics. Med Phys. 2021;48(7):4075-4084. https://doi.org/10.1002/mp.14834
34. Stuschke M, Thames HD. Fractionation sensitivities and dose-control relations of head and neck carcinomas: analysis of the randomized hyperfractionation trials. Radiother Oncol. 1999;51(2):113-121. https://doi.org/10.1016/s0167-8140(99)00042-0
