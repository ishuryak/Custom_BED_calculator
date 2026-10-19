### What this calculator does

For a radiotherapy regimen of $n$ fractions of $d$ Gy given over $T$ days, it computes three biologically effective doses (BEDs):

- the tumor BED for head and neck squamous cell carcinoma (HNSCC), under two models of tumor repopulation
- the BED for early (acute mucosal) toxicity
- the BED for late normal-tissue toxicity.

Each value is compared with a comparison (control) regimen, by default 35 × 2 Gy in 47 days. HYPNO's two arms can also be chosen as the comparison. Gy₁₀ and Gy₃ denote BED computed with α/β = 10 Gy and 3 Gy.

The tab "Candidate regimens (model screen)" lists shorter, hypofractionated regimens whose tumor BED change from the comparison meets a chosen minimum (which may be negative, as a non-inferiority margin), while the early and late toxicity BEDs stay within chosen limits. The "Sensitivity" tab shows how the tumor BED difference of the regimen being evaluated changes across tumor parameter sets.

This is a research and teaching tool, meant for hypothesis generation. It is not a treatment-planning system and must not be used on its own to choose treatment for a patient. Differences of a few Gy are within the uncertainty of the model parameters, and any regimen considered for a trial needs clinical precedent, organ-at-risk planning and toxicity stopping rules.

### Tumor BED, dose-independent (DI) repopulation model

Accelerated repopulation starts at a fixed time $T_k$ and proceeds at a fixed rate $\lambda_{DI}$ (the Withers "hockey stick"). It costs $\lambda_{DI}/\alpha$ Gy of tumor BED per day after $T_k$. An optional background growth rate $g$ (default 0) costs $g/\alpha$ Gy per day over the whole course:

$$BED_{DI} = n d \left(1 + \frac{d}{\alpha/\beta}\right) - \frac{\lambda_{DI}\,\max(0,\ T - T_k)}{\alpha} - \frac{g\,T}{\alpha}$$

### Tumor BED, dose-dependent (DD) repopulation model

Repopulation is a response to cell killing:

- It starts when the cumulative log cell kill reaches a constant $C$.
- Its rate is $\lambda_{DD}$ times the mean fraction of tumor cells killed per day.

With $K = n d (1 + d/(\alpha/\beta))$:

$$BED_{DD} = K - \frac{\lambda_{DD}\left(1 - e^{-\alpha K / T}\right)\max\left(0,\ T - \frac{C\,T}{\alpha K}\right)}{\alpha} - \frac{g\,T}{\alpha}$$

**Calibration at a fixed anchor.** $C$ is not entered by the user. It is solved so that $BED_{DD} = BED_{DI}$ for the **anchor regimen** (default 35 × 2 Gy in 47 days), for whatever tumor parameters are entered. The anchor is separate from the comparison regimen: choosing another comparison, such as HYPNO's control arm, changes the differences shown but not the DD model. Versions 1 and 2.0 calibrated at the comparison regimen, so changing the comparison changed the DD model itself. That behavior is available as a labeled legacy option. The DI parameters ($\lambda_{DI}$, $T_k$) still set the DD calibration, and the shared $g$ cancels in it.

For some parameter combinations no repopulation onset after the start of treatment can make the two models agree at the anchor. This happens, for example, when $\lambda_{DD}$ is too small to produce the DI repopulation. The DD model is then marked as unavailable and its values are hidden. Searches that require or rank by DD return no regimens with an explanation, and the default "DI or DD" criterion with ranking by the larger gain uses DI alone. With $\lambda_{DI}$ = 0, neither model has accelerated repopulation.

**How this DD model relates to the published one.** The equations follow Shuryak, Hall and Brenner [1] (Eq. 1–5). The parameters do not. In [1], each model was fitted separately to 33 trial arms (α/β 10 Gy, $g$ = 0.0154 per day), and the DD fit had its own α (0.224 Gy⁻¹), maximum rate (1.17 per day) and $C$ (14.5). Here, both models share one α and $C$ is calibrated to the DI model. The two versions differ in practice. With the fitted parameters of [1], the DD model places repopulation in the 35 × 2 Gy course from about day 36 and gives small tumor BED gains for the GORTEC (+3.2 Gy₁₀) and TROG (+2.1 Gy₁₀) accelerated arms and about no change for HYPNO (+0.5 Gy₁₀). With the calibrated version and this calculator's presets, the DD model predicts losses for all three. The fitted values are not offered as defaults. DD results here are therefore a calibrated variant of the dose-dependent hypothesis, not a reproduction of the fitted model.

### Early toxicity BED (acute oral and pharyngeal mucosa)

Mucosa repopulates quickly after an onset time $T_{k,e}$, at rate $\lambda_e = \ln 2 / T_{p,e}$:

$$BED_{early} = n d \left(1 + \frac{d}{(\alpha/\beta)_e}\right) - \frac{\lambda_e\,\max(0,\ T_e - T_{k,e})}{\alpha_e}$$

Here $T_e$ is the overall time in **elapsed days** (last treatment day minus first), which is the calendar count minus 1. Fowler's parameters and tolerance zone were derived on that clock [5, 9]: 70 Gy in 35 fractions over 7 weeks counts as 46 days and gives 53.1 Gy₁₀, as in his Table 1. The calculator converts once, inside the calculation.

### Late toxicity BED

Late-responding tissues are assumed not to repopulate during treatment, so $\alpha$ does not enter:

$$BED_{late} = n d \left(1 + \frac{d}{(\alpha/\beta)_l}\right)$$

It is evaluated at the prescription dose per fraction, not at the dose to individual organs at risk. The Results tab also shows it at α/β 2 and 4 Gy, because one α/β can reverse a comparison: 20 × 2.9 Gy is below 35 × 2 Gy at α/β 3 Gy (114.1 vs 116.7 Gy₃) but above it at α/β 2 Gy (142.1 vs 140.0). Equivalent doses in 2-Gy fractions are $EQD2 = BED / (1 + 2/(\alpha/\beta))$.

### How the regimen search works

- **Candidates.** Once-daily treatment, Monday to Friday, starting on a Monday. The number of fractions runs over the chosen range, below the comparison number. The dose per fraction runs over multiples of the chosen step (default 0.05 Gy) above the comparison dose per fraction, up to the chosen maximum (default 5 Gy). Rows above 3.1 Gy per fraction, the largest curative fraction size in the benchmark trials (HYHOPE), are flagged "beyond precedent".
- **Time.** The overall time must be at least the chosen number of days shorter than the comparison (default 7 days).
- **Tumor.** The tumor BED change from the comparison must be at least the chosen minimum (default +1 Gy) under the chosen model: DI or DD (default), both, DI only, or DD only. A negative minimum allows a loss, as in a non-inferiority design. The column "Credited by" says which model meets the criterion.
- **Early toxicity.** Choose one of four limits. The default is the early BED of HYPNO's 55 Gy in 20 fractions (55.9 Gy₁₀ with the default acute parameters), computed with the same acute model as the candidates. That regimen was given to 395 patients in a randomized trial, 75% with weekly cisplatin, with grade 3+ acute mucositis in 50.9%. The limit therefore means "modeled acute BED no higher than a regimen with known randomized acute toxicity", not "safe". The others are the comparison's early BED plus an allowed excess, the lower edge of the Fowler zone (59 Gy₁₀ for radiotherapy alone, or 54 Gy₁₀ with concurrent chemotherapy), and the Fenwick et al. boundary [15] on cell-kill BED (α/β 10 Gy, no repopulation term), $69.5\,(T/32.2)/\sin(T/32.2) - 3.5$ Gy₁₀ with $T$ in inclusive days and the angle in radians, defined only up to 50 days. It is an empirical boundary from historical radiotherapy-alone schedules: all schedules below it were tolerable in that data set, which is not a guarantee for new schedules, and chemotherapy is not included. It does not depend on Fowler's parameters. Meade et al. [16] advise careful phase 1 testing of a regimen that at least two of three acute models predict to be intolerable.
- **Why HYPNO is the default acute limit.** It is anchored to observed acute toxicity in a randomized trial in the target setting (chemoradiation for locally advanced HNSCC), whereas Fowler's 59 Gy₁₀ zone was derived for radiotherapy alone. Its limits: HYPNO's population (no HPV testing, 87% tobacco users, 12 centers), its weekly cisplatin 35 mg/m² (not 3-weekly 100 mg/m²) and its target volumes may differ from a new trial's; regimens much shorter than 26 days rely on the mucosal repopulation model for the time dependence; and the acute model cannot separate regimens 1–2 Gy₁₀ apart (see the HYHOPE example below). A candidate at this limit still needs acute-toxicity monitoring and stopping rules.
- **Late toxicity.** Choose the comparison's late BED plus an allowed excess (default 0), or the late BED of HYPNO's 55 Gy in 20 fractions (105.4 Gy₃ at α/β 3 Gy). The default comparison allows more: 70 Gy in 35 fractions is 116.7 Gy₃.
- **Validity.** A model whose comparison BED is negative cannot credit any regimen, and a regimen with a negative BED cannot qualify on that endpoint.
- **Output.** For each number of fractions, the dose per fraction that ranks highest is listed. The ranking is by the larger change of the two models (default, optimistic), the smaller change (conservative), or one model's change. The ZIP download also holds every qualifying regimen, not only the listed one, and a scenario file that reproduces the search.
- **Next step.** The column "Next dose step would exceed" names what one more dose step would break: early toxicity, late toxicity, or the tumor criterion (or a combination). "End of dose search range" means the maximum dose searched was reached. "None (a higher dose ranks lower)" means a higher dose would still meet every limit but ranks lower, which can happen under the DD model.

With the defaults (scenario D18 and the HYPNO acute limit), three regimens gain at least 1 Gy₁₀ under the DI model, led by 29 × 2.25 Gy in 39 days (ΔDI +1.2, ΔDD −2.4 Gy₁₀). None gains under the DD model. A shorter course at HYPNO's acute BED gives less tumor BED than 35 × 2 Gy under both models, so a non-inferiority margin (a negative minimum change) is the natural setting for hypofractionation questions. The list depends strongly on the tumor repopulation settings. An empty list is a result, not an error.

### Calculation details

- **Smooth threshold.** Every $\max(0, x)$ is computed as $\ln(1 + e^{3x})/3$, as in the original calculator. This differs from the sharp threshold by less than $\ln 2/3$ = 0.231 days. It changes a BED by at most about 0.23 days times the daily repopulation loss ($\lambda/\alpha$ Gy per day), most when the overall time is close to the onset time: with the defaults, at most 0.23 Gy for the DI BED and 0.18 Gy for the early BED. For the DD BED it is up to about 0.3 Gy with the defaults over the default search range (up to 5 Gy per fraction) and 0.4 Gy up to 8 Gy, and larger with a small $\alpha$ or a large $\lambda_{DD}$. The threshold type is recorded in every export.
- **Overall time.** When computed from the calendar, $T$ counts calendar days from the first to the last fraction inclusive, so 35 fractions take 47 days. The tumor models use this count, as in [1]. Tumor parameters from sources that count elapsed days are converted when entered as presets (Hartley: onset day 22 elapsed = day 23 inclusive). A typed overall time cannot be shorter than the calendar allows.
- **Units of repopulation estimates.** Published "Gy per day" values are often in 2-Gy-equivalent (EQD2) units. In tumor BED at α/β 10 Gy they are 1.2 times larger: Withers' 0.6 Gy/day is 0.72 Gy₁₀/day and the HYPNO design's 0.65 Gy/day is 0.78 Gy₁₀/day [2, 10]. The tumor panel shows the current loss in both units.

### Default parameters and sources

**Tumor (DI and DD models).**

- Default values (version 2.3, scenario D18): $\alpha$ = 0.35 Gy⁻¹, $\alpha/\beta$ = 10 Gy, $\lambda_{DI}$ = 0.23 day⁻¹ (doubling time 3.0 days), $\lambda_{DD}$ = 0.40 day⁻¹, $T_k$ = 28 days, $g$ = 0. The DI loss is 0.66 Gy₁₀ per day after day 28, at the low end of the source estimates below (0.66–0.84 Gy₁₀/day), which makes shortening the course worth less. Each value is within published ranges, but the set is not fitted to outcome data, and its DD predictions conflict with HYPNO (see Scenario D18 below). The version 1 defaults ($\alpha$ 0.2, $\lambda_{DI}$ 0.2, $\lambda_{DD}$ 0.5, a loss of 1.0 Gy₁₀ per day) remain available as a preset.
- Presets. Each sets only the values its source gives, so DI-only presets keep the DD maximum rate:
  - **Withers et al. [2]:** about 0.6 Gy per day in 2-Gy fractions after about 4 ± 1 weeks. The preset sets $\lambda_{DI} = \alpha \times 0.6\,(1 + 2/(\alpha/\beta))$ with $\alpha$ = 0.2 Gy⁻¹ and $T_k$ = 28 days, and follows the tumor α/β: 0.72 Gy₁₀/day at α/β 10 Gy, 0.84 Gy BED/day at 5 Gy.
  - **Hartley et al. [7]:** $\alpha$ = 0.3 Gy⁻¹, α/β 10 Gy, doubling time 3 days, onset day 22 elapsed (23 inclusive), 0.77 Gy₁₀/day.
  - **HYHOPE protocol [6]:** $\alpha$ = 0.35 Gy⁻¹, α/β 10 Gy, doubling time 3 days, $T_k$ = 21 days, 0.66 Gy₁₀/day.
  - **Scenario D18 (the default):** $\alpha$ = 0.35 Gy⁻¹, α/β 10 Gy, $\lambda_{DI}$ = 0.23 per day, $\lambda_{DD}$ = 0.40 per day, $T_k$ = 28 days. An analyst scenario chosen on 2026-10-18 so that DD repopulation rates and daily losses fall in ranges taken from DI-model studies. It is not a fitted set, and with it the DD model predicts a tumor BED deficit of 5.8 Gy₁₀ for HYPNO's test arm against its control, which HYPNO's result (non-inferior locoregional control) does not support. Treat its DD values as a pessimistic bound.
- Clinical estimates of the tumor α/β for head and neck cancer are about 10 Gy, with wide uncertainty: 10.5 (6.5–29) Gy from hyperfractionation trials with 1.1–2.55 Gy fractions [3], summarized in [4]. Use at larger fractions is an extrapolation. A lower estimate of 6.4 Gy (2.7–10.2) comes from an unpublished overview cited in the HYPNO rationale.

**Early toxicity (acute mucosa).**

- Default values: $(\alpha/\beta)_e$ = 10 Gy, $\alpha_e$ = 0.35 Gy⁻¹, $T_{k,e}$ = 7 elapsed days, doubling time $T_{p,e}$ = 2.5 days, as proposed by Fowler et al. [5]. The HYHOPE protocol uses the same mucosal values [6].
- With these values, Fowler et al. place the acute tolerance zone at 59–61 Gy₁₀ for radiation alone [5]. It is an empirical gap between schedules reported as tolerable and intolerable, not a probability, and Fowler's own Table 1 marks 70 Gy in 35 fractions (53.1 Gy₁₀) as "near tolerance?". Fowler later described a 59–63 Gy₁₀ "Grey Zone" [9]. For concurrent chemotherapy, Fowler et al. suggest a decrement of possibly 3–5 Gy₁₀ [5]. The zone shown with chemotherapy (54–58 Gy₁₀) spans that range, and the Fowler limit in the search uses its lower edge. These zones apply only with Fowler's parameter values. Fowler et al. also advise that a regimen approaching the zone should first be tested in a small dose-escalation study observing acute reactions [5].
- With chemotherapy, the 53.1 Gy₁₀ of 70 Gy in 35 fractions is about 1 Gy₁₀ below the lowered zone (54 Gy₁₀), not well below it.
- **Acute toxicity observed with 55 Gy in 20 fractions** (55.9 Gy₁₀ modeled). HYPNO, randomized, 75% with weekly cisplatin: grade 3+ acute mucositis 50.9% (190/373) vs 54.7% (208/380) after 66 Gy in 33 fractions [10]. Jacinto et al., 20 patients with cisplatin: grade 3 mucositis 40%, median actual treatment time 29 days [11]. Miller et al., retrospective, older patients, almost all without chemotherapy: grade 3 mucositis (RTOG) in 72 of 113 completers at the end of treatment, 26% four weeks later [14]. Madhava et al., with carboplatin: grade 3 mucositis in 16 of 19, lasting 4 weeks or more in 3 [13]. These cohorts differ in population, grading and time point, so their rates cannot be compared to estimate a chemotherapy effect.
- A modeled acute BED at or below the comparison does not guarantee tolerability. HYHOPE dropped its 12-fraction level (44.4 Gy at 4 fractions/week, with weekly cisplatin) after an apparent increase in grade 2–3 mucositis and pain in a parallel postoperative trial. The protocol lists its modeled acute BED as 52 Gy₁₀ (51.3 Gy₁₀ with its stated parameters), against 51.4 Gy₁₀ for the 15-fraction level that continued, both computed by the protocol on its inclusive clock (19 days) [6]. On this calculator's elapsed clock the 15-fraction level is 52.2 Gy₁₀.
- Other estimates exist: $\alpha_e$ = 0.3 Gy⁻¹ [7], and a later refit on 108 schedules [8].

**Late toxicity.**

- Default value: $(\alpha/\beta)_l$ = 3 Gy, the conventional value for late normal-tissue effects [6, 9].
- Estimates for head and neck late effects vary with the endpoint. For grade 2+ late effects in one hyperfractionation trial the estimate was 3.3–4.0 Gy (confidence limits 2.2–6.7), depending on the assumed dose-response slope [3]. Spinal cord α/β is lower, about 0.87 Gy, so a single late α/β cannot stand in for the cord.
- A small decrease in modeled late BED has not been shown to be a clinical benefit. In randomized altered-fractionation trials, a 15% decrease in late mucosal BED was needed before any sparing of late mucosal toxicity was observed [7].

**Clinical context.**

- HYPNO (792 patients) compared 55 Gy in 20 fractions at 5 fractions/week with 66 Gy in 33 fractions at 6 fractions/week, not with 70 Gy in 35 fractions. Weekly cisplatin 35 mg/m² was given to 75%, and HPV status was not tested. At 3 years, locoregional control was 50.7% vs 51.2% (HR 1.098, 90% CI 0.93–1.30) and grade 3+ late effects 18.8% vs 20.2% (HR 0.926, 0.68–1.25). The test arm passed a 10-percentage-point non-inferiority margin for both [10]. The margin does not support escalation beyond 55 Gy in 20 fractions, and the trial does not select between the DI and DD models.
- HYHOPE (phase 1, 24 patients): favorable-risk HPV-positive oropharynx cancer, reduced elective volumes, 3 weekly doses of cisplatin 40 mg/m², designed against 60 Gy in 30 fractions as a de-escalation. 46.5 Gy in 15 fractions was the maximum tolerated dose. Grade 3 oral mucositis 5/24, any acute grade 3 toxicity 7/24, no late grade 3+ at a median of 26.6 months, 2-year locoregional recurrence 13.8% (95% CI 3.2–32%) [12]. It is not a comparison for unselected HNSCC.
- Placing a predicted ΔBED in Gy beside an observed difference in percentage points does not test a model. A fitted dose-response relation is needed for that.

**References.**

1. Shuryak I, Hall EJ, Brenner DJ. Dose dependence of accelerated repopulation in head and neck cancer: supporting evidence and clinical implications. Radiother Oncol. 2018;127(1):20-26. https://doi.org/10.1016/j.radonc.2018.02.015
2. Withers HR, Taylor JM, Maciejewski B. The hazard of accelerated tumor clonogen repopulation during radiotherapy. Acta Oncol. 1988;27(2):131-146. https://doi.org/10.3109/02841868809090333
3. Stuschke M, Thames HD. Fractionation sensitivities and dose-control relations of head and neck carcinomas: analysis of the randomized hyperfractionation trials. Radiother Oncol. 1999;51(2):113-121. https://doi.org/10.1016/s0167-8140(99)00042-0
4. van Leeuwen CM, Oei AL, Crezee J, et al. The alfa and beta of tumours: a review of parameters of the linear-quadratic model, derived from clinical radiotherapy studies. Radiat Oncol. 2018;13:96. https://doi.org/10.1186/s13014-018-1040-z
5. Fowler JF, Harari PM, Leborgne F, Leborgne JH. Acute radiation reactions in oral and pharyngeal mucosa: tolerable levels in altered fractionation schedules. Radiother Oncol. 2003;69(2):161-168. https://doi.org/10.1016/s0167-8140(03)00231-7
6. Moon D (PI). HYHOPE: Phase I study of de-intensified hypofractionated radiation therapy for human papilloma virus-associated oropharynx cancer. Protocol STU 2020-1079, version 4.1 (2024). UT Southwestern Simmons Comprehensive Cancer Center. ClinicalTrials.gov NCT04580446. https://clinicaltrials.gov/study/NCT04580446
7. Hartley A, Sanghera P, Kazi W, et al. Correlation of currently used radiobiological parameters with local control and acute and late mucosal toxicity in randomised studies of altered fractionation for locally advanced head and neck cancer. Clin Oncol (R Coll Radiol). 2011;23(1):29-33. https://doi.org/10.1016/j.clon.2010.08.007
8. Pardo-Montero J, Parga-Pazos M, Fenwick JD. Classification of tolerable/intolerable mucosal toxicity of head-and-neck radiotherapy schedules with a biomathematical model of cell dynamics. Med Phys. 2021;48(7):4075-4084. https://doi.org/10.1002/mp.14834
9. Fowler JF. 21 years of biologically effective dose. Br J Radiol. 2010;83(991):554-568. https://doi.org/10.1259/bjr/31372149
10. Bentzen SM, Gupta T, Jacinto AA, et al. Hypofractionated versus normofractionated accelerated radiation therapy with or without cisplatin for locally advanced head and neck squamous cell carcinoma (HYPNO): a randomized, open-label, phase 3, noninferiority trial. Int J Radiat Oncol Biol Phys. 2026;125(2):387-398. https://doi.org/10.1016/j.ijrobp.2025.12.056
11. Jacinto AA, Batalha Filho ES, Viana LS, et al. Feasibility of concomitant cisplatin with hypofractionated radiotherapy for locally advanced head and neck squamous cell carcinoma. BMC Cancer. 2018;18:1026. https://doi.org/10.1186/s12885-018-4893-5
12. Moon DH, Avkshtol V, Hughes RS, et al. HYHOPE: a phase 1 study of deintensified hypofractionated radiation therapy for human papillomavirus-associated oropharynx cancer. Int J Radiat Oncol Biol Phys. 2025;123(3):694-702. https://doi.org/10.1016/j.ijrobp.2025.05.009
13. Madhava K, Hartley A, Wake M, et al. Carboplatin and hypofractionated accelerated radiotherapy: a dose escalation study of an outpatient chemoradiation schedule for squamous cell carcinoma of the head and neck. Clin Oncol (R Coll Radiol). 2006;18(1):77-81. https://doi.org/10.1016/j.clon.2005.08.014
14. Miller D, et al. A multicentre retrospective cohort study of hypofractionated radiotherapy in head and neck cancer. Clin Oncol (R Coll Radiol). 2026;59:104332. https://doi.org/10.1016/j.clon.2026.104332
15. Fenwick JD, Lawrence GP, Malik Z, Nahum AE, Mayles WP. Early mucosal reactions during and after head-and-neck radiotherapy: dependence of treatment tolerance on radiation dose and schedule duration. Int J Radiat Oncol Biol Phys. 2008;71(2):625-634. https://doi.org/10.1016/j.ijrobp.2008.02.008
16. Meade S, Sanghera P, Glaholm J, Hartley A. Models of acute mucosal tolerance to radiotherapy alone applied to synchronous chemoradiation schedules in head and neck cancer. Tumour Biol. 2014;35(3):2017-2023. https://doi.org/10.1007/s13277-013-1267-0

### Limitations

- **Concurrent chemotherapy.** The chemotherapy checkbox only lowers the displayed Fowler tolerance zone for acute mucositis, and the Fowler early limit when that option is used in the search. No BED in this calculator, tumor, early or late, includes an effect of chemotherapy, although chemotherapy increases both tumor control and acute toxicity. A limit relative to the comparison assumes the same chemotherapy in both regimens, but shorter courses usually deliver fewer cycles (HYPNO: 4 vs 5 weekly courses, 121 vs 149 mg/m² cumulative).
- **Screening surrogates.** The BED limits are screening surrogates, not validated safety constraints, and the acute model cannot separate regimens at the resolution the search ranks them (see the HYHOPE 12-fraction example).
- **Fowler's model** has no treatment volume or mucosal area, comes from pre-intensity-modulated data, assumes complete repair between fractions and applies to oral and pharyngeal mucosa only [5].
- Human papillomavirus (HPV) status of the tumor, tumor volume and stage are not modeled.
- Dose to individual organs at risk is not modeled. BEDs use the prescription dose per fraction.
- Two or more fractions per day are not supported. Incomplete repair between same-day fractions would raise the BED.
- Treatment gaps are not modeled.
- The tumor parameters are not fitted to outcome data in this calculator (see "How this DD model relates to the published one").

### Version

Version 2.3 of the Custom BED calculator. Version 1 remains available at https://ishuryak.shinyapps.io/custom_bed_calculator/. It still carries a calibration error at tumor α/β other than 10 Gy (at α/β 5 Gy, 35 × 2 Gy gives a DD BED of 73.8 Gy instead of 79.0).

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
