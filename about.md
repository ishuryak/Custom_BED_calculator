### What this calculator does

For a radiotherapy regimen of $n$ fractions of $d$ Gy given over $T$ days, it computes three biologically effective doses (BEDs):

- the tumor BED for head and neck squamous cell carcinoma (HNSCC), under two models of tumor repopulation
- the BED for early (acute mucosal) toxicity
- the BED for late normal-tissue toxicity.

Each value is compared with a reference regimen (default 35 × 2 Gy in 47 days). Gy₁₀ and Gy₃ denote BED computed with α/β = 10 Gy and 3 Gy.

The tab "Regimens meeting the BED criteria" then lists shorter, hypofractionated regimens whose tumor BED gain over the reference meets a chosen minimum, while the early and late toxicity BEDs stay within chosen limits.

This is a research and teaching tool, meant for hypothesis generation. It is not a treatment-planning system and must not be used on its own to choose treatment for a patient. Differences of a few Gy are within the uncertainty of the model parameters, and any regimen considered for a trial needs clinical precedent and toxicity stopping rules.

### Tumor BED, dose-independent (DI) repopulation model

Accelerated repopulation starts at a fixed time $T_k$ and proceeds at a fixed rate $\lambda_{DI}$ (the Withers "hockey stick"). It costs $\lambda_{DI}/\alpha$ Gy of tumor BED per day after $T_k$:

$$BED_{DI} = n d \left(1 + \frac{d}{\alpha/\beta}\right) - \frac{\lambda_{DI}\,\max(0,\ T - T_k)}{\alpha}$$

### Tumor BED, dose-dependent (DD) repopulation model

Repopulation is a response to cell killing:

- It starts when the cumulative log cell kill reaches a constant $C$.
- Its rate is $\lambda_{DD}$ times the mean fraction of tumor cells killed per day.

With $K = n d (1 + d/(\alpha/\beta))$:

$$BED_{DD} = K - \frac{\lambda_{DD}\left(1 - e^{-\alpha K / T}\right)\max\left(0,\ T - \frac{C\,T}{\alpha K}\right)}{\alpha}$$

$C$ is not entered by the user. It is solved so that $BED_{DD} = BED_{DI}$ for the reference regimen, for whatever tumor parameters are entered. The DI parameters ($\lambda_{DI}$, $T_k$) therefore also change the DD results.

The original calculator fixed this calibration at $\alpha/\beta$ = 10 Gy. This version recalculates it for any $\alpha/\beta$ and any reference regimen.

For some parameter combinations no repopulation onset after the start of treatment can make the two models agree at the reference. This happens, for example, when $\lambda_{DD}$ is too small to produce the DI repopulation. The DD model is then marked as unavailable and its values are hidden. Searches that require or rank by DD ("both models", "DD model", or ranking by the smaller gain or by the DD gain) return no regimens with an explanation, and the default "DI or DD" criterion with ranking by the larger gain uses DI alone. With $\lambda_{DI}$ = 0, neither model has repopulation and both equal the plain linear-quadratic (LQ) BED.

Both models follow Shuryak, Hall and Brenner [1], without the slow background tumor-growth term ($\gamma T$ in [1]). Including it would credit shorter regimens with a small additional tumor BED gain. With the growth rate of [1] (0.0154 per day) and $\alpha$ = 0.2 Gy⁻¹, this is about 0.08 Gy per day saved.

### Early toxicity BED (acute oral and pharyngeal mucosa)

Mucosa repopulates quickly after an onset time $T_{k,e}$, at rate $\lambda_e = \ln 2 / T_{p,e}$:

$$BED_{early} = n d \left(1 + \frac{d}{(\alpha/\beta)_e}\right) - \frac{\lambda_e\,\max(0,\ T - T_{k,e})}{\alpha_e}$$

### Late toxicity BED

Late-responding tissues are assumed not to repopulate during treatment, so $\alpha$ does not enter:

$$BED_{late} = n d \left(1 + \frac{d}{(\alpha/\beta)_l}\right)$$

It is evaluated at the prescription dose per fraction, not at the dose to individual organs at risk. Equivalent doses in 2-Gy fractions are $EQD2 = BED / (1 + 2/(\alpha/\beta))$.

### How the regimen search works

- **Candidates.** Once-daily treatment, Monday to Friday, starting on a Monday. The number of fractions runs over the chosen range, below the reference number. The dose per fraction runs over multiples of the chosen step (default 0.05 Gy) that are above the reference dose per fraction, up to the chosen maximum (default 8 Gy).
- **Time.** The overall time must be at least the chosen number of days shorter than the reference (default 7 days).
- **Tumor.** The tumor BED gain over the reference must be at least the chosen minimum (default 1 Gy) under the chosen model: DI or DD (default), both, DI only, or DD only.
- **Early toxicity.** By default, early BED may not exceed the reference's plus an allowed excess (default 0). Alternatively, it may not exceed the lower edge of the Fowler tolerance zone: 59 Gy₁₀ for radiotherapy alone, or 54 Gy₁₀ with concurrent chemotherapy.
- **Late toxicity.** Late BED may not exceed the reference's plus an allowed excess (default 0).
- **Output.** For each number of fractions, the dose per fraction that ranks highest is listed. The ranking is by the larger gain of the two models (default), the smaller gain, or one model's gain.
- **Next step.** The column "Next dose step would exceed" names what one more dose step would break: early toxicity, late toxicity, or the tumor criterion (or a combination). "End of dose search range" means the maximum dose searched was reached. "None (a higher dose ranks lower)" means a higher dose would still meet every limit but ranks lower, which can happen under the DD model.

With the default tumor parameters, no regimen gains tumor BED under the DD model while keeping early BED at or below the reference, so the listed regimens qualify under the DI model. The list depends strongly on the tumor repopulation settings. With the Withers or HYHOPE presets, no regimen qualifies at the default limits.

### Calculation details

- **Smooth threshold.** Every $\max(0, x)$ is computed as $\ln(1 + e^{3x})/3$, as in the original calculator. This differs from the sharp threshold by less than 0.24 days.
- **Effect of the smooth threshold.** It changes a BED by at most about 0.23 days times the daily repopulation loss ($\lambda/\alpha$ Gy per day). With the defaults this is at most 0.23 Gy for the DI BED and 0.18 Gy for the early BED, and largest when the overall time is close to the onset time. For the DD BED it is up to about 0.4 Gy with the defaults over the search range, and larger with a small $\alpha$ or a large $\lambda_{DD}$.
- **Overall time.** When computed from the calendar, $T$ counts calendar days from the first to the last fraction inclusive, so 35 fractions take 47 days. This matches the convention of [1] and of the HYHOPE protocol [6]. A typed overall time cannot be shorter than the calendar allows.

### Default parameters and sources

**Tumor (DI and DD models).**

- Default values: $\alpha$ = 0.2 Gy⁻¹, $\alpha/\beta$ = 10 Gy, $\lambda_{DI}$ = 0.2 day⁻¹, $\lambda_{DD}$ = 0.5 day⁻¹, $T_k$ = 28 days. These are the defaults of version 1 of this calculator and are meant to be changed by the user. They imply a DI repopulation loss of 1 Gy of tumor BED per day after day 28.
- The model structure follows [1]. The fixed-onset idea follows Withers et al. [2].
- Presets:
  - **Withers et al. [2]:** a dose increment of about 0.6 Gy per day, in 2-Gy fractions, is needed to compensate repopulation that starts after about 4 ± 1 weeks. In BED with $\alpha/\beta$ = 10 Gy this is 0.6 × 1.2 = 0.72 Gy₁₀ per day, so the preset sets $\lambda_{DI}/\alpha$ = 0.72 Gy/day with $\alpha$ = 0.2 Gy⁻¹ and $T_k$ = 28 days.
  - **HYHOPE protocol [6]:** tumor $\alpha$ = 0.35 Gy⁻¹, doubling time 3 days ($\lambda_{DI} = \ln 2/3$ per day) and $T_k$ = 21 days.
- Clinical estimates of the tumor $\alpha/\beta$ for head and neck cancer are about 10 Gy, with wide uncertainty: 10.5 (6.5–29) Gy [3], summarized in [4].

**Early toxicity (acute mucosa).**

- Default values: $(\alpha/\beta)_e$ = 10 Gy, $\alpha_e$ = 0.35 Gy⁻¹, $T_{k,e}$ = 7 days, doubling time $T_{p,e}$ = 2.5 days, as proposed by Fowler et al. [5]. The HYHOPE protocol uses the same mucosal values [6].
- With these values, Fowler et al. place the acute tolerance zone at 59–61 Gy₁₀ for radiation alone [5]. Fowler later described a 59–63 Gy₁₀ "Grey Zone" [9]. For concurrent chemotherapy, Fowler et al. state that more information about the decrement, possibly 3–5 Gy₁₀, is required [5]. The zone shown with chemotherapy (54–58 Gy₁₀) spans both ends of that range, and the Fowler limit in the search uses its lower edge. These zones apply only with Fowler's parameter values.
- The standard 35 × 2 Gy regimen gives about 52 Gy₁₀, so the default limit "early BED ≤ reference" is stricter than that tolerance zone. The HYPNO regimen of 55 Gy in 20 fractions gives about 55 Gy₁₀. It was delivered in a randomized phase 3 trial with non-inferior grade 3+ late effects [10]. Its acute burden is substantial with concurrent chemotherapy: grade 3 mucositis in 40% of 20 patients with cisplatin [11], and in 84% of 19 patients with carboplatin in a series cited in the HYHOPE protocol [6].
- A modeled acute BED at or below the reference does not guarantee tolerability. HYHOPE dropped its 12-fraction level (44.4 Gy at 4 fractions/week, with weekly cisplatin) after an apparent increase in grade 2–3 mucositis and pain in a parallel postoperative trial. The protocol lists its modeled acute BED as 52 Gy₁₀ (51.3 Gy₁₀ with its stated parameters), against 51.4 Gy₁₀ for the 15-fraction level that continued, computed the same way [6]. The two levels had essentially the same modeled acute BED but different observed toxicity.
- Other estimates exist: $\alpha_e$ = 0.3 Gy⁻¹ [7], and a later refit on 108 schedules [8].

**Late toxicity.**

- Default value: $(\alpha/\beta)_l$ = 3 Gy, the conventional value for late normal-tissue effects [6, 9].
- Estimates for head and neck late effects vary with the endpoint. For grade 2+ late effects in a hyperfractionation trial the estimate was 4.0 (3.3–5.0) Gy [3].
- A small decrease in modeled late BED has not been shown to be a clinical benefit. In randomized altered-fractionation trials, a 15% decrease in late mucosal BED was needed before any sparing of late mucosal toxicity was observed [7].

**Clinical context.**

- In the randomized HYPNO trial (792 patients), 55 Gy in 20 fractions at 5 fractions/week was non-inferior to 66 Gy in 33 fractions at 6 fractions/week, for both locoregional control and grade 3+ late effects. Concurrent weekly cisplatin (35 mg/m²) was permitted [10]. Both arms are shown on the Results tab.

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

### Not modeled

- the effect of concurrent chemotherapy on tumor control or on any BED (only the Fowler tolerance zone is adjusted)
- human papillomavirus (HPV) status of the tumor
- tumor volume and stage
- dose to individual organs at risk
- two or more fractions per day, where incomplete repair between fractions would raise late-tissue BED
- treatment gaps.

### Version

Version 2 of the Custom BED calculator. Version 1 remains available at https://ishuryak.shinyapps.io/custom_bed_calculator/.

Igor Shuryak, MD, PhD, is144@cumc.columbia.edu
