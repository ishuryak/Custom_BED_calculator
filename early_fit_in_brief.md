### How the early-toxicity model was fitted

This is the research option, **Trial-count fit**. The default early-toxicity model is Fowler 2003.

**Data.** The fit uses the number of patients with grade 3 or worse acute mucositis in each arm of 22 randomized trials of altered fractionation in head and neck cancer. Where a trial's scale names confluent mucositis as its top grade, that count is used instead. The 22 trials have 51 arms and 6,042 evaluable patients. Every count was extracted twice, independently, from the full text, and a fitted count was used only when both extractions agreed. A sample was also checked against the original tables. The list "Trials used for the early-toxicity fit" shows each trial, its reference, arms, schedules, concurrent therapy, endpoint and counts.

**Model.** The model is a Bayesian hierarchical logistic regression, fitted in Stan. It has four parts:
- Each trial has its own baseline risk, drawn from a common distribution, with a shift for trials graded on the symptom-based CTCAE version 4 scale.
- Within a trial, the log-odds of mucositis rise with the difference in early BED between an arm and that trial's control arm. Early BED is E = D + Q/(α/β) − r × max(0, T − Tk), with T in elapsed days, where Q allows incomplete repair between fractions given on the same day.
- A separate term accounts for differences in planned concurrent chemotherapy between arms.
- An arm-level random effect allows more spread between arms than the binomial distribution alone, which did not fit the data.

**What was estimated.** The model estimates the mucosal α/β, the regrowth rate r, the slope and the chemotherapy effect. Some quantities were fixed instead:
- **Onset of regrowth:** fixed at day 7, because every trial arm lasted longer than about two weeks, so the data cannot locate it.
- **Repair half-time:** fixed at 3.2 hours. Fits at 2 hours, 4 hours and complete repair gave nearly the same results.

**Priors (chosen in advance).**
- Mucosal α/β: lognormal centered on 10 Gy (95% about 5 to 20 Gy).
- Regrowth rate: lognormal centered on 0.8 Gy of BED per day.
- Slope: half-normal(1) per 10 Gy of BED.
- Chemotherapy effect: normal(0, 1) on the log-odds scale.
- Between-trial spread of baselines: half-normal(1).
- Arm-level spread: half-normal(0.5).

**Checks.** All parameters passed the sampler diagnostics: no divergent transitions, R-hat below 1.01, and effective sample sizes above 400 with complete schema checks. Simulation-based calibration with 500 simulated data sets showed no detectable departure.

**Results (posterior medians, 90% intervals).**
- **Regrowth rate:** 1.55 Gy of BED per day (1.07 to 2.45), about twice Fowler's 0.79.
- **Mucosal α/β:** 12.9 Gy (7.6 to 21.8). This is set mostly by its prior.
- **Concurrent chemotherapy:** at the same schedule it multiplied the odds of severe mucositis by about 1.7.

**Validation.** Each trial was left out in turn, and the model, refitted without it, predicted it.
- **Overall:** the fitted model did better than Fowler's parameters in 18 of 22 trials, by about 0.28 log-score units per trial (95% interval 0.10 to 0.55).
- **Where the advantage came from:** mostly accelerated trials. For shorter hypofractionated courses the evidence is mixed. In HYPNO (55 Gy in 20 fractions against 66 Gy in 33), the comparison depends on how its 4 vs 5 weekly cisplatin courses are treated. ELAN-RT, in frail older patients with a split course, favored Fowler's parameters.
- **Historical tolerability labels:** these also favored Fowler's parameters, although Fowler's parameters were derived from largely the same schedules.
- **Assessment timing:** few trials scored mucositis after radiotherapy ended, when it often peaks.

**How the calculator uses it.** The point values are the posterior medians, and Uncertainty uses 4,000 joint posterior draws. The +5 limit is in this model's own BED. The model gives no patient risk, and Fowler 2003 remains the default.
