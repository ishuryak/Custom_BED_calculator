<p class="fitted-short"><strong>In short.</strong> The fitted models combine results of randomized trials with stated assumptions about plausible parameter values. The resulting parameter samples describe uncertainty under those assumptions, not the full uncertainty of the biology. Tumor control and acute mucositis were fitted separately. Late toxicity is not fitted: it uses an assumed α/β. The calculator's Refit support and the shares on the Uncertainty page say how often a model criterion is met across the samples. They are not a patient's probability of tumor control or of a side effect.</p>

{{CHECK}}

<h4 id="fitted-bayes">What "Bayesian" means here</h4>

A Bayesian fit combines two ingredients. The **prior** describes what was assumed about each parameter before the trial results were used, for example that tumor α/β most likely lies between {{AB_LO}} and {{AB_HI}} Gy. The **likelihood** describes how likely the observed trial results are under a given set of parameter values. The **posterior** combines both: it is the range of parameter values compatible with the assumptions and the trials together.

The posterior is represented by a large set of samples. Each sample is one complete, internally consistent set of parameter values. The calculator evaluates a schedule once per sample, so the spread of the results carries the uncertainty of all parameters together, including how they depend on each other.

When the trials say little about a parameter, its posterior stays close to its prior. This happens with tumor and mucosal α/β, because the fitted trials span a narrow range of dose per fraction. A Bayesian fit makes this visible instead of hiding it behind one fixed value, but the result for that parameter is then mainly an assumption.

<h4 id="fitted-what">What was fitted</h4>

- **Tumor control:** local or locoregional control at 2 to 5 years in the arms of randomized trials of altered fractionation without concurrent cytotoxic chemotherapy ({{TUMOR_ARMS}} arms or strata in {{TUMOR_UNITS}} fitting units). A fitting unit is a trial, or a stratum within one, with its own baseline.
- **Acute mucositis:** the number of patients with grade 3 or worse (or confluent) acute mucositis in {{EARLY_ARMS}} arms of {{EARLY_FAMILIES}} randomized trial families ({{EARLY_N}} evaluable patients). This is the research option "Trial-count fit". The default early model is Fowler 2003, which is not fitted here.
- **Late toxicity:** not fitted. The late BED uses α/β = {{LATE_AB}} Gy with complete repair and no time factor. On the Uncertainty page its α/β is drawn from an assumed lognormal range (median {{LATE_MED}} Gy, central 95% {{LATE_LO}}–{{LATE_HI}} Gy). This range is an analyst assumption, not an estimate from trial data.

In both fits, the information about dose and time comes mainly from differences between randomized arms of the same trial, because each trial (or unit) has its own baseline.

<h4 id="fitted-tumor">Tumor model</h4>

{{TUMOR_BRIEF}}

**How the model works.** The tumor fit uses a Bayesian binomial model with a Poisson tumor-control response. For arm $i$ of fitting unit $j$, the number of patients with tumor control follows a binomial distribution with probability $\mathrm{TCP}_i = \exp(-e^{\eta_i})$, where

$$\eta_i = c_j - \alpha\,(\mathrm{BED}_i - \overline{\mathrm{BED}}_j) + \lambda\,\max(0, T_i - T_k) + g\,(T_i - \overline{T}_j) + \sigma z_i .$$

Here $\mathrm{BED}_i$ is the linear-quadratic dose at the draw's α/β, $T_i$ the overall time in days, and the bars are the unit's means. Each fitting unit has a separate baseline $c_j$ with a fixed weak normal prior. The arm-level effects $\sigma z_i$, with an estimated common spread σ, allow more variation between arms than binomial sampling alone. This differs from the acute-mucositis model, which also estimates a common distribution for the trial baselines. The control counts are reconstructed by rounding the reported control rate times the number of patients analyzed. Most reported rates are Kaplan–Meier or cumulative-incidence estimates, so these counts are an approximation.

**The two tumor models are not independent fits.** The dose-independent (DI) model is fitted first. The dose-dependent (DD) model used by the calculator is calibrated to it, one DI sample at a time:

1. α, α/β and background growth $g$ are taken from the DI sample, not refitted.
2. The DD onset constant is solved so that DD and DI give the same regrowth for 35 × 2 Gy in 47 days. Agreement there is imposed, not found.
3. Only the DD regrowth rate, the unit baselines and the arm-level effects are then fitted, with a soft penalty that favors regrowth onset at or after day {{DD_ONSET}} in every arm (a preference, not a hard bound).
4. Each of the {{N_DI}} DI samples used has {{PER_DI}} DD samples fitted given it, giving {{TUMOR_ROWS}} paired rows. These are {{N_DI}} distinct DI samples, each repeated {{PER_DI}} times, not {{TUMOR_ROWS}} independent DI samples.

A DD model fitted on its own, without this calibration, did not sample reliably. It is not used by the calculator, and its prior on the onset constant does not apply here.

<details class="fitted-tech"><summary>Technical details: tumor priors, sampling and checks</summary>

**DI fit, priors (the values the fit was run with).**

{{TUMOR_DI_PRIORS}}

The α/β prior puts 95% of its weight between {{AB_LO}} and {{AB_HI}} Gy.

**Calibrated DD fit, priors and fixed quantities.**

{{TUMOR_DD_PRIORS}}

**Sampling.** {{TUMOR_SAMPLER}}

**Checks.** {{TUMOR_DIAG}}

</details>

**BED and EQD2.** Rates of tumor BED lost per day are given either in Gy of BED or in Gy EQD2 (equivalent dose in 2-Gy fractions). At α/β = 10 Gy, EQD2 = BED / 1.2, so 1.2 Gy of BED per day is 1 Gy EQD2 per day. The two are the same quantity on different scales, and a summary on one scale is not a different estimate from one on the other.

<h4 id="fitted-early">Acute mucositis model (research option)</h4>

{{EARLY_BRIEF}}

<details class="fitted-tech"><summary>Technical details: early-toxicity priors, sampling and checks</summary>

For arm $i$ of trial $j$ the model is $\mathrm{logit}\,p_i = \eta_j + \kappa\,(E_i - E_{j,\mathrm{ref}}) + \gamma\,\mathrm{chemo}_i + \sigma z_i$. Here $E$ is the early BED with incomplete repair and the time-loss term, and $E_{j,\mathrm{ref}}$ is that of the trial's reference arm. The trial baseline $\eta_j = m_0 + \delta\,\mathrm{CTCAE}_j + \tau u_j$ is drawn from a common distribution.

{{EARLY_PRIORS}}

**Fixed in this scenario.** Onset day {{EARLY_TK}} (elapsed days), repair half-time {{EARLY_H}} hours, exact time hinge.

**Sampling.** {{EARLY_SAMPLER}}

**Checks.** {{EARLY_DIAG}}

</details>

<h4 id="fitted-checks">How the fitting was checked</h4>

- **Markov chain Monte Carlo (MCMC)** is the method that draws the posterior samples. Several independent chains are run from different starting points.
- **R-hat** compares the chains with each other. Values near 1 (the limit here is {{RHAT_MAX}}) mean they agree.
- **Effective sample size (ESS)** is the number of independent samples the correlated chain samples are worth. The limit here is {{ESS_MIN}} for every checked quantity.
- **Divergences** are numerical warnings from the sampler that the posterior may not have been explored fully. None were allowed.
- **Posterior predictive checks** simulate trial results from the fitted model and compare them with the observed ones.

These checks show that the computation is adequate. They do not show that the model is biologically right or that it applies to shorter schedules.

<h4 id="fitted-validation">What the validation supports</h4>

- **In-sample replication.** {{PPC_IN}} of {{TUMOR_ARMS}} tumor arms fall inside their 90% posterior predictive interval. This uses each unit's own fitted baseline, so it has little power to detect misfit.
- **Simulated-data calibration** (early fit). Refitting simulated data sets checks that the fitting code recovers its own model. It does not test the model against reality.
- **Held-out trials.** Each trial (family) was left out in turn and predicted by a fit without it. For the early fit the results, including the trials that favored Fowler's parameters, are in the summary above. For the tumor fit they are in the full refit explainer.
- **Transport.** The tumor trials reach {{TUMOR_DMAX}} Gy per fraction. Larger fractions, and shorter courses than the trials used, are extrapolations.

<h4 id="fitted-app">How the samples are used in the calculator</h4>

- **Point values.** The main list uses one set of values: the posterior medians. For the tumor, this is the "trial-fit point scenario".
- **Refit support** (main list) is the share of the {{TUMOR_ROWS}} paired tumor rows in which the schedule meets the tumor criterion under the DI or the DD model. It uses the tumor samples only.
- **Uncertainty page.** Its tumor criterion is set on the page. The default requires both tumor models in the same paired row. Toxicity limits there use the early-model samples ({{EARLY_ROWS}} joint samples for the research option) and the assumed late range.
- **Pairing across endpoints.** Tumor and toxicity samples come from separate fits and are paired independently. This is an assumption, not a fitted joint distribution.

<h4 id="fitted-sources">Sources and reproducibility</h4>

{{SOURCES}}
