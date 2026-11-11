# Early-toxicity plan, step 5: first refit to randomized-arm mucositis counts

Written 2026-11-11. This is the first fit of the plan to the real outcomes. It gives research results
only: the app is unchanged, and Fowler 2003 stays the default unless the step 6 validation shows a
clear gain. Nothing is pushed or deployed.

## 1. Approvals and the step 2 spot-check

Igor Shuryak, 2026-11-11: "C1–C3 as recommended, go (and do the spot-check". This is read as
agreement to start step 5, approval of C1–C3 as recommended in the step 4 record, and a request
that Claude do the step 2 spot-check against the PDFs.

**Spot-check by Claude against the PDF text layers (not by Igor).** All five items in section 4
of the step 3 proposal pass:

| Item | PDF result | Stored | Verdict |
|---|---|---|---|
| Fu 2000 (RTOG 9003), Table 3, mucous membrane grade 3 + 4 | SF 67 + 0, HFX 109 + 1, AFX-S 109 + 3, AFX-C 122 + 1 | 67/268, 110/263, 112/274, 123/268 | Pass |
| Ghoshal 2008, Table 4 against the Discussion | Table 4 conventional arm: grade 1 45, grade 2 70, grade 3 27 (19%), grade 4 0, total 142. Discussion: "19% (23/142)" | 27/142 | Pass. **Conflict resolved:** the table grades sum to 142 only with 27, and 27/142 = 19.0% matches the printed percentage, while 23/142 = 16.2%. The Discussion's 23 is a typographical error |
| Skladowski 2000 (CAIR), Table 4 | CAIR 17 + 28 + 3 = 48 of 50 (III–IV), control 22 + 13 + 0 = 35 of 49. Dische > 15: 31 and 13 | 48/50 and 35/49 (primary), 31/50 and 13/49 (secondary) | Pass |
| Olmi 2003 (ORO 93-01), Table 7 mucosa | A: 9, 10, 27, 7, 1 (sum 54), B: 6, 6, 19, 21, 0 (52), C: 6, 3, 17, 22, 2 (50) | 8/54, 21/52, 24/50 | Pass. Grade 3 + 4 per arm, the row sums equal the arm sizes |
| Tandon 2018, Table 2 | Mucositis "Grade 3": 11 (36.7) control, 16 (53.3) study | 11/30, 16/30 | Pass |

## 2. Decisions applied

- **C1, chemotherapy within trials.** One term γ ~ N(0, 1) on the log-odds scale. Each arm is
  coded against its trial's reference arm by **cumulative planned concurrent chemotherapy dose**:
  +1 more, −1 less, 0 the same. The codes are:

  | Arm | Code | Basis |
  |---|---:|---|
  | GORTEC 99-02 accelerated + chemo | −1 | 2 cycles against 3 |
  | GORTEC 99-02 very accelerated | −1 | none against 3 cycles |
  | RTOG 0129 accelerated | −1 | 2 cisplatin cycles against 3 |
  | Vienna V-CHART + MMC | +1 | adds mitomycin C |
  | ORO 93-01 conventional + chemo | +1 | adds carboplatin/5-FU |
  | Bartelink 2002 multiple daily fractions | −1 | about 150 against 210 mg/m² in total (10 mg/m² on 15 days, 6 mg/m² on 35 days). **Ambiguous:** the daily dose is higher |

  A committed test checks that the code is nonzero exactly where an arm's planned therapy differs
  from its reference arm's. Sensitivity analyses: the 5 families dropped (U08), only a whole
  regimen added or removed counted (U09), and no chemotherapy term (U14).
- **C2, onset Tk.** Scenarios of 7 (primary), 14 and 21 days.
- **C3, overdispersion.** Every unit is fitted twice, binomial and with an arm-level random effect
  σ_arm ~ HalfNormal(0.5). The binomial primary fit fails its posterior predictive check
  (chi-square discrepancy, p = 0.0025), so **the overdispersed set is primary**, as C3 specifies.
- **Unchanged from steps 1–4:**
  - D3 priors: L* = 10κ ~ HalfNormal(1), α/β ~ lognormal(log 10, 0.354), r ~ lognormal(log 0.8,
    0.5), m0 ~ N(0, 2.5), τ0 ~ HalfNormal(1)
  - M1 scale-group offset δ ~ N(0, 1), with G2 = CTCAE v4 and "other" for Pinto
  - M2 W1 sensitivity and M4 sensitivity trials
  - 3.2 h repair half-time and the elapsed clock
- **New priors registered** in `sim_provenance.yaml`: `et_prior_gamma`, `et_prior_sigma_arm`,
  `et_prior_delta`, `et_sbc_reps`. The registry validates.

## 3. Model and runs

`tools/early_tox/step5_model.stan`:
y_i ~ Binomial(N_i, p_i), logit p_i = η_trial + κ·ΔE_i + γ·chemo_i + σ_arm·z_i, with
ΔE = ΔD + ΔQ/(α/β) − r·ΔT_eff against the trial's reference arm, and η_j = m0 + δ·group_j + τ0·z_j.
The benchmark (U01) fixes α/β = 10 and r = 0.79217 (Fowler 2003). The candidate (U02 and others)
samples α/β and r.

`tools/early_tox/step5_fit.py`:
- Revision `f47cd4d`, CmdStan 2.40.0, seed 20261111.
- 4 chains × 1,000 warmup + 1,000 draws, adapt_delta 0.95, one chain at a time, one thread,
  `nice -n 19`. A Windows free-memory check runs before each unit. About 4 s per fit.
- 12 units × 2 likelihoods, plus 15 leave-one-family-out units × 2.
- Outputs are in `2026_11_10_early_toxicity/step5/{overdispersed,binomial}/<unit>/`: summary,
  diagnostics, per-arm predictive check, and draws for the primary units. All numbers are
  collected in `step5/results_table.csv`.

**Sampler gate** (no divergences, R-hat < 1.01, bulk and tail ESS > 400, for the shared
parameters, the baselines and the decision quantities):
- All 27 overdispersed fits pass.
- 3 binomial fits miss narrowly (R-hat 1.0103 for U02 and U03, 1.0136 without RTOG 9003). The
  binomial set is not primary.

**Simulation-based calibration** (`tools/early_tox/step5_sbc.py`, seed 20261112): 500 replicates
on the real design, priors drawn independently in numpy.
- Ranks are uniform for every parameter: chi-square p = 0.076 (L*), 0.063 (α/β), 0.600 (r),
  0.546 (γ), 0.265 (σ_arm), 0.350 (m0), 0.387 (τ0).
- 17 of 500 replicates had a divergence, all with prior draws near the extremes. The real-data fits
  had none.

## 4. Results (overdispersed set, median and 90% interval)

| Unit | L* (per 10 BED) | α/β (Gy) | r (BED/day) | γ (chemo) | WAIC |
|---|---|---|---|---|---:|
| U01 benchmark, Fowler fixed | 0.72 (0.44–1.04) | 10 fixed | 0.79 fixed | 0.50 (0.07–0.94) | 222.4 |
| **U02 candidate, Tk 7 (primary)** | **0.47 (0.27–0.75)** | **12.2 (7.4–21.0)** | **1.51 (1.00–2.41)** | **0.53 (0.15–0.90)** | **218.6** |
| U03 Tk 14 | identical to U02 (as step 4 predicted) | | | | 218.6 |
| U04 Tk 21 | 0.46 (0.25–0.74) | 12.2 (7.3–20.5) | 1.60 (1.03–2.71) | 0.55 (0.16–0.93) | 219.2 |
| U05 L* ~ HalfNormal(0.5) | 0.44 (0.25–0.69) | 12.1 (7.1–20.8) | 1.57 (1.06–2.61) | 0.51 (0.15–0.89) | 217.8 |
| U06 L* ~ HalfNormal(2) | 0.48 (0.27–0.77) | 12.3 (7.2–20.5) | 1.49 (0.98–2.43) | 0.52 (0.13–0.91) | 218.3 |
| U07 α/β prior SD doubled | 0.48 (0.28–0.78) | **17.8 (7.1–48.2)** | 1.48 (0.98–2.41) | 0.52 (0.14–0.90) | 218.7 |
| U08 5 chemo families dropped | 0.51 (0.26–0.84) | 11.3 (6.6–19.6) | **1.08 (0.58–1.91)** | n/a | (22 arms) |
| U09 chemo on/off coding | 0.44 (0.25–0.74) | 12.0 (7.1–21.0) | 1.55 (0.99–2.58) | 0.42 (−0.09–0.94) | 218.2 |
| U10 W1 windows only | 0.47 (0.24–0.84) | 11.5 (6.7–19.3) | 1.48 (0.91–2.52) | 0.45 (−0.03–0.88) | (23 arms) |
| U11 plus IAEA-ACC, Pinto, Poddar | 0.48 (0.27–0.76) | 12.5 (7.3–21.8) | 1.50 (1.00–2.42) | 0.54 (0.18–0.89) | (41 arms) |
| U14 no chemo term | 0.43 (0.23–0.74) | 11.8 (7.1–20.5) | 1.49 (0.94–2.51) | n/a | 220.6 |

Arm-level SD σ_arm is 0.29 (0.17–0.47) in U02.

**Leave-one-family-out (U02 settings, 15 fits):**
- r median ranges 1.30–1.58: lowest without Vienna, highest without PMH.
- γ ranges 0.41–0.70, L* 0.44–0.58, α/β 11.4–12.3.
- No single trial drives the result.

**Decision quantities: odds ratio for grade ≥3 mucositis against 70 Gy / 35** (5 days a week, one
fraction a day):

| Schedule | U01 Fowler benchmark | U02 candidate | U08 chemo families dropped |
|---|---|---|---|
| 15 × 3.5 Gy | 1.93 (1.49–2.57) | 3.79 (2.25–6.73) | 2.37 (1.16–5.19) |
| 20 × 2.75 Gy | 1.22 (1.13–1.33) | 2.27 (1.50–3.49) | 1.56 (0.91–2.83) |
| 25 × 2.4 Gy | 1.11 (1.07–1.17) | 1.70 (1.29–2.27) | 1.32 (0.92–1.99) |
| 30 × 2.2 Gy | 1.16 (1.10–1.24) | 1.39 (1.22–1.59) | 1.23 (1.03–1.50) |

**What the results say, and what they do not:**
1. **Slope and chemotherapy are estimated from the data.** L* is about 0.5 log-odds per 10 BED
   units. Planned concurrent chemotherapy adds about 0.5 log-odds (odds ratio about 1.7) at the
   same schedule.
2. **The time effect r is about twice Fowler's value** (1.5 against 0.79 BED/day). It is stable
   across Tk, the slope priors, the windows, the M4 trials and every single-trial deletion. It falls
   to 1.08 (0.58–1.91), with Fowler's value inside the interval, only when all five
   chemotherapy-differing trials are removed together. Those trials include the two shortest arms
   (Vienna V-CHART at 16 elapsed days and GORTEC very accelerated at 23), which carry most of the
   information on time.
3. **α/β is set by the prior**, as step 4 predicted (U07 widens it to 7–48 Gy). No empirical α/β
   claim is made (stopping rule: nonidentification blocks an empirical claim).
4. **Clinical implication, if it holds in step 6:** shortened courses are more mucosally toxic
   relative to 70/35 than the Fowler limit implies. For 20 × 2.75 Gy, the odds ratio is 2.27
   against 1.22.
5. **WAIC favors the candidate by 3.9 over the benchmark.** That is an in-sample comparison, not
   the step 6 validation.

## 5. Historical classifier (step 5.3)

`tools/early_tox/step5_historical.py`:
- Inputs: 103 complete-course Pardo-Montero records (13 intolerable) and the 107-record
  sensitivity set (15). Partial courses are refused.
- Model: a logistic regression on one exposure.
- Cohort structure: leave-one-source-out cross-validation over 61 sources, with Tk chosen inside
  each training fold.

| Exposure | Held-out log score (primary) | Against Fowler | Sources better / worse |
|---|---:|---:|---|
| Intercept only | −40.18 | | |
| Fowler 2003, Tk 7 | −25.90 | | 49 / 12 against intercept only |
| Alternative time boundary (Tk 0–35 chosen in fold, 3 days on full data, 2–11 across folds) | −26.03 | −0.13 | 51 / 10 |
| Step 5 candidate point (α/β 12.2, r 1.51, 3.2 h repair) | −30.08 | −4.18 | 22 / 39 |

The sensitivity set gives the same ordering (−28.44, −28.55, −32.74). On these labels, Fowler's
exposure reaches a 50% predicted probability of intolerability at 64.8 Gy₁₀.

**This conflicts with the trial fit, and the reason is not settled.**
- **Not an independent test.** Fowler 2003 tuned its tolerance zone on largely the same historical
  schedules, so Fowler holds an in-sample advantage on these labels.
- **Different endpoint.** The labels are a binary tolerable or intolerable judgment of a whole
  schedule (mostly radiotherapy alone), not grade ≥3 rates in randomized arms.
- **A point, not the posterior.** Only the candidate's posterior medians were used, without its
  uncertainty.
- **What it does show.** The historical labels do not support a larger time effect.
- **Before step 6:** the historical-label crosswalk (step 2 gate, still open) is needed, so that
  both data sets can be held out together.

## 6. Tests

All mutation checks were confirmed red, then restored.
- **`tests/python/test_early_tox_step5.py`** (8 tests). Each was turned red by a targeted
  mutation:
  - Olmi chemotherapy code removed
  - reference contrast not subtracted
  - "other" scale group dropped
  - decision Q term zeroed
  - W1 filter broken
  - γ sign flipped in Stan
  - repopulation sign flipped in Stan
  - overdispersion switch not passed
  - skip-if-done disabled
  - family drop ignored
  - SBC α/β prior SD doubled
- **`tests/python/test_early_tox_step5_historical.py`** (6 tests). Mutations confirmed red:
  - partial-course check removed
  - Fowler exposure with repair
  - fold fit on all data
  - folds inverted
  - Tk chosen with the held-out source
- **Order:** tests, then the combined smoke (both likelihoods, all units), then tests again, then
  the full run on the committed revision.
- **Suite:** `python3 -m pytest -q tests/python` gave 259 passed (245 before step 5).
- **Stale-constant audit:** `stale_constant_audit.py --project . --project-root . --strict` gave
  0 errors, after registering the Fowler acute α/β (`fowler_acute_ab`). A changed value was
  confirmed to fail it.

## 7. Next (step 6), and what is still open

1. **Crosswalk.** Build the historical-label to trial-family crosswalk (step 2 gate), so CHART
   and other shared evidence is held out together.
2. **Validation.** Leave-family-out validation of U01 against U02 (new-trial and conditional
   contrast tasks), with the Bayesian bootstrap of family weights (two seeds), leave-one-family
   influence, and HYPNO as the R6 benchmark.
3. **Unresolved schedules (from step 2).** EORTC 22851 course 2, and PMH 26 vs 28 days.
4. **Tension.** The conflict between the trial fit (r ≈ 1.5) and the historical labels (Fowler
   better) is reported in step 6, not resolved by choosing one.
