# Plan: early (acute mucosal) toxicity model with data-based uncertainty

Date: 2026-11-10. Status: PLAN ONLY, nothing in the app has changed. Applies to version 4.5.0
(revision `faf6817`, live at https://ishuryak.shinyapps.io/custom_bed_calculator_v4/).
Exploratory scripts behind every number below: `2026_11_10_early_toxicity/exploratory/`
(`binding.R`, `parse.py`, `fit.py`, `implic.py`. Run order and commands in section 9).

## 1. Why this matters

With the default settings (trial-fit tumor profile, comparison 35 × 2 Gy in 47 days, margins
tumor −5, early +5, late +5), the early-toxicity limit sets the highest allowed dose per fraction
for every course from 5 to 29 fractions. It is also the reason no course of 12 or fewer fractions
is listed: at 12 fractions the tumor condition needs at least 3.90 Gy and the early limit allows
at most 3.85 Gy (`binding.R`). Only at 30 fractions does the late limit bind instead.

| Fractions | Days | Tumor needs d ≥ | Early allows d ≤ | Late allows d ≤ |
|---|---|---|---|---|
| 12 | 16 | 3.90 | 3.85 | 4.20 |
| 13 | 17 | 3.65 | 3.65 | 4.00 |
| 15 | 19 | 3.30 | 3.30 | 3.65 |
| 20 | 26 | 2.65 | 2.80 | 3.00 |
| 25 | 33 | 2.30 | 2.45 | 2.60 |
| 30 | 40 | 2.05 | 2.25 | 2.25 |

So the main list's shortest courses, and its upper dose bounds, are statements about the
early-toxicity model more than about the tumor model.

## 2. What the app uses now, and how good the evidence is

**Model.** Early BED = n·d·(1 + d/(α/β)) − (λ/α)·max(0, T_elapsed − Tk), with Fowler et al.
2003 values α/β 10 Gy, α 0.35 Gy⁻¹, Tk 7 days, Tp 2.5 days (λ/α = 0.792 Gy₁₀ per day). The main
page uses these fixed values and the limit "standard's early BED + 5 Gy₁₀" (standard 53.1 Gy₁₀,
HYPNO 55 Gy/20 is 55.9 Gy₁₀). Toxicity uncertainty appears only on the Uncertainty tab.

**Evidence grade by component** (judgment, reasons given):

| Component | Source | Grade | Main weakness |
|---|---|---|---|
| α/β 10 Gy | Denham 1995: 9.3 (95% 5.8–17.9) | Moderate | Wide interval. Matters most for large d, exactly where the app works |
| Tk 7 days (elapsed) | Fowler 2003, chosen by inspection | Weak | Pardo-Montero 2021 refit on 108 schedules gives 19.4 days (inclusive). Our refit leaves Tk broad (section 3) |
| λ/α 0.79 Gy₁₀/day | Fowler 2003 (Tp 2.5 d "as a datum") | Moderate | Bentzen 2001 (CHART) Dprolif 0.80 (0.7–1.1) Gy/day in 2-Gy equivalents, about 0.96 Gy₁₀/day: consistent |
| Tolerance limit | Fowler zone 59–61 Gy₁₀ (RT alone). App uses relative +5 Gy₁₀ | Weak | The +5 Gy₁₀ margin has no stated risk meaning |
| Applicability to 13–20 fractions of 2.6–3.7 Gy | None | Weak (extrapolation) | Of 116 schedules in the Pardo-Montero data, about 6 exceed 2.2 Gy per fraction and none lies in this region |
| Concurrent chemotherapy | Fowler 2003, Meade 2014 (about 3–5 Gy₁₀ offset) | Weak | Model data are mostly radiotherapy alone. The standard of care is chemoradiotherapy |

**Uncertainty tab priors** (`R/uncertainty.R::toxicity_prior`): α/β lognormal (95% 5–20),
α uniform 0.30–0.40, Tp uniform 2–4 days, Tk uniform 7–19 days, all independent and stated as
assumptions. Two problems: (i) only λ/α = ln2/(α·Tp) and Tk enter the model, so drawing α and Tp
separately adds nothing but obscures what is assumed. (ii) the ranges are not fitted, and the
strong data-induced correlation between Tk and λ/α is ignored.

**Clinical check available now.** HYPNO (Bentzen 2026): grade 3+ mucositis 50.9% (190/373) for
55 Gy/20 in 4 weeks against 54.7% (208/380) for 66 Gy/33 at 6 fractions a week. Cisplatin in 75%.
The current model gives about 55.9 and 55.4 Gy₁₀, consistent with no difference. This is one
arm pair, at 2.75 Gy, and does not test 13–16 fractions.

## 3. Exploratory refit (done 2026-11-10, scratch only)

Data: Pardo-Montero 2021 supporting file (116 entries: 98 Tol, 13 Intol, 2 Disp_Tol,
2 Disp_Intol, 1 Disp). The single unresolved entry was dropped. Disputed entries were counted
with their leaning (115 used, 13 intolerable). BED used incomplete repair, T½ 3.2 h, as in
Pardo-Montero. Model: P(intolerable) = logistic((BED₁₀ − (λ/α)·max(0, T_elapsed − Tk) − B0)/s).
Random-walk Metropolis, 4 chains, all split R-hat ≤ 1.005.

| Parameter | Weak prior on λ/α | λ/α prior from Bentzen 2001 |
|---|---|---|
| B0 (Gy₁₀) | 81.3 (61.9–99.3) | 72.0 (57.9–89.2) |
| Tk (elapsed days) | 28.9 (2.8–36.1) | 13.8 (0.8–24.4) |
| λ/α (Gy₁₀/day) | 2.12 (0.76–17.6) | 0.97 (0.80–1.21) |
| s (Gy₁₀) | 5.6 (3.0–10.4) | 4.6 (2.6–8.3) |
| α/β (Gy) | 10.3 (5.6–19.7) | 11.4 (6.3–20.8) |

Medians with 95% intervals. Readings: (a) the binary data alone cannot separate Tk from λ/α.
(b) with the Bentzen prior the posterior is usable, but Tk still spans almost the whole range.
(c) α/β is barely updated from its prior. Treat these as a feasibility check, not results:
duplicates are not removed and endpoint definitions are not harmonized (Phase A).

**Effect on the main list.** Highest dose per fraction meeting "standard + 5 Gy₁₀", per
posterior draw (`implic.py`):

| Fractions | App now | Refit, informed prior: median (95%) | Refit, weak prior: median (95%) |
|---|---|---|---|
| 13 | 3.65 | 3.45 (3.15–3.75) | 3.25 (2.15–3.80) |
| 15 | 3.30 | 3.15 (2.85–3.35) | 2.95 (2.15–3.40) |
| 20 | 2.80 | 2.65 (2.50–2.80) | 2.55 (2.05–2.85) |
| 25 | 2.45 | 2.40 (2.30–2.45) | 2.35 (2.05–2.50) |
| 30 | 2.25 | 2.20 (2.20–2.25) | 2.20 (2.05–2.25) |

The app's current values sit at the upper end of the data-supported range. If the refit holds,
13–15 fractions likely leave the list (the tumor needs ≥ 3.65 and ≥ 3.30 Gy) and HYPNO's
2.75 Gy at 20 fractions lies above the median limit. This direction must be checked against the
clinical arm data (Phase A2) before anything changes in the app: HYPNO did not show more mucositis
with 55 Gy/20.

Implied risk at the app's current upper doses (informed prior): median P(intolerable) about
0.05–0.06 at 13–20 fractions against 0.009 for the standard. The +5 Gy₁₀ margin therefore
corresponds to a several-fold relative increase in modeled intolerability at the short end.

## 4. Goal and success criteria

Replace fixed, partly hand-picked acute parameters with joint posterior draws fitted to two
independent kinds of clinical data, carry those draws into every page that applies an early limit,
and state the limit on a risk scale. Done when:

1. Every early-toxicity parameter used by the app traces to a fitted posterior or a cited source
   (evidence ledger entries, constant registry entries).
2. Out-of-sample checks (Phase C) pass the pre-stated thresholds.
3. Each main-page row shows the share of draws meeting the early limit, and the table states the
   limit as a modeled risk difference as well as in Gy₁₀.
4. Fowler 2003 values remain available as a named preset that reproduces version 4.5.0 exactly.

## 5. Work plan

### Phase A: data (Claude extracts, Igor checks a sample)

- **A1. Schedule-level tolerability data (Pardo-Montero/Fenwick, 116 entries).** Write
  `data/early_tox_schedules.csv` with source label, fraction doses and times, class, and
  provenance (file, row). Identify duplicates (the same schedule reported under several labels,
  such as K22/F9 and F19/K45) and keep one, or model them with a shared-schedule effect. Decide the
  coding of disputed entries (main: drop "Disp", count leanings. Sensitivity: drop all disputed).
  Record which entries used concurrent chemotherapy, if the source says.
  *Gate:* counts reproduce 98/13/2/2/1. 10 random rows checked against Fenwick 2008 tables.
- **A2. Arm-level acute mucositis from randomized trials.** For each trial in
  `refit_trials_2026_10_25c_main.csv` (19 trials, 45 arms), plus HYPNO, KROG 0201, PMH, Yamazaki
  and any other randomized or prospective hypofractionated series found by a PubMed search
  (protocol: search terms and screening log saved), extract per arm: number evaluable, number with
  peak acute mucositis at grade ≥ 3, scoring system (RTOG, CTCAE version, WHO, Dische), chemo use,
  technique (2D/3D/IMRT). Two independent transcriptions, as done for the tumor data.
  Disagreements resolved against the PDF. Every value needs the downloaded full text (no values
  from memory or abstracts). Missing PDFs reported by name.
  *Gate:* two transcriptions agree after reconciliation. Provenance CSV complete.
- **A3. Literature priors**, each with verbatim quote and page in `References/README.md`:
  Denham 1995 (α/β), Bentzen 2001 (Dprolif, and its onset estimate, if reported in the full text),
  Fowler 2003, Fenwick 2008, Pardo-Montero 2021, Meade 2014 (chemotherapy offset).

### Phase B: model (Stan, in `fit/`, the same toolchain as the tumor refit)

Shared biology across both likelihoods: α/β, Tk, r = λ/α.

- **Likelihood 1 (schedules):** Bernoulli(intolerable), logit = (E − B0)/s, E = BED₁₀(α/β, with
  incomplete repair) − r·softplus(T_elapsed − Tk).
- **Likelihood 2 (trial arms):** Binomial(grade ≥ 3 events | n), logit p_ij = μ_j + E_ij / s₂,
  with μ_j a per-trial intercept (partially pooled). The intercept absorbs scoring system,
  chemotherapy and era. The within-trial contrast identifies α/β, Tk and r. This design is why
  the second data set matters: within-trial differences are much less confounded than
  between-study tolerability labels.
- **Priors:** α/β lognormal(log 10, 0.354). Tk uniform(0, 35) elapsed days. R weakly informative
  lognormal(log 0.8, 1.0). S, s₂ half-normal(5). Bentzen 2001 used **either** as an r prior
  **or** as validation, not both, and if used as a prior then the CHART arms are dropped from
  likelihood 2 (they are the data behind that estimate). Decision D3.
- **Alternatives fitted on the same data:** Fenwick 2008 boundary (no Tk), and a Tk-free
  linear time factor. Compared by leave-one-trial-out log score. Main model chosen by D5.
- **Diagnostics gate (as for the tumor refit):** no divergences, R-hat < 1.01, bulk and tail ESS
  > 400 per parameter. Prior and posterior predictive checks. Simulation-based calibration on
  synthetic data with known parameters (registered in `sim_provenance.yaml`). Posterior
  correlation of Tk and r reported.

### Phase C: validation (pre-stated, before looking at the result)

- Leave-one-trial-out prediction of each trial's arm difference in grade ≥ 3 mucositis. Report
  coverage of the 90% predictive intervals (target: 80–97%).
- Held-out check against HYPNO: predicted risk difference 55/20 minus 66/33 (6 per week) against
  the observed −3.8 percentage points (compute its 95% interval from 190/373 and 208/380).
- Predicted early BED and risk for HYHOPE 46.5/15 and 44.4/12 against what those reports
  describe (the 12-fraction level was dropped for mucositis in a parallel trial).
- Fixed fixture plus seed pinned as a golden-value test.

### Phase D: decision rule on the main page (Igor decides, D2)

Options, all evaluated per posterior draw:

1. Keep "+5 Gy₁₀" and add an "early support" share (proportion of draws meeting it), flagged
   below 80% like the tumor refit support.
2. Risk scale: modeled grade ≥ 3 mucositis risk difference against the standard at most
   δ percentage points (for example 5) in at least 80% of draws.
3. Show both during one release, then retire option 1.

A translation table (Gy₁₀ margin ↔ risk difference at each fraction number) is published either
way, so users can see what "+5" means.

### Phase E: app integration

- `draws/early_tox_<variant>.csv` (4000 joint draws: α/β, Tk, r, B0, s, s₂) with SHA-256 in the
  manifest. Loaded and checked like the tumor draws.
- Replace `toxicity_prior()` acute entries with the draws. Change `bed_early_draws()` to use r
  directly (drop separate α and Tp). Keep the late and organ priors unchanged.
- Main page: early-support column. Joint support pairs tumor and early draws by index (they are
  independent fits, so pairing is valid and stated).
- Presets: "Fowler 2003 (version 4.5)" reproduces 4.5.0 exactly. "Early refit point" uses
  posterior medians of identified quantities.
- Constant registry: register Tk_e, λ_e/α_e, α/β_e, margins. Retire replaced literals in both
  R and Python arms. Run `stale_constant_audit.py --strict`.
- About page, user guide and README updated in the same pass (README parity), with the evidence
  table of section 2.

### Phase F: tests, release

- One committed test per change, each seen red without it and green with it (e.g. draws file
  hash mismatch refused. Fowler preset reproduces 4.5.0 list exactly, 69 doses in 18 rows.
  Early support equals the share computed independently in Python. Incomplete repair term zero
  for once-daily schedules).
- Combined smoke run, full suite (`tools/run_tests.R`), live walkthrough, then `tools/deploy.R`
  only after Igor approves the release.
- `docs/FIXES.md` rows. First line states that the main list moves, with before and after counts.

## 6. Decisions for Igor (recommendation first)

- **D1. Endpoint for likelihood 2:** peak acute mucositis grade ≥ 3 (any scale), with scale as a
  trial-level covariate. Recommended. Alternative: confluent mucositis only (fewer trials).
- **D2. Main-page rule:** option 3 (show both, then move to the risk scale). Recommended.
- **D3. Bentzen 2001:** use as validation, not prior. Recommended, because CHART is in the data.
- **D4. Chemotherapy:** absorbed by trial intercepts. Assume no interaction with schedule.
  Recommended for now. Flag as a limitation, test with a chemo × E term if enough arms.
- **D5. Default model:** keep Fowler as the default until Phase C passes. Then the refit becomes
  the default only if its validation is at least as good. Recommended.
- **D6. Fit α/β** (prior from Denham) rather than fix at 10. Recommended.

## 7. Risks

- Tk may stay unidentified even with arm data (few arms below 25 days). Then the result is wide
  uncertainty, which is still more truthful than a fixed value. The main page will show low early
  support for short courses.
- Mucositis scoring differs across eras and scales. Trial intercepts handle level, not slope.
- Modern IMRT spares mucosa unevenly. Dose-area effects (Bentzen 2001) are not modeled.
- 13 intolerable schedules is a small number of events for likelihood 1.
- Publication and reporting bias in which trials report acute grades.

## 8. Effort (judgment, not measured)

A1 half a day. A2 two to four days (largest item, depends on PDF access). B one to two days.
C one day. D–F two days. Compute is small (minutes per fit).

## 9. Reproducing section 3

```
cd docs/development/2026_11_10_early_toxicity/exploratory
python3 parse.py <path>/PardoMontero_2021_SupportingInfo_108_schedules.xls sched.pkl
python3 fit.py sched.pkl weak 1 draws_weak.npy
python3 fit.py sched.pkl informed 1 draws_informed.npy
python3 implic.py
Rscript binding.R   # edit setwd() to the app folder
```
The supporting file is in `Custom_BED_calculator_v2/References/` (CC BY-NC, not committed).
