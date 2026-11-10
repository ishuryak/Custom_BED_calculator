# Audit of the early-toxicity plan (2026-11-10)

Status note (added 2026-11-10 with plan version 2): F4, F5, F6 and F9 were revised or replaced by the independent review `2026_11_10_early_toxicity_independent_review.md` (section 5), and the posterior tables in section 1 below are provisional because they include 8 partial-course records. Plan version 2 (`2026_11_10_early_toxicity_plan.md`) is the current specification. The text below is unchanged.

Audited document: `2026_11_10_early_toxicity_plan.md` (version 1, unchanged by this audit).
Reviewers: two independent subagent reviews run in parallel, one methodological and statistical
(six-factor methods audit, all exploratory numbers re-run, including `binding.R` under Windows R),
one radiation oncology and radiobiology (every factual claim checked against the local PDFs of
Fowler 2003, Fenwick 2008, Pardo-Montero 2021, Meade 2014, Fowler 2010, Bentzen 2001, HYPNO 2026,
HYHOPE protocol 2024 and outcomes 2025). The consolidation and the corrected re-run below are by
the main session.

**Verdict: needs revision before Phase B.** The aim and the overall design (fit to tolerability
labels plus within-trial arm data, propagate joint draws, keep Fowler as a reproducible preset)
are sound. The exploratory evidence in plan section 3 had a coding error, one false coverage
statement and a misencoded prior. With these corrected, the label-only refit contradicts the
best randomized evidence (HYPNO), which changes what the plan must do.

## 1. Corrected exploratory results (replace plan section 3)

Corrected scripts: `2026_11_10_early_toxicity/exploratory_after_audit/` (the version 1 scripts are
kept unchanged in `exploratory/` so the version 1 numbers stay reproducible). Changes: the two
"Disp_Intol" schedules are now coded intolerable (15 events, not 13), the Bentzen-based prior on
λ/α is lognormal(log 1.053, 0.115), whose 95% range is the stated 0.84–1.32 Gy₁₀/day, the weak
prior is lognormal(log 0.8, 0.5) (95% 0.30–2.1), and draws with no qualifying dose are counted
as "none" instead of being silently dropped.

| Parameter | Weak prior | Bentzen-based prior |
|---|---|---|
| Tk (elapsed days) | 18.1 (2.0–28.2) | 16.6 (1.8–23.3) |
| λ/α (Gy₁₀/day) | 1.19 (0.80–2.19) | 1.09 (0.90–1.33) |
| B0 (Gy₁₀) | 72.8 (59.4–88.5) | 71.7 (57.9–86.5) |

| Fractions | App now | Refit, weak | Refit, Bentzen-based |
|---|---|---|---|
| 13 | 3.65 | 3.25 (2.55–3.70) | 3.35 (3.00–3.60) |
| 16 | 3.25 | 2.85 (2.20–3.25) | 2.95 (2.70–3.15) |
| 20 | 2.80 | 2.50 (none–2.80) | 2.60 (2.40–2.75) |
| 30 | 2.25 | 2.20 (none–2.25) | 2.20 (2.15–2.25) |

**HYPNO face-validity check (new, derived).**

| Quantity | App now (Fowler) | Refit, weak | Refit, Bentzen-based |
|---|---|---|---|
| Early BED, 55/20 minus 35 × 2 Gy (Gy₁₀) | +2.8 | +11.0 (+2.7 to +27.8) | +8.9 (+4.8 to +14.2) |
| Share of draws in which 55/20 passes "+5" | yes | 10% | 3% |
| Early BED, 55/20 minus 66/33 at 6 per week | +0.4 | +4.9 (+0.2 to +13.6) | +3.9 (+1.3 to +7.0) |

Observed in HYPNO: grade 3+ mucositis 50.9% against 54.7%, a difference of −3.8 percentage
points (95% CI −10.9 to +3.3, computed from 190/373 and 208/380). So the label-only refit predicts
the 4-week arm to be worse, and the randomized trial does not show that. The trial cannot exclude
a harm of up to 3.3 points, so this is tension, not proof of error. The fixed Fowler values (Tk 7)
agree better with HYPNO than the refit does. The most likely reasons are that the 1980s–2000s
tolerability labels, mostly from radiotherapy alone and conventional or accelerated schedules,
do not constrain Tk in the 3–4 week range, and that the label endpoint (prolonged, intolerable
reactions) is not the same as peak grade 3+ mucositis in modern chemoradiotherapy. Both readings
are inferences, not shown.

**Consequence for the plan.** A label-only refit must NOT replace the Fowler values. The arm-level
clinical data (Phase A2), including hypofractionated arms, are required, and HYPNO must act as a
pre-stated gate (finding 5).

## 2. Findings

Severity: blocker, major or minor. Each finding gives the evidence, fix, test and verification.

**F1. Disputed-intolerable schedules were coded tolerable (major, verified by both reviewers).**
Evidence: `exploratory/parse.py` line 20 tests `startswith('int')`, so "Disp_Intol" (including
F31, 55 Gy in 16 fractions of 3.44 Gy) became tolerable, and the plan's description of the coding
was wrong. Fix: done in `exploratory_after_audit/parse.py` (`'intol' in cls.lower()`). Test: the
event count is 15 under the main coding and 13 when all disputed entries are dropped. Verification:
`python3 fit.py sched.pkl informed 1 x.npy` prints `n=115 intolerable=15` (observed 2026-11-10).

**F2. The data do reach the binding region (major, verified).** Evidence: F31 has 16 × 3.28 Gy
(tolerable) and 16 × 3.44 Gy (disputed, leaning intolerable), both over 22 days. Fowler 2003 p5
also mentions Manchester 16 × 3.3–3.375 Gy in 21 days as tolerated. Plan section 2 says no
schedule lies in this region, which is wrong. Fix: replace with "three entries lie in or next to
it (F31 twice at 16 fractions, K18 at 20 × 2.55 Gy)" and add F31 as a named posterior predictive
check in Phase C. Test: a parse check listing 13–20 fraction schedules above 2.6 Gy returns the
two F31 rows.

**F3. The Bentzen 2001 prior does not match the model it informs (major, verified from source).**
Evidence: Bentzen's Dprolif of 0.80 Gy/day assumed proliferation active from day 12 or earlier and
ignored incomplete repair (Bentzen 2001 p9–10: a 25-day lag would raise it to 1.3–1.4, and
incomplete repair with T½ 3.2 h would lower it by about 25%). The plan uses it as a prior on λ/α
while Tk is free. The version 1 code also encoded a narrower range than stated (0.76–1.21 instead
of 0.84–1.32). In addition, CHART is in both the label data and the planned arm data. Fix: use
Bentzen 2001 as validation only (decision D3 already recommends this) and drop the "informed"
variant from any default. If a prior is ever needed, encode it on the quantity Bentzen measured
(BED recovered between elapsed days 12 and 45). Test: a prior predictive check of recovered dose
over days 12–45 against Bentzen's estimate.

**F4. The Tk prior ignores mucosal biology (major, inferred).** Evidence: Dörr's data summarized
in Bentzen 2001 p10 (about 60% of an 8 Gy dose repopulated by day 10) and the biopsy data cited by
Fowler 2003 p3 put the onset near 7–14 days. The weak-prior Tk of 18–29 days falls after the
usual peak of mucositis. Also, Fowler's Tk 7 was chosen from biopsy data, not "by inspection"
as plan section 2 says. Fix: main prior lognormal with median 8 and 95% range about 4–16 days,
with uniform(0, 35) as a sensitivity analysis. Correct the wording in section 2. Test: Phase C
HYPNO gate under each prior.

**F5. No pre-stated clinical face-validity gate (major, derived, see section 1).** Fix: add to
Phase C, before any fit: "A model variant becomes a default only if (a) the posterior predictive
difference for HYPNO 55/20 against 66/33 at 6 per week has a 90% interval overlapping the observed
−10.9 to +3.3 points, and (b) 55/20 passes the early limit in at least 50% of draws." Add HYHOPE
(46.5/15 with weekly cisplatin, grade 3 mucositis 5/24, reduced elective volumes) and the
Jacinto 2018 phase 1 (to be read in full before use) as single-arm calibration points with a site
intercept. HYPNO must then be removed from the A2 fitting set (F7). Test: the gate fails the
current label-only draws (it does: 3–10% pass).

**F6. Duplicate schedules break independence (major, verified).** Evidence: 115 rows hold 92
unique schedules (15 duplicate groups). Three of the events are one schedule, 72 Gy in 36 fractions
over 36 days. One schedule (54 Gy in 30 fractions over 12 days) is labeled intolerable in one
source and tolerable in another. About 11 unique events support 5 parameters. Fix: deduplication,
or a shared-schedule random effect, becomes the main analysis in A1, with conflicting labels kept
as two observations. Test: `(n, total dose, T)` is unique in `early_tox_schedules.csv`.

**F7. Double counting between the two data sets, and HYPNO both fitted and held out (major,
verified by reading).** Evidence: A2 lists HYPNO for fitting while Phase C holds it out. Trials
such as CHART appear both as tolerability labels and as planned arm data, so the two are not
"independent kinds" of data. Also KROG 0201, PMH and Yamazaki are listed in A2 as additions but
are already among the 19 trials. Fix: build a source-to-trial crosswalk and drop from the label
likelihood any schedule whose trial enters the arm likelihood. Keep HYPNO held out. Correct the
A2 list. Test: crosswalk overlap count is 0.

**F8. The proposed risk-scale rule mixes two different quantities (blocker for Phase D,
verified by reading).** Evidence: P(intolerable) from the label model is the probability that a
curated review labels a schedule intolerable, not a patient's risk of grade 3+ mucositis. Plan
section 3 ("several-fold relative increase") and Phase D option 2 treat it as patient risk.
A risk difference from the arm model also depends on the trial intercept, which the plan does not
fix. Separately, the "+5" relative rule cancels B0 and s, so it uses only Tk, λ/α and α/β. Fix:
define δ only on the arm model at a stated baseline (recommended: posterior predicted grade 3+
rate for 70/35 in a chemoradiotherapy-era IMRT trial, cited), and relabel the section 3 statement
as a ratio of label probabilities. Test: the δ function refuses to run without a baseline.

**F9. Chemotherapy differs between arms and may interact with time (major).** Evidence: HYPNO
cumulative cisplatin was 121 against 149 mg/m² (HYPNO p8). Cisplatin may suppress compensatory
proliferation, which would interact with the time term. Trial intercepts cannot absorb an
arm-level difference. Fix: extract cumulative platinum dose per arm in A2, add an arm-level term,
and fit chemo × time as a sensitivity analysis. Test: the A2 schema rejects an arm row without
chemotherapy dose.

**F10. Endpoint mismatch between the two likelihoods (major, inferred).** Evidence: "intolerable"
in Fenwick 2008 means prolonged reactions with possible consequential late damage. Peak grade 3+
is common in chemoradiotherapy (about 50% in HYPNO) and, after 3–4 week courses, can peak after
treatment ends, while HYPNO scored "during and immediately after" treatment. Fix: separate slopes
per likelihood (already planned), extract the assessment window per arm, and add feeding tube use
and unplanned breaks as secondary endpoints where reported. Test: the A2 schema rejects an arm
without an assessment window.

**F11. Phase C thresholds cannot separate good from bad models (major, derived by reviewer).**
Evidence: with about 20 trials, an 80–97% coverage target is met by a perfectly calibrated model
only about 84% of the time and by a model with 75% true coverage about 41% of the time. Fix: the
primary test becomes the leave-one-trial-out log score of the refit against the fixed Fowler
values, with coverage reported descriptively with a binomial interval, plus the F5 gate. Test:
the Phase C script exits non-zero if any criterion lacks a numeric threshold.

**F12. Minor corrections.** (a) The text writes max(0, ·) but the code uses softplus with
sharpness 3 per day (at most about 0.2 Gy₁₀ at the kink), and the Tk prior is uniform(0, 40) in
code but (0, 35) in Phase B. State both. (b) At 30 fractions early and late tie at 2.25 Gy, so
"only at 30 fractions does the late limit bind" should read "at 30 fractions the two tie". (c) The
planned test "incomplete repair term zero for once-daily schedules" is wrong: with T½ 3.2 h the
24-hour term adds about 0.12 Gy₁₀, so the test should be "below 0.2 Gy₁₀". (d) HYHOPE level 1
(44.4 Gy in 12) was given 4 times a week over 3 weeks, and it was dropped after a preliminary
signal in a parallel postoperative trial, with no dose-limiting toxicity. State that. (e)
Denham 1995's 9.3 (5.8–17.9) comes from secondary slides, and the abstract gives 7–10 Gy. Get the
full text before using it as a prior. (f) `References/README.md` A7 says Bentzen 2001 is not
downloaded, but it is held in `References_refit_trials/`. (g) Bentzen 2001 reports a dose-area
effect (p7) of 1.21 Gy per cm of field length. Add it to the limitations.

## 3. Missing radiation oncology elements (add to plan section 7)

Mucosal dose under IMRT with simultaneous integrated boost (elective mucosa gets 2.3–2.7 Gy, not
the prescription dose), subsite and HPV status, volume de-intensification (HYHOPE), consequential
late effects as the clinical reason the acute limit matters, supportive-care era, and handling of
treatment breaks and actual overall time. Candidate additional data sets, all to be verified
before use: the eight trials listed by Lester in HYPNO's references, Madhava (55/20 with
carboplatin), Sanghera 2007, Jegannathen (Christie), and Bratman phase 2.

## 4. Investigated and dismissed

- MCMC validity: adaptation stops at burn-in, so the stored chain has a fixed kernel. Smallest bulk ESS 956.
- Jacobian for the half-normal prior on s sampled as log s: correct.
- Incomplete-repair formula: matches Thames' formula (standard 84.1 against 84.0 Gy₁₀).
- Clock conversion (elapsed = inclusive − 1): consistent with the app.
- HYPNO BED arithmetic (55.9 against 55.4 Gy₁₀) and the EQD2 to BED₁₀ conversion of Bentzen's estimate: correct.
- Pairing independent tumor and early draws by index: acceptable.
- Relative "+5" criterion in principle: defensible, since an additive chemotherapy offset cancels.

## 5. Implementation order

1. **Igor:** decisions below (blocks step 3).
2. **Claude, about 1 hour:** revise the plan to version 2 with F1–F12 and section 1 of this audit.
   Gate: Igor approves version 2.
3. **Claude, half a day:** Phase A1 with deduplication (F6) and crosswalk (F7). Gate: counts
   15/13, unique keys, overlap 0.
4. **Claude with Igor's spot check, 2–4 days:** Phase A2 with chemotherapy dose, assessment window
   and technique (F9, F10). Every value from a downloaded full text. Missing PDFs reported by name.
   Gate: two transcriptions reconciled.
5. **Claude, 1–2 days:** Phase B in Stan with the F4 prior, simulation-based calibration on the
   real arm design. Gate: diagnostics as in the plan.
6. **Claude, 1 day:** Phase C with the F5 HYPNO gate and the F11 log-score test. Gate: refit beats
   or ties Fowler on log score and passes HYPNO. **If it fails, stop:** keep Fowler as the default,
   add the refit only as an uncertainty band on the Uncertainty tab, and report.
7. Phases D–F as planned, with F8 settled first.

## 6. Decisions for Igor (recommendation first)

1. **Revise the plan to version 2 now?** Recommended yes.
2. **D3, Bentzen 2001:** validation only. Recommended (F3).
3. **Tk prior:** biological lognormal (median 8 days, 95% about 4–16) as main, uniform as
   sensitivity. Recommended (F4).
4. **HYPNO gate (F5):** a variant that excludes 55/20 in more than half of its draws cannot
   become the default. Recommended.
5. **D2, main-page rule:** keep "+5 Gy₁₀" with an early-support share until the arm model passes
   Phase C. Only then add a risk-difference rule at a cited baseline (F8). Recommended.
6. **Interim app change (optional):** none. The audit does not justify changing the live early
   limit before the arm data are fitted. Recommended.

## 7. Stopping condition

When steps 2–6 pass their gates, proceed to Phases D–F. Reopen this audit only if a gate fails or
new evidence appears (for example a further randomized hypofractionation trial with per-arm
mucositis).
