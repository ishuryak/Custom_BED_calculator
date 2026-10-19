---
title: "HNSCC BED calculator v2: review and action plan"
subtitle: "Mathematical correctness of both repopulation models and clinical fitness for exploring HNSCC trial design"
---

Igor Shuryak, MD, PhD — is144@cumc.columbia.edu

Review date: October 18, 2026. Status: **PLAN ONLY. No code, defaults or app text were changed.**

Reviewed: branch `v2-toxicity-hypofractionation` at commit `6f0209c` (app code unchanged since `83a188c`), the handoff, the fix plan (D1–D18), the two DD parameter analyses, the external (Codex) review F1–F12, the reference catalog, and the 2018 DD paper (Shuryak, Hall and Brenner, Radiother Oncol 2018;127:20–26).

# 1. Bottom line

1. **The arithmetic is correct.** The app's DI and DD formulas match equations 1–5 of the 2018 paper (without the growth term g), and the acute formula matches Fowler 2003. I re-ran the test suite today: 1131 expectations, 0 failures. An independent Python implementation written for this review reproduces every number I checked.
2. **The DD model in the app is not the published DD model, and the plan does not yet say so.** The app shares one α between DI and DD and solves the DD onset constant C so that DD equals DI at the reference. The 2018 model fitted α, the maximum rate and C to the trial data and gave DD its own α. The two give opposite answers for the trials that matter most (section 3.1). The handoff's conclusion that the DD contradiction with accelerated trials is "structural" is wrong as a statement about the DD model. It is a property of the app's calibration.
3. **The proposed D18 DD defaults should not be approved as they stand.** They were chosen by plausibility bands for onset, doubling time and daily loss that come from DI-model analyses (Withers, Fowler). Those bands exclude the region where the fitted DD model lives (late onset, fast rate), and the D18 DD prediction for HYPNO (a 5.8 Gy₁₀ deficit) conflicts with the trial's result. The DI half of D18 (α 0.35, λ_DI 0.23, 0.66 Gy₁₀/day) is well supported and can go ahead.
4. **For trial design, the main output should be a trade-off, not a filtered list.** At HYPNO's acute BED, every regimen shorter than 55 Gy/20 fractions loses tumor BED against HYPNO, under both models and both parameter sets tested (section 4.1). A new regimen therefore has to buy shortening with acute toxicity, or give up tumor BED. The app should show that curve, compare against HYPNO as well as 70/35, and allow a non-inferiority margin.
5. **Several app statements and limits are clinically misleading** (section 4.2): the late limit allows more dose in 20 fractions than HYPNO tested, the search reaches 8 Gy per fraction, chemotherapy is assumed equal between arms, and a Fowler 59 Gy₁₀ default would be a radiotherapy-alone limit in a cisplatin setting.
6. **Six correctness fixes come first** (time convention for early and tumor BED, a fixed DD calibration anchor, published-value tests, Withers preset, negative-BED handling, export provenance). They are small, and none needs a scientific decision except the anchor default.

# 2. What is correct (verified)

| Item | How verified | Result |
|---|---|---|
| DI tumor BED = paper Eq. 1 without gT | Term by term against the 2018 paper | Match |
| DD onset = C·T/(αK) (paper Eq. 2), P_kill = 1 − exp(−αK/T) (Eq. 3), rate γ·P_kill (Eq. 4–5) | Term by term; αK = n(αd + βd²) | Match |
| Calibration solves DD = DI at the reference for any α/β (the v1 bug at α/β ≠ 10 is fixed) | Existing test, plus my Python | Match |
| Acute BED, Fowler parameters | Fowler 2003 Table 1 at T = 46 d: 53.11 (paper 53.1); CAIR at 34 d: 62.61 (paper 62.6) | Match on Fowler's elapsed-day convention only (see M2) |
| Calendar time (inclusive) | 35 fx = 47 d, 20 fx = 26 d; HYPNO 6 fx/week 33 fx = 38 d | Match |
| D18 numbers in the handoff | Python: reference DI 71.51, DD reference onset 23.5 d, HYPNO test − control under DD −5.8 | Match |
| Test suite | Windows R 4.5.1, `testthat::test_dir`, today | 1131 passed, 0 failed, 2 package-version warnings |

# 3. Mathematical and model issues

Severity: **A** = changes conclusions a trial designer would draw; **B** = changes numbers or reproducibility; **C** = documentation or robustness.

## 3.1 M1 (A). The app's DD calibration is not the fitted 2018 model

**What differs.** In the paper, DI and DD were fitted separately to 33 trial arms with α/β 10 and g = 0.0154/day: DI α 0.069, γ 0.035/day, Tk 28.6 d; DD α 0.224, γ 1.17/day, C 14.5 (paper Table 2). The app instead uses one α for both, and computes C so that DD = DI at the reference.

**Why it matters.** The fitted DD model places onset for 35 × 2 Gy at day 36.2 and repopulates fast after that (rate 0.386/day, doubling 1.8 days, 1.72 Gy₁₀/day). Short courses then finish before or soon after their own onset. Recomputed with the fitted parameters (Python, this review):

| Experimental arm vs its control | Fitted 2018, with g: ΔDI / ΔDD (Gy₁₀) | Fitted 2018, g = 0: ΔDI / ΔDD | App, current defaults: ΔDI / ΔDD | App, D18: ΔDI / ΔDD | Observed Δ control (points) |
|---|---|---|---|---|---|
| GORTEC very accelerated | +7.2 / **+3.2** | +1.4 / +1.4 | +11.6 / **−2.6** | +4.1 / **−3.7**¹ | +18.2 |
| TROG accelerated | +2.8 / **+2.1** | −3.1 / +0.3 | +8.1 / **−6.6** | −1.4 / **−8.5**¹ | +8.5 |
| HYPNO 55/20 vs 66/33 | −1.6 / **+0.5** | −4.3 / −0.3 | +0.9 / **−4.7** | −2.4 / **−5.8** | within ±1.4 |
| CHART | −1.4 / −0.8 | −8.8 / −3.1 | −0.1 / −7.3 | not computed | −1.2 |

¹ D18 rows for GORTEC and TROG are against 70/35 in 47 days, not the trial's own control. Within-model differences only: ΔBED values from fitted models with different α are not on the same scale (see M8).

**Reading.** The fitted DD model gets the direction right for GORTEC, TROG and HYPNO. The app's calibrated DD gets all three wrong, and D18 makes it worse. The grid search found the right structure (short courses avoid repopulation only if the reference onset is late, about 35 days, and the rate after onset is fast) and then rejected it, because its bands cap the daily loss at 1.0 Gy₁₀ and the doubling time at no less than 2 days. Those bands are estimates made with the DI model, so using them to select DD parameters assumes DI is right. Under the DD hypothesis, "onset" and "rate" mean something different.

**HYPNO as a quantitative check (inference, flagged).** A tumor dose-response slope of about 1–2.5 percentage points of locoregional control per Gy₁₀ near 70 Gy (a common range; to be verified against a source before it is written anywhere) turns D18's 5.8 Gy₁₀ DD deficit into roughly 6–15 points. HYPNO observed at most 1.4 points and passed a 10-point margin. The D18 DD setting is therefore hard to reconcile with the one randomized hypofractionation trial in this population.

**What to do** (decision N1, section 6):

- Keep the current calibration as an option, labeled "DD calibrated to DI at the anchor (pessimistic for short courses)".
- Choose DD parameters by agreement with outcome data, not by DI-derived bands. Concretely: refit γ_DD and C (α fixed at the chosen value, or fitted) to the 2018 Table 1 summary data with a separate clonogen-number nuisance for each trial, so that only within-trial contrasts inform the repopulation parameters. This produces **new** parameters, which respects your instruction not to reuse the old fitted values as defaults. It is a research step (Phase 3) with its own record.
- Use the published 2018 fitted values **in a test only**, to prove the formulas reproduce the paper (M4). A test fixture is not a default.

## 3.2 M2 (B). Time convention: early BED and tumor BED

- **Early BED (WP1 1a, confirmed).** Fowler, Fenwick and Meade count elapsed days (inclusive T − 1). The app uses inclusive T with Tk 7, so early BEDs are 0.79 Gy₁₀ low (52.3 vs 53.1 for 70/35). Differences between regimens are unchanged, but the absolute Fowler limit is 0.79 Gy₁₀ too lenient.
- **Tumor BED (new).** The same issue applies to tumor parameters taken from Fowler 2010 or Hartley 2011, which were stated on the elapsed convention. The shift does not cancel for tumor BED, because short regimens end before Tk and the reference does not. Under D18, ΔDI for 55/20 vs 70/35 is −1.39 Gy₁₀ with inclusive T and −2.05 with elapsed T. The 2018 paper tabulates inclusive-style times (33 fx = 45 d), so its parameters belong with inclusive T.
- **Fix.** Keep T inclusive everywhere in the interface. Give each parameter set a stated convention and convert inside the formula: `T_model = T − offset`, with offset 1 for Fowler, Fenwick, Hartley, Meade and HYHOPE-derived sets and 0 for the 2018 and v1 sets. Show the convention on the Model page.

## 3.3 M3 (A for trial design). The DD model changes when the comparison reference changes (WP3 3a)

The DD constant C is solved at whatever regimen is entered as the reference. Choosing HYPNO 55/20 as the comparator, which is the natural comparator for a new trial, would silently change the DD model itself (55/20 has DD BED 58.0, 59.8 or 70.1 depending on the reference). **Make the calibration anchor a separate input, default 35 × 2 Gy / 47 d, and never let the comparator move it.** This is a prerequisite for the comparator feature in Phase 5.

## 3.4 M4 (B). No test checks the formulas against published numbers

The current golden values come from your own talk table, and the large grid test checks the code against an independent implementation of the same formulas. Neither proves fidelity to the papers. Add tests that reproduce published values:

| Source | Value to reproduce | Tolerance |
|---|---|---|
| Fowler 2003 Table 1 | 70 Gy/7 w at 46 d = 53.1; CAIR = 62.6; at least 5 more rows | ±0.05 (rounding of the table) |
| Fowler 2010 | 60 fractions in 39 d = 62.8; 42 d = 60.4 | ±0.05 |
| Hartley 2011 | 70/35/46 at α_e 0.3 = 48.0 | ±0.05 |
| Fenwick 2008 boundary (if D11 approved) | the four recomputed regimens in fix plan section 5 | ±0.1 |
| HYHOPE protocol | 15-fraction level 51.4 Gy₁₀ with its stated parameters | ±0.1 |
| Shuryak 2018, Table 2 parameters with g | qualitative Fig. 1: DD TCP for 2 Gy/fx plateaus within 30–40 fractions while DI keeps rising (T1/2, K = 2.1 × 10⁵) | shape test |
| Shuryak 2018, consistency | onset ratio (T/K)ᵢ/(T/K)_ref independent of α and γ (already derived; pin it) | 1e-12 |

## 3.5 Smaller issues

| ID | Sev. | Issue | Fix |
|---|---|---|---|
| M5 | C | Every max(0, x) uses a smooth version (up to 0.23 days, 0.23 Gy₁₀ DI, ~0.4 Gy₁₀ DD at defaults, up to 1 Gy₁₀ in some settings). Published formulas use a hard threshold. | Decision N2: switch to a hard threshold (recommended; v1 agreement then holds except within a day of onset), or keep and document. Calibration has a closed form with a hard threshold. |
| M6 | B | The growth term g is omitted. It matters most with small α (0.22 Gy₁₀/day for the 2018 DI α) and changes the sign of ΔDD for 55/20 vs 70/35 in the fitted model (+1.41 with, −0.03 without). | Add g as an optional input, default 0 (D5a). Required for the 2018 fidelity test. |
| M7 | B | Withers preset assumes α/β 10 (WP1 1b); search runs on a negative reference BED (1c). | As planned. |
| M8 | B | If DI and DD ever have different α (fitted mode, N1), ΔBED in Gy is not comparable across models: 1 Gy₁₀ is 0.069 log cell kill under the 2018 DI fit and 0.224 under DD. | Show Δ(log cell kill) = α·ΔBED next to ΔBED whenever the two α differ, and rank "both models" on that scale. |
| M9 | C | The search allows up to 8 Gy per fraction. LQ extrapolation for tumor and mucosa beyond about 4–5 Gy is weakly supported in this setting, and there is no curative precedent above 3.1 Gy per fraction (C5). | Default d_max 3.5 Gy; flag rows above 3.1 Gy; warning at 5 Gy. |
| M10 | C | The analysis notes evaluate 2–3 fractions/day without incomplete repair. Not an app defect (the app is once daily). | Keep the limitation sentence in the notes. |

# 4. Clinical appropriateness for trial design

## 4.1 The design question the app should answer

With Fowler's acute parameters, mucosa repopulates from day 7 at 0.79 Gy₁₀/day, while tumor (DI) repopulation starts near day 28. Shortening a course below about 26 days therefore removes acute compensation faster than it saves tumor BED. At HYPNO's own acute and late BED (55.9 Gy₁₀ elapsed convention and 105.4 Gy₃), the best regimen for each fraction number gives:

| Fractions | Dose (Gy) | Days | ΔDI vs HYPNO, v1 defaults | ΔDD vs HYPNO, v1 defaults | ΔDI vs HYPNO, D18 | ΔDD vs HYPNO, D18 |
|---|---|---|---|---|---|---|
| 19 | 2.80 | 25 | −2.0 | −1.2 | −2.0 | −1.5 |
| 17 | 3.05 | 23 | −2.5 | −1.3 | −2.5 | −1.6 |
| 15 | 3.20 | 19 | −6.8 | −3.6 | −6.8 | −4.5 |
| 12 | 3.75 | 16 | −8.3 | −4.2 | −8.3 | −5.3 |

(Python, this review; DD anchored at 70/35; 0.05 Gy dose steps; the DI columns agree for both sets here because none of these regimens reaches Tk = 28 d.)

So under the app's own models, a regimen that is both shorter than HYPNO and no more toxic acutely gives less tumor BED. That is a clear, testable hypothesis, and the app does not currently show it, because it compares only against 70/35 and filters for a gain. A trial designer needs:

1. **HYPNO 55/20 as a selectable comparator** (with the DD anchor fixed, M3).
2. **A trade-off view:** for each fraction number, tumor ΔBED (both models) against acute ΔBED, with the HYPNO point and the Fowler zone marked. This replaces the gain filter as the main output; the table stays.
3. **A non-inferiority mode:** allow the minimum tumor gain to be negative (a margin in Gy₁₀, e.g. −2), as in HYPNO's design. The current lower limit of 0 cannot express it.
4. **An evidence-anchored acute limit:** "acute BED of HYPNO 55/20" (delivered to 395 patients, 50.9% grade 3+ mucositis with 75% cisplatin) as a cap option beside "reference + excess", Fowler 59 and Fenwick.

## 4.2 Clinical review findings

An independent radiation-oncology reviewer (subagent, read-only) read app.R, about.md, the plans and `References/README.md`. It could not open the PDFs, so its numbers come from the README and the fix plan, which record full-text checks. Items marked "judgment" are clinical opinion, not verified facts. I checked the arithmetic of item C5.

**Verdict:** the arithmetic is sound, but the tool frames the question wrongly. It looks for a tumor BED gain over a radiotherapy-alone 70/35 reference, uses toxicity surrogates at the prescription dose, and has no anchor to the regimen actually tested (HYPNO 55/20 with cisplatin).

**Wrong or misleading statements now in the app**

| ID | Where | Problem | Fix |
|---|---|---|---|
| C1 | app.R:83–87, about.md:117 | "The reference regimen is assumed to have the same chemotherapy." False for shorter courses: in HYPNO the 55/20 arm received 4 vs 5 weekly courses (121 vs 149 mg/m² cumulative). | Chemotherapy per arm (Phase 5.3); state that fewer cycles lowers both tumor sensitization and mucositis. |
| C2 | about.md:87 | Says the 52 Gy₁₀ reference "is stricter than that tolerance zone". That compares a chemoradiation reference with the radiotherapy-alone zone. Fowler 2003 labels 70/35 (53.1) "Near tol?", and the chemotherapy-lowered zone starts at 54. | Reword: with chemotherapy the margin is about 1 Gy₁₀, not 7. |
| C3 | about.md:87 | Mucositis context cites Jacinto (n = 20) and Madhava via HYHOPE, but not HYPNO's randomized 50.9% vs 54.7%. | Lead with HYPNO; cite Madhava 2006 directly. |
| C4 | app.R:117–121 | Default criterion "DI or DD" and ranking by the larger gain list regimens only DI credits, under a tab titled "Regimens meeting the BED criteria". | Default to the smaller gain (or both models with a non-inferiority margin); show the DD shortfall on every row; rename the tab "Candidate regimens (model screen)". |
| C5 | app.R:129, app.R:124 | d_max 8 Gy and 5 fractions reach SBRT-like schedules for large mucosal and nodal volumes. The largest curative fraction size with precedent here is 3.1 Gy (HYHOPE 46.5/15). | Default d_max 3.5 Gy; flag rows above 3.1 Gy "beyond clinical precedent"; hard limit with warning at 5 Gy. |
| C6 | about.md:93, app.R:90–92 | The late limit (α/β 3, prescription dose, ≤ 70/35 = 116.7 Gy₃) allows 20 × 2.94 Gy (116.4 Gy₃, checked), well above HYPNO's tested 105.4 Gy₃. It also says nothing about mandible, constrictors, larynx, carotid, cord or brainstem. | Late-limit option anchored to HYPNO 55/20; α/β 2/3/4 side by side; OAR panel later (cord α/β 0.87). |
| C7 | app.R:331–335 | HYPNO note omits that its control was 66/33 at 6 fractions/week (not 70/35), that HPV was not tested, that 75% received cisplatin, and the 10-point margin. | Add those facts and WP2 2c. |

**Missing for trial design (ranked by the reviewer)**

1. Non-inferiority tumor criterion ΔBED ≥ −δ instead of "gain ≥ 1 Gy". The reviewer's conversion (judgment, uses Stuschke γ50 3.1): a 10-point locoregional-control margin near 70% control is about −3 Gy₁₀. That is smaller than the DI–DD disagreement the app already shows, which is itself a key output. **Verify the slope against the source before quoting it.**
2. Population-specific references: HPV-negative locally advanced with cisplatin; HPV-positive favorable risk (de-escalation references such as 60/30); postoperative; elderly or cisplatin-ineligible (radiotherapy alone, where 55/20 or 65/30 are the real comparators). Each preset sets the comparator, chemotherapy and acute limit.
3. Organ-at-risk fraction-size ceiling: the user gives an OAR's maximum dose as a fraction of prescription (e.g. 1.07 for a hotspot inside the PTV); the tool reports EQD2 at that OAR's α/β against conventional limits and the largest allowed fraction size for each fraction number.
4. Robustness across scenarios: DI/DD × literature presets × tumor α/β 6.4/10/14 × late α/β 2/3/4 × chemotherapy increment 0/1.8/4.4 Gy₁₀. Show the regimens that pass in all or most scenarios. This is the deferred uncertainty analysis in scenario form.
5. Phase I ladder generator: fixed fraction number (15 or 20), starting at HYPNO's acute BED, steps of about 1.5–2 Gy₁₀; HYHOPE's dose-limiting toxicity definition (grade 4, or grade 3 mucositis persisting at 3 months); a design that handles late toxicity (TITE-CRM or BOIN, judgment); Fowler 2003's advice to run an escalation study near 59–61 Gy₁₀.
6. Stopping-rule templates (judgment): grade 3 mucositis beyond 4 weeks, feeding tube at 3 months, consequential late effects, locoregional-control futility look.
7. Treatment-gap and start-day sensitivity (HYPNO's measured time difference was 11.5 days, against 12 assumed).
8. A volume caveat beside the results (no acute model includes volume).

**Benchmark rows to add** (numbers from References/README section F): HYPNO annotations (75% cisplatin 35 mg/m² weekly, HPV not tested, grade 3+ mucositis 50.9%/54.7%); HYHOPE 46.5/15 with cisplatin 40 mg/m² × 3, favorable-risk HPV-positive, de-escalation vs 60/30, median 18 days, grade 3 oral mucositis 21% (5/24); Jacinto 55/20 with cisplatin, grade 3 mucositis 40% (n = 20); Miller 2026 55/20 almost all without chemotherapy, 64% (72/113) at end of treatment, mean age 72; Madhava 55/20 with carboplatin, 84% (16/19). Accelerated anchors (UK 65/30, DAHANCA 6/week, CHART, RTOG concomitant boost) only after source checking.

**On the defaults:** D18's tumor DI values are clinically defensible. A Fowler 59 Gy₁₀ default acute limit is not defensible for a cisplatin trial: it is the radiotherapy-alone zone, about 3 Gy₁₀ above HYPNO's 55.9, which already produced about 51% grade 3+ mucositis. Recommended default: acute BED ≤ HYPNO 55/20, labeled "tolerated in a randomized trial with 75% cisplatin", with the chemotherapy-lowered 54 as the alternative.

**Surprises the reviewer noted:** Miller's real-world 55/20 without chemotherapy had more grade 3 mucositis (64%) than HYPNO with cisplatin (50.9%), so grading time and patient fitness swamp modeled differences of 1–3 Gy₁₀; HYHOPE dropped a level whose modeled acute BED matched the one it kept; and the fix plan has no item making cisplatin cycles an input.

**My reading.** The reviewer and my own analysis agree on the central points: HYPNO should be an anchor for both tumor comparison and toxicity limits, the gain filter should become a non-inferiority and trade-off view, and the acute surrogate is coarse (resolution of a few Gy₁₀ at best). Items 2, 3, 5 and 6 are new features beyond the current plan; I put 1, 2 (as presets), 4 and 7 into Phase 5 and the rest into later rounds.

# 5. Action plan

Rules for every code step (from your global rules): one committed test per change, observed failing without the change and passing with it; a mutation of the shipped line turns it red; the live app path runs the edited code; one `docs/FIXES.md` row per landed change; one commit per cluster citing the step ID; push to `v2-toxicity-hypofractionation` only; `main`, the live v1 app and deployment untouched.

## Phase 0. Before any code

| Step | What | How | Verify |
|---|---|---|---|
| 0.1 | Decision round | Ask the decisions in section 6 in one numbered message | Answers recorded in `docs/development/` with date and name |
| 0.2 | Scaffold | `python3 ~/.claude/install_auditors.py ~/.claude --init-project <repo>`: creates `evidence_ledger.yaml`, `sim_provenance.yaml` (N/A content; no simulated data), `audit_manifest.yaml`, `constant_registry.yaml`. Register the default constants (α, λ, Tk, Fowler zone, d_max, reference regimen) in the constant registry | `stale_constant_audit.py --project . --strict` exits 0 |
| 0.3 | Copy this document into `docs/development/2026_10_18_review_and_action_plan.md` | Plain copy | Linked from `docs/FIXES.md` |

## Phase 1. Correctness fixes (no new science)

Order matters: 1.1 changes absolute early BEDs, which every later acute test depends on; 1.2 must precede any comparator work.

| Step | Change | Files / how | Verify (test, red without the change) |
|---|---|---|---|
| 1.1 | Time convention (M2, D1). Add `t_offset` per tissue parameter set; early default offset 1; tumor offset from the preset's source | `bed_models.R`: `bed_early(..., t_offset = 1)` uses `T - t_offset`; same argument for `bed_di`/`bed_dd`; `default_params()` carries `t_off_e`, `t_off_tumor` | `test_bed_models.R::"early BED reproduces Fowler 2003 Table 1"` (53.1, 62.6 and 5 more rows); red on current code (52.3). Second test: tumor offset changes ΔDI for 55/20 by λ/α·1 |
| 1.2 | Separate DD anchor from comparator (M3, D4) | `compute_all(m, d, T, p, ref, anchor = default_ref())`; `dd_calibration(anchor, ...)`; UI: anchor panel, collapsed, default 35×2/47 | Test: changing `ref` leaves `bed_dd` of 55/20 unchanged; red on current code (58.0 vs 59.8 vs 70.1) |
| 1.3 | Published-value fidelity tests (M4), plus optional g (M6) | New file `tests/testthat/test_published_values.R`; `bed_di`/`bed_dd` gain `g = 0` | Each value in the M4 table; the 2018 shape test needs g, so it is red until g exists |
| 1.4 | Withers preset from α/β (1b) | `tumor_presets()` computes `0.6*(1+2/ab)*alpha`; `preset_matching()` includes α/β | Test at α/β 5: λ_DI/α = 0.84; red on current code (0.72) |
| 1.5 | Negative BED (1c, D3) | `suggest_regimens()` returns reason `negative_reference` per endpoint | Test with α 0.02: search returns that reason; red now |
| 1.6 | Non-inferiority margin | `LIMITS$min_gain` lower bound −15; help text "negative = allowed loss" | Server test: min_gain −2 accepted and reaches `suggest_regimens()`; red now (validation rejects) |
| 1.7 | Dose ceiling (M9, C5) | `d_max` default 3.5; rows above 3.1 Gy flagged "beyond clinical precedent"; warning at 5 Gy | Server test: flag column present for d > 3.1; warning shown at d_max 6 |
| 1.10 | Wording fixes C1–C3, C7 and tab rename (C4) | about.md and app.R help text only | Text diff; server test checks the new tab title |
| 1.8 | Export provenance (1d, D10) | CSV header lines `# key: value` with every parameter, anchor, comparator, app commit hash, date; absolute BEDs added | Round-trip test: parse header, recompute rows, equal to 1e-9 |
| 1.9 | Records (1e, 1f) | References/README Gupta provenance; about.md γ → g; OAR formula in the first plan | Text diff only |

Phase 1 exit check: full suite green; mutation checks (restore each old line in a tempfile copy, confirm red); one combined server smoke run of all tabs on the same commit; reference early BED 53.1 at defaults.

## Phase 2. Smooth threshold (if N2 = hard)

| Step | Change | Verify |
|---|---|---|
| 2.1 | Replace `smooth_max0` by `pmax(0, x)`; closed-form calibration `onset = T0 − y` | Test: at T = Tk + 0.1, DI loss equals λ·0.1/α exactly; the v1-agreement test is relaxed to "equal unless within 1 day of an onset", with the reason in the test |

## Phase 3. DD parameters from outcome data (research step, before any DD default changes)

| Step | What | How | Verify |
|---|---|---|---|
| 3.1 | Data | Transcribe the 33 arms of 2018 Table 1 (d, n, T, patients, observed TCP) into `docs/development/<date>_dd_refit/arms.csv`, double-entered and diffed | Two independent transcriptions agree |
| 3.2 | Fit | Binomial likelihood, TCP = exp(−K_trial·S), one K per trial (nuisance), shared repopulation parameters; fit (a) DD with α fixed at the chosen DI α and C, γ free; (b) DD with α free; (c) DI for comparison; g fixed at 0.0154 and at 0 | Recovers the published Table 2 values when stage-specific K is replaced by the paper's setup as closely as possible (sanity, not equality); profile-likelihood CIs |
| 3.3 | Check against HYPNO | Predicted Δ for 55/20 vs 66/33 under each fit | Reported beside the observed ≤1.4 points |
| 3.4 | Record | `docs/development/<date>_dd_refit.md` with the decision it informs | Reviewed by you before Phase 4 |

This replaces the DD half of D18. If you prefer not to refit, the alternative is to keep the DI-calibrated DD and label it plainly as the pessimistic bound (N1 option b).

## Phase 4. Defaults and presets

| Step | Change | Verify |
|---|---|---|
| 4.1 | Tumor DI defaults α 0.35, λ_DI 0.23 (Tp 3.0 d), Tk 28, α/β 10, offset per M2 (D18 DI part) | Test pins each value; reference DI 71.5 (or 70.8 with elapsed offset; recompute) |
| 4.2 | DD default from Phase 3 (or labeled pessimistic calibration) | Test: DD = DI at the anchor; HYPNO Δ shown on page |
| 4.3 | Presets: v1 defaults, Withers, Fowler 2010, Hartley 2011, HYPNO design (0.65 Gy/day), HYHOPE (D12) | `preset_matching()` covers every preset field incl. λ_DD and offset; test each preset calibrates |
| 4.4 | Default acute limit: HYPNO 55/20 acute BED, with concurrent chemotherapy checked by default (N4) | Server test: default search non-empty or explains itself |

## Phase 5. Trial-design features

| Step | Feature | How | Verify |
|---|---|---|---|
| 5.1 | Comparator selector: 70/35, HYPNO 55/20, HYPNO 66/33 (6/wk), custom | `ref` reactive driven by selector; anchor stays (1.2) | Test: comparator change moves Δ columns, not absolute DD |
| 5.2 | Trade-off plot | For each m, best d at a sweep of acute caps; plot ΔDI and ΔDD vs Δearly; mark HYPNO, Fowler zone | Snapshot of the data frame behind the plot; values equal section 4.1 table at the HYPNO cap |
| 5.3 | Chemotherapy as acute BED increment (D16) | Input in Gy₁₀, default 0; Meade scaling (4.4 Gy₁₀ per 300 mg/m² cisplatin) as help text; applied per arm, so a shorter arm with fewer cycles can carry less | Test: increment adds exactly to early BED of the arm it is set for |
| 5.4 | Late α/β sensitivity (WP4 4b) | Late BED at 2, 3, 4 Gy side by side; flag regimens whose late verdict flips | Test: 20 × 2.9 flips between α/β 2 and 3 |
| 5.5 | Model agreement (3d) and scenario robustness | Count of regimens passing under both models; run the search over the named presets and show for each regimen how many scenarios it passes (deterministic scenarios, not Monte Carlo, honoring the deferral of full uncertainty analysis) | Test: robustness count equals an independent loop over presets |
| 5.6 | Fenwick 2008 boundary option (D11) | Constant 32.2, T inclusive, ≤ 50 d | Fidelity test from M4 table |
| 5.7 | Late limit anchored to HYPNO 55/20 (C6, N8) | Option beside "reference + excess" | Test: 20 × 2.94 Gy fails under the HYPNO late anchor, passes under the 70/35 anchor |
| 5.8 | Default criterion and ranking (C4, N9) | Rank by the smaller gain; ΔDD shortfall column always shown | Server test of the defaults |
| 5.9 | Population presets (reviewer item 2) | Each sets comparator, chemotherapy increment and acute limit | Test each preset loads and searches |
| 5.10 | Gap and start-day sensitivity (item 7) | Add k interruption days to T; show ΔBED per day of gap | Test: one gap day after Tk lowers DI by exactly λ/α |

Later rounds add the organ-at-risk fraction-size ceiling, the phase I ladder generator and stopping-rule templates (reviewer items 3, 5, 6).

## Phase 6. Evidence and framing text (WP2, WP4, D17)

Every number traced to `References/README.md` section F or to a PDF page. Includes: HYPNO acute table and 10-point margin; HYHOPE 2025 outcomes with population caveat; Jacinto, Miller 2026, Madhava 2006; Fowler's advice on a small dose-escalation study near the zone; Stuschke ranges and extrapolation caveat; the statement that the app's DD calibration differs from the fitted 2018 model, with the table in 3.1; limitations (screening surrogates, chemotherapy, volume, OAR dose per fraction).

Verify: run `reader_first_audit.py` and the acronym and US-spelling checks on about.md; every number on the page has a row in `evidence_ledger.yaml`.

## Phase 7. Combined check and review

1. Full test suite, mutation checks, one end-to-end server test, headless-browser walkthrough of every tab, on one commit whose hash is written into the run log.
2. Independent review: mathematics (independent reimplementation incl. new features), radiation oncology, statistics of the refit.
3. Push to the branch. No deployment (D8).

## Later rounds (not in this plan)

Organ-at-risk panel (cord α/β 0.87 Gy, brainstem, mandible, carotid, larynx/constrictors, with dose per fraction from a fraction of prescription); full uncertainty analysis on joint fitted-parameter uncertainty from the Phase 3 fit; deployment under a new app name after `rsconnect::listDeploymentFiles()`.

# 6. Decisions needed (recommendation first)

Previous decisions D1–D17 keep the recommendations in the fix plan, except where changed below.

1. **N1, DD parameters:** (a) refit to outcome data, Phase 3 (recommended); (b) keep DI calibration and label DD as pessimistic for short courses; (c) approve D18 DD values as proposed (not recommended, section 3.1).
2. **N2, threshold:** hard max(0, x) (recommended) or keep the smooth version from v1.
3. **D18 split:** approve the DI part (α 0.35, λ_DI 0.23, Tk 28) now; hold the DD part for N1 (recommended).
4. **N4, default acute limit:** HYPNO 55/20 acute BED with chemotherapy on (recommended); Fowler 59 Gy₁₀ (D18 proposal); reference + 0 (current).
5. **N5, comparator default:** keep 70/35 as default with HYPNO selectable (recommended); or HYPNO as default.
6. **N6, non-inferiority mode:** allow negative minimum gain (recommended).
7. **N7, dose ceiling:** default 3.5 Gy, rows above 3.1 Gy flagged, warning at 5 Gy (recommended); or 4.5 Gy.
8. **N8, late limit:** add a HYPNO 55/20 late anchor as an option (recommended) or make it the default.
9. **N9, default ranking:** rank by the smaller of ΔDI and ΔDD and show the DD shortfall (recommended); or keep the larger gain.
10. **D4 (anchor):** now required, not optional (recommended yes, default 35 × 2 Gy / 47 d).
11. **Scope of this round:** Phases 0–2 and the wording fixes now; Phase 3 refit next; Phases 4–5 after the refit (recommended).

# Appendix A. Auditor triage

| Auditor | Verdict | Reason |
|---|---|---|
| methods-audit | RUN | This document (plan audit of model and calibration) |
| radiation-oncology-reviewer | RUN | Independent subagent, section 4.2 |
| radiation-epidemiology-review | RUN (partial) | Dose-response and repopulation plausibility, sections 3.1, 4.1 |
| statistician-review | RUN (partial) | Calibration logic; the Phase 3 refit needs a full review when built |
| paper-code-audit | RUN | Equations vs 2018 paper and Fowler, section 2 |
| surprising-result-skeptic | RUN | The "structural DD contradiction" claim, section 3.1 |
| stale-constant-audit | RUN later | No registry yet; Phase 0.2 creates it |
| pipeline-audit | N/A | No multi-stage pipeline; a Shiny app with unit tests |
| manuscript-audit, evidence-audit, reader-first-audit | N/A now | No manuscript; reader-first and evidence apply to about.md in Phase 6 |
| simulation-audit | N/A | No simulated data |
| biomedical-data-audit | N/A | No patient-level data |
| ml-practitioner-review | N/A | No learned model |
| tough-reviewer | N/A | No grant |
| subagent-scan | RUN | One clinical subagent with a fixed return brief |

Missing scaffold (not created, per your instruction not to execute): `evidence_ledger.yaml`, `sim_provenance.yaml`, `audit_manifest.yaml`, `constant_registry.yaml` in the repository. Phase 0.2 creates them.

# Appendix B. Reproduction

Scripts used for this review (scratch, not part of the app): `chk.py` (2018 fitted models on the trial pairs), `chk2.py` (Fowler reproduction, D18 numbers, time convention), `chk3.py` (regimens at HYPNO acute BED). They are in `review_scripts_2026_10_18/` next to this document (run with `python3 <script>`), and go to `docs/development/2026_10_18_review_scripts/` with this document in Phase 0.3.
