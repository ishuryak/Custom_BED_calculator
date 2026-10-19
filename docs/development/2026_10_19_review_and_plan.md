# Version 2.4 review, verification of the external review, and the version 2.5 plan

19 October 2026. Reviewed snapshot: branch `v2-toxicity-hypofractionation`, commit `ad16504` (app 2.4.0, equation version 2.2). Reviewer: Claude (Opus 5.5) for Igor Shuryak. External review verified here: `docs/development/2026_10_19_external_review_codex.md` (copy of `HNSCC_BED_v24_review_and_fix_plan_2026-10-19.md`, findings F01–F08).

## 1. Bottom line

The DI, DD, acute, Fenwick and late BED equations are implemented correctly for the calculator's declared model, and every published or source value I rechecked reproduces. I found no arithmetic error in the formulas. Every counterexample in the external review reproduces exactly.

What must change before the agreement output is used to support a trial proposal:

1. **Validity (F01).** The Sensitivity and Model agreement tabs count scenarios in which a comparison BED is negative. In the external review's stress case, a regimen is reported as "15 of 15" when only 3 scenarios are evaluable. The main search already blocks this. The two newer tabs do not.
2. **Claims (F02).** Several sentences turn BED arithmetic into clinical statements the calculator cannot support ("within the parameter uncertainty", "close to HYPNO in tumor effect", "pessimistic bound", "non-inferiority margin" for a BED threshold, "can only lose dose" before the DI onset). I agree with the external review on all but one: the D18 DD deficit against HYPNO can be judged against the trial, conditionally, with a cited dose-response slope (section 4, N1). That sentence is restated, not deleted.
3. **Tab independence, export and replay (F03–F05).** The agreement tab fails when unrelated candidate-tab fields are cleared. The agreement analysis cannot be exported. Replay accepts files from unsupported model versions and the scenario ID depends on key order.
4. **Labels, provenance, records (F06–F07)** and **release checks (F08).**

No equation, default parameter, onset law or default acute policy needs to change. The plan in section 6 keeps the D18 defaults and the HYPNO acute limit that Igor chose on 2026-10-19.

## 2. What I checked

- Read `R/bed_models.R`, `R/export.R`, `app.R`, `about.md`, `References/README.md`, the seven test files and `docs/FIXES.md` in full.
- Baseline suite on `ad16504`: **1504 expectations, 0 failed, 0 errors, 0 skipped** (Windows R 4.5.1. Warnings: testthat, shiny and DT built under R 4.5.3).
- Probe script `docs/development/2026_10_19_review_scripts/probe_v24.R` (output reproduced in section 3).
- Sources reread for the points at issue: Shuryak 2018 (own copy, Publications/Articles/accelerated_repopulation.pdf), HYPNO (Bentzen 2026), Stuschke 1999, Fowler 2010, the HYHOPE protocol title and signature pages, Fenwick 2008 abstract formula.

## 3. Verification of the external review

| ID | External finding | My verdict | Evidence (probe_v24.R on `ad16504`) |
|---|---|---|---|
| F01 | Sensitivity/agreement count scenarios with a negative comparison BED | **Confirmed** | λ_DI 0, 35 × 2 Gy over 150 days: 32 × 2.10/44 d and 33 × 2.05/45 d listed with 15/15. Only 3 scenarios have positive DI and DD comparison BED. Version 1 defaults at α/β 10 give comparison DI −38 Gy and ΔDI +103.31 Gy, counted as a pass |
| F02 | Unsupported outcome and uncertainty claims | **Confirmed, one modified** (N1) | g counterexample: 17 × 3.05 vs HYPNO, ΔDI −2.4602 (g 0), −2.3282 (g 0.0154), −1.6031 (g 0.1), although both end before day 28. α/β 2: ΔDI +0.2968 at acute 55.78 < 55.87 Gy₁₀ |
| F03 | Agreement depends on candidate-tab inputs. Export depends on the Results regimen | **Confirmed** by reading `app.R:489-510, 600-609, 644-654` (agreement calls `search_args()`, which validates `min_gain`, `min_days`, `m_range`; `scenario()` reads `regimen()`). Server tests added in step 4 reproduce it | |
| F04 | The agreement analysis cannot be exported. Filename says 2.1 | **Confirmed** (`app.R:611-617`) | |
| F05 | Replay ignores versions, hinge, clock. ID depends on key order | **Confirmed** | Equation 9.9 + hard hinge + elapsed tumor clock replays to the same 3 rows. Reversing `params` changes the ID (a98024bb355b → 0f3766cc8351) |
| F06 | HYPNO acute cap labeled Gy₁₀ at mucosal α/β 5. Preset provenance overstated. Fenwick comment says elapsed days | **Confirmed** | Cap 70.99 BED Gy at (α/β)_e 5, labeled "Gy₁₀" (`app.R:528`); `R/bed_models.R:68-71` |
| F07 | Stale source catalogue, versions and the 2018 dataset size | **Confirmed** | 2018 paper p. 21: "16 publications ... 7283 patients treated with 35 different treatment regimens" (about.md says 33 arms). HYHOPE local protocol: "Version 4.2" on the signature page (about.md says 4.1). `References/README.md` A1 says Fowler 2003 "Full text not seen" though section F reads it. E1 says g is omitted (it is implemented). Fowler 2010 "0.7 Gy/day" is "0.7 Gy EQD2 per day" in the source text. `DESCRIPTION` and file headers say 2.1 |
| F08 | No dependency lock. Launcher tied to one R path. Browser walkthrough needed | **Confirmed** | `run_app.bat` hard-codes `C:\Program Files\R\R-4.5.1` |

Fixture values in the external review's sections 5 and 10 all reproduce (D18 35/2: 71.5143/71.5143. 20 × 2.75: 70.1245/62.8473. 33 × 2/38 d: 72.6286/68.6403. HYHOPE 15 × 3.1: 60.9150/56.3705. Acute 53.1054, 55.8660, 55.4350, 52.2011 Gy₁₀. Fenwick T 47/26/19: 98.574, 74.169, 70.203). The fitted 2018 DD model reproduces +3.17 (GORTEC), +2.06 (TROG) and +0.54 Gy₁₀ (HYPNO).

## 4. Additional findings

**N1. The D18 DD conflict with HYPNO can be stated conditionally, with a source (modifies F02).** The external review says compatibility "cannot be determined from the BED difference alone". True as a strict statement, but a slope makes it checkable. Stuschke and Thames 1999 (Table 1 and Discussion 4.1) estimate γ₅₀ = 3.2 (95% CI 1.6–4.8) for HNSCC tumor control from the randomized hyperfractionation trials, and 0.8 ± 0.4 from RTOG 83-13. HYPNO's local control (about 51%) is near the middle of the curve, where γ₅₀ applies. Under D18, test minus control is −2.50 Gy₁₀ (DI) and −5.79 Gy₁₀ (DD), which are −3.4% and −8.4% of the control EQD2. With γ₅₀ 1.6–3.2, that predicts a loss of about 5–11 points (DI) and 13–27 points (DD). HYPNO's upper 90% confidence limit of the hazard ratio, 1.30, corresponds to about a 9-point loss from 51.2%. So with the trial-based slope range, the D18 DD prediction lies outside what HYPNO's result allows, and the DI prediction is at or inside it. With the low RTOG 83-13 slope, both are compatible. This is a linearized, illustrative translation (it ignores chemotherapy, heterogeneity between the populations, and HR-to-percentage conversion error). **Action:** restate about.md's D18 sentence in this conditional form with reference [3]. Drop "pessimistic bound". This is a judgment that goes beyond the external review and should be confirmed by Igor (decision A1 in section 7).

**N2. "The margin does not support escalation beyond 55 Gy in 20 fractions" (`app.R:464`, `about.md:119`).** HYPNO does not say this, and the sentence has no clear meaning (escalation of what, and supported how?). Replace with: "HYPNO tested one schedule. Its result does not by itself establish any other schedule."

**N3. DD onset shows "Inf" on the agreement tab when λ_DI = 0** (no repopulation in either model). Probe: three rows with `dd_onset` = Inf, printed by `sprintf("%.1f")` as "Inf". Show "none".

**N4. Stale comment** in `benchmark_regimens()` describing an `own_control` column that does not exist.

**N5. Results table units.** Toxicity rows say "(Gy)" with no α/β. Label them with the α/β used (part of F06).

**N6. Coupled (legacy) anchor in the agreement tab.** In legacy mode the agreement tab is calibrated at the sidebar comparison, not at its own comparison. Decision (A2): in legacy mode each analysis calibrates at its own comparison, which is what the option's label says, and each tab states the anchor it used.

Checked and found correct (no action): HYHOPE tumor onset 21 days is on the protocol's inclusive clock (the protocol counts 15 fractions as 19 days), so the preset needs no conversion. The HYHOPE 12-fraction level gives 51.32 Gy₁₀ on that clock at 19 days, as about.md states. HYPNO's 10-point margin applies to both co-primary endpoints (HRc 1.31 for loco-regional control). The HYPNO control arm's 38 days (6 fractions/week, Monday start) is correct. The default candidate rows and agreement rows match the external review's independent enumeration.

## 5. Clinical plausibility

- **Default candidate list (vs 35 × 2 Gy).** 29 × 2.25 Gy/39 d, 28 × 2.30/38 d and 30 × 2.20/40 d gain about +1.1–1.2 Gy₁₀ under DI and lose 2.0–2.9 Gy₁₀ under DD. These are modest, 1-week shortenings with fraction sizes well inside common practice. They are DI-only candidates: under the DD model shortening a 7-week course costs tumor effect. The list is a statement about the DI model, and the "Credited by" column says so.
- **Agreement rows (vs HYPNO 55/20).** 17 × 3.05 Gy/23 d, 18 × 2.90/24 d and 19 × 2.80/25 d lose 1.5–2.8 Gy₁₀ under both models and save 1–3 days. Clinically, a 1–3-day shortening of a 4-week course is a small gain to set against a modeled loss of tumor effect. The case for a trial would have to come from logistics (fewer visits) rather than from expected efficacy. 17 × 3.05 Gy sits at HYPNO's modeled acute BED, where HYPNO observed 50.9% grade 3+ mucositis, and fraction sizes of 3.05 Gy have curative precedent only in the favorable-risk, reduced-volume HYHOPE setting among the benchmarks. A shorter course also fits fewer weekly cisplatin doses (HYPNO: 4 vs 5 planned).
- **D18 DD.** See N1. With the trial-based slope, D18's DD rate is too pessimistic about HYPNO's 4-week schedule to be fully consistent with that trial. That is a reason to report DD results under D18 as one scenario, not as a constraint. The DI model under D18 is consistent with HYPNO.
  - *Annotation added 2026-10-19 (version 2.6, review R4): superseded. Both sentences assume that Stuschke's γ₅₀ applies unchanged on the calculator's time-corrected EQD2 scale, which neither source establishes. Under that assumption the illustrated DI loss is 5.5–10.9 points (γ₅₀ 1.6–3.2), not shown to be consistent with HYPNO, and the illustrated DD loss exceeds the 9.3-point proportional-hazards conversion only for γ₅₀ above about 1.1. See `2026_10_19_v26_review_response.md`. The text above is kept as written on the day.*
- **Acute model.** HYPNO's two arms have modeled acute BEDs of 55.87 and 55.43 Gy₁₀ and observed grade 3+ mucositis of 50.9% and 54.7%: the model's near-equality is consistent with the trial, but the ordering is reversed by 0.4 Gy₁₀, which shows the resolution limit in practice.
  - *Annotation added 2026-10-19 (version 2.6): "resolution limit" is not established. A 0.4 Gy₁₀ ordering reversal against a 3.8-point difference in observed rates is one comparison of two arms with different schedules and chemotherapy exposure; it shows the model does not order these two arms, not a measured resolution of the model.*
- **Late BED** is a prescription-dose surrogate. It says nothing about cord, constrictors, larynx or mandible. Unchanged limitation.

## 6. Common plan of action (merges F01–F08 and N1–N6)

Order follows dependencies. Every change gets a test that fails on `ad16504` for the stated reason and passes after the change. A mutation of each shipped line must turn a test red.

**Step 1. One validity rule for every analysis (F01).**
- What: add `tumor_endpoint()` in `R/bed_models.R` returning, per model, the raw comparison and candidate BED, `valid` (both finite and > 0, model calibrated), a reason (`dd_unavailable`, `negative_reference`, `negative_candidate`) and the delta (NA when invalid). Use it in `suggest_regimens()`, `sensitivity_table()` and `agreement_regimens()`.
- Agreement returns `n_scenarios`, `n_valid_both`, `n_pass_both`, `n_invalid`. Ranges use evaluable scenarios only (NA when none). The tab shows "3 of 3 evaluable (12 of 15 not evaluable)".
- Sensitivity table shows "not evaluable (negative comparison BED)" and the summary counts evaluable scenarios.
- Verify: the 150-day counterexample gives 15 total, 3 valid in both, 3 passes, 12 not evaluable. The default three agreement rows stay at 15 of 15. All-invalid gives NA ranges.

**Step 2. Clinical text (F02, N1, N2).** Rewrite the listed sentences in `about.md` and `app.R`: "allowed BED loss" instead of "non-inferiority margin". No "within the parameter uncertainty". Conditional pre-onset explanation with g. D18 sentence per N1. HYHOPE example without a "1–2 Gy₁₀" resolution claim. Hartley 15% as an aggregate observation. N2 sentence. Verify: a test greps the rendered About page and tab text for the withdrawn phrases.

**Step 3. Units, provenance and labels (F06, N3–N5).** Early-BED unit follows (α/β)_e everywhere (search note, cap note, Results rows, agreement table). "Beyond precedent" becomes "Above benchmark fraction sizes (> 3.1 Gy)". Fenwick comment fixed. "none" for an infinite DD onset. A provenance table of every preset value (direct estimate, conversion, chosen value, inherited, analyst scenario) in `about.md` and `References/README.md`. The Sensitivity tab states the inherited DD rate, g and anchor.

**Step 4. Independent tabs (F03, N6).** Split `search_args()` into shared toxicity settings and candidate-only settings. Agreement validates only its own inputs and explains its own empty results with the shared message function. The export no longer reads the Results regimen except as optional context. Legacy coupling calibrates each analysis at its own comparison and each tab shows its anchor. Verify: server tests clearing candidate gain, days and range keep three agreement rows. Clearing the Results fraction count keeps the export. Changing the sidebar comparison does not change agreement in fixed-anchor mode.

**Step 5. Export and replay (F04, F05).** Schema 2: `analysis` = candidate or agreement; `selection_version`; `scenario_set` (definitions of all 15 scenarios with inherited values). Code identity (SHA-256 of the three source files). Agreement ZIP adds `scenario_detail.csv` (candidate × scenario, raw BEDs, validity, reasons). Canonical ID: object keys sorted recursively, arrays kept in order, UTF-8, LF. `rerun_scenario()` rejects unsupported schema/equation/selection versions, hinge and clocks, checks types and ranges, recomputes the calibration and reports a mismatch with the saved one. Schema-1 candidate files migrate explicitly. Download names use `APP_VERSION`. Verify: fresh-process round trips of both analyses (rows, all qualifying rows, scenario detail, within 1e-8 Gy). Each unsupported field rejected. Reordered keys keep the ID. A changed scientific input changes it. A corrupted calibration is detected.

**Step 6. Sources and versions (F07).** about.md: 35 regimens from 16 publications, 7283 patients. HYHOPE protocol version 4.2. References/README.md: Fowler 2003 full text read. G implemented. Fowler 2010 rate in EQD2 per day. Parameter-origin ledger. One version source: `APP_VERSION` drives the UI title and download name, and a test checks `DESCRIPTION`. Headers updated. `docs/FIXES.md` row for 2.5.

**Step 7. Release checks (F08).** `renv.lock` (R 4.5.1 package versions as installed), launcher that finds Rscript, full suite, mutation checks, combined smoke (full app via testServer and a headless-browser walkthrough of every tab and both downloads), push to the branch. **No deployment.**

Deferred (separate scientific work, as both reviews agree): DD refit to outcome data, chemotherapy effect, OAR/NTCP module, repair-aware or twice-daily schedules, a robust-search mode that optimizes over scenarios.

## 7. Decisions taken without the PI (confirmed)

**Igor Shuryak confirmed A1–A5 on 2026-10-19, after reviewing the implemented version 2.5.**

- **A1.** D18/HYPNO sentence restated conditionally with Stuschke's γ₅₀ (N1) rather than removed.
- **A2.** Legacy coupling calibrates each analysis at its own comparison (N6). Fixed anchor stays the default, so default results do not move.
- **A3.** A scenario counts as evaluable only when both tumor models are usable for both regimens. The "pass" count is among evaluable scenarios, and the tab always shows how many were not evaluable. Rows are not removed for having non-evaluable scenarios.
- **A4.** The agreement margin stays a BED margin applied at each scenario's α/β. The export states the EQD2 equivalents. No fixed-EQD2 option is added now.
- **A5.** App version 2.5.0. Equation version stays 2.2 (no formula changed). New selection version 2.5 (validity rule extended to sensitivity and agreement).

## 8. Implementation record

All seven steps were implemented on 2026-10-19 in commits `73b54a5` and `725514a` (branch `v2-toxicity-hypofractionation`). The second applies the validity rule to the Results-tab, benchmark and HYPNO-note differences, which still showed raw differences against a non-positive BED; I found it while searching for other copies of the old logic. Nothing was deployed. Evidence rows: `docs/FIXES.md`, version 2.5.

- **Tests first.** `tests/testthat/test_v25_review_fixes.R` (26 tests) was written before the code changes and run against `ad16504` in a separate worktree: 25 failed on the reproduced defects, and one guard (fixed anchor ignores the sidebar comparison) passed, as it should.
- **Suite on `725514a`:** 1651 expectations, 0 failed, 0 errors, 0 skipped (baseline 1504).
- **Mutations:** all 14 mutations (`docs/development/2026_10_19_review_scripts/mutations.R`, output in `mutations_output.txt`) turn at least one test red on `725514a`: sensitivity validity bypassed, invalid scenarios counted as passing, wrong agreement comparator, candidate inputs leaking into the agreement, replay ignoring hinge or equation version, key-order-dependent ID, wrong unit label, sidebar anchor in the agreement, unchecked calibration, export requiring the Results regimen, negative candidate credited, Inf onset printed, and the old display rule.
- **Live in the app:** headless-browser walkthrough on `725514a`, 23 of 23 checks (all tabs, both downloads, cleared candidate input, legacy coupling, unavailable DD, custom mucosal α/β, MathJax). It found one more defect, the agreement table's last column cut off at the card edge, fixed by horizontal scrolling and pinned by a test.
- **Defaults unchanged:** the probe after the fixes (`2026_10_19_review_scripts/probe_after_output.txt`) reproduces the default candidate and agreement rows to every printed digit. The stress case now gives 3 evaluable passes and 12 non-evaluable scenarios.
- **Existing tests changed** only where they pinned a renamed label or message (`test_app_server.R`: the fraction-size column and the agreement empty message) or a renamed column (`test_v24_agreement.R`: `n_both_within` → `n_pass_both`).

Not done, by design: deployment; the deferred scientific work in section 6. A1–A5 confirmed by Igor Shuryak on 2026-10-19 (section 7).
