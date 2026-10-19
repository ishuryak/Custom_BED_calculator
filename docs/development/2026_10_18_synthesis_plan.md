# Synthesis of two reviews and the plan for version 2.1

Date: 2026-10-18 (late). This file is the controlling plan for version 2.1. It supersedes `2026_10_18_fix_plan.md` and `2026_10_18_handoff.md` as the starting point. Those stay as historical records.

## Inputs

1. **Review A (Claude):** `HNSCC_BED_calculator/REVIEW_AND_ACTION_PLAN_2026_10_18.md` (copied here as `2026_10_18_review_A_claude.md`), with its check scripts in `2026_10_18_review_scripts/`.
2. **Review B (Codex, second review):** `C:\Users\ishur\Documents\Codex\2026-10-18\rev-2\outputs\HNSCC_BED_second_review_and_action_plan_2026-10-18.md` (Stages 0–9, findings N1–N8, numerical fixtures). Not copied. Summarized below.

Igor Shuryak's instruction (2026-10-18): synthesize both reviews, implement the plan autonomously, commit and push to the branch, and do not publish the app live.

## Where the reviews agree

- The arithmetic of the DI, DD, acute and late formulas is correct and matches the 2018 paper's equations (without the growth term g) and Fowler 2003.
- The Fowler acute convention is elapsed days. The app's inclusive count makes acute BED 0.79 Gy₁₀ low, so the absolute Fowler cap is too lenient.
- The DD calibration must be anchored to a fixed regimen, separate from the comparison regimen.
- The Withers preset must follow the tumor α/β. A negative reference BED must block the search for that endpoint.
- The exports need full provenance. Late BED should be shown at α/β 2, 3 and 4. HYPNO must be compared with its own control (66/33 at 6 fractions/week).
- D18 is not a validated parameter set. The rejected 7-day DD onset floor stays out.
- The handoff's claim that the DD model is "structurally" contradicted by accelerated trials is too strong. Review A shows that the fitted 2018 DD model predicts gains for GORTEC (+3.2 Gy₁₀) and TROG (+2.1 Gy₁₀) and about no change for HYPNO, while the app's calibrated DD predicts losses. Review B shows that sign counts and Gy-versus-percentage-point comparisons cannot select a model.

## Corrections from Review B that Review A missed (all verified)

| Item | Check |
|---|---|
| HYPNO's 0.65 Gy/day is in EQD2 (2-Gy-equivalent) units: 66 Gy + 9 d × 0.65 = 71.85 Gy (HYPNO p. 389). In BED at α/β 10 it is 0.78 Gy₁₀/day. | Text read in the local PDF, lines "0.65 Gy/day … 71.85 Gy". Review A's "HYPNO design 0.65 Gy/day" wording is corrected. |
| Ignoring incomplete repair understates cell-kill BED. The analysis notes say it "overstates". | Algebra: the repair cross-term is nonnegative. |
| `bed_early()` is called with `T − 1` in two archived scripts. Converting inside the low-level function would double-convert. | grep of `dd_table.R:30`, `dd_table_05.R:30`. The conversion goes in `compute_all()`. |
| Strict-band grid: no combination passes all 13 regimens (maximum 10), not "only degenerate solutions". | Taken from Review B's recount. Recorded in the errata. |
| HYPNO exact outcomes: 3-year LRC 50.7% vs 51.2%, HR 1.098 (90% CI 0.93–1.30). Late grade 3+ 18.8% vs 20.2%, HR 0.926 (0.68–1.25). Mucositis 190/373 vs 208/380. Cisplatin 597/792. | Read in the local PDF (Table 3 and Results). |
| HYHOPE, Jacinto, Miller, Madhava denominators and time points. | Read in the local PDFs (7/24 any acute grade 3, 5/24 oral mucositis, 26.6 months follow-up. Jacinto 40%, 29 days median. Miller 72/113 at end of treatment, 119 referred. Madhava 16/19). |
| Numerical fixtures (legacy and D18 BEDs, onsets, C/α, Withers at α/β 5 = 82.04, Fowler 2010 examples 62.79 and 60.41, matched λ_DD 0.6655, g effects 3.619 and 2.002 Gy, late reversal of 20 × 2.9 Gy). | Reproduced by an independent Python implementation (`2026_10_18_review_scripts/fix.py`). |

## Where the reviews differ, and what version 2.1 does

| Topic | Review A | Review B | Version 2.1 (and why) |
|---|---|---|---|
| Tumor defaults | Approve the DI half of D18 | Keep v1 values as default. D18 as a named scenario | **Defaults unchanged.** D18 and the source sets become presets. A default change alters every number on the page and needs Igor's explicit decision. Presets give the same access with no silent change. |
| Default acute limit | HYPNO 55/20 acute BED | Keep the reference-relative limit. Fowler cap as a separate policy. An empty list is informative | **Reference-relative stays default.** HYPNO 55/20 is added as a selectable policy, beside Fowler. Same reasoning. |
| DD parameters | Refit to the 2018 trial arms | Separate model-fitting project | **Not in 2.1.** Both reviews say a fit is new scientific scope. The page states that the app's DD calibration differs from the fitted 2018 model and shows the consequence. |
| Hinge | Hard max recommended | Keep softplus for compatibility. Record it | **Keep softplus**, record it in exports, test that published examples reproduce. |
| Export | Header block in the CSV | ZIP with rectangular CSV, scenario JSON, README | **ZIP**, as Review B: a CSV with a header block is not a normal CSV. |
| Default ranking | Smaller of ΔDI and ΔDD | Expose disagreement. Best-model sort is an optimistic option | **Ranking default unchanged**, labeled "optimistic", plus a "Credited by" column (DI and DD, DI only, DD only) and the smaller delta on every row. |
| Dose ceiling | 3.5 Gy default | Not discussed | **Default 5 Gy**, and every row above 3.1 Gy (HYHOPE, the highest curative fraction size with precedent here) is flagged. Default results are unchanged (all rows are below 3.1 Gy). |
| Non-inferiority | Allow a negative minimum gain | Not discussed | **Allowed** (down to −20 Gy), labeled as an allowed loss. |
| Chemotherapy increment, Fenwick boundary | Phase 5 | Stage 8, optional, individually gated | **Deferred.** The existing checkbox keeps its meaning. |
| Trade-off plot | Phase 5 | Not discussed | **Deferred.** The HYPNO acute policy plus the new sensitivity tab cover the main use. |

## Scope of version 2.1 (in commit order)

Every step: tests written first and observed failing on the pre-change code (`6f0209c`), then passing. A `docs/FIXES.md` row. The app's own entry point exercised by server tests.

1. **Acute clock.** `compute_all()` converts inclusive T to Fowler's elapsed T once (`acute_model_time()`), recorded as `clock_e = "elapsed"`. The low-level `bed_early()` keeps its time argument on the model clock. Tests: Fowler 2003 53.1 and 62.6. Fowler 2010 62.79 and 60.41. Reference acute BED 53.10544. 22 × 2.75 Gy fails the 59 Gy₁₀ cap (59.71).
2. **Fixed DD anchor.** `compute_all(m, d, T, p, anchor)` and `suggest_regimens(p, ref, anchor = default_anchor())`. The app has an anchor panel (default 35 × 2 Gy, 47 days) and a "legacy: calibrate at the comparison regimen" option. Comparison regimen presets: custom (default 35 × 2 Gy), HYPNO test 55/20 (26 d), HYPNO control 66/33 at 6 fractions/week (38 d). Tests: changing the comparison leaves absolute DD BED unchanged. Legacy mode reproduces 58.0214 and 70.1242.
3. **Presets.** Every preset sets every field it claims. Legacy v1 (all five tumor values), D18 scenario (all five), Withers rule (λ_DI/α = 0.6(1 + 2/(α/β)), follows α/β), Hartley 2011 (α 0.3, α/β 10, Tp 3 d, Tk 22 elapsed = 23 inclusive), HYHOPE protocol (α 0.35, α/β 10, Tp 3 d, Tk 21). DI-only presets leave λ_DD unchanged and say so. Tests: Withers at α/β 5 gives reference BED 82.04. Matching covers every field.
4. **Search validity and trial-design options.** Negative or nonfinite reference BED blocks the endpoint with a reason. Membership columns. Smaller delta. Precedent flag. Negative minimum gain. Acute policy "HYPNO 55/20". Late policy "HYPNO 55/20". All qualifying rows kept for export. Tests: λ_DI 1/day at α 0.2 gives reason `negative_reference_di` instead of a +78 Gy apparent gain. Search counts legacy strict 7, Fowler 59 17, D18 strict 0, D18 Fowler 11, chemotherapy 54 0.
5. **Background growth g** (default 0) in both tumor models, also when DI repopulation is 0. Tests: g 0.0154 at α 0.2 lowers 47-day BED by 3.619 Gy and 26-day BED by 2.002 Gy. Calibration unchanged.
6. **Export ZIP** (`results.csv` full precision, `all_qualifying.csv`, `scenario.json`, `README.txt`), deterministic scenario ID. Round-trip test reruns the search from `scenario.json` and reproduces the rows within 1e-8 Gy.
7. **Results and sensitivity.** Late BED at α/β 2 and 4 on the Results tab. Benchmark rows (HYPNO test vs its own control, HYHOPE) with verified context. A Sensitivity tab: the evaluated regimen under each tumor preset at tumor α/β 7, 10 and 15, with ΔDI, ΔDD, Δacute and Δlate at α/β 2, 3 and 4. Labeled "selected scenarios, not a probability".
8. **Text.** about.md (anchor, clock, presets with sources and units, HYPNO numbers, the fitted-versus-calibrated DD statement, chemotherapy wording, limitations), README (version 2.1, legacy warning), References/README (Gupta provenance, HYPNO 0.65 in EQD2 units), errata for the two DD analysis notes.
9. **Release checks.** Full suite, mutation checks of the shipped lines, server tests, headless-browser walkthrough of every tab with screenshots, push to `v2-toxicity-hypofractionation`. No deployment.

## Not in version 2.1 (recorded so nothing is silently dropped)

- DD refit to outcome data (Review A Phase 3) and any DD model revision.
- Chemotherapy acute increment q, Fenwick boundary, repair-aware schedules, 2–3 fractions/day.
- Organ-at-risk panel. The correct inverse for a later OAR module is d = r/(2f)·[√(1 + 4B/(n r)) − 1], with no f² inside the root (Review B Stage 1.6. The first plan document has the wrong form).
- Portable rebuild of the archived DD analysis scripts (Review B Stage 7). They are frozen records of commit `6f0209c`. The errata say which statements are wrong.
- `renv.lock` and a test-only CI workflow (Review B Stage 9.4).
- Phase I ladder generator, stopping-rule templates, population presets (clinical reviewer items 2, 5, 6).
- Deployment (needs `rsconnect::listDeploymentFiles()` and Igor's approval).

## Decisions taken without Igor (autonomous instruction), to confirm later

| ID | Decision | Reversible by |
|---|---|---|
| A1 | Tumor defaults stay at v1 values. D18 is a preset | Changing `default_params()` and the UI defaults |
| A2 | Default acute and late limits stay reference-relative | Changing the radio-button defaults |
| A3 | Softplus hinge kept | Not reversed in 2.1 |
| A4 | Default d_max 5 Gy (was 8). Rows above 3.1 Gy flagged | One UI default |
| A5 | Export is a ZIP | `downloadHandler` |
| A6 | Default anchor 35 × 2 Gy / 47 d. Legacy coupling available | Anchor panel |
