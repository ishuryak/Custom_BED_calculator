# Plan: act on the 2026-10-25 review and move v3 to the verified-trial refit

Status: APPROVED. Igor Shuryak accepted all 12 recommendations in section 5 on 2026-10-25 ("proceed with your recommendations").

Review: `HNSCC_BED_calculator/HNSCC_BED_v3_comprehensive_review_and_action_plan_2026-10-25.md` (41 pages, external Codex review; evidence archive `HNSCC_BED_v3_review_evidence_2026-10-25.zip`). Repository state at planning time: branch `refit_trials_2026_10_22`, HEAD 677ce21 (local, unpushed). Live v3 runs f951cd3.

## 1. Goal

Igor (2026-10-25): v3 should use the verified-trial refit, not the older fit to the 35 arms of the 2018 Table 1.

Today the app takes three things from that older fit:

| App element | Source now | Source after this plan |
|---|---|---|
| Uncertainty tab draws | `draws/tumor_{di,dd}_{primary,paper_restrictions}.csv`, free DI and free DD, 35 arms | DI from the refit main set; DD = calibrated DD paired row by row with its parent DI draw |
| "Fitted values" preset and fitted-DD mode | `draws/refit_point_values.json` | Medians from the new main fit |
| Report card refit column | the same point values | the same, with a dataset/date/model badge |

## 2. What the review says about this goal

The review agrees with the direction but sets gates first. It asks for corrected source data, a fixed fit cache identity, full diagnostics and a stated likelihood before any new draws ship. It also lists app defects that are independent of the fit and that matter for any draws.

I agree with nearly all of it. Three points change in our favor if the app switches to calibrated DD:

- **M1 (unstated DI/DD coupling) mostly goes away.** Each calibrated DD draw is fitted conditional on one DI draw, so the row pairing is defined by the model rather than by file order. It becomes a declared modular coupling, p_DI(θ|y)·p_DD(φ|θ,y). The toxicity draws are still paired independently, and that pairing must be declared.
- **M7 (shipped DD divergences) goes away for the app.** The free DD fits are divergent (17 in main, 62 in chemo) and not identified. The calibrated DD main run had 0 divergences. 11 of its 200 runs have λ_DD R-hat above 1.01, so it still needs more iterations.
- **α/β is not identified by these data.** Its posterior (10.7, 95% 6.6–19.2) is close to its prior (lognormal around 10 Gy, 95% about 5–20). The point preset should keep 10 Gy, as the current export does.

One point I would qualify. M6/C01 asks for an estimator-aware likelihood. Most arms report only a KM or CIF percentage, with no SE, CI or risk set. A full estimator-aware likelihood is therefore not possible for most of the data. I propose to label the binomial likelihood as an approximation and add an effective-N sensitivity analysis (decision 5).

## 3. Triage of every finding

Verdicts: **accept** (do as the review says), **modify** (do it differently, reason given), **defer** (later batch, not a release blocker for the refit).

| ID | Finding (short) | Verdict | Stage |
|---|---|---|---|
| S9 | Clearing an enabled organ's α/β turns a block into "no check" | accept | A |
| S1 | Old uncertainty rows exported under new inputs | accept | A |
| S2 | Uncertainty export drops applied OAR limits, so replay changes toxicity share 0 to 1 | accept | A |
| S8 | Fixed-toxicity uncertainty ignores organ α/β controls (hidden 2/3) | accept | A |
| S3 | Replay instructions omit dependencies | accept | A |
| S4 | Replay accepts incompatible metadata | accept | A |
| S5 | Skipped-day policy differs between searches | accept, one shared policy | A |
| S6, S7 | Release Git identity, direct calendar API validation | accept | F, A |
| OAR-01 | Physical mean-dose limits compared as EQD2 | accept (decision 12) | A |
| OAR-02 to 05 | HYHOPE PRVs and metrics, no invented D0.03cc, range ends, version metadata | accept; unassessed mandatory metrics listed visibly | A |
| D1 | DVH maximum taken from bin midpoints | accept: keep interval bounds, indeterminate when the top bin straddles the limit | A |
| D2, D3 | DVH validation; quoted identifier headers pass the filter | accept | A |
| M1 | Cross-model coupling | modify: resolved by calibrated DD pairing; toxicity pairing declared | E |
| M2 | Variable-α/β BED labeled Gy₁₀ | accept (decision 8) | E |
| M3 | Conservative quantile drops invalid draws | accept: invalid draws count as failures; result "unavailable" when the quantile is not evaluable | A |
| M4 | Full-draw toxicity check does not confirm eligibility | accept: re-check on all draws; failing choices marked, next candidate evaluated | A |
| M5 | Fit cache reused by file name, provenance overwritten | accept: run identity hash, refuse mismatched resume | C |
| M6, C01 | Rounded KM/CIF rates treated as counts | modify (decision 5) | B, D |
| M7 | Divergences; incomplete diagnostics | accept: full diagnostic gate (decision 9) | C |
| M8 | MCSE "at most" claim | accept: quantity-specific MCSE from chain IDs, or labeled IID benchmark | E |
| M9 | Day-7 onset is a soft penalty | accept: wording plus fitted vs evaluated violation counts | C, E |
| M10 | LOTO is a conditional plug-in score | accept: label it; add grouped study holdout (C08) | C |
| M11 | Intercept prior links trial means to repopulation | modify: intercept-prior sensitivity (N(0,3) vs N(0,10)) plus wording; no re-centering | D |
| C02 | Endpoints pooled across estimands | accept: keep `estimate_type`/`endpoint_definition` in the fit CSV; KM-only and CIF-only sensitivities | B, D |
| C03 | RTOG 9003 rates and denominators from different reports | accept (decision 4) | B |
| C04 | KROG HYPO weighting omits one T2 patient (n 28.219 should be 28.243) | accept | B |
| C05 | Mixed fraction sizes averaged | accept (decision 7) | B |
| C06 | Clock conventions mixed | accept (decision 6) | B, D |
| C07 | "Common chemo policy" is not equal exposure | accept: relabel; drop the α-direction remark | E |
| C08 | Holdout by unit, not study | accept: CAIR, Yamazaki and RTOG 9003 held out whole | C |
| C09 | Agreement and onset wording | accept: the review's replacement sentences | E |
| C10 | Quote checks are not semantic checks | accept: component-sum and cohort-match asserts in the builder | B |
| Sec. 5 | "Eligible" uncertainty rows confused with tumor pass (5×5 Gy example) | accept: tumor-pass status column, new wording | A |
| Sec. 5 | Source page "no dose-response" and "never stored" text | accept | E |
| Sec. 6 | Research-hypothesis cover sheet in exports | defer, after release | later |

Nothing in the review is dismissed. The fixed anchor, D18, the both-model criterion, the +1 Gy threshold and the conservative ranking all stay, as the review recommends.

## 4. Stages

Each change ships with a test that fails before the change and passes after it. The fail is observed against the parent commit.

**Stage A. App fixes that do not depend on the fit.** S9, S1, S2, S8, S3–S5, S7, OAR-01 to 05, D1–D3, M3, M4 and the "eligible" wording. The review's two golden fixtures come first:
- HYPNO against itself with cord 100% and applied limits: toxicity share 0, and replay gives 0.
- D18 against 35×2/47, threshold −3, cord 100%: zero candidates, and clearing cord α/β never yields three unchecked rows.

Stage A touches `R/oar.R`, `R/mod_oar.R`, `R/uncertainty.R`, `R/mod_uncertainty.R`, `R/export.R`, `R/schedule.R` and `app.R`. The export schema version is bumped, and old files replay under their own semantics or are refused with a message.

**Stage B. Data corrections in the builder.**
- KROG weighting fixed.
- RTOG 9003 per decision 4.
- Fraction components stored, and K computed as Σd + Σd²/(α/β) (decision 7).
- Source clock kept beside the normalized T, with a nominal-clock column.
- Endpoint estimator, time point and the exact source percentage carried into the CSV.
- Builder asserts that component counts sum and that rate and N come from matching reports.

Output: `data/refit_trials_2026_10_2X_main.csv` plus sensitivity files. The 10-25 files stay, labeled superseded.

**Stage C. Fit hardening.**
- Run identity: a hash of data bytes, Stan sources, priors, sampler settings, seed and versions. Resume only on a match. The conditional DD run also binds the parent DI file hash and draw IDs.
- Full diagnostics for every parameter and the key contrasts: divergences, R-hat, bulk and tail ESS, E-BFMI, tree depth.
- More outer DI draws for calibrated DD (decision 10), with a replication check.
- Grouped study holdout. LOTO labeled as a conditional score.
- Onset violations reported for fitted and evaluated regimens separately.

**Stage D. Refit.** Main DI and calibrated DD, plus these sensitivities:
- effective N halved;
- intercept prior N(0,10);
- nominal clock;
- KM-only and CIF-only arms;
- the alternative RTOG 9003 source;
- the existing ARTSCAN, new-trial and common-chemotherapy-policy sets.

Each fit must pass the decision-9 gate, or it does not ship.

**Stage E. App integration.**
- `export_draws.py` writes `tumor_di_refit.csv` and `tumor_dd_calibrated.csv`, with a `di_draw` key and chain/draw IDs. The manifest carries data hash, run identity and diagnostics.
- `load_draws` refuses files whose manifest gate fails.
- The free 35-arm draws are kept in the repository as historical and not offered in the app.
- The "paper restrictions" variant is replaced by a fixed-α/β-10 refit variant; that label can say Gy₁₀.
- Point values, the report card column and badges are rebuilt from the new main fit.
- The Model and sources page and guide text follow the review's wording (C07, C09, sec. 5, M2, M8, M9).
- APP_VERSION 3.1.0; equations stay 3.0; uncertainty and selection contract versions are bumped.

**Stage F. Release.**
- Full test suite, the review's acceptance checklist (sec. 8) and the local browser walkthrough.
- Bundle check, then deploy to the existing `custom_bed_calculator_v3` app only. The v2 app is never touched.
- Rollback record pointing at f951cd3 and its draws.
- Live re-check and a fresh-download replay.

**Rough cost.**

| Stage | Estimate |
|---|---|
| A | The largest: about 25 tests |
| B, C | Moderate |
| D | About 1 hour for DI fits; about 1.5–2 hours per set for calibrated DD at 1000 outer draws (200 took 19 minutes) |
| E, F | Moderate |

A and B/C are independent, so they can run in parallel.

**Knock-on.** The Pardee preliminary data quote the 10-25 HYPNO numbers (DI −2.0, DD −3.3). Stage B can move them, so the Pardee text must be re-checked after stage D.

## 5. Decisions for Igor

1. **What the refit replaces.** (a) Only the 35-arm fit: uncertainty draws, fitted preset, report card column, with D18 staying the default. Or (b) also make refit medians the default deterministic preset. **Recommend (a).** The review says to keep D18. The refit's α/β is prior-driven, and the refit extrapolates beyond 2.55 Gy per fraction. (b) can be revisited once the report card shows both side by side.
2. **DD in the uncertainty tab.** Calibrated DD paired with DI, or free DD from the new fit. **Recommend calibrated DD.** The free DD is not identified and has divergences, and the pairing fixes M1.
3. **Sequence.** (a) Full A–F before deploying, or (b) ship the current 10-25 fits now as "provisional" and correct later. **Recommend (a).** Stage B will move the numbers, and (b) would put two different refits in front of users.
4. **RTOG 9003.** Use the internally consistent Fu 2000 pair (N 268/263/274/268, rates 46.0/54.4/47.5/54.5%) in the main fit, with Beitler 2014 updated rates as a sensitivity. **Recommend yes.** Both PDFs are in hand, and the alternative is guessing the 1073 vs 1076 allocation.
5. **Likelihood.** Keep binomial, labeled as an approximation, and add an effective-N-halved sensitivity. **Recommend yes.** An estimator-aware likelihood needs SEs that most papers do not report.
6. **Clock.** Delivered median/mean T stays main, with a nominal-calendar sensitivity. **Recommend yes.**
7. **Mixed fraction sizes.** Use component sums Σd + Σd²/(α/β) in the fit, which needs a small Stan input change. **Recommend yes.** It affects Ghoshal, RTOG 9003, RTOG 88-09 and Vienna.
8. **Uncertainty units.** Relabel as "BED (Gy, α/β varies by draw)" and keep the +1 threshold at each draw's α/β. Or convert to ΔEQD2 per draw. **Recommend relabel.** It keeps the threshold's meaning close to the deterministic tab. EQD2 can be added later as a display option.
9. **Release gate for shipped draws.** 0 post-warmup divergences, R-hat ≤ 1.01, and bulk and tail ESS ≥ 400 for all sampled parameters and the HYPNO contrast. Calibrated DD runs that miss the gate are rerun with longer chains. **Recommend yes.**
10. **Outer DI draws for calibrated DD.** 1000 instead of 200, with 4×100 inner draws each. **Recommend 1000.** The review's outer-cluster SE for HYPNO sign agreement is 1.9 points at 200.
11. **Deployment target.** The same v3 app (`custom_bed_calculator_v3`) as version 3.1.0, with a rollback record to f951cd3. **Recommend yes.** v2 is untouched.
12. **Organ mean-dose limits (OAR-01).** Compare source physical limits on physical dose, with the source fractionation stated. Show EQD2 conversions as illustrative only, never as pass/fail. **Recommend yes.**
