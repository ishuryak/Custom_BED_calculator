# Early-toxicity steps 1–6: response to the independent review (Phases A–D)

Written 2026-11-11. The review (`../HNSCC_BED_early_toxicity_independent_review_and_action_plan_2026-11-11.md`,
outside the repository) found nine issues, F1–F9, in commit `a42cf3f`. Claude checked every
finding before acting: `docs/development/PLAN_independent_review_2026_11_11.md`, section 1.

**Decisions, Igor Shuryak, 2026-11-11** (plan approved: Phases A–D):
- D2: decide on the app research option after A–D.
- D3: adequate ascertainment means assessment reaching at least 4 weeks after RT with no post-RT gap
  over 3 weeks.
- D5: Tandon stays 0 with a −1 sensitivity.
- D1, D4, D6: as the review recommended.

**What changed:**
- The app is unchanged, and Fowler 2003 stays the default.
- No count in the data changed.
- The scientific results of steps 5 and 6 reproduce exactly. They are now backed by complete
  diagnostics, run identities and archived raw chains.
- Several statements in the earlier records were wrong and are corrected in section 4.

## 1. What was done, per finding

| Finding | Change | Test (each seen red when its line is mutated) |
|---|---|---|
| F7 extraction identity | `data/early_tox/source/extraction_map.csv` names the exact row of each extraction that confirms each outcome. The builder checks report, numbers, endpoint class (never a grade 2 row), unique outcome IDs and one confirming use per row. The loader refuses an arm entered twice. Outcomes gain an `extraction_rows` column, all other outcome fields are byte-identical | `test_early_tox_arms.py::test_a_count_moved_to_another_arms_value_is_refused` (O055 8→15 refused), `::test_one_extraction_row_cannot_confirm_two_outcomes`, `::test_grade_2_and_wrong_report_rows_are_refused`, `::test_duplicate_outcome_ids_and_map_rows_are_refused`, `test_early_tox_step5.py::test_an_arm_entering_the_analysis_twice_is_refused` |
| F6 ascertainment, population, totals | `source/ascertainment.csv`: last post-RT assessment, largest post-RT gap, endpoint time point and analysis population per family. Site metadata corrected (Gupta M, Grover, DAHANCA, EORTC 22851). Totals are derived from the loaded rows. U10 is relabeled the nominal-window sensitivity. New units U15 (candidate) and U16 (benchmark) use the D3 rule | `test_early_tox_step3_step4.py::test_ascertainment_rule_fixtures` (DAHANCA's 8-week gap fails), `::test_adequate_set_and_recorded_time_points`, `::test_a_family_without_an_ascertainment_row_is_refused`, `test_early_tox_step5.py::test_reported_totals_come_from_the_loaded_rows` |
| F1 planned chemotherapy | `source/planned_therapy.csv`: agent, courses, cumulative planned dose and status (known or assumed_equal) per arm. A completeness test refuses a known planned-dose difference without a code or declared scenarios. HYPNO (4 vs 5 weekly courses, about 75% intended) is scored under three declared scenarios. Tandon, Sikdar, Shivananjappa, Mani and Tallari are "assumed equal", because their papers state no planned difference | `test_early_tox_step5.py::test_planned_therapy_differences_are_coded_or_declared_as_scenarios` (red with HYPNO's scenarios removed: "planned 140.0 vs 175.0 but no code or scenario"), `::test_hypno_scenarios_and_tandon_sensitivity_are_declared`, `test_early_tox_step6.py::test_hypno_is_scored_under_three_declared_chemotherapy_scenarios`, `::test_probability_mixture_matches_a_pure_code_and_rejects_bad_weights` |
| F4 resume and eligibility | One run identity over data, settings, repair half-time, Stan source, sampler and versions. Units are written to `*.partial-<pid>` and published by atomic rename after their manifest. Reuse is refused on mismatch, legacy output, an unreadable manifest or changed outputs. Step 6 has `--rescore-only`. `promote` is replaced by `score_criterion_met` and `numerical_rule_pass`, subsets and failed or legacy gates are never eligible. `decision` stays with Igor | `test_early_tox_step5.py::test_unit_state_*`, `::test_run_identity_changes_with_each_defining_input`, `::test_overdispersion_switch_reaches_the_fit_and_finished_units_are_reused` (sampling marker), `test_early_tox_step6.py::test_cache_without_identity_or_with_another_identity_is_refused`, `::test_subsets_failed_gates_and_legacy_gates_are_never_eligible` |
| F5 gates and SBC | Exact-schema gate: every sampled quantity, including all `zarm`, must have finite R-hat and ESS. Fixed benchmark α/β and r are checked as constants. E-BFMI ≥ 0.3 in every chain, treedepth saturation reported. Raw CmdStan chains archived in `~/early_tox_chain_archive` (114 fits, 1.2 GB), with their hashes in each manifest. SBC uses the exact discrete bin probabilities (11/101, then 10/101) and lists its divergent replicates | `test_early_tox_step5.py::test_gate_requires_every_quantity_finite_and_checks_constants`, `::test_sbc_bins_cover_101_discrete_ranks_and_divergent_reps_are_listed` |
| F2 conditional intervals | `contrast_check()` weights draws, baselines and arm effects by the observed reference-arm likelihood. It simulates the other arms' counts and compares the same continuity-corrected log OR on replicates and data. The latent log OR, predicted counts and weight ESS are reported separately | `test_early_tox_step6.py::test_conditional_prediction_responds_to_the_reference_count`, `::test_uninformative_reference_gives_the_unconditional_prediction`, `::test_conditional_check_handles_zero_and_all_events_and_shares_the_reference`, `::test_conditional_check_runs_on_fit_shaped_draws_with_a_scale_group_offset` |
| F8 repair | The repair half-time is a setting, `halftime_h`, part of the run identity. Units U17–U22 cover 2 h, 4 h and complete repair. Step 6 `--structural` refits both models on all primary families at each half-time and scores the external checks. The historical comparison is now factorial in repair | `test_early_tox_step6.py::test_repair_half_time_reaches_the_external_contrasts`, `test_early_tox_step5_historical.py::test_matched_repair_variants_differ_only_in_repair` |
| F3 applicability | `tools/early_tox/applicability_table.py` writes `2026_11_10_early_toxicity/applicability_table.csv` (66 arms). The schedule class follows a stated rule. The table carries therapy status, technique, scale, assessment schedule, role and each family's score | `test_early_tox_applicability.py::test_schedule_classes_of_known_arms` (red for both defects of the first draft rule) |
| F9 app integration | Not implemented. It is Phase E, held for decision D2. The fixture is recorded: 526 of 4,000 expanded U02 draws give E ≤ 0 for 70/35, so the legacy positivity filter must not be reused | none yet |

**Two defects found by the combined smoke run, not by the unit tests.** Both are fixed in `04a3082`
with tests:
- The conditional check read fitted draws without their scale-group offset.
- The fixed-constant check demanded agreement to 1e-9, against CmdStan's 6-figure summaries.

## 2. Runs (revision `04a3082`, step 6 records `378d58f`, which differs only in the applicability script)

`step5_review/` holds 42 overdispersed units (U01–U22 and 22 leave-one-family-out units) and 20
binomial units, 62 fits in all.
- **Overdispersed:** all 42 pass the complete gate.
- **Binomial:** 3 of 20 miss R-hat < 1.01 (U01 1.014, U05 1.011, U07 1.012). The binomial set is
  not primary: its predictive check fails, p ≤ 0.016.

`step6_review/`, run with `--structural`: 44 held-out fits, 2 all-primary fits and 6 structural
fits.
- All pass the complete gate.
- The rescore-only run from the same caches completed without sampling.
- **Restart:** the run was interrupted once and restarted. Units U01–U03 were recognized as complete
  by identity and output hashes, and were not refit.

## 3. Results

**Held-out scores reproduce exactly.** The seed and data are unchanged, so the draws are identical.

| Task | Mean difference per family (95% Bayesian bootstrap, both seeds) | Families better | numerical_rule_pass |
|---|---|---:|---|
| New trial | +0.284 (+0.109 to +0.55) | 18/22 | true |
| Conditional contrast | +0.279 (+0.103 to +0.55) | 18/22 | true |

`numerical_rule_pass` is now true *with* complete gates, full scope and archived chains. It is still
not a promotion decision.

**Corrected conditional calibration (F2).** 23 of 29 held-out contrasts fall inside the 90%
observed-scale predictive interval for the benchmark, and 25 of 29 for the candidate. The median
widths are 2.08 and 1.89. These match the review.

**HYPNO under the chemotherapy scenarios (F1).** All-primary fits, 300 inner draws, the conditional
score difference is candidate minus Fowler-parameter benchmark:

| Scenario for the 55/20 arm | Difference |
|---|---:|
| Code 0 (as before) | −0.686 |
| Code −1 (one course fewer, everyone) | +0.595 |
| 75% / 25% probability mixture | +0.376 |

The ranking depends on the chemotherapy assumption. HYPNO does **not** show that the candidate
over-predicts hypofractionation toxicity. These are sensitivity scenarios, not a corrected estimate,
and arm × chemotherapy counts were not found in the paper.

**Other external checks** (difference): ELAN-RT −0.42 (frail older patients, 2-week split), Dwivedi
+0.13, IAEA-ACC +0.22, Poddar +0.13, Pinto +0.07, Sanchiz +1.05.

**Tandon sensitivity (D5)**, held-out fold, conditional score:
- Benchmark: −3.02 at code 0, −4.02 at −1.
- Candidate: −2.53 and −3.25.
- The candidate is ahead in both: +0.49 and +0.77.

**Repair half-time (F8).** Both models refitted on all primary families:

| Half-time | Candidate r | α/β | OR 20 × 2.75 vs 70/35, candidate | OR, benchmark | HYPNO code 0 / −1, difference | ELAN, difference |
|---|---|---|---|---|---|---|
| 2 h | 1.52 | 13.1 | 2.47 | 1.27 | −0.69 / +0.58 | −0.38 |
| 3.2 h (primary) | 1.55 | 13.0 | 2.40 | 1.24 | −0.68 / +0.60 | −0.42 |
| 4 h | 1.57 | 13.0 | 2.35 | 1.22 | −0.66 / +0.62 | −0.44 |
| Complete | 1.51 | 13.5 | 2.57 | 1.30 | −0.68 / +0.53 | −0.34 |

The repair assumption moves nothing material.

**Adequate peak ascertainment (D3).** Only CAIR and ORO 93-01 qualify: 5 arms, 255 patients.
- U15 gives r = 1.11 (0.51–2.27) and an odds ratio for 20 × 2.75 Gy of 3.05 (0.72–33).
- U16 gives 1.44 (1.08–1.98).
- This set cannot settle the question.

Most trials either scored only during RT or first saw patients 4 or more weeks after RT. That is a
limitation of the evidence as a whole, now stated rather than hidden. Nominal-window U10 (13
families) gives r 1.55 (1.00–2.62), the same as U02.

**Historical labels, matched repair (F8).** Leave-one-source-out log scores:
- **3.2 h repair:** Fowler −27.54, candidate −30.29.
- **Complete repair:** Fowler −25.90, candidate −29.28.

Fowler stays ahead on these labels. This is a descriptive compatibility check: Fowler 2003 was
tuned on largely the same schedules.

**SBC, corrected bins (F5).** No detected rank departure, the smallest p is 0.148 (expanded run)
and 0.067 (first run). The divergent replicates are now listed: 16 of 500 with 41 transitions, and
17 with 29.

**Applicability (F3).** Experimental hypofractionated arms in the primary fit, per-family difference:

| Trial | Difference |
|---|---:|
| Tandon (2.4 Gy) | +0.49 |
| Shivananjappa (2.2 Gy SIB) | +0.03 |
| Sikdar (2.2 Gy SIB) | +0.02 |
| Grover (2.2 Gy SIB) | −0.17 |
| PMH (reference arm 51 Gy/20 at 2.55 Gy) | −0.32 |

The evidence within hypofractionation is mixed, not one-sided. The largest gains come from Vienna,
RTOG 9003 and GORTEC 99-02. That alone does not show the time term is specific to acceleration.

## 4. Statements corrected

| Earlier statement | Correct statement |
|---|---|
| "The data now hold … 6,256 patients" (step 6 record), "4,500 patients" (step 5 and earlier messages) | 51 arms, 22 families, **6,042** evaluable, 2,945 events. The 15-family set had **4,599** evaluable |
| "HYPNO and ELAN-RT … both favor Fowler … the candidate over-predicts" | HYPNO's ranking reverses under plausible chemotherapy scenarios. ELAN-RT favors Fowler but is a frail, split-course population |
| "The candidate's advantage comes from accelerated trials (≤ 2 Gy)" | The primary fit already includes five hypofractionated comparisons with mixed results. The largest gains are in accelerated or multiple-fraction trials, which does not establish that the effect is specific to acceleration |
| "23 of 29 within the 90% predictive interval" (both models) | That interval was not conditional on the reference arm. Corrected: 23/29 benchmark, 25/29 candidate |
| SBC "calibrated, smallest p 0.28" | Exact bins: no detected rank departure, smallest p 0.148, 16/500 replicates had divergences |
| `promote: true` | `score_criterion_met` and `numerical_rule_pass` are numerical criteria. The decision is Igor's |
| The step 6 benchmark is "Fowler" | It is a **Fowler-parameter statistical benchmark**: Fowler's α/β and r with fitted slope, chemotherapy term, baselines and arm variation, 3.2 h repair and an exact hinge. The app's Fowler limit uses complete repair, a smooth hinge and a fixed BED threshold |

## 5. Recommendation and what remains

**Keep Fowler as the app default.** The trial-count model improves the average held-out log score
across the 22 recorded families (about +0.28 per family, positive lower bound, 18 of 22). That result
now has complete diagnostics. However:
- The evidence within hypofractionation is mixed.
- HYPNO's ranking depends on its chemotherapy assumption.
- ELAN-RT favors the Fowler-parameter benchmark in a different population.
- Few trials observed the post-treatment mucositis peak.
- The historical labels favor Fowler.

A research scenario remains possible (decision D2, after this report). It would need Phase E as the
review specifies:
- its own evaluator and manifest
- no positivity filter
- relative outputs only, without patient-risk numbers

**Open, not blocking:**
- ELAN-RT second extraction, which needs a clean publisher PDF.
- Arm × chemotherapy counts for HYPNO, if obtainable.
- The historical-label crosswalk, only if independent validation across both evidence streams is
  ever claimed.

## 6. Verification

- `python3 -m pytest -q tests/python` passes. Tests that compile or sample Stan are marked
  `sampling`.
- `stale_constant_audit.py --project . --project-root . --strict`: 0 errors.
- Combined smoke on `04a3082`:
  - step 5 smoke units U02, U15, U16, U19 and U22
  - step 6 subset with `--structural`: gates pass, labeled ineligible
  - rescore-only from the same cache: no sampling
- Full runs as in section 2.

Every number in section 3 is read from `step5_review/`, `step6_review/`,
`step5_review/historical/summary.json`, the two `sbc/summary.json` files and `applicability_table.csv`.
