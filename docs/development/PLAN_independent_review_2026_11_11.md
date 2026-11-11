# Plan: applying the independent review of early-toxicity steps 1–6 (2026-11-11)

Review: `../HNSCC_BED_early_toxicity_independent_review_and_action_plan_2026-11-11.md` (outside the
repository), reviewed commit `a42cf3f`. Claude checked every finding before writing this plan, on
2026-11-11.

**Approved by Igor Shuryak, 2026-11-11:** Phases A–D. D2: decide after A–D. D3: at least 4 weeks after RT with no post-RT gap over 3 weeks. D5: Tandon stays 0 ("assumed equal") with a −1 sensitivity. D1, D4 and D6 as recommended.

## 1. Verification of the findings

| Finding | Claim | Checked by | Verdict |
|---|---|---|---|
| F1 | HYPNO planned 4 vs 5 weekly cisplatin courses; the model coded 0 | Bentzen 2026 PDF p. 5 ("4 courses in the HFX arm and 5 courses in the NFX arm"), p. 8 (prescribed in 75%) | **Confirmed.** Rescoring with the HFX arm coded −1 reverses the HYPNO result: −0.686 → +0.613 (reviewer: +0.586) |
| F2 | `contrast_check()` does not condition on the reference count; `inner` unused; latent interval against an observed statistic | Code reading (`step6_validate.py`, `contrast_check`) | **Confirmed.** The 23/29 coverage and the external "pred" intervals are wrong as labeled. The integrated log scores are unaffected |
| F3 | The primary set already includes hypofractionation (PMH 51/20, Tandon 60/25, Grover, Shivananjappa and Sikdar 66/30). Dwivedi is not accelerated | `data/early_tox/arms.csv`, `step6/per_family.csv` | **Confirmed.** Per family: PMH −0.32, Grover −0.17, Tandon +0.49, Shivananjappa +0.03, Sikdar +0.02. My statement that "the advantage comes from ≤2 Gy acceleration and the hypofractionation tests favor Fowler" was too strong |
| F4 | Step 5 skips a unit when `diagnostics.json` exists; step 6 caches have no identity; `promote` ignores the gates and the subset scope | Code reading | **Confirmed** (`promote=bool(...)` never reads `gates_all_pass`) |
| F5 | Gates skip NaN and unlisted parameters (zarm); SBC bins wrong; SBC divergences not reported | Code reading; `step5_expanded/sbc/ranks.csv` | **Confirmed.** 16/500 replicates with 41 divergences (expanded), 17/500 with 29 (first run). The first bin covers 11 ranks |
| F6 | Totals; EORTC 22851 window; DAHANCA assessment gap; site metadata | Horiot 1997 Table 4 (during treatment 123/245 and 160/240, 6 weeks after 84/245 and 168/240); Mortensen p. 70; `outcomes.csv` join | **Confirmed.** 51 arms, **6,042** evaluable, 2,945 events. The earlier set was **4,599**. I reported 6,256 and 4,500: both wrong |
| F7 | The builder confirms counts by report and numbers only | Negative control: O055 changed to 15/30 | **Confirmed**: accepted as "both extractions" |
| F8 | Repair half-time fixed at 3.2 h, no 2/4 h sensitivity; the historical comparison mixes repair assumptions | Code reading | **Confirmed** |
| F9 | Legacy positivity filters would break the fitted model's same-regimen support | 4,000 expanded U02 draws | **Confirmed**: 526 draws give E ≤ 0 for 70/35 |

## 2. Planned changes (one item per finding; each change gets a test seen red first)

**Phase A: data and builder. No sampling.**
- **A1 (F7).** The builder matches each extraction by report, arm, endpoint class and threshold, as
  well as the numbers.
  - Store the matched extraction-row indices, and use an explicit arm-label map where labels differ.
  - Require unique outcome and arm IDs and one primary outcome per arm.
  - Negative controls that must fail: O055 changed to 15/30, a grade 2 row, a pharyngeal row, a
    duplicate ID, a swapped arm.
  - Positive control: an unchanged rebuild reproduces all six CSVs byte for byte.
- **A2 (F6).** Record source-population metadata without changing any count:
  - endpoint time point (EORTC 22851 "during treatment")
  - assessment schedule and the largest post-RT gap (DAHANCA: last RT week, then 8–9 weeks)
  - as-treated, completer and evaluable status
  - site corrections: Gupta M (no oral cavity, 14 glottic), Grover (no oral cavity), DAHANCA (oral
    cavity included), EORTC 22851 (nasopharynx included)
  - Replace W1/W2 with one declared ascertainment rule (see decision D3). This changes only U10's
    membership.
- **A3 (F1).** Add structured planned-therapy fields: drug, dose, number of courses, the proportion
  intended to receive chemotherapy, and status (known, assumed or unknown).
  - Audit every weekly-cisplatin arm. Tandon 2018 does not state a planned number of courses; it
    reports at least 5 cycles (200 mg/m²) in 83% vs 73% (p. 4). It stays 0 with status "assumed
    equal". Grover caps both arms at 6 courses and stays 0.
  - HYPNO becomes an external scenario set, not a fixed 0.
  - **If the audit changes a primary arm's code, the data version changes and the affected fits are
    redone (Phase D).**
- **A4 (F6).** Every reported total is derived from the loaded analysis set (6,042, 2,945, 51, 22).
  Correct the step 5 and step 6 records.

**Phase B: run identity and gates. No sampling.**
- **B1 (F4).** Add a run-identity hash over:
  - the normalized fitting data, families, chemotherapy codes and fraction times
  - the Stan source, priors, repair half-time, Tk and sampler settings
  - the package versions

  A unit is reused only when its identity, required outputs and gates all agree. The completion
  manifest is written last, by atomic rename.
- **B2 (F4).** A `--rescore-only` mode for step 6 that refuses to sample. Replace `promote` with
  `numerical_rule_pass`, which is false for subsets or any failed gate, plus a separate
  `decision: research_only | default` set by Igor.
- **B3 (F5).** Exact-schema gates: every sampled parameter, including zarm and the decision
  quantities, must be finite. Fixed benchmark α/β and r are checked as constants. Report energy and
  treedepth.
- **B4 (F5).** SBC summary with discrete-rank bin probabilities 11/101 and 10/101, plus the list of
  divergent replicates. "No detected rank departure" replaces "calibrated".
- **Retained artifacts.** The current step 6 NPZ caches (4,000 rows, passing narrow gates) are kept
  as an identified research generation. Raw chains were not archived, so diagnostics beyond the
  stored gate values cannot be recovered without refitting. This is stated rather than inferred.

**Phase C: corrections from saved draws. No new fits.**
- **C1 (F2).** A true conditional predictive check:
  - weight posterior and baseline draws by the reference-arm likelihood
  - simulate the counts of the other arms
  - compute the same continuity-corrected log OR on replicates as on the observed data
  - label the latent interval separately

  Tests: changing the reference count changes the prediction; a no-information reference recovers
  the unconditional prediction; 0 and all events; multiarm trials. Regenerate `contrasts.csv`,
  `external.csv` and the coverage figures (reviewer: 23/29 for Fowler, 25/29 for the candidate).
- **C2 (F1).** HYPNO scenarios from the saved all-primary draws:
  - chemotherapy code 0
  - code −1
  - a 75/25 probability mixture formed before the binomial likelihood

  Reported as sensitivity, not as a corrected estimate. Seek arm × chemotherapy counts (not in the
  paper as read).
- **C3 (F8).** Historical classifier rescored with matched repair: Fowler and candidate both at
  3.2 h, and both with complete repair. The reviewer found the ordering holds: −27.54 vs −30.29 and
  −25.90 vs −29.28. Label it a descriptive compatibility check.
- **C4 (F3).** Applicability table generated from the data. Per family it gives:
  - dose per fraction, fractions per week and per day, gaps
  - chemotherapy strategy, technique, scale, assessment schedule
  - role (primary held-out, external, contextual)
  - per-family score

  Rewrite the step 6 interpretation with the reviewer's replacement text, adapted, and correct my
  earlier statements to Igor.

**Phase D: necessary fits only (about 1 CPU hour in total, one core).**
- **D1.** Refit U10 under the corrected ascertainment rule, together with a matched benchmark.
- **D2 (F8).** Repair sensitivity at h = 2 h, 4 h and complete repair for U01 and U02 on all
  primary data, with the HYPNO and ELAN external checks. Held-out refits only if a variant's held-out
  claim will be reported.
- **D3.** If A3 or decision D5 changes a primary chemotherapy code, refit U01, U02 and all 44
  held-out fits under the new data version, and rerun step 6. Under the recommended D5, none.
- **D4.** An ELAN-RT second extraction, only if a clean publisher PDF arrives.

**Phase E: app integration (step 7). Only after decision D2 below.** Research evaluator, manifest
and domain rule as in review F9 points 1–7. The legacy Fowler path and the 4.5.0 list (69 doses in
18 rows) are unchanged and tested. This needs its own approval before deployment.

## 3. Decisions for Igor (recommended answer first)

- **D1. Default.** Keep Fowler. *Recommended* (the reviewer agrees).
- **D2. Research scenario in the app.** Add it as an explicitly selected scenario after Phases A–D
  pass. *Recommended* if you want it. The alternative is Fowler only, which removes Phase E.
- **D3. Estimand.** Keep the approved mixed-window primary analysis (M2) with honest labels. Refit
  the window-restricted sensitivity (U10) under a declared ascertainment rule: assessment reaching
  at least 4 weeks after RT **with** no post-RT gap longer than 3 weeks. Under this rule DAHANCA and
  EORTC 22851 leave U10. *Recommended.*
- **D4. HYPNO.** Report chemotherapy scenarios, not one corrected estimate. *Recommended.*
- **D5. Chemotherapy coding where the planned courses are not stated (Tandon).**
  - *Recommended:* keep 0, marked "assumed equal", because the paper states no planned difference
    and the delivered amounts were similar (at least 5 cycles in 83% vs 73%).
  - Report a −1 sensitivity from the saved draws for the held-out Tandon fold.
  - This changes no primary input, so no D3 refits are needed.
  - Alternative: code −1 and refit everything.
- **D6. Historical labels.** Keep them as a compatibility check. No crosswalk unless independent
  validation is claimed. *Recommended.*

## 4. Order and stopping condition

A1 → A2 → A3 → A4 → B1 → B2 → B3 → B4 → C1 → C2 → C3 → C4 → D1 → D2 → (D3 if needed) → records,
FIXES rows and push → stop for decision D2 before any Phase E work.

Each step's tests must pass, red when broken, before the next. The full Python suite runs before
every commit.

When Phases A–D pass, stop. Reopen only for a failed criterion or new evidence. No further general
review is planned.

## 5. Status (2026-11-11)

Phases A–D are done. The record is `2026_11_11_early_toxicity_review_response.md`. Phase E waits
for decision D2.

