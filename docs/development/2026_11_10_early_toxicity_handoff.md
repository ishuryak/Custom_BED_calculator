# Handoff: early (acute mucosal) toxicity model with uncertainty

Written 2026-11-10 at the end of the session. Start here, then read the plan.

## 1. Where things stand

- **The plan** is `docs/development/2026_11_10_early_toxicity_plan.md` (version 2). It is the
  current specification. Version 1 is kept as `2026_11_10_early_toxicity_plan_pre_2026_11_10.md`,
  and the diff between them is `2026_11_10_early_toxicity_plan_diff_2026_11_10.patch`.
- **Code.** Repository `github.com/ishuryak/Custom_BED_calculator`. The local branch
  `v45_review_2026_11_09` tracks `origin/v4`. Commit `9c35f08` added the plan and its inputs. The
  commit that adds this file adds records only.
- **The app is unchanged.** The live version is 4.5.0 (revision `faf6817`) at
  https://ishuryak.shinyapps.io/custom_bed_calculator_v4/. Nothing was deployed.
- **No refit has been run under plan version 2.** Igor Shuryak's instruction (2026-11-10): "fix the
  plan, then commit and push to github, do not do the refit yet".
- **Git identity.** This WSL environment has no `user.name` or `user.email`. Commit `9c35f08`
  passed the identity used by earlier commits (`ishuryak`) with `git -c` for that commit only.
  No git configuration was changed.

## 2. The question and the short answer

Igor Shuryak asked how good the evidence for the early-toxicity BED is, since it appears to limit
the shorter courses, and whether it can be improved with realistic uncertainty.

- **It does limit them, conditionally.** Under the default settings, the early limit sets the
  highest allowed dose per fraction for 5–29 fractions. It is also why no course of 12 or fewer
  fractions is listed. This is a model result, and the 12/13 boundary differs by only one 0.05-Gy
  grid step.
- **The evidence is moderate for conventional and accelerated schedules and weak for the short
  courses.** The Fowler 2003 values are a reasonable reference. However, the onset of mucosal
  regrowth is weakly identified, the "+5 Gy₁₀" margin is a policy choice with no patient-risk
  meaning, and the historical data have few schedules above 2.6 Gy per fraction.
- **It can be improved** by fitting patient event counts from randomized trial arms, keeping the
  historical tolerability labels in a separate model, and showing uncertainty in separate layers.
  Plan sections 4–8 give the details.

## 3. Documents, in the order they were written

All are in `docs/development/`.

| File | What it is | Status |
|---|---|---|
| `2026_11_10_early_toxicity_plan_pre_2026_11_10.md` | Plan version 1 | Superseded |
| `2026_11_10_early_toxicity_plan_audit.md` | Audit of version 1 by two subagent reviewers (findings F1–F12) | F4, F5, F6 and F9 revised or replaced by the review. A status note at the top says so |
| `2026_11_10_early_toxicity_independent_review.md` | Independent review by Codex (corrections R1–R10, steps 1–8, tests 1–18) | Adopted in full in version 2 |
| `2026_11_10_early_toxicity_plan.md` | Plan version 2 | **Current** |
| `2026_11_10_early_toxicity/exploratory/` | Version 1 scripts (`binding.R`, `parse.py`, `fit.py`, `implic.py`) | Kept for reproducibility only |
| `2026_11_10_early_toxicity/exploratory_after_audit/` | Scripts corrected after the audit | Kept for reproducibility only. Still include the 8 partial courses |
| `2026_11_10_early_toxicity/verification/verify_review_claims.py` | Recomputes the review's numeric claims | All matched on 2026-11-10 |

## 4. What changed between version 1 and version 2, and why

1. **Eight partial-course records are excluded.** They are spreadsheet rows 541–576 of the
   Pardo-Montero 2021 supporting file. Pardo-Montero excludes them as the first parts of courses
   with a planned gap. The selections are 103 records (primary, all disputed labels dropped) and
   107 (sensitivity, disputed labels counted by their leaning).
2. **Schedule identity uses the full fraction-dose and time vectors.** The audit's
   (fractions, total dose, duration) key wrongly merged four pairs with different timing.
   Evidence identity (report, cohort, trial, arm, endpoint) is kept separate from schedule
   identity.
3. **Randomized contrasts do not identify the onset day Tk** when both arms end after onset,
   because Tk cancels from the difference. Version 1 claimed they did. Identification is now
   checked before fitting, and Tk is fixed within named scenarios if it is not identified.
4. **HYPNO is an external benchmark, not a pass/fail gate.** The audit's gate compared against
   70/35, which HYPNO did not randomize, and used an overlap rule that very uncertain predictions
   pass. The audit's claim that the label-only refit "contradicts" HYPNO is withdrawn: that
   depends on a patient-risk slope the label model does not estimate.
5. **The 8-day onset prior is a labeled sensitivity scenario,** not evidence. Fowler's 7 days stays
   in the reference.
6. **Chemotherapy is part of each randomized strategy.** Planned therapy defines the arm. Delivered
   dose is recorded but not adjusted for automatically.
7. **One dose evaluator** (with incomplete repair) for fitting and prediction, plus a separate
   legacy evaluator that reproduces 4.5.0.
8. **Uncertainty shares are defined precisely.** Candidate and reference are evaluated in the same
   draw. No numeric placeholder is used for "no qualifying dose".
9. **The label-only refit posteriors are provisional and no longer quoted**, because they include
   the partial courses and their draws were not archived.

## 5. Numbers verified this session

Recomputed independently with `verify_review_claims.py` (argument: the `.xls` below):

- Selection counts 116 / 115 / 111 / 107 / 103 records, with 15 / 15 / 13 / 15 / 13 intolerable
  labels.
- 92 coarse against 96 full-vector schedule groups, and 13 event-bearing groups under either key.
- Fowler no-repair early BEDs: 53.10544 (70/35), 55.86597 (55/20) and 55.43495 Gy₁₀ (66/33 at six
  per week). The HYPNO contrast is +0.43 Gy₁₀.
- Repair fixture 84.12426 Gy₁₀ for 70/35 (T½ 3.2 h), and stress fixture +0.28402 Gy₁₀ for 5 × 8 Gy
  at 24 h.
- HYPNO grade ≥ 3 mucositis 190/373 against 208/380, a difference of −3.80 points (Wald 95% CI
  −10.92 to +3.33).

Data file: `Custom_BED_calculator_v2/References/PardoMontero_2021_SupportingInfo_108_schedules.xls`
(CC BY-NC, not committed). Reference PDFs are in `Custom_BED_calculator_v2/References/`, and
Bentzen 2001 is in `References_refit_trials/eric_2026_10_22/Part 3/`. All are git-ignored.

## 6. Open decisions

Plan section 9 gives the full text. Ask them in one numbered round, with the recommendations,
before starting step 1.

| ID | Topic | Needed before | Recommendation |
|---|---|---|---|
| D1 | Primary endpoint group | Step 3 | Maximum grade ≥ 3 acute mucositis, scale version recorded and mapped |
| D2 | Primary clinical context | Step 1 | Broad-field definitive chemoradiotherapy for locally advanced head and neck squamous cell carcinoma (HNSCC) |
| D3 | Prior scales for the clinical benchmark | Step 5 | Adopt the review's proposed baseline and slope priors |
| D4 | "+5" margin when α/β varies | Step 7 | Numerical BED increment at the sampled α/β, with the EQD2 version as a sensitivity |
| D5 | Clinical tolerance δ for a risk rule | Step 7 | No default. Without one, show the risk difference without pass/fail labels |
| D6 | Denham 1995 full text | Step 2 | Obtain it. Until then, no prior uses its 9.3 (5.8–17.9) interval |

Settled by version 2 unless Igor Shuryak objects:
- Fowler stays the default until step 6 shows a conclusive advantage, and the live app is
  unchanged until then.
- Bentzen 2001 is a consistency check only.
- HYPNO is a benchmark.
- The two likelihoods are fitted separately.

**Update 2026-11-10 (later):** D1–D5 were decided as recommended, D6 was replaced (Denham 1995 is
not obtainable, so its interval is cited via Kaanders 1999), and D7 (follow the gates, no pilot)
was added. Step 1 is done: `2026_11_10_early_toxicity_step1_target_definition.md`. Four PDFs were
added to the local git-ignored reference folder: Sanghera 2007, Kaanders 1999, Lester 2024 and
Jegannathen 2011. All four are ScienceDirect web pages printed to PDF, so wide tables are clipped
and about one line is hidden at each page break. Use the publisher PDF for any table that becomes
data.

## 7. Next session

1. ~~Ask D1–D6 (one round).~~ Done. **Start at step 2**, using section 7 of the step 1 record as
   the starting index.
2. **Step 1** (about 2 hours): write the structured target definition. No code change.
3. **Step 2** (3–5 days, Igor Shuryak spot-checks): build `data/early_tox/` (reports, schedules,
   fractions, labels, arms, outcomes, evidence families). Every value comes from a downloaded full
   text, and missing PDFs are reported by name with the missing-reference table. Gate: the
   selection counts in section 5 reproduce.
4. **Steps 3–4:** endpoint mapping, then identification and prior predictive checks.
5. **Step 5 is the first refit.** Do not start it before steps 1–4 pass and Igor Shuryak agrees.

Rules that apply throughout:
- Each code change needs one committed test, seen failing without the change and passing with it,
  plus a `docs/FIXES.md` row.
- Long fits follow the analysis-script convention: checkpoints, seeds, thread limits and logs.
- No deployment without explicit approval.
