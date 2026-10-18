# Code review of BED calculator v2, 2026-10-18

## Scope and method

The code reviewed was `app.R`, `R/bed_models.R`, `about.md`, `README.md` and the tests, after the 6 fractions/week option was removed at the user's request (test `only 5 fractions per week (Mon-Fri) is offered`, observed failing before the change and passing after).

Four independent reviewers ran in parallel. None edited project files.

1. **Engineering (R/Shiny):** probed the server with `testServer` and called the functions directly.
2. **Radiobiology mathematics:**
   - reimplemented all BEDs independently in Python, agreeing within 4e-12 on 2 parameter sets × 290 regimens
   - solved the DD calibration independently (brentq), agreeing within 1e-13 over 1,728 parameter sets
   - checked the formulas against the 2018 paper, equations 2–5.
3. **Radiation oncologist:** clinical meaningfulness, and every claim on the "Model and sources" page checked against `References/README.md`.
4. **Tests, documentation, usability:** 20 new mutants, parity between documentation and code, writing rules.

The findings marked "verified" were reproduced by the main session in `verify1.R` and `verify2.R`.

## Confirmed correct

- **Formulas:** all four BED formulas and EQD2 match an independent implementation.
- **DD calibration:** BED_DD equals BED_DI at the reference for any tumor α/β. The v1 bug is fixed, v1 and v2 agree at α/β = 10 (3e-14).
- **Paper match:** the DD form matches equations 2–5 of the 2018 paper. The γT term is dropped consistently from both models.
- **Acute mucosa:** the formula and parameters match the Fowler 2003 abstract and Fowler 2010 Table 4. 35 × 2 Gy gives 52.3 Gy10.
- **Smooth threshold:** the "< 0.24 days" statement is correct (ln2/3 = 0.231 d).
- **Calendar:** the inclusive day count matches the 2018 paper's Table 1, HYHOPE and Pardo-Montero 2021. For acute BED, a one-day convention shift cancels against the reference.
- **Search:** toxicity BEDs rise monotonically in dose per fraction. The full grid search does not miss a regimen when DD is non-monotone, and the 12 default rows match an independent Python search.
- **Speed:** about 0.007 s per search.
- **Citations:** all ten on the sources page match the verified sources. Reference [1], the user's own 2018 paper, is in Publications, it is not in `References/README.md`.

## Defects

| # | Severity | Defect | Evidence |
|---|---|---|---|
| D1 | critical | λ_DD = 0 is allowed by the input limits but gives BED_DD = NaN. The verdict badges crash and suggestions contain an all-NA row. | verified: NaN, 1 NA row |
| D2 | major | When the DD calibration is impossible (onset before treatment, or C infinite), the warning shows only on the Results tab. DD values are still shown and used in badges and suggestions. This happens for ordinary inputs (λ_DD = 0.2 gives onset at day −16.2) and in 27% of a 1,728-point parameter grid. | verified |
| D3 | major | The warning text is wrong at λ_DI = 0 ("day Inf … before treatment") and at λ_DD = 0 ("day −Inf"). | verified |
| D4 | major | The "Limiting toxicity" label compares room left in Gy10 against Gy3 and is sometimes wrong. Example: tol_early = 10, m = 20, d = 2.90 is labeled "early", but the next step breaks the late limit (+0.35) while early stays within its limit (+9.04 ≤ 10). | verified, 7 of 524 cases in one scan |
| D5 | major | The 8 Gy cap on the dose grid silently truncates results. The best dose can sit at the cap with the wrong label. | engineering probe |
| D6 | major | The reference overall time is never derived from the calendar. A 30-fraction reference keeps T = 47 d, and manual T is not checked against the number of fractions (m = 20, T = 1 is accepted). | engineering + QA probes |
| D7 | minor | m is rounded silently. The fraction slider (5–34) is not tied to the reference. The dose grid is anchored at ref_d + step instead of round doses. Validation messages show internal ids. | probes |
| D8 | major (deployment) | `rsconnect::deployApp` ignores `.gitignore`, so it would upload `References/` (paywalled PDFs) and `tests/` to shinyapps.io. Needs `.rscignore`. | engineering |
| D9 | major (tests) | 11 of 20 new mutants survived. Untested: best dose per fraction number, tol_late, the binding label, eqd2_early, a non-default reference inside the search, the calibration `ok` flag, and the λ = 0 boundary. Two tests partly recompute the code's own answer. There are no server-side tests. | QA mutation table |
| D10 | minor (documentation) | The search space is not documented in `about.md`. Also: the smooth threshold can move BED_DD by up to 1 Gy (it flipped DD pass/fail in 2 grid cells, none at defaults), dropping γT costs short regimens about 1.6 Gy, the DI onset also moves DD through the calibration. Plus README wording, undefined acronyms (HNSCC, HPV, Gy10) and "licence". | math + QA |
| D11 | minor (usability) | The plot separates DI and DD by color only. The CSV "Δ" header may garble in Excel, and the CSV is not rounded. The plot title says "gain" when ΔDD is negative. | QA |

## Scientific and clinical issues (need the user's decision)

| # | Issue | Evidence |
|---|---|---|
| S1 | "Early BED ≤ 35 × 2 Gy" rejects the HYPNO regimens that 792 patients tolerated. 55 Gy/20 fx gives 55.07 against 52.31, HYPNO's own control arm, 66 Gy/33 fx at 6 fx/week, gives 54.64. Every default suggestion delivers less total dose (51–56 Gy in 19–22 fx) than 55/20. Fowler's tolerance zone is 59–61 Gy10 for radiotherapy alone, 3–5 Gy10 lower with concurrent chemotherapy. | verified |
| S2 | The default suggestions depend on the tumor repopulation default λ_DI/α = 1 Gy/day. At Withers' 0.6 Gy/day (reference [2] on the page) no regimen qualifies, and 4 qualify when early BED may exceed the reference by 3 Gy. | verified |
| S3 | Several "suggestions" are trivial: gains of +0.15 to +0.55 Gy, and 32 × 2.10 Gy saves 3 days. | clinical review |
| S4 | Green and red badges read as safety clearance. HYHOPE dropped its 12-fraction arm for mucositis at a modeled acute BED of 52 Gy10, below the reference. | References/README A4 |
| S5 | There is no concurrent-cisplatin flag and no organ-at-risk conversion. Late BED is evaluated at the prescription dose only. | clinical review |
| S6 | The tab name "Suggested" and the disclaimer do not cover use for trial design. | clinical review |


## Decisions (Igor Shuryak, 2026-10-18)

"Fix all bugs and go with your recommendations." Concretely:

- **S1:** keep "early BED ≤ reference" as the default. Show the HYPNO arms as comparison rows and the Fowler zone next to the acute BED, and add an option to use the Fowler lower edge as the early limit.
- **S2:** keep the version 1 tumor defaults. Add Withers and HYHOPE presets and a note that the list depends on them.
- **S3:** by default, hide regimens gaining less than 1 Gy or saving fewer than 7 days. Both thresholds are adjustable.
- **S4/S6:** neutral wording in place of the badges. The tab is renamed and the disclaimer covers trial design. Add the HYHOPE 12-fraction note.
- **S5:** add a concurrent-chemotherapy checkbox. The organ-at-risk panel is deferred.
- **6 fractions/week:** removed at the user's request ("not really practical").

## Round 1 fixes

All defects D1–D11 were fixed.

- **Tests:**
  - The model tests now compare against an independent implementation: a day-by-day calendar, a `uniroot` calibration and a brute-force search.
  - New server tests use `shiny::testServer`.
  - The suite grew from 234 to 1029 expectations.
- **Red/green check:** each fix was undone in a fresh copy of the app, with the suite run in a separate process.
  - 28 of 29 cases were caught.
  - The survivor removed an explicit λ_DD = 0 check that the general calibration path already covers, so behavior was unchanged.
- **Runner bug found on the way:** the first runner reused folder names (the tests call `set.seed`), so planted bugs piled up across runs. It was fixed by using `tempfile()` and `callr`.

## Round 2 review

Four reviewers ran again.

- **First-round fixes:** D1–D6 and D8 were confirmed. D7 was partly fixed, because the slider fix introduced a new bug (below).
- **Calibration status:** independent `brentq` checks over 20,000 random parameter sets found 0 status mismatches.
- **HYHOPE tables:** reproduced with the app's formulas, which also confirms the calendar convention.

Round 2 findings, all fixed:

| # | Finding | Fix |
|---|---|---|
| R1 | The Withers preset used 0.6 Gy/day as BED. Withers gives it in 2-Gy fractions (p. 140), which is 0.72 Gy₁₀/day. | λ_DI/α = 0.72. The text is corrected. |
| R2 | Retyping the reference fractions (35 → 3 → 30) left the search range at [2, 2], and the page then reported "no regimen". | `search_range_for_ref()` resets the range, and it is tested. A live browser check showed 35 → 3 → 30 → 35 giving 5–34, 2–2, 5–29 and 5–34. |
| R3 | The Fowler limit stayed in use after the acute parameters were changed. | The search now shows a message and stops. |
| R4 | "Smaller of ΔDI and ΔDD" silently became ΔDI when DD was unavailable. | Returns `dd_unavailable` with a message. |
| R5 | Empty searches gave a generic message. | The message now gives the reason: empty dose range, or no shorter regimen. |
| R6 | The preset menu kept a stale name after edits, and re-selecting a preset did nothing. | Added a "Custom" entry and `preset_matching()`, checked live. |
| R7 | about.md: "tolerated" for HYPNO had no source. | Reworded. The text now cites Jacinto 2018 (40% grade 3 mucositis with cisplatin) and the 84% carboplatin series quoted in the HYHOPE protocol. |
| R8 | about.md: the smoothing-effect numbers were wrong. The effect is at most 0.23 d × λ/α. | Corrected. |
| R9 | about.md: Fowler's chemotherapy wording was overstated, and the 59–63 Gy₁₀ grey zone (Fowler 2010) was missing. | Corrected. |
| R10 | about.md: the HYHOPE note needed precision: postoperative parallel trial, 52 listed vs 51.3 computed, 51 for the continued 15-fraction level. | Added. |
| R11 | The Results tab showed an HYPNO comparison but no DD-vs-trial context. | Added test-minus-control BED differences: +0.9 Gy DI and −4.7 Gy DD at defaults, although the trial found non-inferiority. |
| R12 | The search tab did not say that HYPNO 55/20 exceeds the default early limit. | Added a dynamic note. |
| R13 | Wording issues: overall time shown as "26.00", "above it" verdicts, a semicolon in a message, undefined DI/DD/EQD2/LQ, the title did not define BED. | All fixed. |
| R14 | Tests: 17 of 28 new mutants survived, mainly in server code. | Gaps closed with new tests. The suite is now 1118 expectations, and a rerun caught 21 of 21 targeted cases (one after tightening the HYPNO note test). |

## Not addressed (scope)

- Accelerated regimens (6 fractions/week, twice daily) cannot be entered. The HYPNO control arm is shown as a fixed comparison row only.
- There is no organ-at-risk panel.
- There is no uncertainty analysis; it is planned next.
- Deployment through `rsconnect` has not been run. rsconnect is not installed, so the `.rscignore` exclusions are untested until `rsconnect::listDeploymentFiles()` is run.
- The project has no git repository. Tests are not committed, and creating a repository needs the user's approval.
