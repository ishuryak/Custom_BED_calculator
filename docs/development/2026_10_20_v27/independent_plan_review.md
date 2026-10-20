# HNSCC BED calculator: review of the proposed v2.7 fixes

**20 October 2026 | Detailed review and execution plan**

**Reviewed plan:** `the plan file (plan_original.md)` in `C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator`.

**Reviewed application:** `Custom_BED_calculator_v2`, version 2.6.0, branch `v2-toxicity-hypofractionation`, local HEAD `0c51ae69a4e080921d7137bfbc59adefad883f8b`. File and line references below refer to this snapshot; use the named function or text as the locator after editing.

**Status:** Review only. No application code, original plan, guide, Git history, repository settings or deployment was changed. The proposed v2.7 release has not been implemented or tested. The accompanying evidence contains tests of the existing v2.6 code and independent reproductions of the plan's numerical claims.

## 1. Recommendation and the four decisions

Proceed with E1-E4 after correcting the instructions in this document. Their purpose is appropriate: make the source branch discoverable, avoid presenting a one-model gain as the default result, describe sensitivity counts accurately, and put modeling limitations beside the results. No equation change is needed to accomplish those four goals.

Do not implement the original plan verbatim. It contains a helper-placement instruction that can prevent startup, a test that cannot establish browser defaults, incomplete walkthrough changes, and several misleading clinical phrases. It also leaves documentation changes until after the proposed verification or deployment.

| Decision | Recommendation | Reason |
|---|---|---|
| 1. Repository entry point | Add the v2 branch link in the app and a short v2 banner on `main`. Keep the default branch unchanged for this release. | This solves discoverability with a small, inspectable change. A repository-wide default-branch switch is unnecessary for E1. |
| 2. Candidate minimum change | Keep **+1 Gy**, with **both models** and **smaller change** selected by default. Explain the empty result. | No regimen passes this particular screen. Do not introduce an allowed loss merely to populate the opening table. Keep negative thresholds available as deliberate exploration. |
| 3. Display vocabulary | Use **selected parameter settings** in the display. Retain machine-readable `scenario` names and export columns. | The settings are a deterministic sensitivity grid, not a sample from an uncertainty distribution. |
| 4. Release | Implement, finish the guide, pass automated checks, inspect locally with `run_app.bat`, then deploy the verified v2 working tree and inspect the live app. | Local inspection is necessary but insufficient by itself. Pushing to GitHub and deploying to shinyapps.io are separate actions. |

These are review recommendations for your execution, not a record that you have approved or performed the changes.

## 2. What was checked, and what the evidence establishes

The review compared the plan with `app.R`, model and export code, `about.md`, release scripts, existing tests, the v2.6 browser walkthrough, README, and the 20 October guide. Independent checks reproduced the Candidate thresholds, ranking order, the distinct Agreement analysis, and replay behavior. The current full test suite was run on an unchanged workspace copy of the application.

**Current baseline: 2,543 expectations; 0 failures; 0 errors; 0 skips; 0 test warnings.** This verifies the existing tests against v2.6, not the proposed changes. The release library was selected explicitly through `BED_LIBRARY`; its 69 installed package DESCRIPTION versions matched the lockfile. The evidence records the paths, commands and package versions. The sandbox prevented Git introspection in the original checkout and affected the environment-report package readability check; the raw log retains a Git-discovery warning outside the test results, and the original HEAD/ref and source hashes were recorded separately. Do not treat those environment-report limitations as application test failures or conceal them as a completely identical deployment environment.

The new v2.7 UI tests, mutations, local walkthrough and live checks described below remain **to be performed after implementation**. This pass did not re-audit every publication from scratch or verify a new live deployment. It checks the targeted plan against the current code and relevant clinical interpretation. The original plan records that GitHub's default branch is `main`; remote repository metadata could not be independently retrieved through the web tool in this pass. Recheck that fact immediately before the README action.

Seven independent Python grid reproductions matched the R calculations with no row/order mismatches and maximum numerical discrepancy `1.83e-14` Gy. Seven candidate/agreement JSON replay and rebuilt-scenario checks also passed in a controlled version-only simulation. These provide specific evidence for preserving the numerical engine; they do not certify an unimplemented v2.7 release.

The local README identifies the previously deployed code as `cd5d835`; current local HEAD is `0c51ae6`. Preserve the distinction between an inspected local revision and a documented deployed revision.

### Numerical reproduction of E2

Use Scenario D18, tumor α/β 10 Gy, comparison 35 × 2 Gy in 47 inclusive days, the fixed 35 × 2/47 DD anchor, fractions 5-34, dose step 0.05 Gy, maximum 5 Gy/fraction, minimum 7 days saved, HYPNO acute cap and comparison late cap. These are the current opening controls apart from the criterion and ranking being deliberately varied.

| Candidate screen | Result | Interpretation |
|---|---|---|
| Current: either model, larger change, +1 Gy | 29 × 2.25; 28 × 2.30; 30 × 2.20 | All three meet the threshold under DI only. |
| Proposed: both models, smaller change, +1 Gy | No rows | Expected and informative. |
| Both, smaller change, 0, -1 or -2 Gy | No rows at each threshold | A modestly negative threshold still does not produce a qualifying regimen in this grid. |
| Both, smaller change, -3 Gy | 30 × 2.20; 29 × 2.25; 28 × 2.30 | All three meet an explicitly allowed modeled loss under DD. |

The nonempty conservative fixture should have these values:

| Regimen and inclusive time | ΔDI, Gy₁₀ | ΔDD, Gy₁₀ |
|---|---|---|
| 30 × 2.20 Gy; 40 days | +1.12000000 | -2.04954215 |
| 29 × 2.25 Gy; 39 days | +1.18839286 | -2.40804467 |
| 28 × 2.30 Gy; 38 days | +1.12628571 | -2.86484061 |

The conservative score is the smaller evaluable model difference. Here it is the DD difference, so the least negative DD change ranks first. Use a tolerance of `1e-6` Gy for these rounded numerical fixtures and check the exact fraction/dose order separately.

**The Agreement tab is a different screen.** Its defaults compare against HYPNO 55 Gy/20 fractions in 26 days, allow -3 Gy, and require only 1 day saved. Its reference-relative late cap is approximately 105.417 Gy₃, versus 116.667 Gy₃ in the default Candidate screen. Its three schedules are 17 × 3.05, 18 × 2.90 and 19 × 2.80 Gy, not the three Candidate schedules above. Do not describe navigating to Agreement as simply changing the Candidate threshold to -3 Gy.

### Mathematical and clinical interpretation

The proposed UI defaults leave `suggest_regimens(..., criterion = "either", rank_by = "best")` unchanged. This is reasonable: the app supplies explicit arguments, while existing direct callers retain their behavior. The equations, DD anchor calibration, validity masks, grid construction, cap calculations and ranking algorithm do not need to change for this release.

Requiring both models is a useful conservative screen **within these model assumptions**. It is not proof that a schedule is clinically conservative, safe or noninferior. Both are approximations; the DD implementation is calibrated to DI at an anchor and shares α within each paired setting. Agreement is therefore correlated and can coexist with common error. The published DD hypothesis was fitted to clinical data; the calculator's calibrated variant should continue to be identified separately. [Shuryak, Hall and Brenner, 2018](https://pubmed.ncbi.nlm.nih.gov/29534828/).

An allowed -3 Gy BED change is not HYPNO's clinical noninferiority margin. HYPNO compared specified treatment strategies, including permitted weekly cisplatin, and defined margins on clinical endpoints. The calculator has no validated conversion from its BED difference to that endpoint margin. [HYPNO primary publication](https://doi.org/10.1016/j.ijrobp.2025.12.056).

## 3. Findings that must change the written plan

### F1. The proposed caveat-helper location can stop application startup

**Priority: release blocker. Files: `app.R:105`, `app.R:653`; original plan E4.1-E4.4.**

`context_strip()` is defined inside `server()`. The plan says to put `caveat_strip()` beside it, then call that helper directly while constructing the top-level `ui`. At that point the server has not run, and the helper is not in the UI's scope. Following those instructions literally produces “could not find function 'caveat_strip'.”

**Fix:** define the static helper at top level before `ui <- page_sidebar(...)`, for example immediately after `reason_text()`. Keep `context_strip()` where it is unless there is a separate reason to move it. Give each instance a unique ID. Add an app-construction test and an actual fresh-page startup check, not just a server-output check.

### F2. The proposed fresh-session unit test does not test defaults

**Priority: release blocker in the verification plan. Files: original E2 tests; `tests/testthat/test_app_server.R:6-16`.**

A no-input `shiny::testServer()` run has `input$criterion == NULL` and `input$rank_by == NULL`. It does not execute the browser input bindings that read the `selectInput()` defaults. The proposed assertions therefore fail both before and after the intended UI changes. The read-only probe in the evidence demonstrates this on the current app. Shiny describes `testServer()` as reactive server testing with a mock session and explicitly set inputs. [Official testServer documentation](https://rstudio.github.io/shiny/reference/testServer.html).

**Fix:** use three complementary checks: inspect the rendered UI's selected options, test server behavior with explicitly supplied current-default inputs, and inspect a genuinely fresh browser session before changing any controls.

Do not reuse the legacy helper named `set_defaults()` as if it represented current page defaults. It explicitly sets the old calculator preset, reference acute cap, maximum dose 8 Gy, and either/best. Preserve it for historical regression coverage; introduce a separately named v2.7 fixture.

### F3. An empty table cannot verify conservative ranking

**Priority: high. Files: original E2 tests; `app.R:197-204`.**

`all(rows$d_di >= threshold & rows$d_dd >= threshold)` is true for zero rows. Checking only the proposed empty screen cannot distinguish worst ranking from best ranking. A ranking regression could pass that test.

**Fix:** assert the selected value of each control separately. Keep the zero-row +1 Gy assertion, and add the three-row -3 Gy fixture above with exact order. Mutate criterion and ranking independently. Restoring `rank_by = "best"` must fail a ranking-specific check even though the +1 Gy table remains empty.

### F4. E2 and E3 contain clinical and logical overstatements

**Priority: high. Files: `app.R:711-716`, `app.R:876-880`, `about.md:59`, original E2.2/E2.6/E3.3.**

- “DI or DD lists regimens that only one model supports” is false: OR is inclusive and can include regimens meeting the threshold under both models.
- “Which repopulation model is true” implies that one of these two parameterized approximations must be the truth.
- “Both models support” can read as clinical endorsement. Say that both modeled differences meet a chosen threshold.
- The 15 settings do not all share one α or one calibrated constant C. α changes across presets and C is recalculated within each setting. The shared α and calibration dependence concern the paired models within a setting.

Use the replacement text in Section 4. Keep the number of settings, the valid denominator and the invalid count visible. Do not convert the count to a confidence percentage.

### F5. E4 should distinguish omitted effects from assumptions of equality

**Priority: high. Files: original E4.1; `app.R:161-166`; guide chemotherapy section.**

No chemotherapy term is calculated in tumor, acute or late BED. “Assumed equal in both regimens” can imply a modeled cancellation that does not exist. A shorter course may change chemotherapy exposure even with the same weekly prescription. The published HYPNO schedules used four versus five planned weekly courses when chemotherapy was given. [HYPNO treatment description](https://doi.org/10.1016/j.ijrobp.2025.12.056).

Likewise, “gaps not modeled” is incomplete. A longer manually entered overall time affects the Results and comparison calculations. The model does not represent the location, sequence or repair consequences of individual interruptions, and the candidate search constructs regular weekday schedules.

“No uncertainty” is also ambiguous. There is substantial uncertainty; the application does not calculate uncertainty intervals. Use that distinction explicitly.

### F6. The documentation update list is incomplete

**Priority: high. Files: `about.md:67`, guide controls and explanatory sections, README historical text.**

The About paragraph beginning “With the defaults … three regimens” would become wrong. The plan updates the surrounding control descriptions but misses this paragraph. Retain the useful three-row example as an explicitly selected either/best example, with its “DI-only hypotheses” qualification, and state the new empty default first.

The guide also contains two substantive descriptions needing correction during this documentation release: absolute statements that α has no DI effect before onset, and assertions that the HYPNO cap already includes chemotherapy so adding a chemotherapy term necessarily double-counts it. The current model does not support those statements. Exact replacements are in Section 4.

### F7. Replay compatibility is narrower than “old files replay exactly”

**Priority: high for release claims. Files: plan Scope/V3; `R/export.R` compatibility checks.**

Compatible v2.6 scenarios with recorded search arguments retain the same calculations if only UI defaults and labels change. A simulation changing `APP_VERSION` to 2.7.0 confirmed this for the checked candidate and agreement fixtures. A synthetic fixture marked as v2.5 is deliberately refused by the current replay validator; explicit migration is a reanalysis. This checks the compatibility guard, not every historical export.

**Fix:** say “Compatible v2.6 scenario files replay with the same numerical results and selection rules.” Recorded schema, equation, selection, scenario-set, hinge and clock compatibility still governs acceptance. Do not promise byte-identical ZIPs or unchanged app/code metadata. A freshly created scenario using new UI defaults is also a different analysis from the old default scenario.

Keep `SELECTION_VERSION`, `EQUATION_VERSION`, schema and scenario-set versions unchanged only while their semantics actually remain unchanged. Update both `APP_VERSION` and `DESCRIPTION` to 2.7.0.

### F8. The copied walkthrough needs structural changes

**Priority: release blocker in the acceptance harness. File: `docs/development/2026_10_19_v26/walkthrough_v26.R`; original V6.**

The current script hardcodes the local target and v2.6 identifiers. It resets criterion/ranking to either/best after edge tests, so a later “default” test is not examining the opening page. Its `check()` records FAIL but does not make the process fail. Simply replacing check 6 and adding a few text searches is insufficient.

**Fix:** parameterize local versus live URL, expected version and output directory; check initial defaults before changing inputs; use explicit nonempty fixtures for table/plot/download checks; restore every changed input or start a new session; collect failures and exit nonzero; close the browser and stop any app process the script started. Keep the v2.6 script as historical evidence.

### F9. The release order tests an unfinished product

**Priority: high. Files: original “Order of work” steps 3-7; `DESCRIPTION:4`; `tools/deploy.R`.**

The plan runs final checks before updating some shipped documentation, and rebuilds the guide after deployment. Finish those changes first. `DESCRIPTION` is omitted from the version list. The deploy script uploads the local directory, uses the first configured shinyapps.io account, and does not itself enforce the test or clean-revision gates.

**Fix:** finish app, About, tests, guide and release descriptions before final validation. Record the verified revision and upload that working tree. Check the account and application target explicitly. Update the README's *live deployment fact* only after successful deployment; before then label 2.7 as the release candidate. Post-deployment timestamp and revision records may be a separate documentation commit.

## 4. Exact changes to make

### E1. Source link and repository banner

In `about.md`, the Version section around line 169, add this Markdown link:

```markdown
Source code for version 2:
[v2-toxicity-hypofractionation branch](https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation).
```

Keep the existing version/equation/selection identifiers and explanation of the historical v1 calibration issue. The branch link follows ongoing v2 development; it is not an immutable release snapshot. The release record and export code identity provide the exact revision. An additional verified commit permalink is useful but optional.

Avoid making the app depend on a perpetually correct claim about GitHub's default branch. If you include “the default branch holds v1,” verify it at release and maintain it when repository settings change.

For the recommended `main` README banner, insert immediately below its title:

```markdown
**Version 2:** use the [live v2 calculator](https://ishuryak.shinyapps.io/custom_bed_calculator_v2/) or its [source branch](https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation). This `main` branch contains version 1, whose DD calibration is incorrect when tumor α/β differs from 10 Gy.
```

Make this a separate README-only commit on `main`. Do not merge the v2 application into main as part of this task. the colleague's route to v1 is unknown; the lack of an app link and the repository entry point explain a discoverability problem, not proof of how he navigated.

### E2. Defaults, empty state and About text

At `app.R:197-202`, add `selected = "both"` to the criterion `selectInput()` and `selected = "worst"` to the ranking `selectInput()`. Keep option values, option order and helper-function defaults unchanged. Keep the `min_gain` initial value at 1.

In `output$suggest_msg`, preserve the existing early return for a specific reason such as unavailable DD or invalid reference. After that return, add a branch for `identical(input$criterion, "both")` with this text:

> No regimen in the searched range meets the selected minimum tumor BED change under both models and all other limits. An empty result is a valid outcome of this screen. The Model agreement tab explores a separate comparison: by default HYPNO 55 Gy/20 fractions, an allowed change of -3 Gy under each model, and at least 1 day saved. A limit defined relative to the comparison also uses that tab's comparison. To explore an allowed loss here, explicitly enter a negative minimum change. “DI or DD model” requires the threshold under at least one evaluable model; it can include regimens that fail the threshold under the other model. Changing these settings changes the hypothesis being screened.

The message should not imply that relaxing toxicity limits is a corrective action for an error. The existing non-both branch may remain, but any revised advice should describe deliberate exploration rather than how to obtain a positive result.

Update `about.md` at these locators:

- **Line 34, unavailable DD:** remove “default” from either/best. Add: “The default both-model criterion and smaller-change ranking then return no regimens, with an explanation.”
- **Line 59, Tumor bullet:** state both as default and retain +1 Gy. Replace the proposed “which model is true” rationale with: “Requiring both avoids listing a regimen by default when it meets the threshold under only one of these model assumptions. It does not establish clinical benefit or account for all model uncertainty.”
- **Line 63, Output bullet:** mark smaller change as default/conservative; larger change as optimistic.
- **Line 67, default example:** replace the opening with: “With the defaults (Scenario D18, HYPNO acute limit, both models, minimum change +1 Gy and smaller-change ranking), no candidate regimen meets all limits. If the criterion is explicitly changed to ‘DI or DD model’ and ranking to the larger change, three regimens meet +1 Gy under DI only: …” Retain the existing three numerical examples and “DI-only hypotheses” explanation. Replace “regimens that both models support” with “regimens whose modeled changes both meet the chosen allowed-loss threshold.”
- **Two-model agreement example around line 84:** replace “No regimen is supported by both models without a loss” with the scoped statement: “With the default parameters, grid and limits, no regimen meets a nonnegative tumor BED change under both models.”

If you later choose a -3 Gy opening Candidate threshold, treat that as a separate documented decision. Update all default text and tests, and note that the Sensitivity summary reads `min_gain` from Candidate: its threshold and counts will also change. Do not silently take that alternative while implementing the recommended plan.

### E3. Counts and sensitivity wording

At `app.R:833-842`, format counts from the existing data columns:

```r
counts <- ifelse(
  a$n_invalid > 0,
  sprintf("%d of %d evaluable settings (%d of %d not evaluable)",
          a$n_pass_both, a$n_valid_both,
          a$n_invalid, a$n_scenarios),
  sprintf("%d of %d settings", a$n_pass_both, a$n_valid_both)
)
```

Use the display header **“Selected parameter settings with both ≥ allowed change”**. The denominator is already in every cell; avoiding a hardcoded “of 15” in the header prevents future inconsistency and keeps the table narrower. Rename the two adjacent range headers to “Δ DI over evaluable settings” and “Δ DD over evaluable settings.” This is display-only.

In the Agreement help text say:

> The ranges and counts use 15 selected parameter settings: five tumor presets × three tumor α/β values. Only settings in which both model differences are evaluable enter the joint count and ranges. Within each setting, the models share α and DD is calibrated to DI at the anchor. These are sensitivity checks, not independent confirmations, confidence intervals or probabilities.

In `output$sens_summary`, keep the current per-model numerators and denominators. Change “evaluable scenarios” to “evaluable parameter settings”; retain the separate not-evaluable count. Append:

> These selected parameter settings are not independent confirmations or a probability sample. Within each setting, DI and DD share α and DD is calibrated to DI at the anchor; α and the calibrated constant can differ between settings. Uncertainty intervals are not calculated.

Use “selected parameter settings” in the Sensitivity introductory help as well. A display column can be titled “Preset” rather than “Scenario” if desired, but do not rename export fields. Keep `n_pass_both`, `n_valid_both`, `n_invalid`, `n_scenarios`, per-model validity fields, scenario IDs and JSON structure unchanged.

### E4. Shared caveat component

Define this **before `ui <-` at top level**, not beside the server-local `context_strip()`:

```r
caveat_strip <- function(id) {
  div(
    id = id, class = "bed-caveat alert alert-light",
    tags$strong("Model limits: "),
    "Chemotherapy effects, HPV status, stage, tumor volume, ",
    "treatment volumes and organ-at-risk doses are not modeled. ",
    "Overall duration is included; individual gap timing and ",
    "interruption effects are not modeled separately. ",
    "Default tumor parameters are exploratory. ",
    "Uncertainty intervals are not calculated. ",
    "See Model and sources."
  )
}
```

Use `tags$strong()` rather than literal Markdown asterisks inside ordinary R text.

Insert `caveat_strip("caveat-results")` immediately beside the Results calculation table, `caveat_strip("caveat-candidates")` beneath the existing Candidate warning, and `caveat_strip("caveat-agreement")` near the Agreement context/controls **before** its potentially long table. Keep the line visible on narrow screens and avoid placing the only caveat below a long result table.

“See Model and sources” may remain a textual navigation instruction. If it is made clickable, implement and test actual tab navigation; do not style inert text as a link.

For the guide, clarify that entering a longer overall time changes the aggregate time penalty but does not simulate a fraction-by-fraction interrupted schedule. The search continues to generate regular once-daily weekday schedules.

### Documentation corrections to include in the same release

In `BED_calculator_v2_references_and_interface_2026_10_20.md`, update the controls table around lines 339-340, Candidate examples, Agreement count labels, screenshots and limitations. Search the whole document for “default,” “DI or DD,” “larger,” “15 of 15,” “scenario,” “before onset,” “chemotherapy,” and “double-count” to find dependent statements rather than editing only two cells.

**Guide around lines 267 and 271, pre-onset DI explanation:** replace absolute claims with:

> In the DI model, α converts both background growth and accelerated repopulation into BED loss. With g = 0, accelerated repopulation is very small for a course ending well before the onset because the implementation uses a smooth threshold. With g > 0, shortening changes DI BED by g/α per day through background growth, including before the accelerated onset. Changing α can therefore affect pre-onset results.

This follows the implemented form `BED_DI = K - [lambda_di * hinge(T - Tk) + g*T] / alpha`. The smooth hinge has a tail; it is not identically zero before `Tk`.

**Guide around lines 405 and 410, HYPNO cap and chemotherapy:** replace the implicit-inclusion/double-count assertion with:

> The HYPNO cap is the radiation-only acute BED calculated for the HYPNO test schedule using the selected mucosal parameters. Its clinical context includes patients receiving chemotherapy, but the formula contains no chemotherapy effect. Comparable modeled acute BED does not ensure comparable toxicity across chemotherapy schedules, populations or treatment volumes. Adding a quantitative chemotherapy adjustment would require a separately justified model and calibration; it is not automatically double-counting an existing term.

The app's chemotherapy checkbox should continue to be described precisely: it changes the Fowler tolerance-zone convention, not any BED value and not the HYPNO acute cap.

**About limitations at lines 158 and 164:** replace “A limit relative to the comparison assumes the same chemotherapy in both regimens” with “A limit relative to the comparison compares radiation-only modeled BEDs; differences in systemic therapy and its effects must be considered separately.” Keep the useful warning about different course lengths and chemotherapy exposure. Replace the unqualified “Treatment gaps are not modeled” with “Overall duration is included, but individual gap timing and interruption effects are not modeled separately.” Apply the same clarification to the guide's comparison-relative early-limit row around line 355 and limitations around lines 415-416.

Update misleading historical README wording that describes v1 as comparator-calibrated. The corrected current About text distinguishes v1's fixed-anchor α/β-dependent calibration error from v2.0's comparator-dependent recalibration. Keep history accurate even in old-version sections.

Rebuild both clean and tracked guide documents from the revised source and preserve a `_pre_` copy. Because the existing tracked-change helper does not mark table-cell edits, provide an explicit table-change list and visually inspect both versions. Do not claim the tracked file captures every edit unless it actually does.

Optional source-description refinements, distinct from the required corrections: qualify the guide's HYHOPE late-toxicity statement by its reported 26.6-month median follow-up; restrict Hartley's 15% late-BED observation to the trials and endpoint analyzed; and distinguish original Fowler defaults from later refitted models. These need no new model implementation. Their source locators are in the companion clinical review.

## 5. Verification specification

Add a dedicated `tests/testthat/test_v27_colleague_review.R`. These are acceptance requirements, not claims that new tests have already passed.

### A. UI selection and startup

Load `app.R` in an isolated environment with the app directory as the working directory, following the existing source pattern. Render `htmltools::renderTags(e$ui)$html`. Locate the individual `<select id="criterion">` and `<select id="rank_by">` elements and their selected options. Assert `both` and `worst`; separately assert the `min_gain` input value is 1. A global search for “both” or “worst” is insufficient because both are already choices.

Require that sourcing/building the UI succeeds. This catches the helper-scope error. In a real fresh browser session, inspect the bound values before setting any inputs and confirm the visible labels agree.

### B. Server and numerical behavior

Build an explicitly named current-opening fixture with D18 values from `default_params()`, g = 0, α/β 10, default anchor, custom 35 × 2/47 comparison, HYPNO acute cap, reference late cap with zero excess, fraction range 5-34, dose maximum 5, step 0.05, 7 days saved, both/worst and +1 Gy. Set the normal Results-regimen inputs and Agreement inputs too, as their outputs need them. Do not inherit omitted values from legacy `set_defaults()`.

Check the following independently:

- At +1, zero Candidate rows and the new explanatory message. Confirm the table does not retain rows from a previous search.
- At -3, exactly three rows in the Section 2 conservative order, with the recorded DI and DD values and both valid differences.
- At explicit either/best/+1, the historical three-row order remains 29, 28, 30 fractions. This protects direct exploration and compatibility.
- At both/worst with unavailable DD, return the DD-specific explanation rather than the generic no-match paragraph.
- At Agreement defaults, retain the separate HYPNO comparison, -3 threshold, 1-day saving and its three schedules. Check context strings as well as rows.
- Confirm that the existing v2.6 validity, ranking, boundary, calibration and export tests still pass.

### C. Counts, denominators and caveats

Assert the default Agreement cell is “15 of 15 settings” and the header names selected parameter settings. Also test a case with some invalid paired settings: the pass denominator must be `n_valid_both`, the invalid count must be separately visible, and ranges must retain the existing joint validity mask. Use an existing v2.6 invalid-setting fixture to avoid inventing a second validity rule.

For Sensitivity, assert the selected-setting qualification and the per-model valid denominators. Do not test only that the word “independent” appears: require the intended negative statement and the within-setting qualification.

Inspect each tab's DOM for its dedicated caveat ID and shared class. Assert one instance on each intended panel, the common text, and no duplicate IDs. Existing panels already contain some chemotherapy or organ-at-risk language, so a broad keyword search is not a reliable failing-before test.

### D. Exports and replay

Freeze representative **compatible v2.6** candidate and agreement scenario files before editing. Include either/best, both/worst with a nonempty -3 fixture, and an empty valid search. Replay after the change and compare ordered numerical tables, all qualifying rows, validity masks, ranking eligibility, resolved limits and meaningful scenario inputs. Do not compare only three displayed rounded values.

Check that a new export records the chosen criterion/rank/minimum, app version 2.7.0 and the correct code identity, while CSV field names and schema/selection/equation/scenario-set identifiers retain their intended values. A fresh empty default search should still export a valid scenario and empty candidate result sets.

Keep a test that incompatible old files are refused and require explicit migration. Do not weaken validation to make a broad release-note promise true. Re-export timestamps, archive bytes and code hashes may legitimately differ.

### E. Red/green and mutation checks

Classify tests correctly. New feature tests should fail against v2.6 for their intended missing behavior; compatibility and unchanged-math tests should already pass. It is wrong to require every newly written regression assertion to be red.

Before committing, temporarily restore each relevant old behavior, one at a time, in a disposable copy or patch, run the specific test, and restore the intended implementation:

- Missing source URL -> E1 link test fails.
- Criterion selected as either -> criterion-default test fails.
- Rank selected as best -> ranking-default test fails, independently of the empty search.
- Old count labels -> E3 display test fails.
- Remove one caveat instance -> the test for that panel fails.
- Put the helper in server scope -> app-construction test fails.

Record the failing assertion and successful restoration. Do not retain mutation changes in the release tree.

### F. Browser walkthrough and local manual acceptance

Create `docs/development/2026_10_20_v27/walkthrough_v27.R`, preserving the v2.6 script. Implement local and live modes: local may start a test app; live navigates to the supplied URL and must not start a local process. Give each run a new evidence directory. Read expected version from the source or pass it explicitly; remove v2.6 filename/title/export assumptions.

A walkthrough must fail the process when any assertion or readiness/idle timeout fails. Accumulate a failure count, capture diagnostics, clean up, then return a nonzero exit status. Put lifecycle code inside `main()` with `on.exit()` or use an equivalent `tryCatch(..., finally = ...)`; do not rely on top-level `on.exit()` to clean up browser/server processes. Deliberately fail one harmless assertion once to verify this mechanism.

Perform these checks in order:

1. Fresh session: title/version 2.7, both/worst selected, +1 Gy, empty Candidate table, correct empty message and correct sidebar versus Agreement contexts.
2. All five tabs: Results, Candidate, Model agreement, Sensitivity, Model and sources open without errors. Expand all sidebar groups. Confirm the three caveats, count wording, source link and rendered equations.
3. Nonempty Candidate fixture: explicitly set -3 Gy and both/worst; check order 30, 29, 28, signs, plot and ZIP contents. Then explicitly test either/best/+1; do not call either state the default.
4. Agreement: restore its defaults and shared toxicity controls, confirm its distinct three schedules, settings counts and ZIP. Vary its comparison once and verify the context updates.
5. Validity: exercise unavailable DD and an invalid input; check the correct explanation and no stale table. Restore every changed control or reload a fresh session.
6. Chemotherapy checkbox: Fowler zone changes as documented; tumor/early/late BED and HYPNO cap do not acquire a chemotherapy term.
7. Sensitivity: sidebar regimen identity is clear, counts respond to its inputs and Candidate threshold, and wording does not imply confidence.
8. Narrow viewport: no hidden default selector, clipped caveat or unreadable table navigation. Check source-link destination, not merely its visible text.
9. Reload fresh: default controls and empty table return. Save logs, screenshots, downloaded scenarios and the tested revision.

Run the corrected automated walkthrough locally, then inspect with `run_app.bat`. After deployment, run the live mode plus the same key manual checks. Local success cannot establish that the live server received the same build.

## 6. Ordered implementation and release procedure

### Phase 0. Establish the baseline and preserve evidence

Open PowerShell in the v2 application directory, not the enclosing project. Copy commands from the editable Markdown companion to preserve long lines:

```powershell
Set-Location -LiteralPath 'C:\Users\ishur\OneDrive\Desktop\Office\HNSCC_BED_calculator\Custom_BED_calculator_v2'
git rev-parse --show-toplevel
git branch --show-current
git status --short
git rev-parse HEAD
```

Confirm the root and branch match this document. If there are pre-existing edits, preserve and inspect them before applying the plan; do not reset or overwrite them. Record the starting revision, current deployment record and baseline scenario fixtures. If HEAD has moved, compare the named locations and numerical fixtures before assuming these line numbers still apply.

Use the documented R 4.5.1 release library:

```powershell
$bedR = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'
& $bedR tools/run_tests.R
if ($LASTEXITCODE -ne 0) { throw 'Baseline tests failed' }
```

If the release library is missing, run `tools/restore_library.R` first and resolve restoration errors. A successful `run_app.bat` using fallback system packages does not establish the locked environment.

### Phase 1. Add acceptance tests and preserve compatibility fixtures

Implement the tests in Section 5. Run the new feature tests against current code and record expected failures. Keep the existing either/best regression fixtures unchanged. Save the v2.6 replay inputs and full comparison outputs before editing.

### Phase 2. Implement E1-E4 and version metadata

Make the source link, explicit UI selections, corrected message, settings labels and global caveat helper. Update `APP_VERSION` in `R/bed_models.R` and `Version` in `DESCRIPTION` to 2.7.0. Keep biological parameter defaults, helper defaults, formulas, caps, schema and selection identifiers unchanged.

Run targeted tests after each item. Inspect the Git diff for accidental model changes, renamed machine-readable fields and edits outside the stated scope. Perform the independent mutations and restore the final implementation.

### Phase 3. Finish all shipped documentation and the walkthrough

Update About, README release-candidate description, `docs/FIXES.md`, and `docs/development/2026_10_20_eric_wang_review.md`. Include each E item, exact decision, numerical baseline, compatibility qualification and validation evidence. State that the fresh page changes from three Candidate rows to zero; helper search behavior for explicit existing arguments is unchanged.

Complete the guide source, clean/tracked renderings, table-change list and screenshots. Finish the walkthrough corrections. Do not claim that v2.7 is live yet.

### Phase 4. Verify the finished release candidate

Run the full locked-library suite, corrected local browser walkthrough, replay comparisons and manual inspection. Require zero failures/errors and investigate any new warning or skip. The expectation count may exceed 2,543; that number is a baseline, not a target to preserve by removing tests.

Inspect the guide renderings and compare the source/README/About default descriptions. Record final source hashes and environment. If an app, About, dependency or functional test change follows these checks, rerun affected checks and the full suite before release.

### Phase 5. Commit and push the reviewed v2 changes

Review and stage the named app/test/document files deliberately; avoid an indiscriminate `git add .` that could capture downloaded ZIPs, temporary outputs or account files. Confirm `git diff --cached` matches the intended scope. Commit on `v2-toxicity-hypofractionation`, record the resulting SHA, and verify the working tree contains no unreviewed changes.

Push only that branch:

```powershell
git push origin v2-toxicity-hypofractionation
if ($LASTEXITCODE -ne 0) { throw 'Push failed' }
git rev-parse HEAD
git ls-remote origin refs/heads/v2-toxicity-hypofractionation
```

Verify the local and remote SHAs agree. If you commit final test logs afterward, distinguish the evidence-only commit from the tested source content and verify unchanged application hashes. Do not claim a specific commit was tested when only a different file state was tested.

### Phase 6. Add the separate main README banner

First verify current remote state:

```powershell
gh api repos/ishuryak/Custom_BED_calculator --jq '.default_branch'
git fetch origin
git show origin/main:README.md
```

Use a separate clean worktree from `origin/main`, or GitHub's README editor, for the one-banner change. For a local worktree, create a uniquely named branch such as `codex/v2-source-banner` from `origin/main`, edit only README, review and commit it, and use the repository's normal PR/merge route to main. Do not switch or reset the verified v2 deployment directory. Confirm the resulting main diff is README-only and that the rendered link targets the v2 branch.

The banner can be published independently of the v2.7 deployment because a v2 app already exists. It must not claim that 2.7 is live before it is.

### Phase 7. Deployment preflight and upload

Return to the verified v2 directory. Confirm the committed source state is still the one tested, with no application changes since local inspection. The existing deploy script names `custom_bed_calculator_v2` and guards against the v1 name; retain that.

Inspect the account selection. `tools/deploy.R` filters to shinyapps.io and takes `acc$name[1]`. Require that the selected account is the intended `ishuryak` account. If more than one account exists, either select `ishuryak` explicitly in the script and validate that exactly one match exists, or resolve account selection before running it. A suitable hardening immediately after obtaining accounts is:

```r
acc <- acc[acc$server == "shinyapps.io" &
           acc$name == "ishuryak", , drop = FALSE]
if (nrow(acc) != 1L) stop("Expected one ishuryak shinyapps.io account")
```

If you make this script change, include it in the reviewed release diff and final checks. Never place tokens or secrets in the plan, logs or repository.

Inspect `rsconnect::listDeploymentFiles(".")` using the release/deploy libraries configured by the script. Verify the application, About and needed R files are present; confirm tests, reference PDFs and development outputs are excluded as intended. Check the deployment record targets the v2 application/account. The current script uploads local files, not whatever happens to be on GitHub. Official deployment documentation confirms the roles of the app directory, account, stored records and bundle exclusions. [rsconnect deployApp documentation](https://rstudio.github.io/rsconnect/reference/deployApp.html).

After the local acceptance gate is satisfied, deploy:

```powershell
& $bedR tools/deploy.R
if ($LASTEXITCODE -ne 0) { throw 'Deployment failed; do not mark 2.7 live' }
```

### Phase 8. Live acceptance, record and rollback

Open a fresh session at `https://ishuryak.shinyapps.io/custom_bed_calculator_v2/`. Run the corrected walkthrough in live mode and inspect the version, defaults, empty-state text, three caveats, settings labels, source link and both downloads. Check a downloaded scenario for app version, explicit selected settings and code identity consistent with the deployed build.

Only after these pass, update the README live record with deployment time and uploaded source revision. Keep this factual post-deployment update separate from runtime changes. Preserve the deployment logs and acceptance artifacts. A browser page that merely loads is insufficient evidence.

If deployment or live acceptance fails, do not mark the release complete. Keep the failed-run evidence. Restore the recorded last known good v2 source in a separate directory, restore its documented dependencies, verify the intended account/app, and redeploy that snapshot to the same **v2** application. Use a known successful deployment record or verified source snapshot; do not assume local HEAD was the last deployed build. Do not force-push, reset unrelated work, redeploy v1, or use `upload = FALSE` as if it selected an arbitrary historical release.

Drafting a factual reply to the colleague is optional after successful deployment. Report only changes and checks actually completed. No message needs to be sent as part of implementing this plan.

## 7. Completion checklist and remaining scope

- E1 link opens the v2 branch; main's banner is accurate and README-only.
- Fresh Candidate page selects both/worst/+1 and explains zero rows.
- The deliberate -3 fixture has the expected three rows and conservative order; either/best remains available.
- Agreement retains its distinct comparison and limits; settings counts expose valid and invalid denominators.
- Shared caveats render on all three intended tabs without startup or layout errors.
- About, README, guide and screenshots agree with actual defaults and accurately describe chemotherapy, gaps and pre-onset DI behavior.
- App and DESCRIPTION versions are 2.7.0; unchanged numerical and file-format semantics retain their identifiers.
- Compatible v2.6 replay fixtures match; incompatible older files remain refused unless explicitly migrated.
- Full tests, failure-sensitive walkthrough, local inspection and live checks pass with recorded source identity.
- Live release metadata is written only after verified deployment; rollback evidence is retained.

This release improves transparency and the opening screen. It does not validate a trial regimen or resolve the uncertainty in HNSCC repopulation parameters. For a trial-design hypothesis, the useful output remains a clearly specified comparison with explicit population, systemic therapy, target volumes, organ-dose constraints, modeled tradeoffs and a separately justified clinical endpoint and margin. The two-model screen helps expose dependence on assumptions; it cannot substitute for those design decisions.

**Deliverables:** this editable review, a PDF reading copy, and a companion evidence archive containing independent review notes, numerical/replay scripts and outputs, baseline tests, package checks and source manifests. No proposed implementation is included as an applied app change.
