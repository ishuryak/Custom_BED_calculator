# Version 4.5.0: implementing the independent review of 4.4.0 (2026-11-09)

## Source and decisions

Igor Shuryak asked on 2026-11-09 to "analyze, check and implement" the independent review `HNSCC_BED_v44_independent_review_and_action_plan_2026-11-09.md` (kept outside the repository, in the project folder). It supersedes parts of an earlier same-day site review (`HNSCC_BED_v44_site_review_2026-11-09.md`). The review's five decisions were adopted as recommended, under that instruction. Igor did not rule on them one by one, so each is listed here for confirmation before deployment:

| Decision | Adopted | Note |
|---|---|---|
| D1. v4 as the default branch, READMEs repaired | READMEs repaired on this branch (v4) and in a separate `main` README change. The default-branch switch and PR #2 need repository permissions and were **not** done. | Needs Igor. |
| D2. Sensitivity and both-model screen | Kept, as sections of Custom search (advanced), not deleted. | The earlier site review recommended deletion; the independent review's retain-under-advanced option keeps the named-assumption and AND-criterion tasks. |
| D3. Merge Report card, Results, Compare; rename Search | Done: six destinations including the main page. | Settings and sources holds Model settings and Model and sources as two tabs. |
| D4. Trial-fit point scenario as default | Done, with complete profiles and the D18 fallback. | Changes the default list: 60 → 69 doses, both-model 17 → 65. |
| D5. EQD2 first | Done for the main tumor column only. Screening stays in BED. | The Uncertainty page still reports BED (see "Not done"). |

## Checks of the review before implementing

- Grid fixtures reproduced exactly with the review's R block: 60 and 69 points, 58 shared, removed 19/2.75 and 21/2.55, both-model 17 and 65, support 0.66025 at 18 × 2.90.
- Preset defect confirmed in code: the matching observer ignored g and the DD rate, and D18 did not reset g.
- Jackson, CAIR and Princess Margaret statements checked against the full texts (`References_refit_trials/fulltext_text/9165135.txt`, `10799721.txt`, `17011446.txt`, `17920715.txt`).
- Windows checkout defect confirmed: a `core.autocrlf=true` clone of `origin/v4` matched 0 of 4 draw hashes.
- One review number differs in the last digit: ΔDD for 35 × 2.1 against 35 × 2 Gy at 47 days is 3.0185, quoted as 3.019. The explainer says 3.02.
- An earlier statement of the site review ("the main list barely changes" under the refit default) was incomplete: the rows per fraction number barely change, but model agreement changes from 17 to 65 of the qualifying doses.

## What was implemented

See `docs/FIXES.md` (2026-11-09) for each fix with its proof, test and red evidence. Commits on `v45_review_2026_11_09`, merged with the evidence track `v45_docs_2026_11_09` (trial list, arm table, explainer, About refit section, `refit_in_brief.md`, Jackson check).

## Jackson dose-encoding check (review 2.1, first half of step 3)

`fit/jackson_exposure_check.py`, record `docs/development/2026_11_09_jackson_sensitivity.md`. At the fitted DI point values the accelerated arm's predicted control is 62.97% with the uniform 66 Gy encoding, 61.56% with the weighted exposure and 61.51% with the prescription mixture. Per draw the mixture minus uniform is −1.43 points (−1.93 to −0.96). No refit was run. A separate observation: at these values DI predicts a Jackson contrast of about +19 points against +4.8 observed; the arm stays inside its in-sample interval because of the arm-level effect.

## Not done, and why

- **Targeted DI refit with the Jackson amendment** (review 2.1 step 3, second half, and step 4): a scientific decision that starts a versioned fit release. Recommended next step only if the encoding change is judged material; the point-value shift above is about 1.4 points in one arm.
- **Uncertainty page in EQD2**: it still reports BED. Converting it needs per-draw conversion before quantiles and a new label set; deferred to keep this release a display change on the main page only.
- **About length**: reordered (purpose, assumptions and fitted-model summary first, release history last) but not cut to the 3,000 words suggested by the earlier site review. It is about 9,200 words.
- **Default branch, PR #2, credential revocation**: need Igor's GitHub and shinyapps.io accounts.
- **Python CI job** (review 6.5): not added. `tests/python` passes locally (208).
- **Deployment**: not done; needs Igor's approval and a resolved credential status.

## Measurements

- Baseline suite on `3b72f85`: 6269 expectations, 0 failed (1 error caused by an edit during the run).
- Local walkthrough on the new layout: 54 of 54 checks.
- Helper agent (evidence track): about 28 minutes.
