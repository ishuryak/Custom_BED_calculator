# Version 4.3.0: fixes from the external review of 4.2.1 (2026-10-31)

Source: `../HNSCC_BED_v4_review_and_action_plan_2026-10-31.md` (12 findings, F01 to F12, with evidence archive `HNSCC_BED_v4_review_evidence_2026-10-31.zip`, outside this repository). Igor Shuryak asked on 2026-10-31 for the review to be checked and acted on autonomously. Base: branch `v4` at 917825b (4.2.1, live). Work branch: `review_v4_2026_10_31`. Not pushed, not deployed.

## Check of the review

Every numerical claim used below was reproduced before any change, on this machine (R 4.5.1, pinned library `library-2026-10-15-R4.5.1`):
- Defaults: 60 qualifying points in 15 rows (30 to 16 fractions), 43 DI-only, 17 both, 0 DD-only.
- λDD = 2: 30-fraction DD changes +2.515833, +3.506194, +2.837178, +1.946087, +1.087580 Gy₁₀ at 2.05 to 2.25 Gy (maximum at 2.10 Gy, minimum at 2.25 Gy).
- Negative-reference fixture (λDI = 2, direct DD α 0.35, λ 0.4, C/α 50, g 0): standard DI −24.57143, DD 73.88970; 30-fraction doses 2.15 to 2.25; 4.2.1 printed "2.10 Gy fails ()".
- Direct DD α 0.005, λ 5, C/α 10, g 0: standard DD BED −284.35 (finite, negative), and 4.2.1 showed no notice.
- Custom anchor 20 × 2.75 Gy/26 days: standard DI/DD 71.51429/81.49437 while 4.2.1 said the models agree at 35 × 2 Gy.
- The transfer observer sent `d = 8` for row 11 (20 fractions) in 4.2.1.

The review's environment failures did not reproduce here. The unmodified full suite on 917825b gave 5691 expectations, 0 failed, 0 errors (the review's run had 1 failure and 2 errors from file-access and DT `setwd` problems and an incomplete copy). pytest was available here (the review had none): 114 passed on 917825b.

## Decisions taken (all as recommended by the review, Section 8)

| Question | Decision |
|---|---|
| Change Main margins, either rule or defaults? | No. Text and validation only. |
| Five BED units as fixed EQD2 across α/β? | No. Keep the BED rule, show units and the EQD2 equivalent 5/(1 + 2/(α/β)). |
| Make the row transfer weaken Search, Both-model or OAR criteria? | No. Scope the message. The optional "use Main policy in Search" button was not built (not requested; can be added). |
| Add an exact Main export? | Yes, new schema `main-1`, separate from Candidate/agreement exports. |
| Edit old identity hashes to add anchors? | No. Legacy identities resolve a fixed historical anchor only after the parent DI Stan-data hash reproduces with it; otherwise refused. |
| Refit to strengthen guards? | No. No fit was run. |
| F12 policy | Bounded retry: a failed draws load is retried at most every 300 s; a success is kept for the process. |

Judgment calls not specified by the review: the EQD2 equivalent is shown to one decimal (4.2 Gy at α/β 10, previously "about 4 Gy"); the transfer refusal shows a warning notification and sends nothing; the Main export stores support as "unavailable" (no column) when the draws cannot be read.

## Not done

- No deployment and no push. The live site still runs 4.2.1.
- No live browser walkthrough of 4.3.0 (needs deployment or a local browser session).
- Optional explicit "use Main policy in Search" action (F04) not built.
- F06/F07: interrupted-repair restart covered only by the existing run-twice test.
- The guide's Word version was built with pandoc and its text checked; it has not been opened in Word for a visual check.

## Verification

See `docs/FIXES.md` (row 2026-10-31) for each finding's test, red-on-4.2.1 result and the combined run.
