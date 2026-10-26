# Version 3.2.1: response to the independent review of version 3.2 (2026-10-26)

Review: `HNSCC_BED_v32_review_and_action_plan_2026-10-26.md` and its evidence ZIP, in the project
folder above this repository. Plan: `PLAN_Codex_v32_review_2026_10_26.md` (same folder).
Branch `fixes_v32_review_2026_10_26` from f436e57 (reviewed HEAD, live 3.2.0 = 4876ae4).

## Decisions (Igor Shuryak, 2026-10-26, "go with recommendations")

1. Scope: F01-F07 and M01-M04 now. M05 (aligning the uncertainty-search default with Search, an
   optional tumor-share filter) and the scientific work package S1-S5 are deferred.
2. Two structures with the same DVH name: refuse, naming both line numbers (alternative: keep both
   as "Cord (1)"/"Cord (2)" and let the user choose).
3. F05: rollback around the existing two renames (alternative: immutable generation folders with a
   pointer file, which would also change how the app loads draws).
4. M02: remove the unused `dvh_n` input (alternative: keep it as a read-only note).
5. After the suites pass: push the branch and stop, and Igor runs the guarded deploy.

## Verification of the review before any change

Every finding was checked against the code before planning. F01 was reproduced by running the
review's two synthetic files through `parse_dvh_eclipse()` on f436e57 (parotid volumes 100/10/0
instead of 100/90/0, and only the second "Cord" block survived). F02, F03, M01 and M02 were confirmed
by reading the code. F04 and F05 were confirmed by their red tests. No finding was disputed.

## What was done, per item

Each test was run against f436e57 (a detached worktree outside the repository) and observed red
for the stated reason, then green on the fix.

| Item | Change | Red on f436e57 |
|---|---|---|
| F01 | Eclipse rows must have the header's field count (bracketed units), volume read from the column the header names, blank or repeated names refused before storing. CSV: per-line field count, blank names, a structure split into separate runs. | Both review files accepted, 4 Eclipse assertions failed. CSV cases were refused already, by later checks with misleading messages. |
| F02 | Per-quantity check before min/max (`gate.invalid_quantities`), constant sampled quantity fails, counts exactly 0, R-hat > 0, R gate requires one finite number per value. Shared table `tests/fixtures/gate_cases.csv`. Policy `gate-2026-10-27`. | Python 8 failures (negative counts, R-hat 0, NaN or Inf in one quantity, constant quantity, absent quantity), R 7 failures including the real loader accepting `1e999`. |
| F03 | Residual divided by the pair's alpha, alpha must be finite, positive and equal within the pair. | All 4 tests: residual 5e-5 log cells (0.0011 Gy) passed the 1e-4 "Gy" gate. |
| F04 | `posterior_target()` and `check_target()` in `fit/repair_dd_ess.py`, the target is recorded in `repair_ess.json`. | 9 failures, all "DID NOT RAISE". |
| F05 | `publish()` with rollback, recovery of an interrupted backup, completeness check, unique staging folder and export lock. | With the f436e57 logic placed in `publish()`, a failed second rename left no `draws/` (2 tests) and an incomplete stage was accepted. |
| F06 | Guide paragraph rewritten, 17 x 3.05 Gy at 12/15 pinned. | Guide assertions failed on the 3.2 text, the 12/15 test passes on both (selection unchanged, as intended). |
| F07 | `picks_sha256()` hashes `<i8` bytes, repair count 752, wording items. | The old encoding gives different hashes for int32 and int64, text assertions failed on the 3.2 text. |
| M01 | `deploy_identity()`, `deploy.R` has no `system("git ...")`. | Function absent, empty identities accepted by preflight. |
| M02, M03 | `dvh_n` removed, `fmt_tox_share()`. | 3548/4000 and 3611/4000 both showed "90%". |
| M04 | Status note in `References/README.md`, section G heading kept (pinned by a 2.6 test) with "historical". | Not tested (text only). |

## Evidence the fixes are live

- Export run for real on the fit folders: `python3 fit/export_draws.py` (first into a scratch
  folder, then into `draws/`). The four draw CSVs and `refit_point_values.json` are byte-identical
  to 3.2.0. `draws_manifest.json` changed only in the gate policy text and version, the new
  `per_quantity_check` note, and the two anchor residuals (3.89e-8 and 4.48e-8, which were log
  cells labelled Gy, now 9.840e-7 and 7.752e-7 Gy, matching the review). No staging, lock or backup
  folder was left behind.
- The new repair check was run against all ten 2026-10-25c DD folders (8 x 200 + 2 x 1000 = 3600
  fits) with their own data files and prior settings: all pass, and the prior-c10 folder passes only
  with sd 10.
- Repair count derived from the `repair_ess.json` records: 752 = 672 bulk only + 76 bulk and tail
  + 4 tail only.

## Things found along the way

- The per-quantity part of F02 cannot be applied retroactively. The fit folders keep one row of
  aggregate diagnostics per fit, not per-quantity summaries, so the 3600 existing fits were
  re-checked only on their aggregates. The manifest says so (`per_quantity_check`). No refit was
  done, because nothing identified an actual failing quantity.
- No model in this project samples a fixed quantity (paper-restricted `ab` and `g` are transformed
  parameters, and `log_ab_raw`/`g_free` still sample from their priors), so a constant sampled
  quantity is now a failure rather than an exclusion. The policy text says so.
- Monte Carlo errors: the "up to 0.25 Gy" for 95% limits had no recorded scenario. Re-measured with
  `2026_10_26_v321/mcse_benchmark.R` (cluster bootstrap over the 1000 DI draws, B = 2000, seed 1)
  for 18 x 2.9 Gy and HYPNO 55 Gy/20 against 35 x 2 Gy: medians 0.06-0.09 Gy, 95% limits up to
  0.28 Gy, shares up to 1.45 points (`mcse_benchmark.csv`). Guide and About now say "up to about
  0.3 Gy" and name the scenarios. The point estimates reproduce the review's section 4.4 exactly.
- The plan proposed a second selection hash (`picks_sha256_v2`). It was not added: a new identity
  field would change `identity_sha256` and make resuming an existing fit folder refuse. Making the
  existing encoding explicit (`<i8`) gives the same portability, and every recorded identity still
  verifies.
- `make_tracked_docx.py` passed the two changed table cells through unmarked (a known limitation).
  They were marked as whole-cell replacements afterwards, and the tracked guide was checked by
  round trip: accept-all equals the clean DOCX and reject-all equals the 3.2 DOCX, paragraph by
  paragraph (198 each).
- The review's full-suite errors (`digest(file = TRUE)` readability) did not reproduce here: the
  baseline run on f436e57 gave 5394 expectations, 0 failed, 0 errors.

## Not done (deferred by decision 1)

M05. S1-S5 (scenario dossier, estimator-aware likelihood, repair, duration and population
sensitivities, Monte Carlo and shortlist stability, external validation). Real de-identified TPS
exports for parser certification. A maximum-grid memory and latency benchmark.
