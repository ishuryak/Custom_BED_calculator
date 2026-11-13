# Software and provenance review — 12 November 2026

This is the software contribution to the comprehensive review. It covers the completed early-toxicity pipeline changes, retained fit and validation artifacts, the new methods-metadata generator, and the latest frequentist comparison plan. It does not replace the separately assigned clinical, mathematical, runtime R, or browser reviews.

## Identity and scope

The reviewed original is `C:/Users/ishur/OneDrive/Desktop/Office/HNSCC_BED_calculator/Custom_BED_calculator_v45`. Local HEAD and local `origin/v45_review_2026_11_09` both resolve to `79e9c812f91482fa57ffa1dddabe89dc329a6303`. This is a statement about the existing local remote-tracking reference, not a fresh fetch of the remote service.

The canonical snapshot is `work/review_2026_11_12/software/source`. It was extracted with Git's `core.autocrlf=false`, then compared against Git objects and the original physical working tree: **1,542 tracked files, zero snapshot blob mismatches, zero original-working-tree byte differences**. Each path, Git blob, SHA-256, and comparison result is in `source_manifest.json`. Ordinary and explicitly configured `git status` continued to fail with “must be run in a work tree”; the object/byte comparison provides the useful tracked-file evidence without claiming a successful status command.

Two visible development documents are untracked and were copied separately to `latest_plans/`:

| Document | SHA-256 |
|---|---|
| `docs/development/2026_11_12_frequentist_comparison_plan.md` | `44d2324f5e7820d9b9155cef43c17f8d5d32bc09e5d25623f3c369c95a0cf292` |
| `docs/development/2026_11_12_frequentist_comparison_plan_audit.md` | `c1b445873cdb7c7026d727b960c247216326c5e1321dbce567f8002defd95deb` |

They are part of this review, but not part of the committed snapshot. The ignored Step6 NPZ caches were audited read-only at their original location. No original application files, caches, Git state, or deployed services were changed. No fits, MCMC, dependency installations, or raw-chain reconstruction were performed.

## Verification actually performed

Native Python was `C:/Users/ishur/AppData/Local/Programs/Python/Python312/python.exe`, version 3.12.10, with NumPy 1.26.4, SciPy 1.15.3, and pandas 2.2.3. `pytest` and `cmdstanpy` were unavailable. A fresh attempted pytest invocation failed with `No module named pytest`; **there is no full-Python-suite pass claim**. The numerical gate/cache probes use an explicit import stub whose sampler raises if called. The methods generator's own data-only fallback permits its actual command to run without CmdStan.

| Executed check | Result and scope |
|---|---|
| `python -B source/tools/build_model_methods.py --root source --check` | Passed on the current canonical snapshot. |
| Gate probes | Valid four-chain gate passes; NaN and positive infinity placed in every EBFMI chain position refuse; wrong chain count refuses; incomplete legacy manifest refuses. |
| Export source selection | Current U02 accepts. U04, U15, U17 and U19 refuse under the exact approved-primary profile check. |
| Synthetic exporter diagnostics | Rehashed diagnostics containing nonfinite gate summaries are still accepted: confirmed residual S1 below. |
| Current cache reuse | One actual 4,000-draw current cache loads under its reconstructed recorded identity; modifying one draw without updating its external cache manifest refuses. |
| Step5 retained outputs | All 62 manifests have required output files, matching file hashes and diagnostic identities/gates, and four recorded chain hashes. Fifty-nine recorded gates pass; the other three are the already identified binomial U01/U05/U07. |
| Step6 retained outputs | All 52 caches have 4,000 draws, finite numeric arrays, recorded passing complete-schema gates, and four recorded chain hashes. All 52 external cache-manifest hashes and identities match actual NPZs. |
| Existing raw-chain recheck record | All 114 rows are marked checked; 111 pass. All 114 identities match current 62 Step5 and 52 Step6 artifacts. Each row has four finite EBFMI values; minimum 0.6525243985397455. This verifies the retained record's consistency, not the inaccessible raw-chain files themselves. |
| Published Step6 reports | All seven files match `report_manifest.json`. |
| Source builder reconstruction | Five CSVs reproduce byte for byte: 66 arms, 68 outcomes, 28 arm schedules, 1,053 fractions, and 42 evidence families. Synthetic grade 1 and grade 2 mucositis rows refuse the grade 3/4 class; the accepted standard grade ≥3 code passes. |
| Methods negative controls | Changed early Stan prior refuses. Altered tumor and early exported draw bytes both incorrectly leave `--check` passing: S2. |

The cache identity positive control uses the recorded fit metadata versions (CmdStan 2.40.0, cmdstanpy 1.3.0, NumPy 2.4.3) only to reconstruct the existing deterministic hash. Actual array operations use installed NumPy 1.26.4; no claim is made that the fitting environment was executed. R tests and runtime validation are reported separately by the parent reviewer.

## Closure of relevant previous findings

| Previous requirement | Current evidence | Status |
|---|---|---|
| R4: Prevent a sensitivity fit being published as the approved primary; bind source outputs and diagnostics | `tools/early_tox/export_fit_for_app.py:57–89` checks manifest schema, gate/smoke flags, every recorded output hash, diagnostic identity, and exact primary profile. Four formerly accepted nonprimary units now refuse. | Substantially resolved; finite-number validation residual S1 remains. |
| R4: Incomplete publication must not appear complete | `export_fit_for_app.py:193–207` verifies staged hashes, removes the old completion manifest, installs data files, and installs the new manifest last. | Fail-closed design is present by source inspection. This is not a transactional old-generation-preserving swap; an interrupted publication can require retry. No crash-injection test was performed here. |
| R6: Every EBFMI chain must be finite and present | `step5_fit.py:330–357` checks all chain values and exact count; fixed constants also require finite mean/standard deviation. | Resolved in direct gate probes, with retained recheck consistency confirmed. |
| R6: Prevent reuse of changed NPZ contents under an unchanged embedded identity | `step6_validate.py:124–145,176–206` adds an external checksum/identity manifest; cache arrays are validated at lines 148–173. | Resolved for the demonstrated prior defect. Unmodified real cache succeeds; edited cache refuses. All existing cache entries match. |
| R6: Preserve fit and report output identity | Complete Step5 manifests, Step6 cache manifest and report manifest are present and consistent. | Resolved for current retained artifacts. The raw-chain archive remains outside this review's accessible environment. |
| R7: Do not infer severe mucositis equivalence from the word mucositis | `build_arms.py:120–156` requires an accepted standard pairing or recorded source-specific equivalence, and duplicate outcome/arm-endpoint guards are at lines 203–222. | Resolved for the demonstrated grade 1/2 collision. All current builder outputs reproduce. Clinical correctness of each curated equivalence remains part of the clinical review. |
| R1/R2/R3/R5/R8 | Clinical descriptions, trial-table content, calendar and runtime-loader checks, and menu behavior | Assessed by the parent/clinical reviewers; no independent blanket closure asserted here. |

## Remaining findings

### S1 — P2: The exporter still accepts nonfinite diagnostic gate values

**Location:** `tools/early_tox/export_fit_for_app.py:74–75`.

`read_unit()` checks inequality thresholds, but does not first require finite numeric scalars. In a disposable copy of actual U02, changing `rhat` to negative infinity and `ess_bulk`, `ess_tail`, and `ebfmi_min` to positive infinity, then updating the diagnostic file's recorded checksum, still passes `read_unit()`. The fit identity, model profile, draws and other output hashes remain unchanged. This is a malformed newly generated/repackaged metadata case; an ordinary unrecorded file modification is correctly caught by the new output hash checks.

**Effect:** The publication boundary does not enforce the finite-diagnostics contract now correctly implemented in the sampling gate. A passing Boolean and favorable infinities can be republished as an approved package. The current U02 diagnostics are finite; this is not evidence that its existing numerical results are wrong.

**Fix:** Validate strict finite numeric scalar gate fields before comparisons, reject Boolean substitutes, require nonnegative integer divergences and a complete four-value finite EBFMI vector, and require the recorded minimum and complete flag to agree with that vector. Validate missing/nonfinite-quantity lists and fixed-constant status consistently with the producer. Prefer a shared diagnostic-record validator used by export and other public readers rather than another independent reduced threshold expression. Preserve the accepted gate policy: tree-depth hits are currently reported, not an automatic failure; do not silently redefine it.

**Acceptance:** Actual U02 passes unchanged; legitimate nonprimary units still refuse; each required scalar with NaN, either infinity, missing values, wrong types, and wrong dimensions refuses before publication. A finite low EBFMI chain cannot be hidden by a favorable recorded minimum. A broken staged input leaves no newly published completion manifest. These tests need only fixtures, not new fits.

### S2 — P2: Methods freshness certifies manifest claims without checking deployed draw-file bytes

**Locations:** `tools/build_model_methods.py:144–166,197–214,232–239`. Generation writes directly at lines 255–256.

The new generator performs useful strong checks: reconstructing the tumor and early design/identity, verifying current scientific source against recorded identities, and checking prior statements. But its reported draw hashes are copied from manifests. The `sources` list contains draw manifests, not actual CSVs. In a minimal disposable current package, changing `tumor_di_primary.csv`'s first `alpha` value to 100 leaves `--check` passing. Restoring that file and changing `early_fit_draws.csv`'s first `r_e` to 100 also leaves it passing. A changed scientific Stan prior correctly refuses, providing a positive test that the identity checks are active.

**Effect:** The release-facing assertion “model_methods.json is current” can be true while the packaged draws disagree with the provenance it displays. Runtime loader rejection may limit downstream use, but does not make this standalone metadata freshness check complete. No mismatch was found in the current package.

**Fix:** Before generating/certifying methods metadata, validate actual referenced draw hashes and declared schema/row counts through the same package validators used for release. Include actual CSV digests among the generator's source dependencies. Check paired tumor files as a pair and early draws as their recorded generation. A failed check must not rewrite the valid methods artifact. When writing a new metadata artifact, use a temporary file and `os.replace` after successful validation to avoid leaving a truncated JSON file after interruption.

**Acceptance:** Current package still passes and reproduces the same scientific metadata. A changed value, removed row, stale counterpart, or mismatched manifest in each referenced draw file refuses. A deliberately changed model or design continues to refuse. Tests should include both primary tumor files and the early file. A failed build leaves the previous `model_methods.json` byte-identical. Repairing these guards does not require changing fit results.

### S3 — P2 implementation dependency: Bootstrap resume must be identity-based, not merely row-present

**Location:** untracked frequentist plan, line 101 (`fit_ml.py`, `bootstrap.csv`, “skip-if-present on resume”). No implementation exists yet, so this is a specification risk rather than a current runtime defect.

The plan gives deterministic per-replicate seeds, which is appropriate, but does not define a generation identity, completion schema, handling of partial CSV writes, duplicated replicate IDs, or changed data/settings on resume. A row-presence test can mix old and new fitting specifications or treat an interrupted row as a completed fit.

**Recommended implementation:** Reuse the established generation-identity approach. Bind the exact design bytes, likelihood/model version, optimizer/bounds/multistart settings, integration settings, bootstrap design, target count, base seed/RNG and relevant software/source version. Validate this manifest before resume. Store each replicate atomically with its ID, seed, identity, status and full diagnostics, then derive the combined CSV deterministically. Keep failed/incomplete fits visible; do not replace failed replicates until a success appears or silently reduce the denominator. Define which retries use the identical seed and which deliberate new run creates a new generation.

**Acceptance before the full bootstrap:** An unchanged completed run resumes with zero reruns; a changed identity refuses mixing; a truncated checkpoint is recoverable and never counted complete; duplicate IDs refuse; deliberate optimizer failure remains a recorded failure with the planned total visible; interruption/restart yields the same ordered successful records as uninterrupted execution. Run these cheap tests and the planned B=20 smoke before costly execution.

### S4 — P2 plan correction: The early T1 compares unlike hash objects

**Location:** untracked frequentist plan, line 111.

T1 requires a design SHA-256 to match the U02 manifest identity. U02's identity is not a pure design hash: `step5_fit.py:305–312` hashes data, settings, sampler, Stan source, environment metadata and extra unit/family fields. A correct frequentist design cannot be directly checked against that composite hash.

**Fix:** Retain two explicit checks: (1) reconstruct the archived Bayesian composite identity using the current methods-generator approach, proving the reference fit lineage; (2) compare canonical design arrays and a separately defined canonical design digest between the Bayesian data builder and the frequentist implementation. Do not use array dimensions alone. Record factor/family order and schedule arrays so a silent ordering change is visible.

**Acceptance:** Current design arrays and their canonical digest match; dropping an arm, changing one count, fraction time, or family order fails. Composite identity reconstruction separately succeeds under recorded metadata and fails if a source-model or source-data input changes. This is a plan/test correction, not a need to refit U02.

## Retain, recompute, refit, and order

1. **Retain the existing valid fits and caches.** The current 62 Step5 output sets, 52 Step6 caches and seven report artifacts have consistent hashes and identities. Retain the distinction between 59/62 passing Step5 fits and the three known failed binomial sensitivities; do not promote failed sensitivities or relabel them as passing. No new MCMC is justified by S1–S4.
2. **Fix boundary validation first (S1), then metadata certification (S2).** Add the focused fixture tests and use the actual current U02 as the positive control. Keep old valid artifacts unchanged during negative tests. Re-exporting the same valid draws or regenerating methods JSON is inexpensive if needed; compare resulting scientific arrays/metadata and report any solely provenance-related difference.
3. **Amend the untracked frequentist plan before implementation (S3/S4).** Preserve the accepted offline-comparison scope. Add identity-bound checkpoint semantics and separate source-lineage/design-equality assertions, together with the mathematical review's likelihood and quadrature requirements.
4. **Execute the installed-environment test suite when available.** This review's pytest attempt was blocked by a missing dependency, not an assertion failure. Do not state that all Python tests passed. Use the project's intended environment rather than altering the production environment during a review. R test outcomes and browser checks are separate evidence in the integrated report.
5. **Only then run the prospective frequentist smoke and full workload under the user's execution authorization.** They are new analyses, not mandatory reruns of the retained Bayesian fits. The proposed T1–T7 tests have not been executed here. A new Bayesian fit becomes necessary only if data, scientific model, or a material fit-generating input changes, or a retained raw-chain check fails for the relevant unit.

Completion requires the new guard positives/negatives, methods byte-binding tests, exact package hashes, corrected plan identities/checkpoints, and the integrated review's scientific/runtime criteria. There is no added requirement for another general review once those explicit criteria pass.

## Evidence map

All evidence is under `work/review_2026_11_12/software/`. `README.md` contains executable commands and limitations; `evidence_inventory.csv` contains artifact hashes and excludes the full snapshot/archive.

- `identity.py`, `identity.log`, `source_manifest.json`: canonical snapshot and tracked original comparison.
- `latest_plans/`, `untracked_development_files.json`: physical uncommitted plan inputs.
- `guard_probes.py`, `guard_probes.log`: finite-chain, endpoint, exporter positive/negative cases and retained recheck summary.
- `cache_probe.py`, `cache_probe.log`: one real cache positive and changed-byte refusal.
- `methods_probe.py`, `methods_probe.log`, `methods_check.log`: actual generator check, changed-draw residual and scientific-model refusal.
- `artifact_audit.py`, `artifact_audit.log`, `step5_manifest_audit.json`, `step6_cache_audit.json`, `rebuilt/`: retained artifact checks and builder reconstruction.
- `retained_checks.py`, `retained_checks.log`, `retained_identity_check.log`, `report_hashes.log`: all 114 identities, 52 cache-manifest entries, and seven report hashes.
- `pytest_attempt.log`: exact missing-pytest blocker; no test-pass inference.
