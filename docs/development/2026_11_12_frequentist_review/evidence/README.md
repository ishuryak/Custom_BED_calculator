# Review evidence, 12 November 2026

The governing consolidated assessment is the root report
HNSCC_BED_frequentist_plan_and_Bayesian_refit_review_2026-11-12.md.
It contains six frequentist-plan findings and four current Bayesian
software/source/documentation findings, with one ordered execution roadmap.
Domain reports provide additional derivation and reproduction details.

Reviewed source: Custom_BED_calculator_v45, version 4.6.1,
HEAD 79e9c812f91482fa57ffa1dddabe89dc329a6303.
The current frequentist plan and audit were physical uncommitted additions;
their exact copies and digests are included under software/latest_plans
and in the mathematical contribution. They were not implemented or fitted
during this review. All 1,542 tracked files and both plans were checked again
at the end and remained unchanged.

## Evidence map

- plan_numeric_checks.py/.json: independently verified tumor binomial constant,
  historical likelihood and all ten selected-posterior contrast medians.
- math/plan_probes.py/.json/.log: quadrature counterexamples and positive
  adaptive controls; finite-cluster illustration; rounded early contrast
  tolerances; exact unchanged current primary Stan input and artifact checks.
- runtime_checks.R/.json/.log: current calendar and replay corrections,
  legitimate package positive, malformed package negatives that still pass.
- run_r_review.R, r_full.log, r_full_results.csv: complete testthat run,
  6,781 expectations, zero failures, one environment-level error. All 285
  expectations in the v4.6/v4.6.1 test files passed. Not a release-pass marker.
- readability_probe.R/.log: direct dataset read/hash succeeds although
  file.access reports unreadable; four recorded gate flags also checked.
- final_checks.py/.json: final tracked-file/plan equality and all 368 copied
  DT dependency files byte-identical to the pinned release installation.
- software/: canonical identities, cache/report/output inventories,
  exporter and metadata negative probes, builder reconstruction and limits.
- clinical/: clinical report and selected pages 391/394 of the user's local
  HYPNO PDF used for source verification. The full publication is not bundled.

The live main/methods check observed version 4.6.1, deployed revision 7ec2faf,
69 main points, correct early-section routing and the rendered tumor formula.
It was a targeted visual/accessible-UI check, not a new 73-item walkthrough.

## Reproduction and limits

Python numerical probes used 3.12.10, NumPy 1.26.4, SciPy 1.15.3, pandas 2.2.3.
pytest/CmdStanPy were unavailable in that interpreter. Stubbed imports in
specific probes explicitly prohibit sampling; methods-builder --check uses
its existing data-only fallback. No new MCMC/bootstrap was run.

R used 4.5.1, Windows UTF-8 and the project's release library, with a verified
byte-identical existing workspace DT copy to avoid resource-directory access
restrictions. Full-suite error remains recorded, not waived or marked passed.
readability_probe independently checks the blocked digest/gate assertions.

See software/README.md and math/math_review.md for probe commands. Supply
an exact source checkout and, for cache checks, the original retained caches.
The archive intentionally excludes copied application source, installed
libraries, full publications, posterior NPZ files and external raw chains.
Root R scripts and some root Python scripts record absolute review paths;
edit their source/output/library path literals for an isolated workspace
before rerunning elsewhere. Retention comparisons use the earlier reviewed
snapshot; restore that exact revision if reproducing those checks elsewhere.

Proposed acceptance tests and future runner invocations are not claimed to
have run. A zero-exit probe can demonstrate an invalid input was accepted;
read its labeled result instead of interpreting execution success as a gate
pass. Deliberately invalid fixtures were never installed in the project.

EVIDENCE_MANIFEST.csv records the SHA-256 and byte count of every other
archive member. No production edits, Git changes or deployment were made.
