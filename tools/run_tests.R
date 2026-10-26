# Runs the full test suite in the release library only (restore it first with
# tools/restore_library.R). Exit status 1 on any failure or error.
source("tools/env.R")
invisible(use_bed_library(required = TRUE))
environment_report()
source("tools/deploy_preflight.R")
rr <- release_revision()
cat("revision:", if (rr$ok) rr$rev else paste("not certifiable,", rr$why), "\n")
r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE))
cat("\nexpectations", sum(r$nb), "failed", sum(r$failed), "errors", sum(r$error),
    "skipped", sum(r$skipped), "warnings", sum(r$warning), "\n")
if (sum(r$failed) + sum(r$error) > 0) { unlink(".suite_passed"); quit(status = 1) }
# Record the revision the full suite passed on (tools/deploy.R requires it). Only a committed,
# resolvable revision counts; otherwise this run is a numerical check only (review R09).
rr_end <- release_revision()
if (rr$ok && rr_end$ok && identical(rr$rev, rr_end$rev)) {
  writeLines(rr$rev, ".suite_passed"); cat("suite marker written for", rr$rev, "\n")
} else {
  unlink(".suite_passed")
  cat("numerical check only, not release-certified:", if (!rr$ok) rr$why else if (!rr_end$ok) rr_end$why else "HEAD moved during the run", "\n")
}
