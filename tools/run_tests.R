# Runs the full test suite in the release library only (restore it first with
# tools/restore_library.R). Exit status 1 on any failure or error.
source("tools/env.R")
invisible(use_bed_library(required = TRUE))
environment_report()
cat("revision:", tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "unknown"), "\n")
r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE))
cat("\nexpectations", sum(r$nb), "failed", sum(r$failed), "errors", sum(r$error),
    "skipped", sum(r$skipped), "warnings", sum(r$warning), "\n")
if (sum(r$failed) + sum(r$error) > 0) { unlink(".suite_passed"); quit(status = 1) }
# Record the revision the full suite passed on (tools/deploy.R requires it).
# Only a clean tree counts: uncommitted changes are not that revision.
if (length(system("git status --porcelain --untracked-files=no", intern = TRUE)) == 0)
  writeLines(system("git rev-parse HEAD", intern = TRUE), ".suite_passed") else unlink(".suite_passed")
