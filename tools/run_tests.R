# Runs the full test suite in the release library only (restore it first with
# tools/restore_library.R). Exit status 1 on any failure or error.
source("tools/env.R")
invisible(use_bed_library(required = TRUE))
environment_report()
cat("revision:", tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "unknown"), "\n")
r <- as.data.frame(testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE))
cat("\nexpectations", sum(r$nb), "failed", sum(r$failed), "errors", sum(r$error),
    "skipped", sum(r$skipped), "warnings", sum(r$warning), "\n")
if (sum(r$failed) + sum(r$error) > 0) quit(status = 1)
