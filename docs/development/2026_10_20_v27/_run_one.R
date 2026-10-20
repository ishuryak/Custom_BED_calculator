source("tools/env.R"); invisible(use_bed_library(required = TRUE))
r <- as.data.frame(testthat::test_file("tests/testthat/test_v27_colleague_review.R", reporter = "summary", stop_on_failure = FALSE))
print(r[, c("test", "nb", "failed", "error")], right = FALSE)
