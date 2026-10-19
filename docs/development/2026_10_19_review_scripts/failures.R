# Print only failing expectations, one line each: Rscript failures.R [file]
args <- commandArgs(TRUE)
files <- if (length(args)) args else list.files("tests/testthat", "^test.*\\.R$", full.names = TRUE)
for (f in files) {
  res <- testthat::test_file(f, reporter = testthat::ListReporter$new(), stop_on_failure = FALSE)
  for (t in res) for (e in t$results) if (!inherits(e, "expectation_success"))
    cat(sprintf("[%s] %s :: %s\n", basename(f), t$test, substr(gsub("\\s+", " ", conditionMessage(e)), 1, 300)))
}
