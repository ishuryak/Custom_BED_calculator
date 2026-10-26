# Version 3.2.1, review F02 (2026-10-26 v3.2 review): the R draws gate rejects every malformed
# diagnostic, using the same fixture table as tests/python/test_v321_gate.py.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "uncertainty.R"), local = TRUE)
keys <- c("divergences", "treedepth_hits", "max_rhat", "min_ess_bulk", "min_ess_tail", "min_ebfmi")
cases <- utils::read.csv(file.path(app_abs, "tests", "fixtures", "gate_cases.csv"),
                         colClasses = "character", na.strings = character())

test_that("F02: the shared gate fixture table gives the same verdicts in R as in Python", {
  for (i in seq_len(nrow(cases))) {
    dg <- list()
    for (k in keys) if (nzchar(cases[[k]][i])) dg[[k]] <- suppressWarnings(as.numeric(cases[[k]][i]))
    if (cases$max_rhat[i] == "abc") dg$max_rhat <- "abc"
    expect_identical(length(draws_gate_failures(dg)) == 0, cases$expect_pass[i] == "TRUE", info = cases$case[i])
  }
})

test_that("F02: vector, logical and NULL diagnostics fail", {
  good <- list(divergences = 0, treedepth_hits = 0, max_rhat = 1.001, min_ess_bulk = 900, min_ess_tail = 800, min_ebfmi = 0.9)
  expect_length(draws_gate_failures(good), 0)
  expect_gt(length(draws_gate_failures(modifyList(good, list(min_ess_bulk = c(900, 100))))), 0)
  expect_gt(length(draws_gate_failures(modifyList(good, list(divergences = FALSE)))), 0)
  expect_gt(length(draws_gate_failures(modifyList(good, list(min_ebfmi = NULL)))), 0)
})

test_that("F02: the real loader refuses a manifest whose ESS is the JSON literal 1e999", {
  draws_dir <- file.path(app_abs, "draws")
  tmp <- file.path(tempdir(), "draws_gate_1e999"); unlink(tmp, recursive = TRUE); dir.create(tmp)
  file.copy(list.files(draws_dir, full.names = TRUE), tmp)
  mf <- file.path(tmp, "draws_manifest.json")
  txt <- readLines(mf, warn = FALSE)
  # replace the first per-fit bulk ESS of the DD primary file with 1e999
  i <- grep('"min_ess_bulk_per_fit"', txt)[1]
  txt[i] <- sub('("min_ess_bulk_per_fit"\\s*:\\s*)[-0-9.eE+]+', "\\11e999", txt[i])
  expect_match(txt[i], "1e999")
  writeLines(txt, mf)
  expect_error(load_draws("dd", "primary", tmp), "ess_bulk")
  expect_s3_class(load_draws("dd", "primary", draws_dir), "data.frame")
})
