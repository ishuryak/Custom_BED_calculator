# Version 2.3 defaults (Igor Shuryak, 2026-10-19): tumor scenario D18 and the
# HYPNO 55 Gy/20 fractions acute limit. Expected values from the independent
# Python implementation (docs/development/2026_10_18_review_scripts/fix.py).

models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)
std <- list(m = 35, d = 2, T = 47)

test_that("default tumor parameters are scenario D18", {
  p <- default_params()
  expect_equal(unlist(p[c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g")]),
               c(alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.40, Tk = 28, g = 0))
  a <- compute_all(c(35, 20), c(2, 2.75), c(47, 26), p, std)
  expect_equal(a$bed_di, c(71.51428571, 70.12445771), tolerance = 1e-8)
  expect_equal(a$bed_dd, c(71.51428571, 62.84730123), tolerance = 1e-8)
  expect_equal(DEFAULT_EARLY_MODE, "hypno")
  # version 1 values stay available
  expect_equal(unlist(legacy_params()[c("alpha", "lambda_di", "lambda_dd", "Tk")]),
               c(alpha = 0.2, lambda_di = 0.2, lambda_dd = 0.5, Tk = 28))
})

test_that("the default search uses the HYPNO acute limit and is not empty", {
  p <- default_params()
  cap <- compute_all(20, 2.75, 26, p, std)$bed_early
  expect_equal(cap, 55.86597229, tolerance = 1e-8)
  s <- suggest_regimens(p, std, anchor = std, early_cap = cap, min_gain = 1,
                        min_days_saved = 7, d_max = 5)
  expect_equal(nrow(s), 3)
  expect_equal(c(s$m[1], s$d[1], s$T[1]), c(29, 2.25, 39))
  expect_true(all(s$bed_early <= cap + 1e-9))
})

test_that("the interface starts from these defaults", {
  env <- new.env()
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old))
  source("app.R", local = env)
  html <- as.character(env$ui)
  val <- function(id) as.numeric(sub(".*value=\"([^\"]+)\".*", "\\1",
           regmatches(html, regexpr(sprintf("id=\"%s\"[^>]*value=\"[^\"]+\"", id), html))))
  expect_equal(c(val("alpha"), val("ab"), val("lambda_di"), val("lambda_dd"), val("Tk")),
               c(0.35, 10, 0.23, 0.40, 28))
  expect_true(grepl("value=\"hypno\" checked", html, fixed = TRUE))
  expect_true(grepl("<option value=\"d18\" selected>", html, fixed = TRUE))
})
