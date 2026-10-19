# Regression tests for the independent review of version 2.1 (2026-10-18).
models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)
source(testthat::test_path("..", "..", "R", "export.R"), local = TRUE)
p <- legacy_params()   # version 1 tumor values (pinned results)
std <- list(m = 35, d = 2, T = 47)
a20 <- list(m = 20, d = 2, T = 26)

test_that("an unusable model shows no difference and does not enter the score", {
  q <- modifyList(p, list(lambda_di = 1, lambda_dd = 2))    # comparison DI BED = -11
  s <- suggest_regimens(q, std, anchor = a20, tol_early = 10, d_max = 4)
  expect_gt(nrow(s), 0)                                      # DD is usable here
  expect_false(any(s$qual_di))
  expect_true(all(is.na(s$d_di)))                            # was +57.8 to +63.2 in 2.1
  expect_equal(s$score, s$d_dd, tolerance = 1e-12)
  expect_equal(s$d_min, s$d_dd, tolerance = 1e-12)
  expect_true(all(is.na(attr(s, "feasible")$d_di)))
})

test_that("both models negative gives its own reason", {
  q <- modifyList(p, list(lambda_di = 1, lambda_dd = 5))
  expect_equal(attr(suggest_regimens(q, std, anchor = a20, tol_early = 10), "reason"),
               "negative_reference_both")
})

test_that("a candidate with a negative DD BED cannot qualify on DD", {
  q <- modifyList(p, list(alpha = 0.1, lambda_di = 0.2, lambda_dd = 2, Tk = 28))
  expect_gt(compute_all(35, 2, 47, q, std)$bed_dd, 0)        # comparison +46 Gy
  expect_lt(compute_all(34, 5, 47, q, std)$bed_dd, 0)        # 34 x 5 Gy: -39 Gy
  s <- suggest_regimens(q, std, anchor = std, criterion = "dd", rank_by = "dd",
                        early_cap = 1e6, late_cap = 1e6, min_gain = -1000, d_max = 5)
  f <- attr(s, "feasible")
  expect_gt(nrow(f), 0)
  expect_true(all(f$bed_dd > 0))
})

test_that("sensitivity rows use their own tumor alpha/beta", {
  st <- sensitivity_table(list(m = 20, d = 2.65, T = 26), std, std, p)
  v1 <- st[st$preset == "calculator", ]
  # independent: candidate ends before Tk (no repopulation); comparison loses 19 Gy
  hand <- function(ab) 53 * (1 + 2.65 / ab) - (70 * (1 + 2 / ab) - 19)
  expect_equal(v1$d_di, hand(c(7, 10, 15)), tolerance = 1e-3)
  expect_equal(v1$ab, c(7, 10, 15))
})

test_that("exports with absolute caps and the Fenwick option rerun exactly", {
  for (srch in list(list(early_cap = 55.86597229, late_cap = 105.41666667, fenwick = FALSE),
                    list(early_cap = 59, late_cap = NULL, fenwick = FALSE),
                    list(early_cap = NULL, late_cap = NULL, fenwick = TRUE))) {
    args <- c(list(m_range = 5:34, d_step = 0.05, d_max = 4, criterion = "either",
                   rank_by = "best", tol_early = 0, tol_late = 0, min_gain = -10,
                   min_days_saved = 7), srch)
    sc <- build_scenario(p, std, std, args, context = list(comparator = "custom"))
    sc2 <- jsonlite::fromJSON(scenario_json(sc), simplifyVector = TRUE)
    r1 <- do.call(suggest_regimens, c(list(p = p, ref = std, anchor = std), args))
    r2 <- rerun_scenario(sc2)
    expect_gt(nrow(r1), 0)
    expect_equal(r2$m, r1$m)
    expect_equal(r2$d, r1$d, tolerance = 1e-12)
    expect_equal(r2$bed_early, r1$bed_early, tolerance = 1e-8)
  }
})

test_that("scenario ID ignores interface context but not search settings", {
  args <- list(m_range = 5:34, d_step = 0.05, d_max = 5, criterion = "either", rank_by = "best",
               tol_early = 0, tol_late = 0, early_cap = NULL, late_cap = NULL,
               min_gain = 1, min_days_saved = 7)
  s1 <- build_scenario(p, std, std, args, context = list(preset = "calculator", regimen_evaluated = list(m = 20)))
  s2 <- build_scenario(p, std, std, args, context = list(preset = "custom", regimen_evaluated = list(m = 21)))
  expect_identical(scenario_id(s1), scenario_id(s2))
  args$min_gain <- 2
  expect_false(identical(scenario_id(s1), scenario_id(build_scenario(p, std, std, args, context = list()))))
})
