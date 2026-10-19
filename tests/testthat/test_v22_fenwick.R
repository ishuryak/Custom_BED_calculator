# Fenwick et al. 2008 acute boundary (optional early-toxicity limit).
# Expected values from the published formula, computed independently:
# B(T) = 69.5 (T/32.2) / sin(T/32.2) - 3.5 Gy10, T inclusive days, T <= 50.

models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)
p <- legacy_params()   # version 1 tumor values (pinned results)
std <- list(m = 35, d = 2, T = 47)

test_that("Fenwick boundary reproduces the published formula and its domain", {
  expect_equal(fenwick_limit(c(47, 26, 19)), c(98.57419, 74.16888, 70.20307), tolerance = 1e-6)
  expect_true(is.na(fenwick_limit(51)))                  # outside the source domain
  expect_equal(fenwick_cell_kill(35, 2), 84)              # alpha/beta 10, complete repair
  expect_lt(fenwick_cell_kill(20, 2.75), fenwick_limit(26))   # HYPNO 55/20: 70.13 below 74.17
  expect_gt(fenwick_cell_kill(20, 2.9), fenwick_limit(26))    # 20 x 2.9 Gy: 74.8 above
})

test_that("the Fenwick option replaces the early cap in the search", {
  s <- suggest_regimens(p, std, anchor = std, fenwick = TRUE, min_gain = -30,
                        tol_late = 30, d_max = 4)
  expect_gt(nrow(s), 0)
  expect_true(all(fenwick_cell_kill(s$m, s$d) <= fenwick_limit(s$T) + 1e-9))
  f <- attr(s, "feasible")
  expect_true(all(fenwick_cell_kill(f$m, f$d) <= fenwick_limit(f$T) + 1e-9))
  # the next dose step of a row limited by early toxicity crosses the boundary
  lim <- s[grepl("early toxicity", s$next_step), ]
  expect_gt(nrow(lim), 0)
  expect_true(all(fenwick_cell_kill(lim$m, lim$d + 0.05) > fenwick_limit(lim$T)))
})
