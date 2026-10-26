# Version 3.2.1 golden values: the review of 2026-10-26 (section 4.3, 4.4) recomputed these
# independently. A change to any of them is a change of results and must be stated as one.
# Deterministic: D18 defaults, calibrated DD, 18 x 2.9 Gy in 24 days against 35 x 2 Gy in 47 days.
# Uncertainty: the shipped primary draws (draw files unchanged since 3.2.0; tolerance 1e-4 Gy).

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "bed_models.R"), local = TRUE)
source(file.path(app_abs, "R", "uncertainty.R"), local = TRUE)

test_that("golden: deterministic default values (review section 4.3)", {
  p <- default_params(); an <- list(m = 35, d = 2, T = 47)
  a <- compute_all(c(18, 35), c(2.9, 2), c(24, 47), p, an)
  expect_equal(a$bed_di, c(67.337999, 71.514286), tolerance = 1e-6)
  expect_equal(a$bed_dd, c(60.885876, 71.514286), tolerance = 1e-6)
  expect_equal(a$bed_early, c(54.663309, 53.105440), tolerance = 1e-6)
  expect_equal(a$bed_late, c(102.66, 116.666667), tolerance = 1e-6)
})

test_that("golden: uncertainty medians and shares on the shipped primary draws (review section 4.4)", {
  old <- setwd(app_abs); on.exit(setwd(old))
  di <- load_draws("di", "primary"); dd <- load_draws("dd", "primary")
  td <- tumor_deltas(list(m = 18, d = 2.9, T = 24), list(m = 35, d = 2, T = 47), di, dd)
  x_di <- as.numeric(td$delta_di[1, ]); x_dd <- as.numeric(td$delta_dd[1, ])
  expect_equal(median(x_di), -4.1513, tolerance = 1e-4)
  expect_equal(median(x_dd), -5.0844, tolerance = 1e-4)
  expect_equal(mean(x_di >= 1), 0.03); expect_equal(mean(x_dd >= 1), 0.00825)
  expect_equal(mean(x_di >= 1 & x_dd >= 1), 0.0075)
})
