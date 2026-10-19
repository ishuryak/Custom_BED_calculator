# Two-model agreement (version 2.4). Expected values from an independent
# Python implementation of the D18 models (see FIXES.md, version 2.4).

models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)
p <- default_params()                       # scenario D18
an <- list(m = 35, d = 2, T = 47)
hypno <- list(m = 20, d = 2.75, T = 26)

test_that("against HYPNO, both models support 17-19 fractions within -3 Gy10", {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  a <- agreement_regimens(p, hypno, an, margin = -3, early_cap = cap,
                          min_days_saved = 1, d_max = 5, m_range = 5:19)
  expect_equal(a$m, c(17, 18, 19))
  expect_equal(a$d, c(3.05, 2.90, 2.80), tolerance = 1e-12)
  expect_equal(a$T, c(23, 24, 25))
  expect_equal(a$d_di, c(-2.460208, -2.786459, -2.028485), tolerance = 1e-5)
  expect_equal(a$d_dd, c(-1.589516, -1.961425, -1.477145), tolerance = 1e-5)
  expect_equal(a$dd_onset, c(14.2804, 14.9734, 15.4237), tolerance = 1e-3)
  expect_true(all(a$before_di_onset))                # all end before DI day 28
  expect_equal(a$n_scenarios, rep(15, 3))
  expect_equal(a$n_pass_both, rep(15, 3))
  expect_true(all(a$d_di_min <= a$d_di & a$d_di <= a$d_di_max))
  expect_true(all(a$d_dd_min <= a$d_dd & a$d_dd <= a$d_dd_max))
  expect_true(all(a$d_di >= -3 - 1e-9 & a$d_dd >= -3 - 1e-9))
})

test_that("no regimen is supported by both models without a loss", {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  a <- agreement_regimens(p, hypno, an, margin = 0, early_cap = cap,
                          min_days_saved = 1, d_max = 5, m_range = 5:19)
  expect_equal(nrow(a), 0)
})

test_that("scenario count drops when the margin is tighter than some scenarios allow", {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  a <- agreement_regimens(p, hypno, an, margin = -2.5, early_cap = cap,
                          min_days_saved = 1, d_max = 5, m_range = 5:19)
  r17 <- a[a$m == 17, ]
  expect_equal(nrow(r17), 1)                          # D18: -2.46 / -1.59 passes
  expect_lt(r17$n_pass_both, 15)                    # but the Withers rule at alpha/beta 7 gives -2.69
})

test_that("against 70/35, scenario counts need both models and long courses pass the DI onset", {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  a <- agreement_regimens(p, an, an, margin = -3, early_cap = cap,
                          min_days_saved = 7, d_max = 5, m_range = 5:34)
  expect_equal(a$m, c(28, 29, 30))
  expect_false(any(a$before_di_onset))                # 38-40 days, DI onset day 28
  di_only <- 0
  for (i in seq_len(nrow(a))) {
    st <- sensitivity_table(list(m = a$m[i], d = a$d[i], T = a$T[i]), an, an, p)
    expect_equal(a$n_pass_both[i], sum(st$d_di >= -3 - 1e-9 & st$d_dd >= -3 - 1e-9, na.rm = TRUE))
    di_only <- di_only + sum(st$d_di >= -3 & !(st$d_dd >= -3), na.rm = TRUE)
  }
  expect_gt(di_only, 0)                               # scenarios where only DI meets the margin exist
})
