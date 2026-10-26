# Version 3.2.1, review F06-F07 (2026-10-26 v3.2 review): the guide and About describe what the
# code does. The Both-model screen lists a regimen at the current settings; 15/15 is not required.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "bed_models.R"), local = TRUE)
guide <- paste(readLines(file.path(app_abs, "docs", "user_guide", "HNSCC_BED_v321_user_guide.md"), encoding = "UTF-8"), collapse = "\n")
about <- paste(readLines(file.path(app_abs, "about.md"), encoding = "UTF-8"), collapse = "\n")

test_that("F06: a regimen passing 12 of 15 settings is listed by the Both-model screen (review example)", {
  p <- default_params(); an <- list(m = 35, d = 2, T = 47); hy <- list(m = 20, d = 2.75, T = 26)
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  a <- agreement_regimens(p, hy, an, margin = -2.5, early_cap = cap, min_days_saved = 1, d_max = 5, m_range = 5:19)
  r <- a[a$m == 17, ]
  expect_equal(nrow(r), 1)
  expect_equal(r$d, 3.05); expect_equal(r$T, 23)
  expect_equal(r$d_di, -2.460208, tolerance = 1e-6); expect_equal(r$d_dd, -1.589516, tolerance = 1e-6)
  expect_equal(r$n_scenarios, 15); expect_equal(r$n_pass_both, 12)
})

test_that("F06: the guide no longer says the screen requires every one of the 15 settings", {
  expect_false(grepl("in each of 15 parameter settings", guide, fixed = TRUE))
  expect_true(grepl("Inclusion does not require 15 of 15", guide, fixed = TRUE))
  expect_true(grepl("not a joint clinical posterior", guide, fixed = TRUE))
})

test_that("F07: the repair count matches repair_ess.json (752 = 748 bulk + 4 tail only)", {
  recs <- Sys.glob(file.path(app_abs, "fit", "output_refit_2026_10_25c_*", "dd_calibrated_*", "repair_ess.json"))
  n <- sum(vapply(recs, function(f) length(jsonlite::fromJSON(f, simplifyVector = FALSE)$repairs), integer(1)))
  expect_equal(n, 752)
  for (txt in list(guide, about)) {
    expect_true(grepl("752 conditional DD fits", txt, fixed = TRUE))
    expect_false(grepl("the 748 conditional DD fits", txt, fixed = TRUE))
  }
  expect_true(grepl("2 Gy twice daily in the Jackson trial", about, fixed = TRUE))
  expect_false(grepl("up to 0.25 Gy", about, fixed = TRUE))
})
