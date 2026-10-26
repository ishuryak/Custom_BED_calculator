# Version 3.2.1, review M02-M03 (2026-10-26 v3.2 review).

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "uncertainty.R"), local = TRUE)

test_that("M03: toxicity shares near the required share show one decimal, the count and the side", {
  a <- fmt_tox_share(3548 / 4000, 4000, 0.9); b <- fmt_tox_share(3611 / 4000, 4000, 0.9)
  expect_identical(a, "88.7% (3548 of 4000; below the required 90%)")
  expect_identical(b, "90.3% (3611 of 4000; meets the required 90%)")
  expect_false(identical(a, b))
  expect_identical(fmt_tox_share(0.9, 4000, 0.9), "90.0% (3600 of 4000; meets the required 90%)")
  expect_identical(fmt_tox_share(0.70, 4000, 0.9), "70%")
  expect_identical(fmt_tox_share(0.99, 4000, 0.9), ">95% (3960 of 4000)")
  expect_identical(fmt_tox_share(0.02, 4000, 0.6), "<5% (80 of 4000)")
})

test_that("M02: the unused DVH fraction-count input is gone and the help says what is used", {
  src <- paste(readLines(file.path(app_abs, "R", "mod_oar.R"), encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl("dvh_n", src, fixed = TRUE))
  expect_true(grepl("the plan's own number of fractions is not used", src, fixed = TRUE))
})
