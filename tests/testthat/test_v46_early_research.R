# Version 4.6: research mode for early toxicity (the trial-count fit), switchable on the main page.
# Fowler 2003 stays the default. docs/development/2026_11_11_early_toxicity_review_response.md
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
fit <- load_early_fit(file.path(app_abs, "draws"))
p_fowler <- modifyList(default_params(), startup_profile(refit_profile(file.path(app_abs, "draws")))$values)
p_fit <- fitted_early_params(p_fowler, fit)
std <- V4_STANDARD

test_that("the fitted model loads, with its point values, and refuses a changed draws file", {
  expect_true(fit$ok)
  expect_equal(nrow(fit$draws), 4000)
  expect_equal(fit$point$ab_e, unname(round(median(fit$draws$ab_e), 4)), tolerance = 1e-3)
  expect_equal(fit$point$Tk_e, 7)
  expect_equal(fit$point$h_e, 3.2)
  tmp <- tempfile(); dir.create(tmp)
  file.copy(file.path(app_abs, "draws", c(EARLY_FIT_MANIFEST, EARLY_FIT_FILE)), tmp)
  cat("1,1,1,1,1\n", file = file.path(tmp, EARLY_FIT_FILE), append = TRUE)
  bad <- load_early_fit(tmp)
  expect_false(bad$ok)
  expect_match(bad$reason, "do not match their manifest")
  expect_false(load_early_fit(tempdir())$ok)
})

test_that("the R evaluator agrees with the Python evaluator used for fitting", {
  fx <- jsonlite::fromJSON(file.path(app_abs, "tests", "testthat", "fixtures", "early_fit_evaluator.json"),
                           simplifyVector = FALSE)
  for (cs in fx$cases) {
    cal <- list(fx_per_week = cs$calendar$fx_per_week, start_weekday = cs$calendar$start_weekday,
                skip_days = as.integer(unlist(cs$calendar$skip_days)))
    expect_equal(as.integer(treatment_days(cs$m, cal$fx_per_week, cal$start_weekday, cal$skip_days)),
                 as.integer(unlist(cs$days)), info = cs$name)
    got <- bed_early_fit(cs$m, cs$d, cs$T, cs$ab, cs$r, fx$Tk, fx$halftime_h, cal)
    expect_equal(got, cs$E, tolerance = 1e-9, info = paste(cs$name, cs$ab, cs$r))
  }
})

test_that("Fowler mode is unchanged: 69 doses in 18 fraction counts", {
  r <- v4_regimens(p_fowler)$rows
  expect_equal(nrow(r), 69)
  expect_equal(length(unique(r$m)), 18)
})

test_that("research mode lists 9 doses in 26 to 30 fractions with the +5 margin", {
  r <- v4_regimens(p_fit)$rows
  expect_equal(nrow(r), 9)
  expect_equal(sort(unique(r$m)), 26:30)
  st <- compute_all(std$m, std$d, std$T, p_fit, default_anchor())
  expect_true(all(r$bed_early <= st$bed_early + 5 + 1e-9))
})

test_that("same regimen gives a zero difference and full support in every draw, even where E <= 0", {
  tx <- early_fit_toxicity_draws(4000, fit)
  e <- bed_early_draws(std$m, std$d, std$T, tx)[1, ]
  expect_gt(sum(e <= 0), 100)                          # about 13% of draws: the legacy filter would fail these
  regs <- data.frame(m = std$m, d = std$d, T = std$T)
  ok <- toxicity_ok_draws(regs, std, tx, list(early = list(mode = "reference", tol = 0),
                                              late = list(mode = "reference", tol = 0)))
  expect_true(all(ok$early_ok))
})

test_that("swapping candidate and reference negates the difference in every draw", {
  tx <- early_fit_toxicity_draws(4000, fit)
  a <- bed_early_draws(20, 2.75, overall_time(20), tx)[1, ]
  b <- bed_early_draws(35, 2, 47, tx)[1, ]
  expect_equal(a - b, -(b - a))
  expect_equal(median(a - b), 18.5, tolerance = 0.05)    # review F9 fixture: median 18.51 for 55/20 vs 70/35
})

test_that("research mode drops the positivity requirement but Fowler mode keeps it", {
  hot <- modifyList(p_fit, list(r_e = 3))               # E(70/35) = 70 + Q/ab - 3 x 39 < 0
  expect_lt(compute_all(std$m, std$d, std$T, hot, default_anchor())$bed_early, 0)
  s <- suggest_regimens(hot, std, m_range = 25:30, tol_early = 5, min_gain = -5)
  expect_false(identical(attr(s, "reason"), "negative_reference_early"))
  cold <- modifyList(p_fowler, list(lambda_e = 3 * p_fowler$alpha_e))
  s2 <- suggest_regimens(cold, std, m_range = 25:30, tol_early = 5, min_gain = -5)
  expect_identical(attr(s2, "reason"), "negative_reference_early")
})

test_that("fixed toxicity in research mode reproduces the point model", {
  tx <- fixed_toxicity_draws(p_fit, 3)
  e <- bed_early_draws(20, 2.75, overall_time(20), tx)[1, ]
  expect_equal(e, rep(compute_all(20, 2.75, overall_time(20), p_fit, default_anchor())$bed_early, 3), tolerance = 1e-9)
})
