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

# Review R5 (v4.6): the runtime package is validated as a whole, and the exact posterior is pinned.
pkg_copy <- function(edit = NULL, draws_edit = NULL) {
  tmp <- tempfile(); dir.create(tmp)
  file.copy(file.path(app_abs, "draws", c(EARLY_FIT_MANIFEST, EARLY_FIT_FILE, "early_fit_trials.csv")), tmp)
  m <- jsonlite::fromJSON(file.path(tmp, EARLY_FIT_MANIFEST), simplifyVector = FALSE)
  if (!is.null(draws_edit)) {
    d <- utils::read.csv(file.path(tmp, EARLY_FIT_FILE))
    utils::write.csv(draws_edit(d), file.path(tmp, EARLY_FIT_FILE), row.names = FALSE)
    m$files[[EARLY_FIT_FILE]]$sha256 <- digest::digest(file = file.path(tmp, EARLY_FIT_FILE), algo = "sha256")
  }
  if (!is.null(edit)) m <- edit(m)
  jsonlite::write_json(m, file.path(tmp, EARLY_FIT_MANIFEST), auto_unbox = TRUE, digits = NA, pretty = TRUE)
  tmp
}

test_that("a package with failed diagnostics, wrong point values or unknown definitions is unavailable", {
  expect_true(load_early_fit(pkg_copy())$ok)
  cases <- list(
    list(function(m) { m$gate$passed <- FALSE; m }, "diagnostics do not meet the gate"),
    list(function(m) { m$gate$divergences <- 99; m }, "diagnostics do not meet the gate"),
    list(function(m) { m$point_values$r_e <- 100; m }, "not the median of the draws"),
    list(function(m) { m$evaluator <- "other-evaluator"; m }, "unsupported evaluator"),
    list(function(m) { m$schema <- NULL; m }, "unknown package schema"),
    list(function(m) { m$fixed$repair_halftime_h <- 2; m }, "unsupported evaluator settings"))
  for (cs in cases) {
    f <- load_early_fit(pkg_copy(cs[[1]]))
    expect_false(f$ok); expect_match(f$reason, cs[[2]], fixed = TRUE)
  }
  bad <- pkg_copy(); writeLines("{not json", file.path(bad, EARLY_FIT_MANIFEST))
  expect_false(load_early_fit(bad)$ok)                      # no uncaught parse error
})

test_that("reordered draws with a matching hash load, but an old exact scenario refuses to replay", {
  rev_dir <- pkg_copy(draws_edit = function(d) { d <- d[rev(seq_len(nrow(d))), ]; d$draw <- seq_len(nrow(d)); d })
  f2 <- load_early_fit(rev_dir)
  expect_true(f2$ok)                                         # same medians, a different exact posterior
  expect_false(identical(f2$draws_sha256, fit$draws_sha256))
  eb <- list(model = "fitted", model_id = fit$model_id, ab_e = fit$point$ab_e, r_e = fit$point$r_e,
             Tk_e = fit$point$Tk_e, h_e = fit$point$h_e, draws_sha256 = fit$draws_sha256,
             package_sha256 = fit$package_sha256, evaluator = EARLY_FIT_EVALUATOR,
             calendar = list(fx_per_week = 5, start_weekday = 1, skip_days = list()))
  for (f in c("oar.R", "export.R")) source(file.path(app_abs, "R", f), local = TRUE)
  expect_true(early_is_relative(apply_early_block(list(), eb, fit)))
  expect_error(apply_early_block(list(), eb, f2), "not the installed one")
})

# Review R3 (v4.6): each regimen's repair term uses its own treatment dates; no inferred calendar.
test_that("custom calendars use their own dates, matching the independent Python evaluator", {
  fx <- jsonlite::fromJSON(file.path(app_abs, "tests", "testthat", "fixtures", "early_fit_evaluator.json"), simplifyVector = FALSE)
  pick <- function(nm) Filter(function(cs) cs$name == nm && abs(cs$ab - fit$point$ab_e) < 1e-9, fx$cases)[[1]]
  a <- pick("20x2.75 skip 4,5,8"); b <- pick("33x2 6/wk Mon")
  cal_a <- list(fx_per_week = 5, start_weekday = 1, skip_days = c(4L, 5L, 8L))
  cal_b <- list(fx_per_week = 6, start_weekday = 1, skip_days = integer())
  ea <- compute_all(20, 2.75, 31, p_fit, default_anchor(), cal_a)$bed_early
  eb <- compute_all(33, 2, 38, p_fit, default_anchor(), cal_b)$bed_early
  # the review's numbers use the app's point values (4 decimals); the fixture uses the unrounded medians
  expect_equal(ea, 31.11404545, tolerance = 1e-6); expect_equal(eb, 29.73698842, tolerance = 1e-6)
  expect_equal(bed_early_fit(20, 2.75, 31, a$ab, a$r, fx$Tk, fx$halftime_h, cal_a), a$E, tolerance = 1e-6)
  expect_equal(bed_early_fit(33, 2, 38, b$ab, b$r, fx$Tk, fx$halftime_h, cal_b), b$E, tolerance = 1e-6)
  # the same m and T on different calendars give different repair terms; never silently the default
  expect_false(isTRUE(all.equal(ea, compute_all(20, 2.75, 31, modifyList(p_fit, list(calendar_e = NULL)), default_anchor(),
                                                list(fx_per_week = 5, start_weekday = 3, skip_days = c(4L, 5L)))$bed_early)))
})

test_that("an overall time that the calendar cannot produce gives no repair-aware result", {
  expect_true(is.na(compute_all(20, 2.75, 40, p_fit, default_anchor(), default_calendar())$bed_early))
  expect_false(is.na(compute_all(20, 2.75, 40, p_fowler, default_anchor())$bed_early))     # Fowler mode unaffected
})

test_that("search candidates and the comparator use their own calendars; equal regimens still give zero", {
  ref6 <- list(m = 33, d = 2, T = 38, calendar = list(fx_per_week = 6, start_weekday = 1, skip_days = integer()))
  s <- suggest_regimens(p_fit, ref6, m_range = 15:30, tol_early = 5, min_gain = -10)
  r0 <- compute_all(33, 2, 38, p_fit, default_anchor(), ref6$calendar)$bed_early
  f <- attr(s, "feasible")
  expect_true(nrow(f) > 0)
  expect_equal(f$d_early, f$bed_early - r0, tolerance = 1e-9)
  tx <- early_fit_toxicity_draws(4000, fit)
  e1 <- bed_early_draws(33, 2, 38, tx, ref6$calendar)[1, ]
  regs <- data.frame(m = 33, d = 2, T = 38); attr(regs, "calendar") <- ref6$calendar
  ok <- toxicity_ok_draws(regs, ref6, tx, list(early = list(mode = "reference", tol = 0), late = list(mode = "reference", tol = 0)))
  expect_true(all(ok$early_ok))                                              # equal schedule: support 1 in every draw
  e2 <- bed_early_draws(33, 2, 38, tx, default_calendar())[1, ]
  expect_true(all(is.na(e2)))                                                # wrong calendar for T: unavailable, not guessed
})
