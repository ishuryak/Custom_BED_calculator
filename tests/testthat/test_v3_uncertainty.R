# Version 3 uncertainty engine (R/uncertainty.R): the vectorized engine must
# agree with the scalar one, a single draw must reproduce the deterministic
# result, unevaluable draws never count as meeting a criterion, and the
# shipped draws must match their manifest.

app_abs <- normalizePath(testthat::test_path("..", ".."))
source(file.path(app_abs, "R", "bed_models.R"), local = TRUE)
source(file.path(app_abs, "R", "schedule.R"), local = TRUE)
source(file.path(app_abs, "R", "uncertainty.R"), local = TRUE)
source(file.path(app_abs, "R", "oar.R"), local = TRUE)
draws_dir <- file.path(app_abs, "draws")
ref <- list(m = 35, d = 2, T = 47)

rand_params <- function(n, seed) {
  set.seed(seed)
  data.frame(alpha = runif(n, 0.01, 0.5), ab = runif(n, 3, 25), lambda = runif(n, 0, 1.5),
             Tk = runif(n, 0, 45), g = runif(n, 0, 0.03), c_over_alpha = runif(n, 5, 600))
}

test_that("vectorized DI and DD agree with the scalar functions on 1000 random cases", {
  pr <- rand_params(1000, 1)
  set.seed(2)
  m <- sample(5:40, 1000, TRUE); d <- runif(1000, 1, 6); T <- overall_time(m) + sample(0:10, 1000, TRUE)
  for (i in seq_len(1000)) {
    dr <- pr[i, ]
    expect_equal(bed_di_draws(m[i], d[i], T[i], dr)[1, 1],
                 bed_di(m[i], d[i], T[i], dr$alpha, dr$ab, dr$lambda, dr$Tk, dr$g), tolerance = 1e-12)
    expect_equal(bed_dd_draws(m[i], d[i], T[i], dr)[1, 1],
                 bed_dd(m[i], d[i], T[i], dr$alpha, dr$ab, dr$lambda, dr$c_over_alpha, dr$g), tolerance = 1e-12)
  }
})

test_that("one draw reproduces the deterministic fitted-DD and DI results of compute_all() exactly", {
  p <- modifyList(default_params(), list(alpha = 0.03, ab = 10, lambda_di = 0.015, Tk = 30, g = 0.012,
                                         dd_mode = "fitted", alpha_dd = 0.028, lambda_dd = 0.47,
                                         c_over_alpha_dd = 390, g_dd = 0.017))
  regs <- data.frame(m = c(20, 17, 30), d = c(2.75, 3.05, 2), T = c(26, 23, 40))
  di <- data.frame(alpha = p$alpha, ab = p$ab, lambda = p$lambda_di, Tk = p$Tk, g = p$g)
  dd <- data.frame(alpha = p$alpha_dd, ab = p$ab, lambda = p$lambda_dd, c_over_alpha = p$c_over_alpha_dd, g = p$g_dd)
  det <- compute_all(regs$m, regs$d, regs$T, p, default_anchor())
  expect_identical(bed_di_draws(regs$m, regs$d, regs$T, di)[, 1], det$bed_di)
  expect_equal(bed_dd_draws(regs$m, regs$d, regs$T, dd)[, 1], det$bed_dd, tolerance = 1e-14)
  tx <- fixed_toxicity_draws(p, 1)
  expect_equal(bed_early_draws(regs$m, regs$d, regs$T, tx)[, 1], det$bed_early, tolerance = 1e-14)
  expect_equal(bed_late_draws(regs$m, regs$d, tx$ab_l)[, 1], det$bed_late)
})

test_that("unevaluable draws are counted and never meet a criterion", {
  # a draw with a huge repopulation loss makes the comparison's DI BED negative
  di <- data.frame(alpha = c(0.03, 0.001), ab = 10, lambda = c(0.015, 1), Tk = 20, g = 0.01)
  dd <- data.frame(alpha = 0.03, ab = 10, lambda = 0.5, c_over_alpha = 400, g = 0.01)[c(1, 1), ]
  regs <- data.frame(m = 20, d = 2.75, T = 26)
  td <- tumor_deltas(regs, ref, di, dd)
  expect_equal(td$valid_di[1, ], c(TRUE, FALSE))
  s <- summarize_draws(regs, td, threshold = -1000)
  expect_equal(s$n_invalid_di, 1)
  expect_equal(s$share_di, 0.5)        # the invalid draw does not pass even a -1000 Gy threshold
  expect_equal(s$share_both, 0.5)
})

test_that("shared draws: a regimen compared with itself has zero difference in every draw", {
  di <- load_draws("di", dir = draws_dir); dd <- load_draws("dd", dir = draws_dir)
  td <- tumor_deltas(data.frame(m = 35, d = 2, T = 47), ref, di, dd)
  expect_true(all(td$delta_di == 0, na.rm = TRUE))
  expect_true(all(td$delta_dd == 0, na.rm = TRUE))
})

test_that("draws files match the manifest, and a changed file is refused", {
  di <- load_draws("di", dir = draws_dir)
  expect_equal(nrow(di), 4000)
  man <- jsonlite::fromJSON(file.path(draws_dir, "draws_manifest.json"))
  expect_equal(attr(di, "sha256"), man$files$tumor_di_primary.csv$sha256)   # attributes survive loading
  # 3.1 (review M8): the Monte Carlo error counts distinct DI draws, not the fit's ESS
  expect_equal(attr(di, "ess"), man$files$tumor_di_primary.csv$n_distinct_di_draws)
  expect_true(all(c("alpha", "ab", "lambda", "Tk", "g") %in% names(di)))
  tmp <- file.path(tempdir(), "draws_tamper"); dir.create(tmp, showWarnings = FALSE)
  file.copy(list.files(draws_dir, full.names = TRUE), tmp, overwrite = TRUE)
  x <- readLines(file.path(tmp, "tumor_di_primary.csv")); x[2] <- sub("^[0-9.]+", "0.5", x[2])
  writeLines(x, file.path(tmp, "tumor_di_primary.csv"))
  expect_error(load_draws("di", dir = tmp), "does not match")
})

test_that("golden values: 50 shipped DI draws, HYPNO 55/20 vs 35 x 2 Gy (pinned 2026-10-26, draws 3.1)", {
  di <- load_draws("di", dir = draws_dir)[1:50, ]
  dd <- load_draws("dd", dir = draws_dir)[1:50, ]
  td <- tumor_deltas(data.frame(m = 20, d = 2.75, T = 26), ref, di, dd)
  g <- c(mean_di = mean(td$delta_di), sd_di = sd(td$delta_di[1, ]))
  # R 4.5.1 (Windows, ucrt); draws SHA-256 in draws/draws_manifest.json. Exact arithmetic, so 1e-10.
  # Repinned 2026-10-26: draws replaced by the 3.1 verified-trial refit (was 1.749, 5.260 on the 10-21 draws).
  expect_equal(unname(g), c(-2.28476945801911, 1.44023467546532), tolerance = 1e-10)
})

test_that("Monte Carlo: toxicity shares from two seeds agree within 4 standard errors", {
  regs <- data.frame(m = 20, d = 2.75, T = 26)
  lim <- list(early = list(mode = "reference", tol = 0), late = list(mode = "reference", tol = 0))
  s <- vapply(c(11, 12), function(sd) mean(toxicity_ok_draws(regs, ref, toxicity_draws(4000, sd), lim)$ok), 0)
  se <- sqrt(2) * share_mcse(mean(s), 4000)
  expect_lt(abs(diff(s)), 4 * se + 1e-12)
})

test_that("toxicity draws are reproducible and leave the caller's random stream untouched", {
  set.seed(5); a <- runif(1)
  set.seed(5); x <- toxicity_draws(100, 7); b <- runif(1)
  expect_equal(a, b)
  expect_identical(toxicity_draws(100, 7), x)
  expect_true(all(x$Tk_e >= 7 & x$Tk_e <= 19 & x$alpha_e >= 0.3 & x$alpha_e <= 0.4))
})

test_that("conditional TCP: zero difference returns the control rate; more BED raises TCP", {
  delta <- matrix(c(0, 5, -5), 3, 1)
  t <- tcp_given_control(delta, 0.03, 0.5)
  expect_equal(t[1, 1], 0.5)
  expect_gt(t[2, 1], 0.5); expect_lt(t[3, 1], 0.5)
  expect_error(tcp_given_control(delta, 0.03, 1), "tcp0 < 1")
})

test_that("display rounding: shares to 5 points with open ends, BEDs to 0.5 Gy", {
  expect_equal(fmt_share(c(0.01, 0.049, 0.05, 0.512, 0.951, 1, NA)),
               c("<5%", "<5%", "5%", "50%", ">95%", ">95%", "not available"))
  expect_equal(fmt_half(c(1.24, 1.26, -0.74)), c("+1.0", "+1.5", "-0.5"))
})

test_that("the dose-per-fraction support of the tumor data is read from the data file", {
  s <- tumor_data_support(file.path(app_abs, "data", "shuryak2018_table1.csv"))
  expect_equal(c(s$d_min, s$d_max), c(1.1, 2.5))
})

test_that("the uncertainty search lists one regimen per fraction number by the stated rule and flags extrapolation", {
  old <- setwd(app_abs); on.exit(setwd(old))
  di <- load_draws("di")[1:400, ]; dd <- load_draws("dd")[1:400, ]
  tx <- toxicity_draws(400)
  lim <- list(early = list(mode = "hypno", tol = 0), late = list(mode = "reference", tol = 0))
  s <- uncertain_search(ref, di, dd, tx, lim, threshold = -5, m_range = 15:30, d_step = 0.1, d_max = 3.5,
                        tox_share_min = 0.5)
  expect_true(nrow(s) > 0)
  expect_false(anyDuplicated(s$m) > 0)
  grid <- attr(s, "grid")
  for (i in seq_len(nrow(s))) {          # the listed dose has the best score among eligible ones
    g <- grid[grid$m == s$m[i] & grid$eligible, ]
    expect_equal(s$score[i], max(g$score))
  }
  expect_true(all(s$outside_data == (s$d > 2.5)))
  expect_error(uncertain_search(ref, di, dd, tx, list(early = list(mode = "fowler"), late = lim$late)),
               "must be 'hypno' or 'reference'")
})

test_that("OAR draw checker: a cord dose within the constraint at every alpha/beta passes every draw", {
  lib <- load_constraints(file.path(app_abs, "data", "oar_constraints.csv"))
  chk <- oar_draw_checker(lib, list(cord = 50))
  tx <- toxicity_draws(200)
  ok <- chk(data.frame(m = 20, d = 2.75, T = 26), tx)        # cord 27.5 Gy: far below
  expect_true(all(ok))
  ok2 <- chk(data.frame(m = 20, d = 2.75, T = 26), tx[1, ])
  expect_equal(dim(ok2), c(1L, 1L))
  hi <- oar_draw_checker(lib, list(cord = 180))(data.frame(m = 20, d = 2.75, T = 26), tx)  # 99 Gy
  expect_false(any(hi))
  expect_null(oar_draw_checker(lib, list(cord = NA)))
})

test_that("the refit numbers quoted on the About page match the shipped draws", {
  # 3.1: alpha/beta varies by draw, so lambda/alpha and g/alpha are in Gy at each draw's alpha/beta, not Gy10
  di <- load_draws("di", dir = draws_dir)
  q <- function(x) stats::quantile(x, c(.025, .5, .975), names = FALSE)
  la <- q(di$lambda / di$alpha); ga <- q(di$g / di$alpha); tk <- q(di$Tk)
  ab <- paste(readLines(file.path(app_abs, "about.md"), encoding = "UTF-8"), collapse = "\n")
  expect_true(grepl(sprintf("%.2f Gy per day (95%% interval %.2f–%.2f)", la[2], la[1], la[3]), ab, fixed = TRUE))
  expect_true(grepl(sprintf("onset near day %.0f (%.0f–%.0f)", tk[2], tk[1], tk[3]), ab, fixed = TRUE))
  expect_true(grepl(sprintf("%.2f Gy per day (%.2f–%.2f) throughout", ga[2], ga[1], ga[3]), ab, fixed = TRUE))
})
