# Tests for R/bed_models.R. Run from the app folder:
#   Rscript -e "testthat::test_dir('tests/testthat')"
# BED_MODELS_FILE can point at a deliberately broken copy (mutation check).

models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)

p <- legacy_params()   # version 1 tumor values (pinned results)
ref <- default_ref()

# ---------------------------------------------------------------------------
# Independent reference implementation ("oracle"), written separately from
# R/bed_models.R: day-by-day calendar, numerical root finding for the DD
# calibration, and brute-force search with explicit loops.
# ---------------------------------------------------------------------------
o_smax <- function(x) ifelse(3 * x > 30, x, log(1 + exp(3 * x)) / 3)
o_lq <- function(n, d, ab) n * d + n * d^2 / ab
o_calendar <- function(n, days_per_week = 5) {
  count <- 0; day <- 0
  while (count < n) {
    day <- day + 1
    if (((day - 1) %% 7) + 1 <= days_per_week) count <- count + 1
  }
  day
}
o_di <- function(n, d, T, q) o_lq(n, d, q$ab) - q$lambda_di / q$alpha * o_smax(T - q$Tk)
o_dd_given <- function(n, d, T, q, Ca) {
  K <- o_lq(n, d, q$ab)
  K - q$lambda_dd / q$alpha * (1 - exp(-q$alpha * K / T)) * o_smax(T - Ca * T / K)
}
o_calibrate <- function(q, rf) {
  target <- o_di(rf$m, rf$d, rf$T, q)
  f <- function(Ca) o_dd_given(rf$m, rf$d, rf$T, q, Ca) - target
  K0 <- o_lq(rf$m, rf$d, q$ab)
  lo <- 1e-9 * K0
  hi <- 1e3 * K0
  if (f(lo) > 0) return(NA_real_)  # even the earliest onset loses too little
  uniroot(f, c(lo, hi), tol = 1e-12)$root
}
# Acute BED on Fowler's clock: elapsed days = inclusive days - 1.
o_early <- function(n, d, T, q) o_lq(n, d, q$ab_e) - q$lambda_e / q$alpha_e * o_smax((T - 1) - q$Tk_e)
o_late <- function(n, d, q) o_lq(n, d, q$ab_l)

# Brute-force search with explicit loops, same rules as documented.
# `anchor` is the regimen at which DD is calibrated (fixed, version 2.1).
o_search <- function(q, rf, m_range, d_step, d_max, criterion, rank_by,
                     tol_early = 0, tol_late = 0, early_cap = NULL,
                     min_gain = 0, min_days_saved = 1, anchor = rf) {
  Ca <- o_calibrate(q, anchor)
  dd <- function(n, d, T) if (is.na(Ca)) NA_real_ else o_dd_given(n, d, T, q, Ca)
  di0 <- o_di(rf$m, rf$d, rf$T, q); dd0 <- dd(rf$m, rf$d, rf$T)
  cap_e <- if (is.null(early_cap)) o_early(rf$m, rf$d, rf$T, q) + tol_early else early_cap
  cap_l <- o_late(rf$m, rf$d, q) + tol_late
  rows <- list()
  for (n in m_range) {
    T <- o_calendar(n)
    if (T > rf$T - min_days_saved) next
    best <- NULL
    k <- 1
    repeat {
      d <- round(k * d_step, 6); k <- k + 1
      if (d > d_max + 1e-9) break
      if (d <= rf$d + 1e-9) next
      gdi <- o_di(n, d, T, q) - di0
      gdd <- dd(n, d, T) - dd0
      di_ok <- gdi >= min_gain - 1e-9
      dd_ok <- !is.na(gdd) && gdd >= min_gain - 1e-9
      tum <- switch(criterion, either = di_ok || dd_ok, both = di_ok && dd_ok,
                    di = di_ok, dd = dd_ok)
      if (!tum) next
      if (o_early(n, d, T, q) > cap_e + 1e-9) next
      if (o_late(n, d, q) > cap_l + 1e-9) next
      sc <- switch(rank_by, best = max(c(gdi, gdd), na.rm = TRUE),
                   worst = min(c(gdi, gdd), na.rm = TRUE), di = gdi, dd = gdd)
      if (is.null(best) || sc > best$score + 1e-12) best <- list(m = n, d = d, T = T, score = sc)
    }
    if (!is.null(best)) rows[[length(rows) + 1]] <- as.data.frame(best)
  }
  if (length(rows) == 0) return(data.frame(m = integer(), d = numeric(), T = numeric(), score = numeric()))
  out <- do.call(rbind, rows)
  out[order(-out$score), ]
}

# The v1 website's DD formula, copied verbatim from v1 app.R (deployed copy).
v1_dd <- function(m, d, T, alpha, r, L_DI, L_DD, Tk) {
  (r*L_DD*(exp(-m*alpha*d*(d + r)/(T*r)) - 1)*log(1 + exp(3*(47*m*d^2 + 47*d*m*r
    + 28*r*log((1 + exp(141 - 3*Tk))^(-L_DI/((exp(-(84*alpha)/47) - 1)*L_DD)) - 1)
    - 3948*r)*T/(47*m*d*(d + r)))) + 3*m*alpha*d*(d + r))/(3*r*alpha)
}

# ---------------------------------------------------------------------------

test_that("golden values: talk table (HNSCC_modeling_update.pptx slide 6)", {
  tab <- data.frame(
    d = c(2, 2.75, 1.8, 4, 2.2, 2.5, 2.9, 3.5, 4, 4, 4.3, 4.5, 4.6),
    m = c(35, 20, 40, 14, 30, 25, 20, 15, 12, 10, 10, 10, 10),
    T = c(47, 26, 26, 40, 40, 33, 26, 19, 16, 12, 12, 12, 12),
    di = c(65.00, 70.12, 84.96, 66.40, 68.52, 73.12, 74.82, 70.87, 67.20, 56.00, 61.49, 65.25, 67.16),
    dd = c(65.00, 58.02, 68.02, 62.02, 63.35, 62.46, 61.13, 59.57, 57.59, 50.41, 54.39, 57.19, 58.63),
    tox3 = c(116.67, 105.42, 115.20, 130.67, 114.40, 114.58, 114.07, 113.75, 112.00, 93.33, 104.63, 112.50, 116.53),
    tox2 = c(140.00, 130.63, 136.80, 168.00, 138.60, 140.63, 142.10, 144.38, 144.00, 120.00, 135.45, 146.25, 151.80))
  a <- compute_all(tab$m, tab$d, tab$T, p, ref)
  expect_equal(round(a$bed_di, 2), tab$di)
  expect_equal(round(a$bed_dd, 2), tab$dd)
  expect_equal(round(a$bed_late, 2), tab$tox3)
  # the talk rounded x.xx5 up; compare within half a cent
  expect_equal(bed_late(tab$m, tab$d, 2), tab$tox2, tolerance = 0.0051, scale = 1)
})

test_that("DD equals v1 website at tumor alpha/beta = 10 wherever DD is calibratable", {
  g <- expand.grid(alpha = c(0.1, 0.2, 0.4), L_DI = c(0.1, 0.2, 0.4),
                   L_DD = c(0.3, 0.5, 0.9), Tk = c(21, 28, 35))
  n_ok <- 0
  for (i in seq_len(nrow(g))) {
    q <- p; q$alpha <- g$alpha[i]; q$lambda_di <- g$L_DI[i]
    q$lambda_dd <- g$L_DD[i]; q$Tk <- g$Tk[i]
    if (dd_calibration(ref, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)$status != "ok") next
    n_ok <- n_ok + 1
    a <- compute_all(c(20, 30, 35), c(2.75, 2.2, 2), c(26, 40, 47), q, ref)
    expect_equal(a$bed_dd, v1_dd(c(20, 30, 35), c(2.75, 2.2, 2), c(26, 40, 47),
                                 q$alpha, 10, q$lambda_di, q$lambda_dd, q$Tk),
                 tolerance = 1e-8)
  }
  expect_gt(n_ok, 30)
})

test_that("all BEDs match the independent implementation over a parameter grid", {
  set.seed(1)
  for (k in 1:40) {
    q <- p
    q$alpha <- runif(1, 0.1, 0.5); q$ab <- runif(1, 5, 15)
    q$lambda_di <- runif(1, 0.05, 0.4); q$lambda_dd <- runif(1, 0.3, 2)
    q$Tk <- runif(1, 18, 35); q$ab_e <- runif(1, 7, 12); q$alpha_e <- runif(1, 0.25, 0.4)
    q$lambda_e <- log(2) / runif(1, 2, 4); q$Tk_e <- runif(1, 5, 12); q$ab_l <- runif(1, 1.5, 4)
    rf <- list(m = sample(30:40, 1), d = runif(1, 1.8, 2.2))
    rf$T <- o_calendar(rf$m)
    n <- sample(5:34, 5); d <- runif(5, 2, 4); T <- sapply(n, o_calendar)
    a <- compute_all(n, d, T, q, rf)
    expect_equal(a$bed_di, o_di(n, d, T, q), tolerance = 1e-10)
    expect_equal(a$bed_early, o_early(n, d, T, q), tolerance = 1e-10)
    expect_equal(a$eqd2_early, o_early(n, d, T, q) / (1 + 2 / q$ab_e), tolerance = 1e-10)
    expect_equal(a$bed_late, o_late(n, d, q), tolerance = 1e-10)
    expect_equal(a$eqd2_late, o_late(n, d, q) / (1 + 2 / q$ab_l), tolerance = 1e-10)
    Ca <- o_calibrate(q, rf)
    if (is.na(Ca)) {
      expect_true(all(is.na(a$bed_dd)))
    } else {
      expect_equal(a$bed_dd, o_dd_given(n, d, T, q, Ca), tolerance = 1e-7)
    }
  }
})

test_that("DD is calibrated to DI at the reference for every tumor alpha/beta (v1 is not)", {
  n_ok <- 0
  for (ab in c(3, 5, 8, 10, 15, 25)) for (Tk in c(21, 28, 35)) for (alpha in c(0.1, 0.3)) {
    q <- p; q$ab <- ab; q$Tk <- Tk; q$alpha <- alpha
    a <- compute_all(ref$m, ref$d, ref$T, q, ref)
    st <- dd_calibration(ref, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)$status
    if (st == "ok") {
      n_ok <- n_ok + 1
      expect_equal(a$bed_dd, a$bed_di, tolerance = 1e-9)
    } else {
      expect_true(is.na(a$bed_dd))                # flagged, never a wrong number
      expect_true(is.na(o_calibrate(q, ref)))     # and the oracle finds no root either
    }
  }
  expect_gt(n_ok, 20)
  r2 <- list(m = 33, d = 2, T = 45)
  a <- compute_all(r2$m, r2$d, r2$T, p, r2)
  expect_equal(a$bed_dd, a$bed_di, tolerance = 1e-9)
  # the v1 formula is wrong at alpha/beta = 5 (documented bug)
  expect_gt(abs(v1_dd(35, 2, 47, 0.2, 5, 0.2, 0.5, 28) - (70 * 1.4 - 19)), 1)
})

test_that("overall time matches a day-by-day Mon-Fri calendar", {
  for (n in 1:80) expect_equal(overall_time(n), o_calendar(n))
  expect_equal(overall_time(c(35, 30, 25, 20, 15, 12, 10, 5, 1)),
               c(47, 40, 33, 26, 19, 16, 12, 5, 1))
})

test_that("early and late toxicity BEDs match hand calculations", {
  hand_early <- 70 * (1 + 2 / 10) - (log(2) / 2.5) * (46 - 7) / 0.35   # 47 inclusive = 46 elapsed
  a <- compute_all(35, 2, 47, p, ref)
  expect_equal(a$bed_early, hand_early, tolerance = 1e-9)
  expect_equal(a$eqd2_early, hand_early / 1.2, tolerance = 1e-9)
  expect_equal(a$bed_late, 70 * (1 + 2 / 3))
  for (ab in c(2, 3, 4)) expect_equal(eqd2(bed_late(35, 2, ab), ab), 70)
  expect_equal(bed_early(5, 6, 5, 10, 0.35, log(2) / 2.5, 7),
               30 * 1.6 - (log(2) / 2.5) * smooth_max0(-2) / 0.35)
  expect_lt(abs(smooth_max0(-2)), 0.001)
})

test_that("DD calibration status: ok, no repopulation, infeasible", {
  cal <- function(q, rf = ref) dd_calibration(rf, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)
  expect_equal(cal(p)$status, "ok")
  expect_gt(cal(p)$onset_ref, 0)
  # no DI repopulation: DD has none either, both equal plain LQ BED
  q <- p; q$lambda_di <- 0
  expect_equal(cal(q)$status, "no_repopulation")
  a <- compute_all(20, 2.65, 26, q, ref)
  expect_equal(a$bed_dd, bed_lq(20, 2.65, 10))
  expect_equal(a$bed_dd, a$bed_di)
  # DD rate 0 cannot reproduce DI repopulation
  q <- p; q$lambda_dd <- 0
  expect_equal(cal(q)$status, "infeasible")
  expect_true(is.na(compute_all(20, 2.65, 26, q, ref)$bed_dd))
  # DD rate too low: onset would have to be before treatment
  q <- p; q$lambda_dd <- 0.2
  expect_equal(cal(q)$status, "infeasible")
  expect_lt(cal(q)$onset_ref, 0)
  expect_true(is.na(o_calibrate(q, ref)))  # the oracle agrees no root exists
})

test_that("suggestions equal an independent brute-force search", {
  settings <- list(
    list(criterion = "either", rank_by = "best"),
    list(criterion = "either", rank_by = "worst", tol_early = 3),
    list(criterion = "di", rank_by = "di", tol_early = 2, tol_late = 1),
    list(criterion = "either", rank_by = "best", min_gain = 1, min_days_saved = 7),
    list(criterion = "either", rank_by = "best", tol_early = 10, tol_late = 0),
    list(criterion = "both", rank_by = "worst", tol_early = 8, tol_late = 5),
    list(criterion = "dd", rank_by = "dd", tol_early = 8, tol_late = 5),
    list(criterion = "either", rank_by = "best", early_cap = 59),
    list(criterion = "either", rank_by = "best", early_cap = 54, d_step = 0.1),
    list(criterion = "dd", rank_by = "dd", tol_early = 8, tol_late = 5, min_gain = 1),
    list(criterion = "either", rank_by = "best", tol_early = 3, min_days_saved = 14))
  refs <- list(ref, list(m = 30, d = 2, T = 40), list(m = 33, d = 2.03, T = 45))
  n_rows <- 0
  for (rf in refs) for (an in list(rf, ref)) for (st in settings) {   # legacy coupling, then fixed anchor
    args <- modifyList(list(m_range = 5:(rf$m - 1), d_step = 0.05, d_max = 6,
                            criterion = "either", rank_by = "best",
                            tol_early = 0, tol_late = 0, early_cap = NULL,
                            min_gain = 0, min_days_saved = 1), st)
    if (is.null(st$early_cap)) args["early_cap"] <- list(NULL)
    s <- do.call(suggest_regimens, c(list(p = p, ref = rf, anchor = an), args))
    o <- do.call(o_search, c(list(q = p, rf = rf, anchor = an), args))
    expect_equal(nrow(s), nrow(o))
    n_rows <- n_rows + nrow(o)
    if (nrow(o) > 0 && nrow(s) == nrow(o)) {
      expect_equal(s$m, o$m)
      expect_equal(s$d, o$d, tolerance = 1e-9)
      expect_equal(s$T, o$T)
      expect_equal(s$score, o$score, tolerance = 1e-6)
      expect_equal(s$weeks_saved, (rf$T - o$T) / 7)
    }
  }
  expect_gt(n_rows, 50)  # the comparison is not vacuous
})

test_that("default search returns non-trivial results, and none when impossible", {
  s <- suggest_regimens(p, ref, min_gain = 1, min_days_saved = 7)
  expect_gt(nrow(s), 0)
  expect_true(all(s$d_di >= 1 - 1e-9 | s$d_dd >= 1 - 1e-9))
  expect_true(all(s$T <= ref$T - 7))
  expect_equal(nrow(suggest_regimens(p, ref, early_cap = 10)), 0)
  expect_true(all(names(empty_suggestions()) %in% names(s)))
})

test_that("next_step names the limit the next dose step actually breaks", {
  checks <- list(list(tol_early = 0, tol_late = 0), list(tol_early = 10, tol_late = 0),
                 list(tol_early = 10, tol_late = 3), list(tol_early = 2, tol_late = 0))
  q <- p; q$ab_l <- 1
  n_checked <- 0
  for (qq in list(p, q)) for (ck in checks) {
    s <- suggest_regimens(qq, ref, tol_early = ck$tol_early, tol_late = ck$tol_late, d_max = 6)
    r0 <- compute_all(ref$m, ref$d, ref$T, qq, ref)
    for (i in seq_len(nrow(s))) {
      nx <- s$d[i] + 0.05
      e_bad <- o_early(s$m[i], nx, s$T[i], qq) > r0$bed_early + ck$tol_early + 1e-9
      l_bad <- o_late(s$m[i], nx, qq) > r0$bed_late + ck$tol_late + 1e-9
      expect_equal(grepl("early toxicity", s$next_step[i]), e_bad)
      expect_equal(grepl("late toxicity", s$next_step[i]), l_bad)
      n_checked <- n_checked + 1
    }
  }
  expect_gt(n_checked, 50)
  # the case found in review: late is the limit although early has less room
  s <- suggest_regimens(p, ref, tol_early = 10, tol_late = 0, d_max = 6)
  expect_equal(s$next_step[s$m == 20], "late toxicity")
  # a row at the top of the dose range says so
  s <- suggest_regimens(p, ref, tol_early = 30, tol_late = 30, d_max = 3)
  expect_true(any(s$d == 3))
  expect_true(all(s$next_step[s$d == 3] == "end of dose search range"))
})

test_that("an infeasible DD model never leaks NA rows or DD criteria", {
  for (ldd in c(0, 0.2)) {
    q <- p; q$lambda_dd <- ldd
    s <- suggest_regimens(q, ref, tol_early = 5)
    expect_gt(nrow(s), 0)
    expect_false(any(is.na(s$m)))
    expect_true(all(is.na(s$d_dd)))
    expect_true(all(s$d_di >= -1e-9))
    for (cr in c("dd", "both")) {
      e <- suggest_regimens(q, ref, criterion = cr)
      expect_equal(nrow(e), 0)
      expect_equal(attr(e, "reason"), "dd_unavailable")
    }
    expect_equal(attr(suggest_regimens(q, ref, rank_by = "dd"), "reason"), "dd_unavailable")
    expect_equal(attr(suggest_regimens(q, ref, rank_by = "worst"), "reason"), "dd_unavailable")
  }
})

test_that("dose grid uses round multiples of the step", {
  s <- suggest_regimens(p, list(m = 33, d = 2.03, T = 45), tol_early = 5, tol_late = 5)
  expect_gt(nrow(s), 0)
  expect_true(all(abs(s$d / 0.05 - round(s$d / 0.05)) < 1e-9))
  expect_true(all(s$d > 2.03))
})

test_that("results stay finite over the allowed input box, including zero rates", {
  box <- expand.grid(alpha = c(0.01, 2), ab = c(0.5, 100), lambda_di = c(0, 2),
                     lambda_dd = c(0, 5), Tk = c(0, 100))
  for (i in seq_len(nrow(box))) {
    q <- p; for (k in names(box)) q[[k]] <- box[[k]][i]
    a <- compute_all(c(1, 10, 35, 60), c(30, 4, 2, 1), c(1, 12, 47, 80), q, ref)
    status <- dd_calibration(ref, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)$status
    cols <- setdiff(names(a), if (status == "infeasible") "bed_dd" else character())
    expect_true(all(is.finite(as.matrix(a[cols]))))
    s <- suggest_regimens(q, ref, d_max = 4)
    expect_false(any(is.na(s$m)))
  }
})

test_that("smooth max helpers are inverse and do not overflow", {
  y <- c(1e-6, 0.01, 1, 10, 100, 1e4)
  expect_equal(smooth_max0(inv_smooth_max0(y)), y, tolerance = 1e-9)
  expect_equal(smooth_max0(1e5), 1e5)
  expect_true(is.finite(inv_smooth_max0(1e5)))
})

test_that("Fowler zone, presets and benchmark regimens match their sources", {
  expect_equal(fowler_zone(FALSE), c(59, 61))
  expect_equal(fowler_zone(TRUE), c(54, 58))
  pr <- tumor_presets()
  expect_equal(pr$withers$lambda_di / pr$withers$alpha, 0.6 * (1 + 2 / 10))  # 0.6 Gy/day in 2-Gy fx = 0.72 Gy10/day
  expect_equal(pr$withers$Tk, 28)                              # about 4 weeks
  expect_equal(pr$hyhope$alpha, 0.35)
  expect_equal(pr$hyhope$lambda_di, log(2) / 3)               # Tp = 3 days
  expect_equal(pr$hyhope$Tk, 21)
  expect_equal(pr$calculator, list(alpha = 0.2, ab = 10, lambda_di = 0.2, lambda_dd = 0.5, Tk = 28))
  for (nm in names(pr)) {
    q <- modifyList(p, resolve_preset(nm, p$ab))
    expect_equal(dd_calibration(ref, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)$status, "ok")
  }
  b <- benchmark_regimens()
  expect_equal(b$m * b$d, c(55, 66, 46.5))
  expect_equal(b$T, c(o_calendar(20, 5), o_calendar(33, 6), o_calendar(15, 5)))
})

test_that("only 5 fractions per week (Mon-Fri) is offered", {
  expect_false("fx_per_week" %in% names(formals(overall_time)))
  expect_false("fx_per_week" %in% names(formals(suggest_regimens)))
  app_src <- readLines(testthat::test_path("..", "..", "app.R"))
  expect_false(any(grepl("Mon-Sat|fpw", app_src)))
})

test_that("empty searches say why", {
  expect_equal(attr(suggest_regimens(p, list(m = 35, d = 9, T = 47), d_max = 8), "reason"),
               "empty_dose_range")
  expect_equal(attr(suggest_regimens(p, ref, m_range = 34, min_days_saved = 14), "reason"),
               "no_shorter_regimen")
  expect_null(attr(suggest_regimens(p, ref, early_cap = 10), "reason"))
})

test_that("search range follows the reference number of fractions", {
  expect_equal(search_range_for_ref(35), c(5, 34))
  expect_equal(search_range_for_ref(30), c(5, 29))
  expect_equal(search_range_for_ref(3), c(2, 2))
  expect_equal(search_range_for_ref(2), c(1, 1))
  # retyping 35 -> 3 -> 30 ends at the full range again (bug found in review)
  r <- search_range_for_ref(3); r <- search_range_for_ref(30)
  expect_equal(r, c(5, 29))
})

test_that("preset matching recognizes presets and edited values", {
  pr <- tumor_presets()
  for (nm in names(pr)) expect_equal(preset_matching(pr[[nm]]$alpha, pr[[nm]]$lambda_di, pr[[nm]]$Tk), nm)
  expect_equal(preset_matching(0.2, 0.2, 27), "custom")
  expect_equal(preset_matching(0.21, 0.2, 28), "custom")
  expect_equal(preset_matching(NA, 0.2, 28), "custom")
  expect_equal(preset_matching(0.35, round(log(2) / 3, 7), 21), "hyhope")  # value echoed by the browser
})
