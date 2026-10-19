# Tests for the version 2.1 changes (docs/development/2026_10_18_synthesis_plan.md).
# Expected numbers come from published tables or from the independent Python
# implementation in docs/development/2026_10_18_review_scripts/fix.py, not
# from R/bed_models.R.

models_file <- Sys.getenv("BED_MODELS_FILE",
                          testthat::test_path("..", "..", "R", "bed_models.R"))
source(models_file, local = TRUE)
export_file <- testthat::test_path("..", "..", "R", "export.R")
if (file.exists(export_file)) source(export_file, local = TRUE)

p <- default_params()
std <- list(m = 35, d = 2, T = 47)
anchor_fn <- function() if (exists("default_anchor")) default_anchor() else std

# --- 1. Acute clock: Fowler counts elapsed days ------------------------------

test_that("acute BED follows Fowler's elapsed-day convention (Fowler 2003, 2010)", {
  a <- compute_all(35, 2, 47, p, std)                 # 35 x 2 Gy, 47 inclusive days
  expect_equal(a$bed_early, 53.10543995, tolerance = 1e-8)
  expect_equal(round(a$bed_early, 1), 53.1)           # Fowler 2003 Table 1, "46 d"
  cair <- compute_all(35, 2, 35, p, std)              # CAIR, 7 fx/week, "34 d"
  expect_equal(round(cair$bed_early, 1), 62.6)
  f10 <- compute_all(c(60, 60), c(1.3, 1.3), c(40, 43), p, std)   # Fowler 2010: 39 and 42 elapsed
  expect_equal(f10$bed_early, c(62.79061740, 60.41411278), tolerance = 1e-8)
  expect_equal(compute_all(20, 2.75, 26, p, std)$bed_early, 55.86597229, tolerance = 1e-8)
  expect_equal(compute_all(22, 2.75, 30, p, std)$bed_early, 59.70979946, tolerance = 1e-8)
  expect_gt(compute_all(22, 2.75, 30, p, std)$bed_early, 59)   # fails the corrected 59 cap
})

test_that("the low-level acute function keeps its time on the model clock", {
  # archived scripts call bed_early(T - 1) themselves; it must not subtract again
  expect_equal(bed_early(35, 2, 46, 10, 0.35, log(2) / 2.5, 7), 53.10543995, tolerance = 1e-8)
})

# --- 2. Fixed DD calibration anchor ------------------------------------------

test_that("the DD anchor is separate from the comparison regimen", {
  expect_true("anchor" %in% names(formals(suggest_regimens)))
  expect_equal(anchor_fn(), std)
  dd_fixed <- compute_all(20, 2.75, 26, p, anchor_fn())$bed_dd
  expect_equal(dd_fixed, 58.02143210, tolerance = 1e-8)
  # comparing against HYPNO's control must not move the DD model
  s <- suggest_regimens(p, list(m = 33, d = 2, T = 38), anchor = anchor_fn(),
                        tol_early = 10, tol_late = 10, min_gain = -30, d_max = 4)
  expect_gt(nrow(s), 0)
  chk <- compute_all(s$m, s$d, s$T, p, anchor_fn())
  expect_equal(s$bed_dd, chk$bed_dd, tolerance = 1e-10)
  # legacy coupling (anchor = comparison) reproduces the historical value
  expect_equal(compute_all(20, 2.75, 26, p, list(m = 20, d = 2.75, T = 26))$bed_dd,
               70.12417477, tolerance = 1e-6)
})

# --- 3. Presets ---------------------------------------------------------------

test_that("Withers preset follows the tumor alpha/beta", {
  expect_equal(withers_lambda(0.2, 5) / 0.2, 0.84, tolerance = 1e-12)
  expect_equal(withers_lambda(0.2, 10) / 0.2, 0.72, tolerance = 1e-12)
  q <- modifyList(p, list(ab = 5, lambda_di = withers_lambda(0.2, 5)))
  expect_equal(compute_all(35, 2, 47, q, std)$bed_di, 82.04, tolerance = 1e-6)
  expect_equal(preset_matching(0.2, withers_lambda(0.2, 5), 28, ab = 5), "withers")
  expect_equal(preset_matching(0.2, 0.144, 28, ab = 5), "custom")
})

test_that("source presets carry their published values and clocks", {
  pr <- tumor_presets()
  expect_equal(pr$hartley$alpha, 0.3)
  expect_equal(pr$hartley$lambda_di, log(2) / 3)
  expect_equal(pr$hartley$Tk, 23)                 # 22 elapsed days = 23 inclusive
  expect_equal(pr$hartley$ab, 10)
  expect_equal(pr$hyhope$ab, 10)
  d18 <- pr$d18
  expect_equal(unlist(d18[c("alpha", "ab", "lambda_di", "lambda_dd", "Tk")]),
               c(alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.40, Tk = 28))
  q <- modifyList(p, d18[c("alpha", "ab", "lambda_di", "lambda_dd", "Tk")])
  a <- compute_all(c(35, 20, 33), c(2, 2.75, 2), c(47, 26, 38), q, std)
  expect_equal(a$bed_di, c(71.51428571, 70.12445771, 72.62857143), tolerance = 1e-8)
  expect_equal(a$bed_dd, c(71.51428571, 62.84730123, 68.64033914), tolerance = 1e-8)
  expect_equal(preset_matching(0.35, 0.23, 28, ab = 10, lambda_dd = 0.4), "d18")
  expect_equal(preset_matching(0.35, 0.23, 28, ab = 10, lambda_dd = 0.5), "custom")
  expect_equal(preset_matching(0.2, 0.2, 28, ab = 10, lambda_dd = 0.5), "calculator")
})

# --- 4. Search validity and trial-design options -----------------------------

test_that("a negative reference BED blocks that endpoint", {
  q <- p; q$lambda_di <- 1        # reference DI BED = 84 - 5 * 19.0 < 0
  expect_lt(compute_all(35, 2, 47, q, std)$bed_di, 0)
  e <- suggest_regimens(q, std, anchor = std, criterion = "di")
  expect_equal(nrow(e), 0)
  expect_equal(attr(e, "reason"), "negative_reference_di")
  s <- suggest_regimens(q, std, anchor = std, tol_early = 10)
  expect_false(any(s$qual_di))
})

test_that("search rows say which model credits them, and flag fraction sizes beyond precedent", {
  s <- suggest_regimens(p, std, anchor = std, tol_early = 8, tol_late = 8, d_max = 4)
  expect_true(all(c("qual_di", "qual_dd", "credited_by", "d_min", "beyond_precedent") %in% names(s)))
  expect_equal(s$d_min, pmin(s$d_di, s$d_dd), tolerance = 1e-12)
  expect_equal(s$beyond_precedent, s$d > 3.1 + 1e-9)
  expect_true(all(s$credited_by %in% c("DI and DD", "DI only", "DD only")))
  f <- attr(s, "feasible")
  expect_true(is.data.frame(f) && nrow(f) >= nrow(s))
})

test_that("search counts under the corrected clock (Review B fixtures)", {
  cnt <- function(q, ...) nrow(suggest_regimens(q, std, anchor = std, min_gain = 1,
                                                min_days_saved = 7, d_step = 0.05, d_max = 8, ...))
  q18 <- modifyList(p, tumor_presets()$d18[c("alpha", "ab", "lambda_di", "lambda_dd", "Tk")])
  expect_equal(cnt(p), 7)
  expect_equal(cnt(p, early_cap = 59), 17)
  expect_equal(cnt(q18), 0)
  expect_equal(cnt(q18, early_cap = 59), 11)
  expect_equal(cnt(q18, early_cap = 54), 0)
  s <- suggest_regimens(q18, std, anchor = std, min_gain = 1, min_days_saved = 7, early_cap = 59)
  expect_false(any(s$qual_dd))   # none credited by DD
})

test_that("absolute late cap (HYPNO 55/20) is honored", {
  s <- suggest_regimens(p, std, anchor = std, early_cap = 55.86597229,
                        late_cap = 105.41666667, min_gain = -10, d_max = 4)
  expect_gt(nrow(s), 0)
  expect_true(all(s$bed_late <= 105.41666667 + 1e-9))
  expect_true(all(s$bed_early <= 55.86597229 + 1e-9))
})

# --- 5. Background growth g --------------------------------------------------

test_that("background growth g subtracts gT/alpha in both models and keeps the calibration", {
  q <- p; q$g <- 0.0154
  a0 <- compute_all(c(35, 20), c(2, 2.75), c(47, 26), p, std)
  a1 <- compute_all(c(35, 20), c(2, 2.75), c(47, 26), q, std)
  expect_equal(a0$bed_di - a1$bed_di, c(3.619, 2.002), tolerance = 1e-9)
  expect_equal(a0$bed_dd - a1$bed_dd, c(3.619, 2.002), tolerance = 1e-9)
  expect_equal(a1$bed_dd[1], a1$bed_di[1], tolerance = 1e-9)
  q0 <- q; q0$lambda_di <- 0            # no accelerated repopulation, growth remains
  expect_equal(compute_all(20, 2.75, 26, q0, std)$bed_dd, 70.125 - 2.002, tolerance = 1e-9)
})

# --- 6. Late BED sensitivity --------------------------------------------------

test_that("late comparison of 20 x 2.9 Gy reverses between alpha/beta 3 and 2", {
  expect_equal(bed_late(20, 2.9, 3), 114.0666667, tolerance = 1e-7)
  expect_equal(bed_late(35, 2, 3), 116.6666667, tolerance = 1e-7)
  expect_equal(bed_late(20, 2.9, 2), 142.1, tolerance = 1e-9)
  expect_equal(bed_late(35, 2, 2), 140, tolerance = 1e-9)
})

# --- 7. Export round trip -----------------------------------------------------

test_that("an exported scenario reruns to the same rows", {
  expect_true(exists("build_scenario") && exists("rerun_scenario"))
  sc <- build_scenario(p = p, ref = std, anchor = std,
                       search = list(m_range = 5:34, d_step = 0.05, d_max = 5,
                                     criterion = "either", rank_by = "best",
                                     tol_early = 0, tol_late = 0, early_cap = NULL,
                                     late_cap = NULL, min_gain = 1, min_days_saved = 7),
                       context = list(comparator = "custom", early_mode = "reference",
                                      late_mode = "reference", chemo = FALSE, coupled = FALSE))
  js <- scenario_json(sc)
  sc2 <- jsonlite::fromJSON(js, simplifyVector = TRUE)
  r1 <- suggest_regimens(p, std, anchor = std, m_range = 5:34, d_step = 0.05, d_max = 5,
                         min_gain = 1, min_days_saved = 7)
  r2 <- rerun_scenario(sc2)
  expect_equal(r2$m, r1$m)
  expect_equal(r2$bed_dd, r1$bed_dd, tolerance = 1e-8)
  expect_equal(r2$bed_early, r1$bed_early, tolerance = 1e-8)
  sc3 <- sc; sc3$search$rank_by <- "worst"
  expect_false(identical(scenario_id(sc), scenario_id(sc3)))
  expect_identical(scenario_id(sc), scenario_id(sc))
})
