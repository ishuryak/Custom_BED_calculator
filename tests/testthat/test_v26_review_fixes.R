# Version 2.6: fixes from the external review of version 2.5.1
# (Documents/Codex/2026-10-18/rev-2/outputs/HNSCC_BED_v251_review_and_action_plan_2026-10-19.md,
# findings R1-R9; record in docs/development/2026_10_19_v26_review_response.md).
# Expected values come from the review's independent oracle (its evidence
# package: math/results/ranking_*.csv, clinical/slope_illustration.json) and
# from a schema-1 file written by the version 2.4 engine (tests/testthat/fixtures).

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
models_file <- Sys.getenv("BED_MODELS_FILE", file.path(app_abs, "R", "bed_models.R"))
source(models_file, local = TRUE)
source(file.path(app_abs, "R", "export.R"), local = TRUE)
read_text <- function(f) paste(readLines(file.path(app_abs, f), encoding = "UTF-8", warn = FALSE), collapse = "\n")
fixture <- function(f) file.path(app_abs, "tests", "testthat", "fixtures", f)

std <- list(m = 35, d = 2, T = 47)
hypno <- list(m = 20, d = 2.75, T = 26)
an <- std
p <- default_params()
# The review's minimal reproduction of R1 (all values inside the UI's limits).
pr1 <- modifyList(p, list(alpha = 0.2, ab = 10, lambda_di = 0.4, lambda_dd = 5, Tk = 28, g = 0))
r1_search <- function(rank_by, criterion = "either", d_step = 0.05, d_max = 2.5, m_range = 32, ...)
  suggest_regimens(pr1, std, anchor = std, m_range = m_range, d_step = d_step, d_max = d_max,
                   criterion = criterion, rank_by = rank_by, tol_early = 30, tol_late = 30,
                   min_gain = 0, min_days_saved = 1, ...)

# --- R1: ranking eligibility ------------------------------------------------

test_that("R1: the conservative ranking never lists a regimen whose DD BED is not positive", {
  s <- r1_search("worst")
  expect_equal(nrow(s), 1)
  expect_equal(s$d, 2.10)                                  # 2.5.1 listed 2.50 (DD BED -0.604)
  expect_equal(s$bed_di, 49.312, tolerance = 1e-9)
  expect_equal(s$bed_dd, 54.76422319, tolerance = 1e-8)
  expect_equal(s$score, 3.312, tolerance = 1e-9)
  f <- attr(s, "feasible")
  bad <- f[f$d == 2.5, ]
  expect_equal(bad$bed_dd, -0.60436688, tolerance = 1e-8)
  expect_false(bad$rank_eligible)
  expect_true(is.na(bad$score) && is.na(bad$d_min))
})

test_that("R1: a DD ranking with no evaluable DD returns an explicit reason, never an NA score", {
  s <- r1_search("dd", d_step = 0.5)
  expect_equal(nrow(s), 0)                                 # 2.5.1 listed 32 x 2.50 with score NA
  expect_identical(attr(s, "reason"), "no_rankable_row")
  f <- attr(s, "feasible")
  expect_equal(nrow(f), 1)
  expect_false(f$rank_eligible)
})

test_that("R1: all 16 criterion and ranking pairs list only rankable rows, best per fraction number", {
  states <- list(
    valid = list(p = p, ref = std, anchor = an, m = 20:34, d_max = 3.5, early_cap = 55.86597),
    dd_invalid_candidate = list(p = pr1, ref = std, anchor = std, m = 28:34, d_max = 3, early_cap = NULL),
    dd_unavailable = list(p = modifyList(p, list(lambda_dd = 0)), ref = std, anchor = an, m = 20:34,
                          d_max = 3.5, early_cap = 55.86597),
    di_reference_negative = list(p = modifyList(legacy_params(), list(lambda_di = 1, lambda_dd = 2)),
                                 ref = std, anchor = list(m = 20, d = 2, T = 26), m = 5:34, d_max = 4,
                                 early_cap = NULL))
  for (st in names(states)) for (cr in c("either", "both", "di", "dd")) for (rk in c("best", "worst", "di", "dd")) {
    z <- states[[st]]
    s <- suggest_regimens(z$p, z$ref, anchor = z$anchor, m_range = z$m, d_max = z$d_max, criterion = cr,
                          rank_by = rk, early_cap = z$early_cap, tol_early = 30, tol_late = 30,
                          min_gain = -50, min_days_saved = 1)
    info <- paste(st, cr, rk)
    expect_false(any(is.na(s$score)), info = info)
    if (rk == "worst") expect_true(all(is.finite(s$d_di) & is.finite(s$d_dd)), info = info)
    if (rk == "di") expect_true(all(is.finite(s$d_di)), info = info)
    if (rk == "dd") expect_true(all(is.finite(s$d_dd)), info = info)
    expect_true(all(is.na(s$d_min) | (is.finite(s$d_di) & is.finite(s$d_dd))), info = info)
    f <- attr(s, "feasible")
    if (nrow(s) == 0) {
      if (!is.null(f) && nrow(f) > 0) expect_identical(attr(s, "reason"), "no_rankable_row", label = info)
      next
    }
    for (i in seq_len(nrow(s))) {
      fm <- f[f$m == s$m[i] & f$rank_eligible, ]
      expect_equal(s$score[i], max(fm$score), tolerance = 1e-12, info = info)
    }
    expect_setequal(s$m, unique(f$m[f$rank_eligible]))
  }
})

test_that("R1: the next-dose column names a ranking that becomes non-evaluable", {
  # The largest d at which DD is still positive for 32 fractions.
  d_grid <- round(seq(2.05, 2.5, by = 0.05), 6)
  x <- compute_all(32, d_grid, overall_time(32), pr1, std)
  last_ok <- max(d_grid[x$bed_dd > 0])
  ref_di <- compute_all(35, 2, 47, pr1, std)$bed_di
  gain <- x$bed_di[d_grid == last_ok] - ref_di - 1e-6       # DI qualifies from last_ok upward only
  s <- suggest_regimens(pr1, std, anchor = std, m_range = 32, d_step = 0.05, d_max = 2.5,
                        criterion = "di", rank_by = "worst", tol_early = 30, tol_late = 30,
                        min_gain = gain, min_days_saved = 1)
  expect_equal(s$d, last_ok)
  expect_identical(s$next_step, "ranking not evaluable (a tumor BED is not positive)")
})

test_that("R1: default candidate rows are unchanged", {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  s <- suggest_regimens(p, std, anchor = an, m_range = 5:34, d_max = 5, early_cap = cap,
                        min_gain = 1, min_days_saved = 7)
  expect_equal(s$m, c(29, 28, 30))
  expect_equal(s$d, c(2.25, 2.30, 2.20))
  expect_equal(s$d_di, c(1.1884, 1.1263, 1.1200), tolerance = 1e-4)
  expect_equal(s$d_dd, c(-2.4080, -2.8648, -2.0495), tolerance = 1e-4)
  # Slack against the HYPNO acute limit (55.865972 Gy10, review table 7.1),
  # with the acute BED written out from Fowler's formula on the elapsed clock.
  acute <- function(n, d, T) n * d * (1 + d / 10) - log(2) / 2.5 * log1p(exp(3 * (T - 1 - 7))) / 3 / 0.35
  expect_equal(s$early_slack, 55.865972 - acute(s$m, s$d, s$T), tolerance = 1e-6)
  expect_equal(s$late_slack, 116.666667 - s$m * s$d * (1 + s$d / 3), tolerance = 1e-6)
  expect_equal(s$fractions_saved, 35 - s$m)
})

test_that("R1: the app explains a search with no rankable regimen", {
  shiny::testServer(app_dir, {
    session$setInputs(m = 20, d = 2.65, auto_T = TRUE, comparator = "custom", ref_m = 35, ref_d = 2,
                      ref_auto_T = TRUE, preset = "custom", alpha = 0.2, ab = 10, lambda_di = 0.4,
                      lambda_dd = 5, Tk = 28, g = 0, anchor_m = 35, anchor_d = 2, ab_e = 10,
                      alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
                      criterion = "either", rank_by = "worst", m_range = c(32, 32), min_gain = 0,
                      min_days = 1, d_step = "0.05", d_max = 2.5, early_mode = "reference",
                      late_mode = "reference", tol_early = 30, tol_late = 30)
    expect_equal(suggestions()$d, 2.10)
    session$setInputs(rank_by = "dd", d_step = "0.5")
    expect_equal(nrow(suggestions()), 0)
    expect_true(grepl("none can be ranked", output$suggest_msg$html, fixed = TRUE))
  })
})

# --- R5: agreement ranges use the evaluable-scenario population ---------------

test_that("R5: both ranges use the scenarios in which both models are evaluable", {
  q0 <- modifyList(p, list(lambda_di = 0))
  ref150 <- list(m = 35, d = 2, T = 150)
  cap <- compute_all(20, 2.75, 26, q0, an)$bed_early
  a <- agreement_regimens(q0, ref150, an, margin = -3, early_cap = cap, min_days_saved = 1,
                          d_max = 5, m_range = 1:34)
  r <- a[a$m == 32, ]
  expect_equal(r$n_valid_both, 3)
  expect_equal(r$n_valid_di, 3)
  expect_equal(r$n_valid_dd, 15)
  expect_equal(r$d_dd_min, 0.98809, tolerance = 1e-5)     # 2.5.1 gave -0.21437 (all 15 DD-valid)
  expect_equal(r$d_dd_max, 2.38883, tolerance = 1e-5)
  det <- attr(a, "scenario_detail"); det <- det[det$m == 32 & det$valid_di & det$valid_dd, ]
  expect_equal(c(r$d_di_min, r$d_di_max), range(det$d_di), tolerance = 1e-12)
  expect_equal(c(r$d_dd_min, r$d_dd_max), range(det$d_dd), tolerance = 1e-12)
})

test_that("R5: all-invalid scenarios give NA ranges", {
  st <- data.frame(valid_di = c(FALSE, TRUE), valid_dd = c(TRUE, FALSE), d_di = c(NA, 1), d_dd = c(2, NA))
  z <- scenario_summary(st, -3)
  expect_true(all(is.na(unlist(z[c("d_di_min", "d_di_max", "d_dd_min", "d_dd_max")]))))
  expect_equal(z$n_valid_both, 0); expect_equal(z$n_pass_both, 0)
})

# --- R2: strict replay boundary ---------------------------------------------

fenwick_scenario <- function() {
  build_scenario(p, std, an,
                 search = list(m_range = 5:34, d_step = 0.05, d_max = 5, criterion = "either",
                               rank_by = "best", tol_early = 0, tol_late = 0, early_cap = NULL,
                               late_cap = NULL, min_gain = 1, min_days_saved = 7, fenwick = TRUE),
                 context = list())
}

test_that("R2: the Fenwick flag must be one true or false value", {
  sc <- jsonlite::fromJSON(scenario_json(fenwick_scenario()))
  direct <- suggest_regimens(p, std, anchor = an, m_range = 5:34, d_max = 5, min_gain = 1,
                             min_days_saved = 7, fenwick = TRUE)
  expect_gt(nrow(direct), 0)
  expect_equal(nrow(rerun_scenario(sc)), nrow(direct))
  for (bad in list("true", 1, NA, c(TRUE, FALSE))) {
    x <- sc; x$search$fenwick <- bad
    expect_error(rerun_scenario(x), "fenwick must be true or false", info = deparse(bad))
  }
  x <- sc; x$search$fenwick <- NULL
  expect_error(rerun_scenario(x), "fenwick must be true or false")
})

test_that("R2: enumerated fields, calibration and the grid size are validated", {
  sc <- jsonlite::fromJSON(scenario_json(fenwick_scenario()))
  bad <- function(f) { x <- sc; f(x) }
  expect_error(rerun_scenario(bad(function(x) { x$search$criterion <- c("both", "either"); x })), "criterion")
  expect_error(rerun_scenario(bad(function(x) { x$search$rank_by <- 1; x })), "rank_by")
  expect_error(rerun_scenario(bad(function(x) { x$analysis <- "both"; x })), "analysis")
  expect_error(rerun_scenario(bad(function(x) { x$calibration$onset_anchor <- 999; x })), "onset_anchor")
  expect_error(rerun_scenario(bad(function(x) { x$search$m_range <- c(5, 5, 6); x })), "m_range")
  expect_error(rerun_scenario(bad(function(x) { x$search$m_range <- 1:101; x })), "m_range")
  expect_error(rerun_scenario(bad(function(x) { x$search$d_step <- 0.001; x$search$d_max <- 30;
                                                x$search$m_range <- 1:34; x })), "exceeds")
})

agree_scenario <- function() {
  cap <- compute_all(20, 2.75, 26, p, an)$bed_early
  build_scenario(p, hypno, an,
                 search = list(m_range = 1:19, d_step = 0.05, d_max = 5, criterion = "both",
                               rank_by = "worst", tol_early = 0, tol_late = 0, early_cap = cap,
                               late_cap = NULL, min_gain = -3, min_days_saved = 1, fenwick = FALSE),
                 context = list(), analysis = "agreement",
                 agreement = list(margin = -3, ab_values = SCENARIO_AB_VALUES))
}

test_that("R2: agreement files cannot carry search fields that contradict the agreement", {
  sc <- jsonlite::fromJSON(scenario_json(agree_scenario()))
  expect_equal(rerun_scenario(sc)$m, c(17, 18, 19))
  bad <- function(f) { x <- sc; f(x) }
  expect_error(rerun_scenario(bad(function(x) { x$search$criterion <- "either"; x })), "agreement file")
  expect_error(rerun_scenario(bad(function(x) { x$search$rank_by <- "best"; x })), "agreement file")
  expect_error(rerun_scenario(bad(function(x) { x$search$min_gain <- 1; x })), "agreement file")
  expect_error(rerun_scenario(bad(function(x) { x$agreement$margin_eqd2$ab_10 <- -9; x })), "margin_eqd2")
  rev_sc <- agree_scenario(); rev_sc$search <- rev(rev_sc$search); rev_sc <- rev_sc[rev(names(rev_sc))]
  expect_identical(scenario_id(rev_sc), scenario_id(agree_scenario()))
})

test_that("R2: source hashes are read from the bytes and flag code changed since loading", {
  ci <- code_identity(app_abs)
  f <- file.path(app_abs, "R", "bed_models.R")
  expect_identical(ci$files[["R/bed_models.R"]],
                   digest::digest(readBin(f, "raw", file.size(f)), algo = "sha256", serialize = FALSE))
  expect_match(ci$files[["R/bed_models.R"]], "^[0-9a-f]{64}$")
  sc <- build_scenario(p, std, an, search = list(), context = list(), code_root = app_abs, code = ci)
  expect_false(sc$code$changed_since_load)
  old <- ci; old$files[["app.R"]] <- "0000"
  sc2 <- build_scenario(p, std, an, search = list(), context = list(), code_root = app_abs, code = old)
  expect_true(sc2$code$changed_since_load)
})

# --- R3: strict replay versus migration -------------------------------------

test_that("R3: a version 2.4 file is refused by strict replay and reanalyzed only on request", {
  sc <- jsonlite::fromJSON(fixture("schema1_v24_r1case.json"))
  orig <- read.csv(fixture("schema1_v24_r1case_results.csv"))
  expect_equal(sc$app_version, "2.4.0")
  expect_error(rerun_scenario(sc), "Exact replay needs that version's engine")
  m <- migrate_scenario(sc)
  expect_equal(m$migration$from$app_version, "2.4.0")
  expect_equal(m$migration$from$schema_version, "1")
  expect_equal(m$migration$from$scenario_id, sc$scenario_id)
  expect_false(identical(m$scenario_id, sc$scenario_id))
  expect_true("search$fenwick = false" %in% m$migration$defaults_filled ||
                identical(sc$search$fenwick, FALSE))
  expect_match(m$migration$note, "not an exact replay")
  r <- rerun_scenario(jsonlite::fromJSON(scenario_json(m)))
  expect_equal(r$d, orig$d)                                  # 2.10 Gy, as the 2.4 engine chose
  expect_equal(r$score, orig$score, tolerance = 1e-8)        # 2.5.1's migration gave 2.50 Gy
  expect_equal(r$bed_dd, orig$bed_dd, tolerance = 1e-8)
  expect_equal(attr(r, "migration")$from$app_version, "2.4.0")
})

test_that("R3: version 2.5 files need migration; unknown and Fenwick-2.1 files are refused", {
  sc <- jsonlite::fromJSON(scenario_json(agree_scenario()))
  # a version 2.5 file is schema 2 with equations 2.2 (version 3 writes schema 3)
  as_v25 <- function(x) { x$selection_version <- "2.5"; x$schema_version <- "2"; x$equation_version <- "2.2"; x$dd <- NULL; x }
  old <- as_v25(sc); old$app_version <- "2.5.1"
  expect_error(rerun_scenario(old), "migrate_scenario")
  expect_equal(rerun_scenario(migrate_scenario(old))$m, c(17, 18, 19))
  expect_error(migrate_scenario(sc), "already of the current versions")
  fut <- sc; fut$selection_version <- "9.9"
  expect_error(rerun_scenario(fut), "not supported")
  expect_error(migrate_scenario(fut), "cannot be migrated")
  f21 <- jsonlite::fromJSON(scenario_json(fenwick_scenario()))
  f21$schema_version <- "1"; f21$equation_version <- "2.1"; f21$selection_version <- NULL
  expect_error(migrate_scenario(f21), "no Fenwick boundary")
  h <- as_v25(sc); h$hinge$type <- "hard"
  expect_error(migrate_scenario(h), "hinge")
})

# --- R4: conditional dose-response illustration -----------------------------

test_that("R4: the illustration reproduces the review's independent values", {
  source(file.path(app_abs, "docs", "development", "2026_10_19_v26", "dose_response_illustration.R"), local = TRUE)
  t <- illustration_table()
  v <- function(model, g, col) t[[col]][t$model == model & t$gamma50 == g]
  expect_equal(v("DI", 1.6, "control_bed"), 72.62857142857143, tolerance = 1e-10)
  expect_equal(v("DD", 1.6, "control_bed"), 68.64033914, tolerance = 1e-8)
  expect_equal(v("DI", 1.6, "delta_bed"), -2.50411372, tolerance = 1e-8)
  expect_equal(v("DD", 1.6, "delta_bed"), -5.79303791, tolerance = 1e-8)
  expect_equal(v("DD", 1.6, "loss_pp_logistic"), 13.356198372098266, tolerance = 1e-6)
  expect_equal(v("DD", 3.2, "loss_pp_logistic"), 25.01401097609651, tolerance = 1e-6)
  expect_equal(v("DD", 1.2, "loss_pp_logistic"), 10.131834251428252, tolerance = 1e-6)
  expect_equal(v("DI", 1.6, "loss_pp_logistic"), 5.546688479946981, tolerance = 1e-6)
  expect_equal(v("DI", 3.2, "loss_pp_logistic"), 10.947831348172082, tolerance = 1e-6)
  expect_equal(attr(t, "ph_upper_loss_pp"), 9.315730122196053, tolerance = 1e-9)
})

test_that("R4: the About page states the illustration's assumptions and no unconditional verdict", {
  ab <- read_text("about.md")
  expect_true(grepl("Conditional illustration, not a prediction", ab, fixed = TRUE))
  expect_true(grepl("a transfer, not a derivation", ab, fixed = TRUE))
  expect_true(grepl("13.4–25.0 percentage points", ab, fixed = TRUE))
  expect_false(grepl("upper 90% confidence limit (a hazard ratio of 1.30, about 9 points) allows", ab, fixed = TRUE))
  expect_false(grepl("both would be compatible", ab, fixed = TRUE))
  expect_false(grepl("roughly 13–27", ab, fixed = TRUE))
})

# --- R6 and R7: comparator scope, history, wording, sources -----------------

test_that("R6: each analysis tab names its own comparison, anchor and numerical limits", {
  shiny::testServer(app_dir, {
    session$setInputs(m = 20, d = 2.65, auto_T = TRUE, comparator = "hypno_control", ref_m = 35, ref_d = 2,
                      ref_auto_T = TRUE, preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23,
                      lambda_dd = 0.4, Tk = 28, g = 0, anchor_m = 35, anchor_d = 2, ab_e = 10,
                      alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
                      criterion = "either", rank_by = "best", m_range = c(5, 32), min_gain = 1,
                      min_days = 7, d_step = "0.05", d_max = 5, early_mode = "hypno",
                      late_mode = "reference", tol_early = 0, tol_late = 0, agree_cmp = "hypno_test",
                      agree_threshold_mode = "own", agree_margin = -3, agree_days = 1)
    sc <- output$search_context$html
    expect_true(grepl("Comparison: 33 × 2 Gy in 38 days (Report card page)", sc, fixed = TRUE))
    expect_true(grepl("DD anchor: 35 × 2 Gy in 47 days", sc, fixed = TRUE))
    expect_true(grepl("55.87 Gy₁₀ (HYPNO 55 Gy/20, fixed)", sc, fixed = TRUE))
    expect_true(grepl("110.00 Gy₃ (comparison 110.00 + 0)", sc, fixed = TRUE))
    ac <- output$agree_context$html
    expect_true(grepl("Comparison: 20 × 2.75 Gy in 26 days (override: HYPNO test arm)", ac, fixed = TRUE))
    expect_true(grepl("105.42 Gy₃ (comparison 105.42 + 0)", ac, fixed = TRUE))
    expect_true(grepl("at least 1 day saved", ac, fixed = TRUE))
    expect_true(grepl("at least 7 days saved", sc, fixed = TRUE))
    expect_true(grepl("the Report card regimen, not a candidate row", output$sens_context$html, fixed = TRUE))
    dd_fixed <- compute_all(20, 2.75, 26, params(), anchor())$bed_dd
    session$setInputs(comparator = "custom")
    expect_equal(compute_all(20, 2.75, 26, params(), anchor())$bed_dd, dd_fixed)
    h <- output$anchor_note$html
    expect_true(grepl("Version 1 used a fixed 35 × 2 Gy, 47-day anchor", h, fixed = TRUE))
    expect_true(grepl("45 days here", h, fixed = TRUE))
    expect_false(grepl("To reproduce version 1, which calibrated at the comparison", h, fixed = TRUE))
  })
  app_txt <- read_text("app.R")
  expect_false(grepl("All differences and toxicity limits are relative to this regimen", app_txt, fixed = TRUE))
  # version 3 (decision 7a), worded for the version 4 pages: one comparison and one threshold
  expect_true(grepl("One comparison and one threshold for every page except the main list", app_txt, fixed = TRUE))
  ab <- read_text("about.md")
  expect_false(grepl("Versions 1 and 2.0 calibrated at the comparison regimen", ab, fixed = TRUE))
})

test_that("R7: default narrative, preset wording, version footer and source records are current", {
  ab <- read_text("about.md")
  expect_false(grepl("under both models, so hypofractionation questions need an allowed loss", ab, fixed = TRUE))
  expect_true(grepl("DI-only hypotheses", ab, fixed = TRUE))
  expect_true(grepl("test minus control −0.5 percentage points", ab, fixed = TRUE))
  expect_false(grepl("Version 2.5 of the Custom", ab, fixed = TRUE))
  shiny::testServer(app_dir, {
    session$setInputs(preset = "withers", ab = 10, alpha = 0.2, lambda_di = 0.144, Tk = 28, lambda_dd = 0.4)
    h <- output$preset_note$html
    expect_true(grepl("calculator's", h, fixed = TRUE) && grepl("choice to complete the preset", h, fixed = TRUE))
    expect_true(grepl("A preset fills only the values it sets", h, fixed = TRUE))
    about <- as.character(output$about$html)
    expect_true(grepl(sprintf("Version %s of the Custom BED calculator (equations %s, selection rules %s",
                              APP_VERSION, EQUATION_VERSION, SELECTION_VERSION), about, fixed = TRUE))
  })
  rr <- read_text("References/README.md")
  for (sec in c("### A2.", "### A3.", "### B1.", "### C1.")) {
    block <- sub("\n### .*", "", substring(rr, regexpr(sec, rr, fixed = TRUE)))
    expect_true(grepl("Current status (2026-10-19)", block, fixed = TRUE), info = sec)
  }
  expect_true(grepl("test minus control −0.5 percentage points", rr, fixed = TRUE))
  expect_false(grepl("3-year LRC difference ≤ 1.4 points", rr, fixed = TRUE))
  expect_true(grepl("not a measured drug conversion", rr, fixed = TRUE))
  expect_true(grepl("## G. Current parameter ledger (version 2.6)", rr, fixed = TRUE))
})

test_that("R1/R7: versions agree across the code and DESCRIPTION", {
  expect_identical(SELECTION_VERSION, "2.6")
  desc <- read.dcf(file.path(app_abs, "DESCRIPTION"))
  expect_equal(unname(desc[1, "Version"]), APP_VERSION)
})
