# Version 3 interface (app.R with the Organs at risk and Uncertainty modules),
# driven with shiny::testServer.

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
ui_html <- function() {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  htmltools::renderTags(env$ui)$html
}

v3_defaults <- list(
  m = 20, d = 2.75, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 26,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  cmp1_m = 17, cmp1_d = 3.05, cmp2_m = 15, cmp2_d = 3.1, cmp3_m = 30, cmp3_d = 2,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent",
  `oar-use_in_search` = FALSE, `unc-variant` = "primary", `unc-tumor_crit` = "di", `unc-tox_source` = "prior",
  `unc-tcp0` = 70, `unc-rule` = "min_q10", `unc-tox_share` = 90)
set_v3 <- function(session, ...) do.call(session$setInputs, modifyList(v3_defaults, list(...)))

test_that("one comparison and one threshold: the both-model screen follows the sidebar unless overridden", {
  shiny::testServer(app_dir, {
    set_v3(session, comparator = "hypno_test", min_gain = -3, agree_margin = 5)   # distinct own value, not in use
    expect_equal(agree_ref()[c("m", "d", "T")], list(m = 20, d = 2.75, T = 26))
    expect_equal(agree_threshold(), -3)
    expect_null(output$agree_override)
    expect_equal(agreement()$m, c(17, 18, 19))             # the v2.4 result, now reached through the global settings
    session$setInputs(agree_threshold_mode = "own", agree_margin = 0)
    expect_equal(agree_threshold(), 0)
    expect_true(grepl("Override on this page", output$agree_override$html))
    expect_true(grepl("threshold 0 Gy instead of the global -3 Gy", output$agree_override$html, fixed = TRUE))
  })
})

test_that("the banner names the evaluated regimen, the comparison, the threshold and the DD set-up", {
  shiny::testServer(app_dir, {
    set_v3(session)
    b <- output$banner$html
    expect_true(grepl("20 × 2.75 Gy in 26 days", b, fixed = TRUE))
    expect_true(grepl("35 × 2 Gy in 47 days", b, fixed = TRUE))
    expect_true(grepl("tumor threshold +1 Gy", b, fixed = TRUE))
    expect_true(grepl("DD calibrated at 35 × 2 Gy in 47 days", b, fixed = TRUE))
  })
})

test_that("report card: values, named differences, limit status and the refit column", {
  shiny::testServer(app_dir, {
    set_v3(session)
    r <- output$report
    expect_true(grepl("Δ vs 35 × 2 Gy in 47 days", r, fixed = TRUE))
    # Version 3.1: the column shows the verified-trial refit with calibrated DD
    # (was "Bayesian refit (for review)" ... "prior-dominated" for the 2018-table fit).
    expect_true(grepl("Verified-trial refit 2026-10-25", r, fixed = TRUE))
    expect_true(grepl("(calibrated DD)", r, fixed = TRUE))
    b <- both()
    expect_true(grepl(sprintf("%.2f", b$bed_di[1]), r, fixed = TRUE))
    expect_true(grepl("within (limit 55.9)", r, fixed = TRUE))  # HYPNO acute limit, HYPNO itself is within
    expect_true(grepl("outside the tumor data", output$report_notes$html))
    # comparison equal to the anchor: DD's comparison value is marked as calibrated
    expect_true(grepl("equal to DI by calibration", r, fixed = TRUE))
  })
})

test_that("organ rows reach the report card, and blank organs are not assessed", {
  shiny::testServer(app_dir, {
    set_v3(session)
    expect_false(grepl("Spinal cord", output$report, fixed = TRUE))
    session$setInputs(`oar-v_cord` = 80, `oar-v_parotid` = 40)
    r <- output$report
    expect_true(grepl("Spinal cord, maximum EQD2", r, fixed = TRUE))
    expect_true(grepl("Parotid (spared side), mean physical dose", r, fixed = TRUE))   # OAR-01 (stage A): a mean limit is physical
  })
})

test_that("fitted DD: the DD BED uses the fitted values, the anchor is not used", {
  shiny::testServer(app_dir, {
    set_v3(session, dd_mode = "fitted", alpha_dd = 0.0272, c_over_alpha_dd = 392.8, g_dd = 0.0171, lambda_dd = 0.467)
    p <- params()
    b <- both()
    expect_equal(b$bed_dd[1], bed_dd(20, 2.75, 26, 0.0272, 10, 0.467, 392.8, 0.0171))
    expect_equal(calibration()$status, "fitted")
    # 3.1 (stage E): the values in this mode are typed by the user (they start at the refit's), so the banner names the mode
    expect_true(grepl("DD values entered directly (not calibrated)", output$banner$html, fixed = TRUE))
    expect_true(grepl("does not use the anchor|instead, so the anchor is not used", output$anchor_note$html))
    expect_false(grepl("equal to DI by calibration", output$report, fixed = TRUE))
  })
})

test_that("calendars: 6 fractions/week, a Thursday start and skipped days change the overall time", {
  shiny::testServer(app_dir, {
    set_v3(session, m = 33, d = 2, fx_week = "6")
    expect_equal(regimen()$T, 38)
    set_v3(session, m = 5, start_day = "4")
    expect_equal(regimen()$T, 7)
    set_v3(session, skip_days = "3, 4")
    expect_equal(regimen()$T, schedule_time(22))
    set_v3(session, start_day = "6")                  # Saturday with 5 per week
    expect_error(regimen(), "Monday to Friday")
    set_v3(session, skip_days = "1")
    expect_error(regimen(), "2 to 400")
    set_v3(session, auto_T = FALSE, T = 20)           # 20 fractions cannot fit in 20 days on Mon-Fri
    expect_error(regimen(), "need at least 26 days")
  })
})

test_that("6 fractions/week also shapes the search grid", {
  shiny::testServer(app_dir, {
    set_v3(session, fx_week = "6", criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 10)
    s <- suggestions()
    expect_true(nrow(s) > 0)
    expect_equal(unname(s$T), schedule_time(s$m, 6, 1))
  })
})

test_that("difference headers name the comparison and the unit; days, not weeks", {
  shiny::testServer(app_dir, {
    set_v3(session, criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 10)
    nm <- names(pretty_suggestions())
    expect_true("Δ BED DI vs 35 × 2 Gy in 47 days (Gy₁₀)" %in% nm)
    expect_true("Δ late BED vs 35 × 2 Gy in 47 days (Gy₃)" %in% nm)
    expect_true("Days saved" %in% nm)
    expect_false("Weeks saved" %in% nm)
    expect_true(sprintf("Above benchmark fraction sizes (> %g Gy)", PRECEDENT_MAX_D) %in% nm)
    expect_true(grepl("Δ BED DI vs 35 × 2 Gy in 47 days (Gy₁₀)", output$benchmarks, fixed = TRUE))
  })
})

test_that("each toxicity limit shows its value next to its control", {
  shiny::testServer(app_dir, {
    set_v3(session)
    expect_true(grepl("Value now: 55.87 Gy₁₀", output$early_limit_value$html, fixed = TRUE))
    expect_true(grepl("Value now: 116.67 Gy₃", output$late_limit_value$html, fixed = TRUE))
    session$setInputs(early_mode = "fenwick")
    expect_true(grepl("depends on each regimen's overall time", output$early_limit_value$html))
  })
})

test_that("organ limits entered as % of prescription restrict the search when chosen", {
  shiny::testServer(app_dir, {
    set_v3(session, criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 30, tol_late = 30)
    free <- suggestions()
    # cord at 75% of the prescription: QUANTEC 50 Gy EQD2 at alpha/beta 2
    session$setInputs(`oar-v_cord` = 75, `oar-use_in_search` = TRUE)
    lim <- suggestions()
    expect_true(nrow(lim) > 0)
    cord_eqd2 <- eqd2_physical(0.75 * lim$m * lim$d, lim$m, 2)
    expect_true(all(cord_eqd2 <= 50 + 1e-9))
    expect_false(isTRUE(all.equal(free$d, lim$d)))
    expect_true(any(grepl("organ-at-risk limit", lim$next_step)))
    # the export records the organ settings and replays exactly
    sc <- jsonlite::fromJSON(scenario_json(scenario()))
    expect_equal(sc$search$oar$library, "QUANTEC")
    r <- rerun_scenario(sc)
    expect_equal(r$d, lim$d)
  })
})

test_that("schema 3 export: calendar and fitted DD replay exactly; a v2.6 file still replays", {
  shiny::testServer(app_dir, {
    set_v3(session, fx_week = "6", dd_mode = "fitted", alpha_dd = 0.03, c_over_alpha_dd = 100, g_dd = 0.01,
           criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 10)
    sc <- jsonlite::fromJSON(scenario_json(scenario()))
    expect_equal(sc$schema_version, "3")
    expect_equal(sc$dd$mode, "fitted")
    expect_equal(sc$search$calendar$fx_per_week, 6)
    r <- rerun_scenario(sc)
    s <- suggestions()
    expect_equal(r$d, s$d); expect_equal(r$T, s$T); expect_equal(r$bed_dd, s$bed_dd, tolerance = 1e-12)
    other <- sc; other$dd$alpha_dd <- 0.031          # a different fitted value is a different scenario
    expect_false(identical(scenario_id(other), scenario_id(sc)))
    bad <- sc; bad$calibration$status <- "ok"         # a saved status that the inputs do not give is refused
    expect_error(rerun_scenario(bad), "calibration status")
  })
  old <- jsonlite::fromJSON(file.path(app_abs, "tests", "testthat", "fixtures", "v26_candidate_both_worst_p1.json"))
  expect_equal(old$schema_version, "2")
})

test_that("uncertainty module: summary, conditional TCP, search, and a ZIP that replays", {
  shiny::testServer(app_dir, {
    set_v3(session, comparator = "hypno_test", min_gain = -3, d_max = 3.5, m_range = c(15, 19))
    tab <- output$`unc-summary`
    expect_true(grepl("Δ tumor BED, DI, vs 20 × 2.75 Gy in 26 days", tab, fixed = TRUE))
    expect_true(grepl("Monte Carlo standard error", tab, fixed = TRUE))
    expect_true(grepl("Change in tumor control", output$`unc-tcp`, fixed = TRUE))
    session$setInputs(`unc-run` = 1)
    st <- output$`unc-search_table`
    expect_true(grepl("listed", st, fixed = TRUE))
    expect_true(grepl("neighbor", st, fixed = TRUE))
    f <- output$`unc-download`
    td <- tempfile(); dir.create(td); zip::unzip(f, exdir = td)
    expect_true(all(c("scenario.json", "evaluated_regimen_summary.csv", "search_with_neighbors.csv", "README.txt") %in% list.files(td)))
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    expect_equal(sc$uncertainty$wording, "shares of simulated parameter draws under the stated models and priors")
    old <- setwd(app_abs); on.exit(setwd(old))
    again <- rerun_uncertainty(sc)
    saved <- read.csv(file.path(td, "evaluated_regimen_summary.csv"))
    expect_equal(again$di_median, saved$di_median, tolerance = 1e-12)
    expect_equal(again$share_tox, saved$share_tox)
  })
})

test_that("under uncertainty the Fowler and Fenwick early limits are refused with a reason", {
  shiny::testServer(app_dir, {
    set_v3(session, early_mode = "fowler")
    expect_error(output$`unc-summary`, "tied to fixed parameter values")
  })
})

test_that("the chemotherapy box says what it changes, and cisplatin weeks are context only", {
  html <- ui_html()
  expect_true(grepl("lower the Fowler zone only (no BED changes)", html, fixed = TRUE))
  shiny::testServer(app_dir, {
    set_v3(session)
    expect_true(grepl("spans 4 treatment weeks (7 for the comparison)", output$cisplatin_context$html, fixed = TRUE))
  })
})

test_that("compare regimens lists the evaluated regimen, three entries and the benchmarks", {
  shiny::testServer(app_dir, {
    set_v3(session)
    t <- output$compare_table
    expect_true(grepl("Evaluated: 20 × 2.75 Gy in 26 days", t, fixed = TRUE))
    expect_true(grepl("A: 17 × 3.05 Gy in 23 days", t, fixed = TRUE))
    expect_true(grepl("HYHOPE", t, fixed = TRUE))
    expect_true(grepl("Dose per fraction outside the tumor data", t, fixed = TRUE))
  })
})

test_that("verdict words use a 0.5 Gy tolerance", {
  shiny::testServer(app_dir, {
    set_v3(session, m = 35, d = 2.005)                 # +0.37 Gy10 tumor BED
    expect_true(grepl("about the same (within 0.5 Gy) under DI", output$verdicts$html, fixed = TRUE))
  })
})

test_that("the refit preset values come from the draws, as medians of identified quantities", {
  env <- new.env(); old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  di <- read.csv(file.path("draws", "tumor_di_primary.csv"))
  rp <- env$refit_preset()
  expect_equal(rp$lambda_di / rp$alpha, median(di$lambda / di$alpha), tolerance = 1e-9)
  expect_equal(rp$Tk, median(di$Tk), tolerance = 1e-6)
  expect_equal(rp$ab, 10)
})

test_that("README claim: at the version 3 defaults the both-model screen opens empty; HYPNO and -3 Gy restore 17-19 fractions", {
  shiny::testServer(app_dir, {
    set_v3(session)
    expect_equal(nrow(agreement()), 0)
    expect_equal(nrow(suggestions()), 0)
    set_v3(session, comparator = "hypno_test", min_gain = -3)
    expect_equal(agreement()$m, c(17, 18, 19))
  })
})
