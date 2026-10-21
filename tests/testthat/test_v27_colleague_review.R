# Version 2.7: interface changes from a colleague's comments of 2026-10-20 and
# the independent review of the plan (docs/development/2026_10_20_v27/).
# E1 source link, E2 two-model default screen, E3 "parameter settings" wording,
# E4 a shared caveat line. Equations, parameters, limits and selection rules
# are unchanged; the replay tests below check that against files frozen from
# version 2.6 (docs/development/2026_10_20_v27/make_v26_replay_fixtures.R).

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
source(file.path(app_abs, "R", "bed_models.R"), local = TRUE)
source(file.path(app_abs, "R", "export.R"), local = TRUE)
fixture <- function(f) file.path(app_abs, "tests", "testthat", "fixtures", f)

# The page as built, without a browser: source app.R and render its UI.
ui_html <- function() {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  htmltools::renderTags(env$ui)$html
}
select_block <- function(html, id) {
  m <- regmatches(html, regexpr(sprintf("(?s)<select id=\"%s\"[^>]*>.*?</select>", id), html, perl = TRUE))
  if (length(m) == 0) "" else m
}
selected_value <- function(html, id) {
  sel <- regmatches(select_block(html, id), regexpr("<option value=\"[^\"]+\" selected>", select_block(html, id)))
  if (length(sel) == 0) NA_character_ else sub("<option value=\"([^\"]+)\" selected>", "\\1", sel)
}

# The current opening page, stated in full (testServer has no browser, so the
# UI's selected values must be supplied; the UI test above checks they agree).
# Deliberately separate from the legacy set_defaults() of test_app_server.R.
opening <- list(
  m = 20, d = 2.65, auto_T = TRUE, T = 26, comparator = "custom",
  ref_m = 35, ref_d = 2, ref_auto_T = TRUE, ref_T = 47,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0,
  anchor_m = 35, anchor_d = 2,
  ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_gain = 1, min_days = 7,
  d_step = "0.05", d_max = 5, early_mode = "hypno", late_mode = "reference",
  tol_early = 0, tol_late = 0, agree_cmp = "hypno_test", agree_margin = -3, agree_days = 1)
set_open <- function(session, ...) do.call(session$setInputs, modifyList(opening, list(...)))

# --- A. UI defaults and startup --------------------------------------------

test_that("A: app.R builds its UI (a helper used by the UI must be defined at top level)", {
  expect_error(html <- ui_html(), NA)
  expect_true(nchar(html) > 10000)
})

test_that("E2: the page opens with both models, the smaller change and +1 Gy", {
  html <- ui_html()
  expect_equal(selected_value(html, "criterion"), "both")
  expect_equal(selected_value(html, "rank_by"), "worst")
  mg <- regmatches(html, regexpr("id=\"min_gain\"[^>]*value=\"[^\"]+\"", html))
  expect_equal(as.numeric(sub(".*value=\"([^\"]+)\".*", "\\1", mg)), 1)
  # every choice is still offered, in the same order
  opts <- regmatches(select_block(html, "rank_by"), gregexpr("value=\"[a-z]+\"", select_block(html, "rank_by")))[[1]]
  expect_equal(opts, sprintf("value=\"%s\"", c("best", "worst", "di", "dd")))
})

test_that("E4: one caveat line on each of the Results, Candidate and Model agreement tabs", {
  html <- ui_html()
  ids <- c("caveat-results", "caveat-candidates", "caveat-agreement")
  for (id in ids) expect_equal(lengths(regmatches(html, gregexpr(sprintf("id=\"%s\"", id), html))), 1, info = id)
  expect_equal(lengths(regmatches(html, gregexpr("class=\"bed-caveat", html))), 3)
  pos <- function(pat) regexpr(pat, html, fixed = TRUE)[1]
  # each tab name appears in its nav link and again on its panel: use the panel
  pane <- function(v) max(gregexpr(sprintf("data-value=\"%s\"", v), html, fixed = TRUE)[[1]])
  tabs <- c(pane("Results"), pane("Candidate regimens (model screen)"), pane("Model agreement"), pane("Sensitivity"))
  expect_true(all(tabs > 0))
  cav <- vapply(ids, function(id) pos(sprintf("id=\"%s\"", id)), numeric(1))
  expect_true(cav[1] > tabs[1] && cav[1] < tabs[2])
  expect_true(cav[2] > tabs[2] && cav[2] < tabs[3])
  expect_true(cav[3] > tabs[3] && cav[3] < tabs[4])
  # placed before the long tables, not under them
  expect_true(cav[1] < pos("id=\"benchmarks\""))
  expect_true(cav[3] < pos("id=\"agree_table\""))
  line <- gsub("\\s+", " ", regmatches(html, regexpr("(?s)<div id=\"caveat-results\"[^>]*>.*?</div>", html, perl = TRUE)))
  for (w in c("Chemotherapy effects", "organ-at-risk doses", "individual gap timing",
              "Uncertainty intervals are not calculated"))
    expect_true(grepl(w, line, fixed = TRUE), info = w)
  expect_false(grepl("assumed equal", line, fixed = TRUE))
})

# --- E1. Source link ---------------------------------------------------------

test_that("E1: the Model and sources page links the version 2 source branch", {
  shiny::testServer(app_dir, {
    set_open(session)
    expect_true(grepl("href=\"https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation\"",
                      output$about$html, fixed = TRUE))
  })
})

# --- B. Server behavior at the opening settings and deliberate variants -------

test_that("E2: at the opening settings no candidate qualifies and the message explains it", {
  shiny::testServer(app_dir, {
    set_open(session)
    expect_equal(nrow(suggestions()), 0)
    msg <- output$suggest_msg$html
    expect_true(grepl("under both models", msg, fixed = TRUE))
    expect_true(grepl("valid outcome", msg, fixed = TRUE))
    expect_true(grepl("Model agreement tab explores a separate comparison", msg, fixed = TRUE))
    expect_true(grepl("can include regimens that fail the threshold under the other model", msg, fixed = TRUE))
    expect_false(grepl("only one model supports", msg, fixed = TRUE))
    # no stale rows after a nonempty search
    set_open(session, min_gain = -3); expect_equal(nrow(suggestions()), 3)
    set_open(session); expect_equal(nrow(suggestions()), 0)
  })
})

test_that("E2: an allowed loss of 3 Gy lists three regimens in conservative order", {
  shiny::testServer(app_dir, {
    set_open(session, min_gain = -3)
    s <- suggestions()
    expect_equal(s$m, c(30, 29, 28))
    expect_equal(s$d, c(2.20, 2.25, 2.30), tolerance = 1e-9)
    expect_equal(s$T, c(40, 39, 38))
    expect_equal(s$d_di, c(1.12000000, 1.18839286, 1.12628571), tolerance = 1e-6)
    expect_equal(s$d_dd, c(-2.04954215, -2.40804467, -2.86484061), tolerance = 1e-6)
    expect_true(all(s$d_di >= -3 & s$d_dd >= -3))
    # the ranking matters here: by the larger change the order is by DI
    set_open(session, min_gain = -3, rank_by = "best")
    expect_equal(suggestions()$m, c(29, 28, 30))
  })
})

test_that("E2: the either/larger screen is still available and unchanged", {
  shiny::testServer(app_dir, {
    set_open(session, criterion = "either", rank_by = "best")
    s <- suggestions()
    expect_equal(s$m, c(29, 28, 30))
    expect_equal(s$d_dd, c(-2.40804467, -2.86484061, -2.04954215), tolerance = 1e-6)
  })
})

test_that("E2: with DD unavailable the DD-specific explanation is shown, not the generic one", {
  shiny::testServer(app_dir, {
    set_open(session, lambda_dd = 0)
    msg <- output$suggest_msg$html
    expect_true(grepl("DD model is unavailable", msg, fixed = TRUE))
    expect_false(grepl("valid outcome", msg, fixed = TRUE))
  })
})

test_that("E2: the Model agreement tab keeps its own comparison, threshold and days", {
  shiny::testServer(app_dir, {
    set_open(session)
    a <- agreement()
    expect_equal(a$m, c(17, 18, 19))
    expect_equal(a$d, c(3.05, 2.90, 2.80), tolerance = 1e-9)
    ctx <- output$agree_context$html
    expect_true(grepl("55 Gy", ctx) || grepl("20 × 2.75", ctx))
    expect_true(grepl("at least -3 Gy under both models", ctx, fixed = TRUE))
  })
})

# --- C. Counts, denominators and wording --------------------------------------

test_that("E3: agreement counts are labeled as selected parameter settings", {
  shiny::testServer(app_dir, {
    set_open(session)
    tab <- output$agree_table
    expect_true(grepl("Selected parameter settings with both", tab, fixed = TRUE))
    expect_true(grepl("15 of 15 settings", tab, fixed = TRUE))
    expect_true(grepl("over evaluable settings", tab, fixed = TRUE))
    expect_false(grepl("Scenarios with both", tab, fixed = TRUE))
  })
})

test_that("E3: with invalid settings the denominator is the evaluable count, shown with the invalid count", {
  shiny::testServer(app_dir, {
    set_open(session, lambda_di = 0, ref_auto_T = FALSE, ref_T = 150, agree_cmp = "sidebar")
    a <- agreement()
    skip_if(nrow(a) == 0, "fixture lists no agreement rows")
    tab <- output$agree_table
    i <- which(a$n_invalid > 0)[1]
    expect_true(grepl(sprintf("%d of %d evaluable settings (%d of %d not evaluable)",
                              a$n_pass_both[i], a$n_valid_both[i], a$n_invalid[i], a$n_scenarios[i]),
                      tab, fixed = TRUE))
  })
})

test_that("E3: the Sensitivity summary qualifies what the settings are", {
  shiny::testServer(app_dir, {
    set_open(session)
    txt <- output$sens_summary$html
    expect_true(grepl("evaluable parameter settings under DI", txt, fixed = TRUE))
    expect_true(grepl("not independent confirmations or a probability sample", txt, fixed = TRUE))
    expect_true(grepl("Within each setting, DI and DD share α", txt, fixed = TRUE))
    expect_true(grepl("can differ between settings", txt, fixed = TRUE))
    expect_true(grepl("Uncertainty intervals are not calculated", txt, fixed = TRUE))
    set_open(session, m = 32, d = 2.1, lambda_di = 0, ref_auto_T = FALSE, ref_T = 150)
    expect_true(grepl("12 of 15 settings are not evaluable", output$sens_summary$html, fixed = TRUE))
  })
})

# --- D. Exports and replay -----------------------------------------------------

replay_cases <- c("candidate_either_best_p1", "candidate_both_worst_m3", "candidate_both_worst_p1",
                  "agreement_default")

test_that("D: compatible version 2.6 files replay to the same results", {
  for (nm in replay_cases) {
    sc <- jsonlite::fromJSON(fixture(sprintf("v26_%s.json", nm)))
    expect_equal(sc$app_version, "2.6.0", info = nm)
    old <- jsonlite::fromJSON(fixture(sprintf("v26_%s_results.json", nm)))
    new <- rerun_scenario(sc)
    rows <- as.data.frame(new)
    expect_equal(nrow(rows), if (length(old$rows)) nrow(old$rows) else 0L, info = nm)
    if (nrow(rows) > 0) {
      for (col in names(old$rows)) {
        o <- old$rows[[col]]; n <- rows[[col]]
        if (is.numeric(n)) { o <- suppressWarnings(as.numeric(o)); expect_equal(n, o, tolerance = 1e-10, info = paste(nm, col)) }
        else expect_equal(as.character(n), as.character(o), info = paste(nm, col))
      }
    }
    feas <- as.data.frame(attr(new, "feasible"))
    if (length(old$feasible)) {
      expect_equal(nrow(feas), nrow(old$feasible), info = paste(nm, "feasible"))
      for (col in intersect(c("m", "d", "T", "d_di", "d_dd", "rank_eligible"), names(old$feasible)))
        expect_equal(feas[[col]], old$feasible[[col]], tolerance = 1e-10, info = paste(nm, "feasible", col))
    }
  }
  expect_equal(nrow(rerun_scenario(jsonlite::fromJSON(fixture("v26_candidate_both_worst_p1.json")))), 0)
})

test_that("D: a new export records the current version and the chosen settings (updated for version 3)", {
  sc <- jsonlite::fromJSON(fixture("v26_candidate_both_worst_p1.json"))
  new <- build_scenario(default_params(), list(m = 35, d = 2, T = 47), list(m = 35, d = 2, T = 47),
                        search = sc$search, context = list())
  expect_equal(new$app_version, APP_VERSION)
  expect_equal(new$selection_version, sc$selection_version)
  # version 3 writes schema 3 / equations 3.0 and still replays the 2.6 file exactly
  expect_equal(c(new$schema_version, new$equation_version), c("3", "3.0"))
  expect_equal(c(sc$schema_version, sc$equation_version), c("2", "2.2"))
  expect_equal(new$dd$mode, "calibrated")
  expect_equal(new$search$criterion, "both")
  expect_equal(new$search$rank_by, "worst")
  expect_equal(new$search$min_gain, 1)
  ag <- jsonlite::fromJSON(fixture("v26_agreement_default.json"))
  expect_equal(SCENARIO_SET_VERSION, ag$scenario_set_version)
  # an incompatible older file is still refused, and migration is explicit
  old <- sc; old$selection_version <- "2.5"
  expect_error(rerun_scenario(old), "migrate_scenario")
})

test_that("D: app and DESCRIPTION carry the same version", {
  desc <- read.dcf(file.path(app_abs, "DESCRIPTION"))
  expect_equal(unname(desc[1, "Version"]), APP_VERSION)
  expect_equal(APP_VERSION, "3.0.0")
})
