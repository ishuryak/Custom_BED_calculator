# Version 4.6: the early-toxicity model switch on the main page, its exports and their replay.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R", "v4_home.R", "v4_export.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

v46_inputs <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = -5,
  preset = "refit", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  home_m_tumor = 5, home_m_early = 5, home_m_late = 5)
pt <- withr::with_dir(app_abs, jsonlite::fromJSON("draws/refit_point_values.json")$values)
v46_inputs <- modifyList(v46_inputs, list(alpha = pt$di$alpha, ab = pt$di$ab, lambda_di = pt$di$lambda,
                                          lambda_dd = pt$dd$lambda, Tk = pt$di$Tk, g = pt$di$g))

FIT6 <- withr::with_dir(app_abs, load_early_fit())

unzip_main6 <- function(f) {
  out <- tempfile(); utils::unzip(f, exdir = out)
  list(sc = jsonlite::fromJSON(file.path(out, "scenario.json"), simplifyVector = TRUE),
       points = utils::read.csv(file.path(out, "points.csv")))
}

test_that("Fowler is the default; the switch changes the main list and shows the research note", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v46_inputs)
    expect_null(input$early_model)                         # nothing chosen yet: params() stays Fowler
    expect_false(early_is_relative(params()))
    session$setInputs(early_model = "fowler")
    expect_equal(nrow(home()$rows), 69)
    expect_null(output$early_model_note)
    session$setInputs(early_model = "fitted")
    expect_true(early_is_relative(params()))
    expect_equal(params()$ab_e, FIT6$point$ab_e)
    expect_equal(nrow(home()$rows), 9)
    expect_equal(sort(unique(home()$rows$m)), 26:30)
    expect_match(output$early_model_note$html, "Research model in use for early toxicity")
    expect_match(output$early_model_note$html, FIT6$model_id, fixed = TRUE)
    session$setInputs(early_model = "fowler")
    expect_equal(nrow(home()$rows), 69)                     # switching back restores the default list
  })
})

test_that("research mode gives the regimen and the comparator their own treatment dates (review R3)", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, list(
      early_model = "fitted", m = 20, d = 2.75, skip_days = "4,5,8", ref_m = 33, ref_d = 2, ref_fx_week = "6")))
    b <- both()
    cal_a <- list(fx_per_week = 5L, start_weekday = 1L, skip_days = c(4L, 5L, 8L))
    cal_b <- list(fx_per_week = 6L, start_weekday = 1L, skip_days = integer())
    expect_equal(regimen()$T, 31); expect_equal(ref()$T, 38)
    expect_equal(b$bed_early[1], compute_all(20, 2.75, 31, params(), anchor(), cal_a)$bed_early, tolerance = 1e-12)
    expect_equal(b$bed_early[2], compute_all(33, 2, 38, params(), anchor(), cal_b)$bed_early, tolerance = 1e-12)
    expect_equal(b$bed_early, c(31.11404545, 29.73698842), tolerance = 1e-6)      # the review's worked values
  })
})

test_that("a research-mode Uncertainty export records each regimen's dates and replays exactly (review R3)", {
  unc <- list(`oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE,
              `unc-variant` = "primary", `unc-tumor_crit` = "di", `unc-tox_source` = "prior", `unc-tcp0` = 70,
              `unc-rule` = "min_q10", `unc-tox_share` = 90, early_mode = "reference", tol_early = 5)
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, c(list(
      early_model = "fitted", m = 20, d = 2.75, skip_days = "4,5,8", ref_m = 33, ref_d = 2, ref_fx_week = "6"), unc)))
    td <- tempfile(); utils::unzip(output$`unc-download`, exdir = td)
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    expect_equal(unlist(sc$early$regimen_calendar$skip_days), c(4, 5, 8))
    expect_equal(sc$early$comparison_calendar$fx_per_week, 6)
    saved <- utils::read.csv(file.path(td, "evaluated_regimen_summary.csv"))
    again <- withr::with_dir(app_abs, rerun_uncertainty(sc))
    for (k in names(saved)) expect_equal(again[[k]], saved[[k]], tolerance = 1e-12, info = k)
    bad <- sc; bad$early$comparison_calendar$fx_per_week <- 5        # dates that cannot give the recorded times
    expect_error(withr::with_dir(app_abs, rerun_uncertainty(bad)), "comparison's recorded calendar cannot give 33 fractions in 38 days")
    bad <- sc; bad$early$regimen_calendar$skip_days <- list()
    expect_error(withr::with_dir(app_abs, rerun_uncertainty(bad)), "regimen's recorded calendar cannot give 20 fractions in 31 days")
  })
})

test_that("a typed overall time the calendar cannot give is explained, not spread over invented gaps (review R3)", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, list(early_model = "fitted", m = 20, d = 2.75, auto_T = FALSE, T = 31)))
    expect_true(is.na(both()$bed_early[1]))
    expect_match(as.character(output$T_shown$html), "end on day 26,\\s+not day 31", perl = TRUE)
    session$setInputs(skip_days = "4,5,8")                       # the gaps entered: dates now give day 31
    expect_false(is.na(both()$bed_early[1]))
    expect_false(grepl("not available", as.character(output$T_shown$html), fixed = TRUE))
    session$setInputs(early_model = "fowler", skip_days = "")    # Fowler mode needs no dates
    expect_false(is.na(both()$bed_early[1]))
    expect_false(grepl("not available", as.character(output$T_shown$html), fixed = TRUE))
  })
})

test_that("no text claims a risk meaning for +5, a universal chemotherapy effect or few post-RT assessments (review R1)", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, list(early_model = "fitted")))
    texts <- c(banner = as.character(output$early_model_note$html),
               brief = paste(readLines(file.path(app_abs, "early_fit_in_brief.md"), encoding = "UTF-8"), collapse = " "),
               about = paste(readLines(file.path(app_abs, "about.md"), encoding = "UTF-8"), collapse = " "))
    for (k in names(texts)) {
      expect_false(grepl("smaller increase in risk|multiplied .*odds|few trials scored", texts[[k]]), info = k)
      expect_match(texts[[k]], "clinical\\s+risk\\s+meaning\\s+is\\s+not\\s+established", info = k)
    }
    expect_match(texts[["brief"]], "not a treatment-specific causal effect", fixed = TRUE)
    expect_match(texts[["brief"]], "only 2 of 22 primary trial families", ignore.case = TRUE)
  })
})

# Review R8: every page and export names the early model and the values actually used.
test_that("the shared banner and Model settings name the active early model and its values", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, list(early_model = "fowler")))
    b <- as.character(output$banner$html)
    expect_match(b, "Early toxicity: Fowler 2003 values (mucosal α/β 10 Gy", fixed = TRUE)
    expect_null(output$early_fit_values)
    session$setInputs(early_model = "fitted", ab_e = 50, Tk_e = 30, Tp_e = 10)    # hidden Fowler fields, edited anyway
    b <- as.character(output$banner$html)
    expect_match(b, sprintf("trial-count fit (research model %s; mucosal α/β %.2f Gy", FIT6$model_id, FIT6$point$ab_e), fixed = TRUE)
    expect_false(grepl("Fowler 2003 values", b, fixed = TRUE))
    v <- as.character(output$early_fit_values$html)
    expect_match(v, sprintf("Mucosal α/β: %.4f Gy", FIT6$point$ab_e), fixed = TRUE)
    expect_match(v, "cannot be edited here", fixed = TRUE)
    expect_equal(params()$ab_e, FIT6$point$ab_e); expect_equal(params()$Tk_e, FIT6$point$Tk_e)
    expect_equal(nrow(home()$rows), 9)                       # the edited Fowler fields change nothing
  })
})

test_that("the Fowler fields are shown only for Fowler; the fitted values only for the fit", {
  e <- new.env(); withr::with_dir(app_abs, sys.source("app.R", envir = e))
  html <- as.character(htmltools::renderTags(e$ui)$html)
  fowler_block <- regmatches(html, regexpr("data-display-if=\"input.early_model != &#39;fitted&#39;\"[^>]*>.*?id=\"Tk_e\"", html))
  expect_length(fowler_block, 1)
  for (id in c("ab_e", "alpha_e", "Tp_e")) expect_match(fowler_block, sprintf("id=\"%s\"", id), fixed = TRUE)
  expect_match(html, "data-display-if=\"input.early_model == &#39;fitted&#39;\"", fixed = TRUE)
})

test_that("Uncertainty names the early source; a research export records which prior entries were used", {
  unc <- list(`oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE,
              `unc-variant` = "primary", `unc-tumor_crit` = "di", `unc-tcp0` = 70,
              `unc-rule` = "min_q10", `unc-tox_share` = 90, early_mode = "reference", tol_early = 5)
  for (src in c("prior", "fixed")) shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, c(list(early_model = "fitted", `unc-tox_source` = src), unc)))
    ctx <- as.character(output$`unc-context`$html)
    if (src == "prior") expect_match(ctx, "its 4,000 joint posterior draws in file order; the Fowler priors are not used", fixed = TRUE)
    else expect_match(ctx, sprintf("mucosal α/β %.2f Gy", FIT6$point$ab_e), fixed = TRUE)
    td <- tempfile(); utils::unzip(output$`unc-download`, exdir = td)
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    pu <- sc$uncertainty$toxicity$prior_use
    expect_match(pu$ab_e, "^not used", info = src)
    expect_match(pu$ab_l, if (src == "prior") "^used" else "^not used", info = src)
    expect_silent(withr::with_dir(app_abs, rerun_uncertainty(sc)))
    bad <- sc; bad$uncertainty$toxicity$prior_use$ab_e <- "used"
    expect_error(withr::with_dir(app_abs, rerun_uncertainty(bad)), "which toxicity prior entries were used")
  })
  shiny::testServer(app_dir, {                                  # Fowler exports are unchanged: no such record
    do.call(session$setInputs, modifyList(v46_inputs, c(list(early_model = "fowler", `unc-tox_source` = "prior"), unc)))
    td <- tempfile(); utils::unzip(output$`unc-download`, exdir = td)
    expect_null(jsonlite::fromJSON(file.path(td, "scenario.json"))$uncertainty$toxicity$prior_use)
    expect_match(as.character(output$`unc-context`$html), "Fowler 2003 model with its parameters drawn from the stated priors", fixed = TRUE)
  })
})

test_that("the Fowler zone limit is refused with the fitted model", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v46_inputs, list(early_model = "fitted", early_mode = "fowler")))
    msg <- tryCatch({ tox_args(); "" }, error = function(e) conditionMessage(e))
    expect_match(msg, "Fowler 2003 zone")
  })
})

test_that("a research-mode main export records the model and replays exactly; another model is refused", {
  withr::with_dir(app_abs, {
    fit <- load_early_fit()
    p <- fitted_early_params(modifyList(default_params(), startup_profile()$values), fit)
    a <- default_anchor(); dom <- tumor_data_support()
    sc <- main_scenario(p, a, NULL, dom); r <- main_compute(p, a, NULL, dom$d_max)
    expect_identical(sc$early$model_id, fit$model_id)
    zf <- tempfile(fileext = ".zip"); id <- write_main_zip(zf, sc, r)
    back <- unzip_main6(zf)$sc
    rr <- main_rerun(back)
    expect_equal(nrow(rr$points), nrow(r$points))
    expect_equal(rr$points$bed_early, r$points$bed_early, tolerance = 1e-9)
    bad <- back; bad$early$model_id <- "early-fit-000000000000"
    expect_error(main_recompute(bad), "not this calculator's")
    # Fowler-mode exports carry no early block and keep their old IDs
    pf <- modifyList(default_params(), startup_profile()$values)
    scf <- main_scenario(pf, a, NULL, dom)
    expect_null(scf$early)
    expect_false(identical(main_id(scf), main_id(sc)))
  })
})

test_that("the main page offers visible links to how each model was fitted and to its trials", {
  e <- new.env(); withr::with_dir(app_abs, sys.source("app.R", envir = e))
  html <- as.character(htmltools::renderTags(e$ui)$html)
  for (id in c("early_how", "early_trials")) expect_match(html, sprintf("id=\"%s\"", id), fixed = TRUE)
  expect_match(html, "How was the early-toxicity model fitted?", fixed = TRUE)
  expect_match(html, "Trials used for the early-toxicity fit", fixed = TRUE)
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v46_inputs)
    hs <- as.character(output$home_settings$html)
    expect_match(hs, "id=\"how_fitted_home\"", fixed = TRUE)                  # tumor method, visible
    expect_match(hs, "id=\"tumor_trials_link\"", fixed = TRUE)                # tumor trials, visible without Details
    expect_silent(session$setInputs(early_how = 1))
    expect_silent(session$setInputs(early_trials = 1))
    expect_silent(session$setInputs(tumor_trials_link = 1))
  })
})

test_that("the early-toxicity trial list is complete, cited and hash-checked; the summary agrees with it", {
  withr::with_dir(app_abs, {
    t <- early_fit_trials(FIT6)
    fitted <- t[t$role == "fitted", ]
    expect_equal(length(unique(fitted$trial)), 22)
    expect_equal(nrow(fitted), 51)
    expect_equal(sum(as.integer(fitted$patients)), 6042)
    expect_true(all(grepl("^[0-9]{6,8}$", t$pmid)))
    expect_equal(sum(t$control_arm == "1"), length(unique(t$trial)))       # one control arm per trial
    dlg <- early_fit_trials_dialog(FIT6)
    expect_true(dlg$ok)
    expect_match(dlg$title, "22 trials, 51 arms", fixed = TRUE)
    brief <- paste(readLines("early_fit_in_brief.md", encoding = "UTF-8"), collapse = " ")
    expect_match(brief, "22 randomized trials", fixed = TRUE)
    expect_match(brief, "51 arms and 6,042 evaluable patients", fixed = TRUE)
    expect_false(grepl("$", brief, fixed = TRUE))                            # no unrendered LaTeX in the modal
    tmp <- tempfile(); dir.create(tmp)
    file.copy(file.path("draws", c(EARLY_FIT_MANIFEST, EARLY_FIT_FILE, EARLY_FIT_TRIALS_FILE)), tmp)
    cat("x\n", file = file.path(tmp, EARLY_FIT_TRIALS_FILE), append = TRUE)
    bad <- early_fit_trials_dialog(FIT6, tmp)
    expect_false(bad$ok)
    expect_match(as.character(bad$body), "does not match", fixed = TRUE)
  })
})
