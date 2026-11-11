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
