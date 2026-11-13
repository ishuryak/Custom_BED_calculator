# Version 4.6.1 (review of 2026-11-12, M2-M5): one page "How the models were fitted (Bayesian)",
# reached from the menu, Details and every "How was ... fitted?" link, built from how_fitted.md and
# the methods record, and failing closed when the record does not describe the deployed draws.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R", "v4_home.R", "v4_export.R",
            "methods_page.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

v461_inputs <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = -5,
  preset = "refit", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  home_m_tumor = 5, home_m_early = 5, home_m_late = 5)

rec <- withr::with_dir(app_abs, load_model_methods())

test_that("the methods record loads and describes the draws the calculator uses", {
  expect_true(rec$ok)
  man <- jsonlite::fromJSON(file.path(app_abs, "draws", "draws_manifest.json"), simplifyVector = FALSE)
  expect_identical(rec$m$tumor$di$draws_sha256, man$files$tumor_di_primary.csv$sha256)
  expect_identical(rec$m$late$status, "not fitted")
})

test_that("every route reaches the one methods page, at its own section", {
  e <- new.env(); withr::with_dir(app_abs, sys.source("app.R", envir = e))
  html <- as.character(htmltools::renderTags(e$ui)$html)
  expect_identical(e$V4_MENU[["How the models were fitted (Bayesian)"]], "fitted")
  expect_match(html, "id=\"go_fitted\"", fixed = TRUE)
  expect_match(html, sprintf("data-value=\"%s\"", FITTED_PAGE), fixed = TRUE)   # the value go_to() selects
  expect_equal(lengths(regmatches(html, gregexpr("id=\"fitted_page\"", html))), 1)
  for (id in c("home_how_fitted", "how_fitted_unc", "how_fitted_settings", "early_how"))
    expect_match(html, sprintf("id=\"%s\"", id), fixed = TRUE, info = id)
  # Stated here, not read from the code: a tumor link must open the tumor section, the early link
  # the early section, and the Details button the page's start.
  want <- c(how_fitted_home = "fitted-tumor", how_fitted_unc = "fitted-tumor", how_fitted_settings = "fitted-tumor",
            early_how = "fitted-early", home_how_fitted = "fitted-bayes")
  expect_setequal(names(FITTED_LINKS), names(want))
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v461_inputs)
    page <- as.character(output$fitted_page$html)
    for (id in names(want)) {
      fitted_route(NULL)
      do.call(session$setInputs, setNames(list(1), id))
      expect_identical(fitted_route(), want[[id]], info = id)
      expect_match(page, sprintf("id=\"%s\"", want[[id]]), fixed = TRUE, info = id)   # the target exists
    }
  })
})

test_that("the page shows the active priors and settings from the record, with no placeholder left", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v461_inputs)
    page <- as.character(output$fitted_page$html)
    expect_false(grepl("{{", page, fixed = TRUE))
    for (p in c(rec$m$tumor$di$priors, rec$m$tumor$dd$priors, rec$m$early$priors))
      expect_match(page, htmltools::htmlEscape(p$prior), fixed = TRUE, info = p$parameter)
    expect_match(page, "lognormal(log 10, 0.354)", fixed = TRUE)
    expect_match(page, "soft penalty favoring onset at or after day 7", fixed = TRUE)
    expect_match(page, "1,000 distinct DI samples, each repeated 4 times", fixed = TRUE)
    expect_match(page, "fixed weak normal prior", fixed = TRUE)
    expect_match(page, "estimated common spread", fixed = TRUE)
    expect_match(page, "Late toxicity:</strong> not fitted", fixed = TRUE)
    expect_match(page, "central 95% 2–4.5 Gy", fixed = TRUE)
    expect_match(page, "Differences between randomized arms", fixed = TRUE)        # tumor brief embedded
    expect_match(page, "Bayesian hierarchical logistic regression", fixed = TRUE)  # early brief embedded
    expect_match(page, "id=\"download_methods\"", fixed = TRUE)
    expect_false(grepl("<h3>How the tumor model was fitted</h3>", page, fixed = TRUE))   # no nested duplicate heading
  })
})

test_that("the early summary's stated priors agree with the record", {
  brief <- paste(readLines(file.path(app_abs, "early_fit_in_brief.md"), encoding = "UTF-8"), collapse = " ")
  pr <- setNames(vapply(rec$m$early$priors, `[[`, "", "prior"), vapply(rec$m$early$priors, `[[`, "", "parameter"))
  expect_identical(pr[["Mucosal α/β"]], "lognormal(log 10, 0.354)"); expect_match(brief, "lognormal centered on 10 Gy", fixed = TRUE)
  expect_identical(pr[["Time-loss rate r"]], "lognormal(log 0.8, 0.5)"); expect_match(brief, "lognormal centered on 0.8", fixed = TRUE)
  expect_identical(pr[["Slope per 10 Gy of early BED"]], "half-normal(1)"); expect_match(brief, "Slope: half-normal(1)", fixed = TRUE)
  expect_identical(pr[["Planned-chemotherapy coefficient"]], "normal(0, 1)"); expect_match(brief, "Chemotherapy adjustment: normal(0, 1)", fixed = TRUE)
  expect_identical(pr[["Spread of trial baselines"]], "half-normal(1)"); expect_match(brief, "spread of baselines: half-normal(1)", fixed = TRUE)
  expect_identical(pr[["Spread of arm-level effects"]], "half-normal(0.5)"); expect_match(brief, "Arm-level spread: half-normal(0.5)", fixed = TRUE)
})

test_that("a record that does not match the deployed draws is refused, and the page says so", {
  tmp <- withr::local_tempdir()
  file.copy(list.files(file.path(app_abs, "draws"), full.names = TRUE, pattern = "\\.(json|csv)$"), tmp)
  expect_true(load_model_methods(tmp)$ok)
  m <- jsonlite::fromJSON(file.path(tmp, "model_methods.json"), simplifyVector = FALSE)
  m$tumor$di$draws_sha256 <- strrep("0", 64)
  jsonlite::write_json(m, file.path(tmp, "model_methods.json"), auto_unbox = TRUE)
  bad <- load_model_methods(tmp)
  expect_false(bad$ok); expect_match(bad$why, "other fitted draws")
  writeBin(raw(0), file.path(tmp, "model_methods.json"))
  expect_match(load_model_methods(tmp)$why, "missing")
  txt <- withr::with_dir(app_abs, methods_page_text(bad))
  expect_match(txt, "technical details of the fits cannot be shown", fixed = TRUE)
  expect_false(grepl("{{", txt, fixed = TRUE))
})

test_that("a placeholder the record cannot fill stops the page instead of showing braces", {
  page <- withr::local_tempfile(fileext = ".md")
  writeLines(c("Text {{TUMOR_ARMS}} and {{NOT_A_FIELD}}"), page)
  expect_error(withr::with_dir(app_abs, methods_page_text(rec, page = page)), "unfilled placeholders.*NOT_A_FIELD")
})
