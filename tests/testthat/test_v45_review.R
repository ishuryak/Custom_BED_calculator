# Version 4.5.0 (2026-11-09): fixes from the independent review of 4.4.0
# (HNSCC_BED_v44_independent_review_and_action_plan_2026-11-09.md).
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

v45_inputs <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = -5,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE)
pt <- withr::with_dir(app_abs, jsonlite::fromJSON("draws/refit_point_values.json")$values)
refit_vals <- list(alpha = pt$di$alpha, ab = pt$di$ab, lambda_di = pt$di$lambda, lambda_dd = pt$dd$lambda,
                   Tk = pt$di$Tk, g = pt$di$g)
banner_text <- function(output) as.character(output$banner$html)

# --- 4.1 Named tumor profiles -----------------------------------------------------------------

test_that("a named profile label never covers values that differ from it (review 4.1)", {
  shiny::testServer(app_dir, {
    # The refit values with an edited DD maximum rate are not the refit.
    do.call(session$setInputs, modifyList(modifyList(v45_inputs, refit_vals), list(preset = "refit", lambda_dd = 0.6)))
    expect_false(grepl("trial-fit point scenario", banner_text(output), fixed = TRUE))
    expect_match(banner_text(output), "tumor values: custom values", fixed = TRUE)
    # The refit values with an edited g are not the refit either.
    do.call(session$setInputs, modifyList(modifyList(v45_inputs, refit_vals), list(preset = "refit", g = 0.02)))
    expect_match(banner_text(output), "tumor values: custom values", fixed = TRUE)
    # D18 values with the refit's g (what choosing D18 after the refit used to leave) are not D18.
    do.call(session$setInputs, modifyList(v45_inputs, list(g = pt$di$g)))
    expect_match(banner_text(output), "tumor values: custom values", fixed = TRUE)
    # Exact values are named.
    do.call(session$setInputs, v45_inputs)
    expect_match(banner_text(output), "tumor values: scenario D18", fixed = TRUE)
    do.call(session$setInputs, modifyList(modifyList(v45_inputs, refit_vals), list(preset = "refit")))
    expect_match(banner_text(output), "tumor values: trial-fit point scenario", fixed = TRUE)
  })
})

test_that("profile matching uses all six fields for complete profiles and the documented fields for partial ones", {
  d18 <- d18_profile()
  expect_identical(sort(names(d18)), sort(TUMOR_FIELDS))
  expect_identical(profile_matching(d18, refit_vals), "d18")
  expect_identical(profile_matching(refit_vals, refit_vals), "refit")
  # Each field on its own: a small edit of a small value is not averaged away by the large ones.
  v <- refit_vals; v$g <- v$g * 1.001
  expect_identical(profile_matching(v, refit_vals), "custom")
  for (k in TUMOR_FIELDS) {
    v <- d18; v[[k]] <- v[[k]] + 0.01
    expect_identical(profile_matching(v, refit_vals), "custom", info = k)
    v <- refit_vals; v[[k]] <- v[[k]] * 1.1 + 0.001
    expect_identical(profile_matching(v, refit_vals), "custom", info = k)
  }
  # Partial presets keep their documented semantics: g is not part of them.
  hy <- modifyList(d18, resolve_preset("hyhope")); hy$g <- 0.01
  expect_identical(profile_matching(hy, refit_vals), "hyhope")
  # An empty field never matches.
  v <- d18; v$g <- NA_real_
  expect_identical(profile_matching(v, refit_vals), "custom")
})

test_that("the app starts from one complete profile, the trial-fit point scenario", {
  e <- new.env(); withr::with_dir(app_abs, sys.source("app.R", envir = e))
  expect_identical(e$STARTUP$name, "refit")
  for (k in TUMOR_FIELDS) expect_equal(e$DEF[[k]], refit_vals[[k]], tolerance = 1e-12, info = k)
  html <- as.character(htmltools::renderTags(e$ui)$html)
  expect_match(html, "<option value=\"refit\" selected>", fixed = TRUE)
  # The six numeric inputs start at the profile values (no observer needed).
  for (k in TUMOR_FIELDS) {
    tag <- regmatches(html, regexpr(sprintf("<input id=\"%s\" type=\"number\"[^>]*>", k), html))
    v <- as.numeric(sub(".*value=\"([^\"]+)\".*", "\\1", tag))
    expect_equal(v, refit_vals[[k]], tolerance = 1e-9, info = k)
  }
})

test_that("an unreadable or malformed point-profile file starts on D18 and says so", {
  tmp <- withr::local_tempdir()
  expect_false(refit_profile(tmp)$ok)
  st <- startup_profile(refit_profile(tmp))
  expect_identical(st$name, "d18"); expect_match(st$why, "cannot read")
  for (bad in list(list(di = list(alpha = 0.04, ab = 10, lambda = 0.03, Tk = 37), dd = list(lambda = 0.4)),
                   list(di = list(alpha = -1, ab = 10, lambda = 0.03, Tk = 37, g = 0.01),
                        dd = list(alpha = 0.04, lambda = 0.4, c_over_alpha = 60, g = 0.01)))) {
    jsonlite::write_json(list(values = bad), file.path(tmp, "refit_point_values.json"), auto_unbox = TRUE)
    rp <- refit_profile(tmp)
    expect_false(rp$ok); expect_match(rp$why, "missing or invalid")
    expect_identical(startup_profile(rp)$values, d18_profile())
  }
  ok <- withr::with_dir(app_abs, refit_profile())
  expect_true(ok$ok)
  expect_equal(unlist(ok$values[TUMOR_FIELDS]), unlist(refit_vals[TUMOR_FIELDS]), tolerance = 1e-12)
})
