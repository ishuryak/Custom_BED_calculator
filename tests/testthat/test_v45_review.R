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

# --- 6.2 / 6.3 Main export consistency and exact replay --------------------------------------

for (f in c("oar.R", "export.R", "v4_export.R")) source(file.path(app_abs, "R", f), local = TRUE)
unzip_main <- function(f) {
  out <- tempfile(); utils::unzip(f, exdir = out)
  list(sc = jsonlite::fromJSON(file.path(out, "scenario.json"), simplifyVector = TRUE),
       points = utils::read.csv(file.path(out, "points.csv")), name = basename(f))
}

test_that("a Main ZIP agrees with its own scenario after a draw-loading failure recovers (review 6.2)", {
  shiny::testServer(app_dir, {
    cache <- get("V4_DRAWS_CACHE", envir = environment(go_to))   # the app's cache, not this file's copy
    rm(list = ls(cache), envir = cache)
    cache$failed_at <- Sys.time()            # a failure just now: no support on the page
    do.call(session$setInputs, v45_inputs)
    tb <- output$home_table
    cache$failed_at <- Sys.time() - 3600     # the retry interval has passed; no input is edited
    z <- unzip_main(output$home_download)
    has_support <- "support" %in% names(z$points)
    expect_identical(z$sc$support$status, "unavailable")   # the page state at the time of the failure
    expect_identical(has_support, identical(z$sc$support$status, "available"))
    expect_equal(nrow(z$points), 60)
    r <- withr::with_dir(app_abs, main_rerun(z$sc))
    expect_identical("support" %in% names(r$points), has_support)
    expect_equal(nrow(r$points), nrow(z$points))
    # After the retry interval the page and the export recover together, with support.
    session$elapse((get("V4_DRAWS_RETRY_SEC", envir = environment(go_to)) + 2) * 1000)
    z2 <- unzip_main(output$home_download)
    expect_identical(z2$sc$support$status, "available")
    expect_true("support" %in% names(z2$points))
    expect_false(identical(z2$sc$analysis_id, z$sc$analysis_id))   # a new, identified analysis
    r2 <- withr::with_dir(app_abs, main_rerun(z2$sc))
    expect_equal(r2$points$support, z2$points$support, tolerance = 1e-12)
  })
})

test_that("exact replay needs the exported ID; an intentional edit is a new analysis (review 6.3)", {
  withr::with_dir(app_abs, {
    p <- default_params(); a <- default_anchor(); dom <- tumor_data_support()
    dr <- list(di = load_draws("di"), dd = load_draws("dd"))
    sc <- main_scenario(p, a, dr, dom); r <- main_compute(p, a, dr, dom$d_max)
    zf <- tempfile(fileext = ".zip"); id <- write_main_zip(zf, sc, r)
    back <- unzip_main(zf)$sc
    expect_equal(nrow(main_rerun(back)$points), 60)                     # intact file replays
    bad <- back; bad$margins$early <- 0; bad$support$margin <- bad$margins$tumor
    expect_error(main_rerun(bad), "does not match its contents")       # edited, old ID
    bad$analysis_id <- NULL
    expect_error(main_rerun(bad), "no analysis ID")                    # edited, ID removed
    bad2 <- back; bad2$analysis_id <- ""
    expect_error(main_rerun(bad2), "no analysis ID")                   # malformed ID
    rr <- main_recompute(bad)                                          # intentional edit
    expect_equal(nrow(rr$points), 18)
    expect_false(identical(attr(rr, "analysis_id"), id))
  })
})

# --- 5: messages follow the actual DD mode and model validity --------------------------------

test_that("Sensitivity never calls a directly entered DD calibrated; Main does not invent a DD loss", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, modifyList(v45_inputs, list(dd_mode = "fitted", alpha_dd = 0.3, lambda_dd = 0.5,
                                                            c_over_alpha_dd = 30, g_dd = 0.01)))
    ctx <- as.character(output$sens_context$html); sm <- as.character(output$sens_summary$html)
    expect_false(grepl("calibrated at", ctx, fixed = TRUE))
    expect_match(ctx, "entered directly, not calibrated", fixed = TRUE)
    expect_false(grepl("DD is calibrated to DI", sm, fixed = TRUE))
    expect_match(sm, "is not calibrated to DI", fixed = TRUE)
    do.call(session$setInputs, v45_inputs)
    expect_match(as.character(output$sens_summary$html), "DD is calibrated to DI at the anchor", fixed = TRUE)
    expect_match(as.character(output$sens_context$html), "DD calibrated at the anchor", fixed = TRUE)
    # DD unavailable (zero DD rate): no claim that the other model estimates a larger loss.
    do.call(session$setInputs, modifyList(v45_inputs, list(lambda_dd = 0)))
    fn <- as.character(output$home_footnotes$html)
    expect_false(grepl("the other estimates a larger loss", fn, fixed = TRUE))
    do.call(session$setInputs, v45_inputs)
    expect_match(as.character(output$home_footnotes$html), "the other estimates a larger loss", fixed = TRUE)
  })
})

# --- 5: six destinations; every former page's inputs and outputs still mounted once -----------

test_that("the menu has six destinations including the main page, and every former task keeps its ids", {
  e <- new.env(); withr::with_dir(app_abs, sys.source("app.R", envir = e))
  html <- as.character(htmltools::renderTags(e$ui)$html)
  expect_identical(names(e$V4_MENU), c("Check a regimen", "Custom search (advanced)", "Uncertainty", "Organs at risk",
                                       "Settings and sources"))
  expect_equal(length(e$V4_MENU) + 1, 6)            # plus "Suggested regimens" (home)
  links <- regmatches(html, gregexpr("id=\"go_[a-z]+\"", html))[[1]]
  expect_setequal(links, sprintf("id=\"go_%s\"", c("home", unname(e$V4_MENU))))
  ids <- c("report", "results", "verdicts", "acute_zone", "benchmarks", "compare_table", "cmp1_m", "cmp2_d", "cmp3_m",
           "suggest_table", "download", "agree_table", "download_agree", "agree_cmp", "agree_margin", "sens_table",
           "sens_tox", "sens_summary", "preset", "alpha", "g", "dd_mode", "anchor_m", "about", "home_table")
  for (id in ids) expect_equal(lengths(regmatches(html, gregexpr(sprintf("id=\"%s\"", id), html))), 1, info = id)
})
