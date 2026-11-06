# Review of version 4.2.1 (HNSCC_BED_v4_review_and_action_plan_2026-10-31), findings
# F01-F05 and F12. Each test was run against 4.2.1 (917825b) and failed there.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

p0 <- default_params()
a0 <- default_anchor()
# BED unit labels as the page writes them (app.R bed_unit), written out here.
bed_unit <- function(ab) if (ab == 10) "Gy₁₀" else if (ab == 3) "Gy₃" else sprintf("Gy (α/β = %g Gy)", ab)
row_of <- function(h, m) h$rows[h$rows$m == m, , drop = FALSE]

# --- F01: range explanation without assumed DD monotonicity -----------------------
test_that("F01: the default list and its worked example are unchanged", {
  h <- v4_regimens(p0)
  expect_equal(nrow(h$rows), 60)
  expect_equal(sum(h$rows$qual_di & h$rows$qual_dd), 17)
  expect_equal(sum(h$rows$qual_di & !h$rows$qual_dd), 43)
  expect_equal(sum(!h$rows$qual_di & h$rows$qual_dd), 0)
  expect_equal(v4_summary(h$rows)$m, 30:16)
  ex <- v4_range_example(row_of(h, 30), p0, a0, h$standard)
  expect_match(ex, "with 30 fractions the doses 2.10–2.25 Gy", fixed = TRUE)
  expect_match(ex, "2.05 Gy fails (tumor BED", fixed = TRUE)
  expect_match(ex, "so the row shows −3.2 to +3.3", fixed = TRUE)
})

test_that("F01: a non-monotone DD row names the doses of its extremes", {
  p <- modifyList(p0, list(lambda_dd = 2))
  h <- v4_regimens(p); x <- row_of(h, 30)
  expect_equal(x$d, c(2.05, 2.10, 2.15, 2.20, 2.25))
  expect_equal(x$d_dd, c(2.515833, 3.506194, 2.837178, 1.946087, 1.087580), tolerance = 1e-6)
  ex <- v4_range_example(x, p, a0, h$standard)
  expect_match(ex, "tumor DD change runs from +1.1 Gy (at 2.25 Gy) to +3.5 Gy (at 2.10 Gy)", fixed = TRUE)
  # The summary column still shows the true extremes.
  expect_identical(v4_summary(h$rows)$d_dd[1], "+1.1 to +3.5")
})

test_that("F01: an unusable reference never credits a schedule in the explanation", {
  p <- modifyList(p0, list(lambda_di = 2, dd_mode = "fitted", alpha_dd = 0.35, lambda_dd = 0.4,
                           c_over_alpha_dd = 50, g_dd = 0))
  h <- v4_regimens(p)
  expect_equal(h$standard$bed_di, -24.57143, tolerance = 1e-5)
  expect_equal(h$standard$bed_dd, 73.88970, tolerance = 1e-5)
  x <- row_of(h, 30)
  expect_equal(x$d, c(2.15, 2.20, 2.25))
  f <- v4_failed_conditions(30, 2.10, p, a0, h$standard)
  expect_length(f, 1)
  expect_match(f, "DI BED of the standard not positive", fixed = TRUE)
  expect_match(f, "DD change −5.3 Gy", fixed = TRUE)
  ex <- v4_range_example(x, p, a0, h$standard)
  expect_false(grepl("fails ()", ex, fixed = TRUE))
  expect_false(grepl("n/a", ex, fixed = TRUE))
  expect_match(ex, "tumor DI change is not available", fixed = TRUE)
})

test_that("F01: unavailable DD and a negative candidate are named, not left blank", {
  p <- modifyList(p0, list(lambda_dd = 0))  # calibrated DD infeasible
  h <- v4_regimens(p)
  f <- v4_failed_conditions(30, 2.05, p, a0, h$standard)
  expect_match(f[1], "DD model not available", fixed = TRUE)
  e <- tumor_endpoint(-1, 70)
  expect_identical(v4_endpoint_words(e, "DI"), "DI BED not positive")
  # Direct DD and a custom anchor still give a full explanation.
  pf <- modifyList(p0, list(dd_mode = "fitted", alpha_dd = 0.3, lambda_dd = 0.5, c_over_alpha_dd = 30, g_dd = 0))
  hf <- v4_regimens(pf)
  if (nrow(hf$rows)) expect_false(grepl("()", v4_range_example(row_of(hf, max(hf$rows$m)), pf, a0, hf$standard), fixed = TRUE))
  an <- list(m = 20, d = 2.75, T = 26)
  ha <- v4_regimens(p0, anchor = an)
  expect_gt(nrow(ha$rows), 0)
  expect_false(grepl("()", v4_range_example(row_of(ha, max(ha$rows$m)), p0, an, ha$standard), fixed = TRUE))
})

test_that("F01: a dose range with a gap is shown as separate runs", {
  expect_identical(v4_range(c(2.10, 2.15, 2.20, 2.40, 2.45), step = 0.05), "2.10–2.20, 2.40–2.45")
  expect_identical(v4_range(c(2.10, 2.15, 2.20), step = 0.05), "2.10–2.20")
  expect_identical(v4_range(3.3, step = 0.05), "3.30")
  x <- data.frame(m = 20, d = c(2.65, 2.70, 2.90), T = 26, d_di = 0, d_dd = 0, d_early = 0, d_late = 0,
                  bed_early = 50, qual_di = TRUE, qual_dd = FALSE)
  expect_identical(v4_summary(x)$dose, "2.65–2.70, 2.90")
  expect_identical(v4_summary(x)$total, "53.0–54.0, 58.0")
})

# --- F02: units, anchor and availability follow the settings -------------------------
test_that("F02: the tumor margin's EQD2 equivalent follows alpha/beta", {
  expect_equal(v4_margin_eqd2(5, 10), 4.1667, tolerance = 1e-4)
  expect_equal(v4_margin_eqd2(5, 2.5), 2.7778, tolerance = 1e-4)
})

test_that("F02: the DD model description names the actual anchor and mode", {
  rel <- v4_model_relation(p0, a0)
  expect_match(rel, "35 × 2 Gy in 47 days (the standard)", fixed = TRUE)
  an <- list(m = 20, d = 2.75, T = 26)
  r0 <- compute_all(35, 2, 47, p0, an)
  expect_equal(c(r0$bed_di, r0$bed_dd), c(71.51429, 81.49437), tolerance = 1e-5)
  rel <- v4_model_relation(p0, an)
  expect_match(rel, "20 × 2.75 Gy in 26 days (not the standard", fixed = TRUE)
  pf <- modifyList(p0, list(dd_mode = "fitted", alpha_dd = 0.3, lambda_dd = 0.5, c_over_alpha_dd = 30, g_dd = 0.01))
  rel <- v4_model_relation(pf, a0)
  expect_match(rel, "is not calibrated to DI", fixed = TRUE)
  expect_false(grepl("shares", rel, fixed = TRUE))
  expect_match(v4_model_relation(modifyList(p0, list(lambda_di = 0)), a0), "rate is zero", fixed = TRUE)
  expect_match(v4_model_relation(modifyList(p0, list(lambda_dd = 0)), a0), "not available", fixed = TRUE)
})

test_that("F02: a finite but negative standard DD BED is unusable", {
  p <- modifyList(p0, list(dd_mode = "fitted", alpha_dd = 0.005, lambda_dd = 5, c_over_alpha_dd = 10, g_dd = 0))
  r0 <- compute_all(35, 2, 47, p, a0)
  expect_true(is.finite(r0$bed_dd) && r0$bed_dd < 0)
  us <- v4_usable_models(p, a0, r0)
  expect_true(us$di); expect_false(us$dd); expect_true(us$dd_available)
  h <- v4_regimens(p)
  expect_false(any(h$rows$qual_dd))
})

# --- F05: transfer validation -------------------------------------------------------
test_that("F05: only a listed point of the current list can be transferred", {
  h <- v4_regimens(p0); s <- v4_summary(h$rows); x <- h$rows
  i20 <- which(s$m == 20)
  expect_equal(i20, 11)
  expect_equal(v4_validate_transfer(i20, "2.75", s, x), list(m = 20, d = 2.75))
  expect_equal(v4_validate_transfer(11L, 2.65, s, x), list(m = 20, d = 2.65))
  for (i in list(NULL, integer(), 0, -1, 1.5, NA, Inf, 16, "abc", c(1, 2)))
    expect_null(v4_validate_transfer(i, "2.75", s, x), info = deparse(i))
  for (d in list(NULL, "", "abc", NA, Inf, 8, 2.85, 2.60, c(2.65, 2.70)))
    expect_null(v4_validate_transfer(i20, d, s, x), info = deparse(d))
  expect_null(v4_validate_transfer(1, "2.75", s[0, ], x[0, ]))
})

# --- F12: support draws are not stuck on a first failure ---------------------------
test_that("F12: a failed draw load is retried after a bounded wait; a success is cached", {
  cache <- new.env(); calls <- 0
  bad <- function() { calls <<- calls + 1; stop("draws unreadable") }
  good <- function() { calls <<- calls + 1; list(di = 1, dd = 2) }
  t0 <- as.POSIXct("2026-10-31 12:00:00", tz = "UTC")
  expect_null(v4_draws(bad, now = t0, retry_sec = 300, cache = cache))
  expect_identical(cache$error, "draws unreadable")
  expect_null(v4_draws(good, now = t0 + 10, retry_sec = 300, cache = cache))  # inside the wait: no retry
  expect_equal(calls, 1)
  expect_identical(v4_draws(good, now = t0 + 301, retry_sec = 300, cache = cache), list(di = 1, dd = 2))
  expect_equal(calls, 2)
  expect_identical(v4_draws(bad, now = t0 + 400, retry_sec = 300, cache = cache), list(di = 1, dd = 2))
  expect_equal(calls, 2)  # success kept, loader not called again
})

# --- page: F01-F05 together ---------------------------------------------------------
v4_defaults <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE)

test_that("F02/F03: the main page text follows the settings and quotes HYPNO correctly", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    intro <- as.character(output$home_intro$html)
    expect_match(intro, "at most 5 Gy₁₀ below the standard 35 × 2 Gy (about 4.2 Gy in 2-Gy fractions", fixed = TRUE)
    expect_match(intro, "BED at most 5 Gy₃ above the standard", fixed = TRUE)
    for (b in c(2.5, 5, 7, 10, 15)) {
      session$setInputs(ab = b)
      intro <- as.character(output$home_intro$html)
      expect_match(intro, sprintf("at most 5 %s below", bed_unit(b)), fixed = TRUE, info = b)
      expect_match(intro, sprintf("about %.1f Gy in 2-Gy fractions", 5 / (1 + 2 / b)), fixed = TRUE, info = b)
      if (nrow(home_rows()) > 0) expect_match(output$home_table, sprintf("Tumor BED change, DI / DD (%s)", bed_unit(b)), fixed = TRUE)
    }
    session$setInputs(ab = 10)
    for (b in c(5, 8)) { session$setInputs(ab_e = b)
      expect_match(as.character(output$home_intro$html), sprintf("(acute mucositis): BED at most 5 %s above", bed_unit(b)), fixed = TRUE) }
    session$setInputs(ab_e = 10)
    for (b in c(2, 4)) { session$setInputs(ab_l = b)
      expect_match(as.character(output$home_intro$html), sprintf("BED at most 5 %s above the standard, at the prescription", bed_unit(b)), fixed = TRUE) }
    session$setInputs(ab_l = 3, home_details = TRUE)
    ft <- as.character(output$home_footnotes$html)
    expect_false(grepl("every BED rises with dose", ft, fixed = TRUE))
    expect_match(ft, "The DD tumor BED need not", fixed = TRUE)
    expect_match(ft, "35 × 2 Gy in 47 days (the standard)", fixed = TRUE)
    expect_match(ft, "190 of the 373 evaluable", fixed = TRUE)
    expect_false(grepl("of 395 patients", ft, fixed = TRUE))
    expect_match(ft, "66 Gy in 33 fractions given 6 a week", fixed = TRUE)
    # Custom anchor: the footnote names it.
    session$setInputs(anchor_m = 20, anchor_d = 2.75)
    expect_match(as.character(output$home_footnotes$html), "20 × 2.75 Gy in 26 days (not the standard", fixed = TRUE)
    session$setInputs(anchor_m = 35, anchor_d = 2)
    # Direct DD with a negative standard BED: announced as unusable.
    session$setInputs(dd_mode = "fitted", alpha_dd = 0.005, lambda_dd = 5, c_over_alpha_dd = 10, g_dd = 0)
    expect_match(as.character(output$home_msg$html), "DD tumor BED is not positive", fixed = TRUE)
  })
})

test_that("F04/F05: the transfer states its scope and refuses a dose outside the row", {
  shiny::testServer(app_dir, {
    sent <- character()
    session$sendInputMessage <- function(inputId, message) sent <<- c(sent, inputId)
    do.call(session$setInputs, v4_defaults)
    session$setInputs(home_table_rows_selected = 11)
    det <- as.character(output$home_detail$html)
    expect_match(det, "keep their own additional criteria", fixed = TRUE)
    expect_false(grepl("and the other menu pages", det, fixed = TRUE))
    # A dose not offered by the row: nothing is sent.
    session$setInputs(home_open_dose = "8", home_open_report = 1)
    expect_false(any(c("m", "d", "ref_m", "min_gain") %in% sent))
    # A listed dose: the regimen and the main-page rules are sent, nothing else.
    session$setInputs(home_open_dose = "2.75", home_open_report = 2)
    expect_true(all(c("m", "d", "ref_m", "ref_d", "min_gain", "tol_early", "tol_late") %in% sent))
    expect_false(any(grepl("^agree_|^oar-|^criterion$|^rank_by$|^m_range$|^d_max$", sent)))
    # A stale row index after the list shrinks: refused.
    sent <- character()
    session$setInputs(home_table_rows_selected = 15, home_open_dose = "3.20")
    session$setInputs(ab_l = 2)  # the list changes
    if (nrow(home_rows()) < 15 || !any(abs(home_points()$d[home_points()$m == home_rows()$m[15]] - 3.2) < 1e-6)) {
      session$setInputs(home_open_report = 3)
      expect_false(any(c("m", "d") %in% sent))
    }
  })
})

# --- F08: exact export and replay of the main list -----------------------------------
for (f in c("oar.R", "export.R", "v4_export.R")) source(file.path(app_abs, "R", f), local = TRUE)

main_roundtrip <- function(p, anchor = a0, draws, context = list()) {
  old <- setwd(app_abs); on.exit(setwd(old))
  dom <- tumor_data_support()
  sc <- main_scenario(p, anchor, draws, dom, context = context)
  r <- main_compute(p, anchor, draws, dom$d_max)
  z <- tempfile(fileext = ".zip"); out <- tempfile("main_unzip_")
  write_main_zip(z, sc, r)
  utils::unzip(z, exdir = out)
  back <- jsonlite::fromJSON(file.path(out, "scenario.json"), simplifyVector = TRUE)
  list(sc = sc, r = r, dir = out, back = back,
       points = utils::read.csv(file.path(out, "points.csv")), summary = utils::read.csv(file.path(out, "summary.csv")))
}
csv_back <- function(x) { f <- tempfile(fileext = ".csv"); utils::write.csv(x, f, row.names = FALSE, na = "NA"); utils::read.csv(f) }
expect_same_table <- function(a, b) {
  expect_identical(names(a), names(b)); expect_equal(nrow(a), nrow(b))
  for (k in names(a)) if (is.numeric(a[[k]])) expect_equal(a[[k]], b[[k]], tolerance = 1e-12, info = k) else
    expect_identical(as.character(a[[k]]), as.character(b[[k]]), info = k)
}
draws0 <- local({ old <- setwd(app_abs); on.exit(setwd(old)); list(di = load_draws("di"), dd = load_draws("dd")) })

test_that("F08: the default main list is exported in full and replays in a fresh R process", {
  rt <- main_roundtrip(p0, draws = draws0)
  expect_equal(nrow(rt$points), 60); expect_equal(nrow(rt$summary), 15)
  expect_true("support" %in% names(rt$points))
  # Same support as the page: 20 x 2.70 Gy is exactly 80%.
  expect_equal(rt$points$support[rt$points$m == 20 & abs(rt$points$d - 2.70) < 1e-9], 0.8)
  expect_identical(rt$back$support$status, "available")
  script <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               sprintf("setwd(%s)", deparse(app_abs)),
               "for (f in c('bed_models','schedule','uncertainty','oar','export','v4_home','v4_export')) source(sprintf('R/%s.R', f))",
               sprintf("sc <- jsonlite::fromJSON(%s, simplifyVector = TRUE)", deparse(file.path(rt$dir, "scenario.json"))),
               "r <- main_rerun(sc)",
               sprintf("write.csv(r$points, %s, row.names = FALSE, na = 'NA')", deparse(file.path(rt$dir, "re_points.csv"))),
               sprintf("write.csv(r$summary, %s, row.names = FALSE, na = 'NA')", deparse(file.path(rt$dir, "re_summary.csv")))),
             script)
  rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  st <- suppressWarnings(system2(rs, shQuote(script), stdout = TRUE, stderr = TRUE))
  expect_null(attr(st, "status"), info = paste(st, collapse = "\n"))
  expect_same_table(rt$points, utils::read.csv(file.path(rt$dir, "re_points.csv")))
  expect_same_table(rt$summary, utils::read.csv(file.path(rt$dir, "re_summary.csv")))
})

test_that("F08: nondefault, empty and unavailable cases replay exactly", {
  cases <- list(
    custom_anchor = list(p = p0, anchor = list(m = 20, d = 2.75, T = 26)),
    direct_dd = list(p = modifyList(p0, list(dd_mode = "fitted", alpha_dd = 0.3, lambda_dd = 0.5, c_over_alpha_dd = 30, g_dd = 0.01)), anchor = a0),
    alpha_beta_clock = list(p = modifyList(p0, list(ab = 7, ab_e = 8, ab_l = 4, clock_e = "inclusive")), anchor = a0),
    empty = list(p = modifyList(p0, list(alpha = 0.01, lambda_di = 0.5, lambda_dd = 0.5)), anchor = a0),
    dd_unavailable = list(p = modifyList(p0, list(lambda_dd = 0)), anchor = a0))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    rt <- main_roundtrip(cs$p, cs$anchor, draws0)
    old <- setwd(app_abs); r <- main_rerun(rt$back); setwd(old)
    expect_same_table(rt$points, csv_back(r$points))
    expect_same_table(rt$summary, csv_back(r$summary))
    if (nm == "empty") { expect_equal(nrow(rt$points), 0); expect_identical(rt$back$result$empty_reason, "negative_reference_both") }
  }
})

test_that("F08: unavailable support is recorded and replayed as unavailable, never recomputed", {
  rt <- main_roundtrip(p0, draws = NULL)
  expect_identical(rt$back$support$status, "unavailable")
  expect_false("support" %in% names(rt$points))
  old <- setwd(app_abs); on.exit(setwd(old))
  r <- main_rerun(rt$back)
  expect_false("support" %in% names(r$points))
  expect_same_table(rt$summary, csv_back(r$summary))
})

test_that("F08: the analysis ID follows scientific inputs only; mismatches are refused", {
  a <- main_roundtrip(p0, draws = draws0, context = list(details = FALSE))
  b <- main_roundtrip(p0, draws = draws0, context = list(details = TRUE))
  expect_identical(main_id(a$sc), main_id(b$sc))
  c2 <- main_roundtrip(modifyList(p0, list(ab_l = 2)), draws = draws0)
  expect_false(identical(main_id(a$sc), main_id(c2$sc)))
  old <- setwd(app_abs); on.exit(setwd(old))
  bad <- a$back; bad$analysis_id <- NULL; bad$support$files$di$sha256 <- strrep("0", 64)
  expect_error(main_rerun(bad), "support draws")
  bad <- a$back; bad$analysis_id <- NULL; bad$main_version <- "0.9"
  expect_error(main_rerun(bad), "main_version")
  # From 4.4.0 the margins are settings; a tumor margin edited without its support margin is refused.
  bad <- a$back; bad$analysis_id <- NULL; bad$margins$tumor <- 3
  expect_error(main_rerun(bad), "support margin differs from its tumor margin")
  bad <- a$back; bad$analysis_id <- NULL; bad$params$alpha <- NA
  expect_error(main_rerun(bad), "alpha")
  bad <- a$back; bad$analysis_id <- NULL; bad$calibration$c_over_alpha <- 1
  expect_error(main_rerun(bad), "calibration")
  expect_error(main_rerun(a$back, draws_loader = function() stop("disk error")), "cannot be read")
  bad <- a$back; bad$analysis_id <- NULL; bad$analysis <- "candidate"
  expect_error(main_rerun(bad), "not a main-list export")
  # Existing Candidate scenarios are untouched by this schema.
  expect_identical(SCHEMA_VERSION, "3")
})

test_that("F08: the page offers the main-list download", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    f <- output$home_download
    expect_true(file.exists(f))
    files <- utils::unzip(f, list = TRUE)$Name
    expect_setequal(files, c("points.csv", "summary.csv", "scenario.json", "README.txt"))
    out <- tempfile(); utils::unzip(f, exdir = out)
    expect_equal(nrow(utils::read.csv(file.path(out, "points.csv"))), 60)
    expect_equal(nrow(utils::read.csv(file.path(out, "summary.csv"))), 15)
  })
})

# --- F10: the test workflow runs on pushes to v4 --------------------------------------
test_that("F10: CI runs on pushes to the v4 branch", {
  wf <- readLines(file.path(app_abs, ".github", "workflows", "tests.yml"))
  br <- grep("^\\s*branches:", wf, value = TRUE)[1]
  expect_match(br, "\\bv4\\b")
})
