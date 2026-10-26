# Stage A of the 2026-10-25 review (docs/development/PLAN_codex_review_2026_10_25.md):
# app correctness fixes that do not depend on the fit. One block per finding ID
# (S1-S9, OAR-01 to OAR-05, D1-D3, M3, M4, "eligible" wording). Every DVH file
# here is synthetic and written to a temporary folder.

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
lib <- load_constraints(file.path(app_abs, "data", "oar_constraints.csv"))
ui_html_a <- function() {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  htmltools::renderTags(env$ui)$html
}

defaults_a <- list(
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
set_a <- function(session, ...) do.call(session$setInputs, modifyList(defaults_a, list(...)))

unzip_to <- function(f) { td <- tempfile("stageA_"); dir.create(td); zip::unzip(f, exdir = td); td }
html_of <- function(x) tryCatch(if (is.list(x)) as.character(x$html) else as.character(x), error = function(e) "")
in_app <- function(expr) { old <- setwd(app_abs); on.exit(setwd(old)); force(expr) }
write_tmp <- function(lines, ext = ".csv") { f <- tempfile("synthetic_dvh_", fileext = ext); writeLines(lines, f, useBytes = TRUE); f }

# --- golden fixtures (review section 4) -----------------------------------------

test_that("golden: D18 vs 35 x 2 Gy/47 d at -3 Gy lists 30/29/28 without organ limits and with a non-binding cord", {
  shiny::testServer(app_dir, {
    set_a(session, min_gain = -3)
    expect_equal(suggestions()$m, c(30, 29, 28))
    session$setInputs(`oar-v_cord` = 50, `oar-use_in_search` = TRUE)   # cord 33 Gy at most: far within 50 Gy EQD2
    expect_equal(suggestions()$m, c(30, 29, 28))
  })
})

# --- S9: an enabled organ check never turns into "no check" ---------------------

test_that("S9: clearing an applied organ's alpha/beta blocks the search instead of dropping the cord limit", {
  shiny::testServer(app_dir, {
    set_a(session, min_gain = -3, `oar-v_cord` = 100, `oar-use_in_search` = TRUE)
    expect_equal(nrow(suggestions()), 0)                      # golden fixture 2: the cord limit binds
    session$setInputs(`oar-ab_cord` = NA_real_)
    expect_error(suggestions(), "organ")                       # blocked, never three unchecked rows
    expect_error(agreement(), "organ")
    expect_error(output$`unc-summary`, "organ")
    session$setInputs(`oar-ab_cord` = 2)
    expect_equal(nrow(suggestions()), 0)
    session$setInputs(`oar-ab_oar` = NA_real_)                 # the other organ alpha/beta too
    expect_error(suggestions(), "organ")
    session$setInputs(`oar-ab_oar` = 3, `oar-use_in_search` = FALSE, `oar-ab_cord` = NA_real_)
    expect_equal(suggestions()$m, c(30, 29, 28))               # not applied: invalid organ inputs do not matter
  })
})

test_that("S9: an entered organ with no verified constraint in the library blocks an applied search", {
  shiny::testServer(app_dir, {
    set_a(session, min_gain = -3, `oar-v_mandible` = 200, `oar-use_in_search` = TRUE)
    expect_error(suggestions(), "no verified constraint")
    expect_error(output$`unc-summary`, "no verified constraint")
    expect_false(grepl("No organ dose is entered", html_of(output$`oar-binding`), fixed = TRUE))
    session$setInputs(`oar-library` = "HYHOPE")               # HYHOPE has a mandible limit
    expect_true(is.data.frame(suggestions()))
  })
})

test_that("S9: the checkers refuse an organ without a constraint instead of passing it", {
  expect_error(build_oar_check(lib, list(mandible = 200), "QUANTEC")(20, 2.75), "no verified constraint")
  chk <- oar_draw_checker(lib, list(mandible = 200), "QUANTEC")
  expect_error(chk(data.frame(m = 20, d = 2.75, T = 26), toxicity_draws(5)), "no verified constraint")
})

# --- S1: stale uncertainty results are never exported under new inputs ---------

test_that("S1: changing any effective input after the search invalidates the uncertainty results", {
  changes <- list(min_gain = 20, `unc-variant` = "paper_restrictions", `unc-tox_source` = "fixed",
                  `unc-rule` = "share_joint", `unc-tox_share` = 95, d_max = 3.3, m_range = c(16, 19),
                  late_mode = "hypno", tol_late = 5, `oar-v_cord` = 40)
  shiny::testServer(app_dir, {
    run <- 0
    for (nm in names(changes)) {
      set_a(session, comparator = "hypno_test", min_gain = -3, d_max = 3.5, m_range = c(15, 19),
            `oar-use_in_search` = TRUE)
      run <- run + 1
      session$setInputs(`unc-run` = run)
      expect_true(grepl("listed", html_of(output$`unc-search_table`), fixed = TRUE), info = nm)
      do.call(session$setInputs, changes[nm])
      expect_false(grepl("listed", html_of(output$`unc-search_table`), fixed = TRUE), info = nm)
      expect_true(grepl("settings changed", html_of(output$`unc-search_msg`)), info = nm)
      td <- unzip_to(output$`unc-download`)
      expect_false(file.exists(file.path(td, "search_with_neighbors.csv")), info = nm)
      sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
      expect_null(sc$search, info = nm)
    }
  })
})

# --- S2 and S8: the applied organ policy is exported, replayed and uses the tab's alpha/beta

test_that("S2 golden: HYPNO against itself with cord 100% applied has toxicity share 0, and the file replays as 0", {
  for (src in c("prior", "fixed")) shiny::testServer(app_dir, {
    set_a(session, comparator = "hypno_test", min_gain = -3, `oar-v_cord` = 100, `oar-use_in_search` = TRUE,
          `unc-tox_source` = src)
    td <- unzip_to(output$`unc-download`)
    saved <- read.csv(file.path(td, "evaluated_regimen_summary.csv"))
    expect_equal(saved$share_tox, 0, info = src)
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    expect_true(isTRUE(sc$oar$applied), info = src)
    again <- in_app(rerun_uncertainty(sc))
    expect_equal(again$share_tox, 0, info = src)
    for (k in names(saved)) expect_equal(again[[k]], saved[[k]], tolerance = 1e-12, info = paste(src, k))
  })
})

test_that("S8: fixed-toxicity uncertainty uses the organ alpha/beta set on the Organs at risk tab", {
  shiny::testServer(app_dir, {
    # cord 83% of 55 Gy = 45.65 Gy in 20 fractions: EQD2 50.81 at alpha/beta 0.5 (exceeds 50), 48.87 at 2 (within)
    set_a(session, comparator = "hypno_test", min_gain = -3, `unc-tox_source` = "fixed",
          `oar-v_cord` = 83, `oar-use_in_search` = TRUE, `oar-ab_cord` = 0.5)
    a <- oar()$table; expect_equal(a$status[a$organ == "Spinal cord"], "exceeds the converted constraint")
    td <- unzip_to(output$`unc-download`)
    x <- read.csv(file.path(td, "evaluated_regimen_summary.csv"))
    expect_equal(x$share_tox, 0)
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    expect_equal(in_app(rerun_uncertainty(sc))$share_tox, 0)
    session$setInputs(`oar-ab_cord` = 2)
    y <- read.csv(file.path(unzip_to(output$`unc-download`), "evaluated_regimen_summary.csv"))
    expect_equal(y$share_tox, 1)
  })
})

test_that("S8: one fixed draw gives the deterministic organ check at several alpha/beta values", {
  p <- default_params()
  regs <- data.frame(m = c(15, 20, 30), d = c(3.2, 2.75, 2.2), T = c(19, 26, 40))
  for (ab in list(c(cord = 0.5, oar = 3), c(cord = 2, oar = 3), c(cord = 1.2, oar = 10))) {
    det <- build_oar_check(lib, list(cord = 80, brainstem = 85), "QUANTEC", ab)(regs$m, regs$d)
    unc <- oar_draw_checker(lib, list(cord = 80, brainstem = 85), "QUANTEC")(regs, fixed_toxicity_draws(p, 1, ab))
    expect_equal(as.vector(unc), det)
  }
})

# --- S3: emitted replay instructions run in a fresh R process --------------------

readme_code <- function(readme) {
  i <- grep("^To rerun in R from the app folder", readme)
  out <- character()
  for (l in readme[-(1:i)]) { if (!grepl("^  ", l)) break; out <- c(out, l) }
  out
}
fresh_replay <- function(td) {
  code <- readme_code(readLines(file.path(td, "README.txt")))
  code <- gsub("'scenario.json'", sprintf("'%s'", normalizePath(file.path(td, "scenario.json"), winslash = "/")), code, fixed = TRUE)
  callr::r(function(dir, code) { setwd(dir); out <- NULL; for (l in code) out <- eval(parse(text = l), envir = globalenv()); as.data.frame(out) },
           args = list(dir = app_abs, code = code), libpath = .libPaths())
}

test_that("S3: the deterministic README replays six-day, Thursday-start and organ scenarios in a clean R process", {
  shiny::testServer(app_dir, {
    set_a(session, fx_week = "6", criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 10)
    s6 <- suggestions(); td6 <- unzip_to(output$download)
    set_a(session, start_day = "4", criterion = "either", min_gain = -10, early_mode = "reference", tol_early = 10)
    s4 <- suggestions(); td4 <- unzip_to(output$download)
    set_a(session, min_gain = -3, `oar-v_cord` = 75, `oar-use_in_search` = TRUE)
    so <- suggestions(); tdo <- unzip_to(output$download)
    r6 <- fresh_replay(td6); expect_equal(r6$d, s6$d); expect_equal(r6$T, s6$T)
    r4 <- fresh_replay(td4); expect_equal(r4$T, s4$T)
    ro <- fresh_replay(tdo); expect_equal(ro$m, so$m)
  })
})

# --- S4: strict uncertainty replay ------------------------------------------------

test_that("S4: replay refuses incompatible or malformed uncertainty metadata, and old 3.0 files", {
  shiny::testServer(app_dir, {
    set_a(session, comparator = "hypno_test", min_gain = -3, d_max = 3.5, m_range = c(15, 19))
    session$setInputs(`unc-run` = 1)
    sc <- jsonlite::fromJSON(file.path(unzip_to(output$`unc-download`), "scenario.json"))
    expect_equal(nrow(in_app(rerun_uncertainty(sc))), 1)
    mut <- list(
      schema = function(x) { x$schema_version <- "999"; x },
      equation = function(x) { x$equation_version <- "999"; x },
      selection = function(x) { x$selection_version <- "999"; x },
      old_uncertainty = function(x) { x$uncertainty_version <- "3.0"; x },
      tox_source = function(x) { x$uncertainty$toxicity$source <- "TYPO"; x },
      rng = function(x) { x$uncertainty$toxicity$rng <- "Wichmann-Hill/Inversion/Rejection"; x },
      prior = function(x) { x$uncertainty$toxicity$prior$ab_e$median <- 11; x },
      seed = function(x) { x$uncertainty$toxicity$seed <- 2.5; x },
      count = function(x) { x$uncertainty$draws$di$n_draws <- 17; x },
      file = function(x) { x$uncertainty$draws$di$file <- "tumor_di_other.csv"; x },
      rule = function(x) { x$search$rule <- "unknown"; x },
      calendar = function(x) { x$search$candidate_calendar$fx_per_week <- 5.5; x },
      skip = function(x) { x$search$candidate_calendar$skip_days <- 3; x },
      share = function(x) { x$search$tox_share_min <- 3; x },
      early = function(x) { x$limits$early$mode <- "TYPO"; x },
      oar_lib = function(x) { x$oar$applied <- TRUE; x$oar$library <- "TYPO"; x },
      oar_flag = function(x) { x$oar$applied <- "yes"; x })
    for (nm in names(mut)) expect_error(in_app(rerun_uncertainty(mut[[nm]](sc))), "Cannot rerun", info = nm)
    expect_error(in_app(rerun_uncertainty(mut$old_uncertainty(sc))), "3.0")
  })
})

test_that("S4: the search part of an uncertainty file replays its listed rows and neighbors exactly", {
  shiny::testServer(app_dir, {
    set_a(session, comparator = "hypno_test", min_gain = -3, d_max = 3.5, m_range = c(15, 19),
          `oar-v_cord` = 60, `oar-use_in_search` = TRUE)
    session$setInputs(`unc-run` = 1)
    td <- unzip_to(output$`unc-download`)
    saved <- read.csv(file.path(td, "search_with_neighbors.csv"))
    sc <- jsonlite::fromJSON(file.path(td, "scenario.json"))
    again <- in_app(rerun_uncertainty(sc, part = "search"))
    expect_equal(nrow(again), nrow(saved))
    for (k in intersect(names(saved), names(again)))
      if (is.numeric(saved[[k]])) expect_equal(again[[k]], saved[[k]], tolerance = 1e-12, info = k) else
        expect_equal(as.character(again[[k]]), as.character(saved[[k]]), info = k)
  })
})

# --- S5: one candidate-calendar policy -------------------------------------------

test_that("S5: both searches build candidates on the regular calendar; skipped days apply to the evaluated regimen", {
  shiny::testServer(app_dir, {
    set_a(session, skip_days = "2", comparator = "hypno_test", min_gain = -3, d_max = 3.5, m_range = c(15, 19), min_days = 1)
    expect_equal(regimen()$T, 29)                       # the evaluated regimen keeps its gap
    session$setInputs(`unc-run` = 1)
    s <- read.csv(file.path(unzip_to(output$`unc-download`), "search_with_neighbors.csv"))
    expect_equal(s$T, schedule_time(s$m))               # 20 fractions: 26 days, not 29
    sc <- jsonlite::fromJSON(file.path(unzip_to(output$`unc-download`), "scenario.json"))
    expect_length(unlist(sc$search$candidate_calendar$skip_days), 0)
  })
  expect_true(grepl("Skipped days apply to the evaluated regimen only", ui_html_a(), fixed = TRUE))
})

# --- S7: the calendar API validates like the interface ----------------------------

test_that("S7: the direct calendar API rejects fractional, missing, negative and out-of-range skipped days", {
  for (bad in list(2.5, NA, -3, 0, 401, c(2, NA), "x"))
    expect_error(schedule_time(10, skip_days = bad), info = paste(bad, collapse = ","))
  expect_error(schedule_time(10, fx_per_week = 5.5))
  expect_error(schedule_time(10, start_weekday = 2.5))
  expect_error(schedule_time(NA))
  expect_error(schedule_time(10.5))
  expect_equal(schedule_time(20, skip_days = c(3, 3, 4)), schedule_time(20, skip_days = c(3, 4)))
  expect_equal(schedule_time(20, skip_days = list()), schedule_time(20))   # an empty JSON array
})

# --- OAR-01: source physical mean limits are checked on physical dose -------------

test_that("OAR-01: a parotid physical mean of 25 Gy in 35 fractions does not pass the 20 Gy physical goal", {
  r <- oar_row(lib, "parotid", list(mode = "gy", value = 25), 35, 2, "QUANTEC", ab = 3)
  expect_equal(r$status, "exceeds the source limit")
  expect_equal(r$limit_physical, 20)
  expect_match(r$basis, "physical")
  expect_equal(r$planned_eqd2, 25 * (25 / 35 + 3) / 5)            # 18.5714: illustration only
  expect_true(isTRUE(grepl("illustrat", r$flags)))
  # the status does not depend on alpha/beta or on DVH shape when the source metric is a physical mean
  expect_equal(oar_row(lib, "parotid", list(mode = "gy", value = 25), 35, 2, "QUANTEC", ab = 10)$status, r$status)
  h <- oar_row(lib, "parotid", list(mode = "gy", value = 21.9), 15, 4, "HYHOPE", ab = 3)
  expect_equal(h$status, "within 5% of the source limit")            # 21.9 < 22 physical, at HYHOPE's own 15 fractions
  # the cord's explicit 2-Gy anchor stays an EQD2 comparison
  expect_match(oar_row(lib, "cord", list(mode = "gy", value = 45), 20, 2.25, "QUANTEC", ab = 2)$basis, "EQD2")
})

test_that("OAR-01: the search checks use the same physical mean comparison", {
  chk <- build_oar_check(lib, list(parotid = 30), "QUANTEC", c(cord = 2, oar = 3))
  # 30% of 70 Gy = 21 Gy physical mean > 20: fails, although its uniform EQD2 (15.96 Gy at 35 fx) is below 20
  expect_false(chk(35, 2))
  expect_true(chk(35, 1.8))            # 18.9 Gy
})

# --- OAR-02 to OAR-05 ---------------------------------------------------------------

test_that("OAR-02: HYHOPE PRVs are organs of their own; a cord PRV above its limit is not hidden by the cord", {
  t <- oar_table(lib, list(cord = list(mode = "gy", value = 36), cord_prv = list(mode = "gy", value = 44)),
                 15, 3.6, "HYHOPE", ab = c(cord = 2, oar = 3))
  expect_equal(t$status[t$organ == "Spinal cord"], "within 5% of the converted constraint")
  expect_equal(t$status[t$organ == "Spinal cord + 5 mm"], "exceeds the converted constraint")
  expect_true(all(c("Brainstem + 5 mm") %in% t$organ))
})

test_that("OAR-02: HYHOPE's strict '<' limit: cord exactly 37 Gy in 15 fractions is not per protocol", {
  expect_equal(oar_row(lib, "cord", list(mode = "gy", value = 37), 15, 37 / 15, "HYHOPE", ab = 2)$status,
               "exceeds the converted constraint")
  expect_false(build_oar_check(lib, list(cord = 100), "HYHOPE")(15, 37 / 15))
})

test_that("OAR-02: unassessed mandatory metrics are listed and block an 'all met' reading", {
  shiny::testServer(app_dir, {
    set_a(session, `oar-library` = "HYHOPE", `oar-v_cord` = 60)
    b <- html_of(output$`oar-binding`)
    expect_true(grepl("Spinal cord + 5 mm", b, fixed = TRUE))
    expect_true(grepl("not assessed", b, fixed = TRUE))
    expect_true(grepl("selected constraint metrics", b, fixed = TRUE))
  })
  cov <- oar_coverage(lib, list(cord = 60), "HYHOPE")
  expect_true(all(c("Spinal cord + 5 mm", "Brainstem", "Brainstem + 5 mm", "Mandible") %in% cov$mandatory_unassessed))
  expect_false(cov$complete)
})

test_that("OAR-03: HYHOPE 'Max' keeps its unspecified volume definition; no D0.03cc is claimed", {
  t <- oar_table(lib, list(cord = list(mode = "gy", value = 30)), 15, 2, "HYHOPE")
  expect_true(isTRUE(grepl("volume definition not specified", t$source_metric[t$organ == "Spinal cord"])))
  src <- c(unlist(lapply(list.files(file.path(app_abs, "R"), full.names = TRUE), readLines, warn = FALSE)),
           readLines(file.path(app_abs, "data", "oar_constraints.csv"), warn = FALSE))
  expect_false(any(grepl("D0\\.03|0\\.03 ?cc", src)))
})

test_that("OAR-04: the upper end of a source range is labeled as the calculator's choice", {
  r <- oar_row(lib, "parotid", list(mode = "gy", value = 18), 15, 2, "HYHOPE")
  expect_true(isTRUE(grepl("upper end is this calculator's choice", r$flags)))
  q <- oar_row(lib, "cochlea", list(mode = "gy", value = 30), 35, 2, "QUANTEC")
  expect_true(isTRUE(grepl("35 Gy", q$flags)))
})

test_that("OAR-05: HYHOPE rows name protocol version 4.2 and the section, and 12 fractions is historical", {
  h <- lib[lib$library == "HYHOPE", ]
  expect_true(all(grepl("v4.2", h$source_short, fixed = TRUE)))
  expect_false(any(h$source_short == "HYHOPE protocol 4.1.4.2"))
  r <- oar_row(lib, "cord", list(mode = "gy", value = 30), 12, 2.5, "HYHOPE")
  expect_true(isTRUE(grepl("removed by amendment", r$flags)))
  expect_true(nzchar(OAR_LIBRARY_META$HYHOPE$version))
})

# --- D1-D3: DVH maxima, validation, identifiers -------------------------------------

test_that("D1: the DVH maximum is an interval; a top bin straddling the limit is indeterminate", {
  dvh <- list(mode = "dvh", value = data.frame(dose = c(0, 40, 50), volume = c(100, 1, 0)))
  r <- oar_row(lib, "cord", dvh, 20, 2.5, "QUANTEC", ab = 2)   # EQD2 at 45: 47.81 (midpoint); at 50: 56.25
  expect_equal(r$status, "indeterminate at this DVH resolution")
  dd <- cumulative_to_differential(c(0, 40, 50), c(100, 1, 0))
  mm <- dvh_means(dd, 20, 2)
  expect_equal(c(mm$max_lo, mm$max_hi), c(40, 50))
  low <- list(mode = "dvh", value = data.frame(dose = c(0, 30, 35), volume = c(100, 1, 0)))
  expect_equal(oar_row(lib, "cord", low, 20, 2.5, "QUANTEC", ab = 2)$status, "within")
  high <- list(mode = "dvh", value = data.frame(dose = c(0, 50, 60), volume = c(100, 1, 0)))
  expect_equal(oar_row(lib, "cord", high, 20, 3, "QUANTEC", ab = 2)$status, "exceeds the converted constraint")
})

test_that("D2: malformed DVHs are refused (NaN, Inf, negative, duplicate, unordered, truncated, short, partial %)", {
  hdr <- "Dose (Gy),Volume (%)"
  bad <- list(
    nan = c(hdr, "0,100", "NaN,50", "50,0"), inf = c(hdr, "0,100", "Inf,50", "60,0"),
    vinf = c(hdr, "0,Inf", "40,50", "60,0"), neg_dose = c(hdr, "-10,100", "40,50", "60,0"),
    neg_vol = c(hdr, "0,100", "40,1", "50,-10"), duplicate = c(hdr, "0,100", "40,50", "40,20", "60,0"),
    unordered = c(hdr, "0,100", "60,0", "40,50"), truncated = c(hdr, "0,100", "40,50", "50,5"),
    one_row = c(hdr, "0,100"), partial_start = c(hdr, "10,80", "40,50", "60,0"),
    over_100 = c(hdr, "0,200", "40,50", "60,0"), all_zero = c(hdr, "0,0", "40,0", "60,0"),
    cgy_as_gy = c(hdr, "0,100", "4000,50", "7000,0"))
  for (nm in names(bad)) {
    f <- write_tmp(bad[[nm]])
    expect_error(read_dvh_upload(f, file.size(f)), info = nm)
  }
  ok <- write_tmp(c(hdr, "0,100", "40,50", "60,0"))
  expect_equal(nrow(read_dvh_upload(ok, file.size(ok))[[1]]), 3)
})

test_that("D3: quoted, BOM-prefixed and padded identifier headers are refused", {
  rows <- c('"SYNTHETIC_REVIEW_ONLY",0,100', '"SYNTHETIC_REVIEW_ONLY",40,1', '"SYNTHETIC_REVIEW_ONLY",50,0')
  for (h in c('"Patient ID","Dose (Gy)","Volume (%)"', '﻿"Patient ID","Dose (Gy)","Volume (%)"',
              '" patient_id ","Dose (Gy)","Volume (%)"', 'MRN,Dose (Gy),Volume (%)',
              '"Dose (Gy)","Volume (%)","Date of Birth"')) {
    lines <- if (grepl("Date of Birth", h)) c(h, "0,100,X", "40,1,X", "50,0,X") else c(h, rows)
    f <- write_tmp(lines)
    expect_error(read_dvh_upload(f, file.size(f)), "identifier", info = h)
  }
  # a blank identifier column is not a value
  f <- write_tmp(c('"Patient ID","Dose (Gy)","Volume (%)"', '"",0,100', '"",40,1', '"",50,0'))
  expect_equal(nrow(read_dvh_upload(f, file.size(f))[[1]]), 3)
  html <- ui_html_a()
  expect_false(grepl("never stored", html, fixed = TRUE))
  expect_true(grepl("temporary file on the server", html, fixed = TRUE))
})

# --- M3: invalid draws in the conservative score -------------------------------------

m3_fixture <- function(k, n = 100) {
  b <- matrix(c(rep(NA_real_, k), rep(10, n - k)), nrow = 1); v <- !is.na(b)
  summarize_draws(data.frame(m = 20, d = 2.75, T = 26), list(delta_di = b, delta_dd = b, valid_di = v, valid_dd = v),
                  list(ok = matrix(TRUE, 1, n)), threshold = 0)
}
test_that("M3: one valid +10 draw among 99 invalid draws is not a conservative +10", {
  s <- m3_fixture(99)
  expect_true(is.na(s$min_q10))
  expect_equal(s$n_valid_both, 1)
})
test_that("M3: the 10th percentile is evaluable only while the invalid draws stay below it", {
  for (k in c(0, 1, 9)) expect_equal(m3_fixture(k)$min_q10, 10, info = as.character(k))
  for (k in c(10, 90, 99, 100)) expect_true(is.na(m3_fixture(k)$min_q10), info = as.character(k))
  # invalid draws count as failures below the valid ones: 5 invalid of 100, values 1..95 -> type-7 q10 of the full ensemble
  b <- matrix(c(rep(NA_real_, 5), 1:95), nrow = 1); v <- !is.na(b)
  s <- summarize_draws(data.frame(m = 1, d = 1, T = 1), list(delta_di = b, delta_dd = b, valid_di = v, valid_dd = v), threshold = 0)
  expect_equal(s$min_q10, 5.9)                       # sorted x[10] = 5, x[11] = 6; h = 10.9
  # negative raw values are kept, never floored
  b2 <- matrix(c(-50, -40, seq(1, 98)), nrow = 1); v2 <- !is.na(b2)
  s2 <- summarize_draws(data.frame(m = 1, d = 1, T = 1), list(delta_di = b2, delta_dd = b2, valid_di = v2, valid_dd = v2), threshold = 0)
  expect_equal(s2$min_q10, unname(stats::quantile(as.vector(b2), 0.1, type = 7)))
})

# --- M4 and the "eligible" wording ----------------------------------------------------

test_that("M4: at a required toxicity share of 100% every listed regimen passes on all draws", {
  shiny::testServer(app_dir, {
    set_a(session, `unc-tox_share` = 100)
    session$setInputs(`unc-run` = 1)
    s <- read.csv(file.path(unzip_to(output$`unc-download`), "search_with_neighbors.csv"))
    listed <- s[s$role == "listed", ]
    expect_true(nrow(listed) > 0)
    expect_true(all(listed$share_tox >= 1 - 1e-12))
    expect_false(any(listed$m == 8 & abs(listed$d - 4.5) < 1e-9))
    expect_false(any(listed$m == 11 & abs(listed$d - 3.8) < 1e-9))
    failed <- s[s$role == "failed_full_check", ]
    expect_true(any(failed$m == 8 & abs(failed$d - 4.5) < 1e-9))
    expect_true(grepl("all draws", html_of(output$`unc-search_msg`)))
  })
})

test_that("M4 and wording: at the default 90% no choice fails the full check; 5 x 5 Gy is not shown as meeting +1", {
  shiny::testServer(app_dir, {
    set_a(session)
    session$setInputs(`unc-run` = 1)
    s <- read.csv(file.path(unzip_to(output$`unc-download`), "search_with_neighbors.csv"))
    expect_equal(sum(s$role == "failed_full_check"), 0)
    msg <- html_of(output$`unc-search_msg`)
    expect_true(grepl("tumor criterion may fail", msg, fixed = TRUE))
    expect_false(grepl("eligible regimen", msg, fixed = TRUE))
    tab <- html_of(output$`unc-search_table`)
    expect_true(grepl("Tumor threshold met", tab, fixed = TRUE))
    five <- s[s$role == "listed" & s$m == 5, ]
    expect_equal(nrow(five), 1)
    expect_lt(five$share_both, 0.05)
    expect_true(all(c("n_valid_both", "n_invalid_di", "n_invalid_dd") %in% names(s)))
  })
})
