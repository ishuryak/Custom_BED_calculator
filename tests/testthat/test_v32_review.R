# Version 3.2: the 2026-10-26 review (docs/development/PLAN_codex_review_2026_10_26.md) and the
# new default evaluated regimen. One block per item. Every DVH file here is synthetic.

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
source(file.path(app_abs, "tools", "deploy_preflight.R"), local = TRUE)
lib <- load_constraints(file.path(app_abs, "data", "oar_constraints.csv"))
ui_html_32 <- function() {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  htmltools::renderTags(env$ui)$html
}
defaults_32 <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
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
set_32 <- function(session, ...) do.call(session$setInputs, modifyList(defaults_32, list(...)))
unzip_32 <- function(f) { td <- tempfile("v32_"); dir.create(td); zip::unzip(f, exdir = td); td }
in_app_32 <- function(expr) { old <- setwd(app_abs); on.exit(setwd(old)); force(expr) }
tmp_lines <- function(lines, ext = ".txt") { f <- tempfile("synthetic_dvh_", fileext = ext); writeLines(lines, f, useBytes = TRUE); f }

# --- default evaluated regimen (Igor Shuryak, 2026-10-26): 18 x 2.9 Gy ----------------------

test_that("the default evaluated regimen is 18 x 2.9 Gy (52.2 Gy in 24 days)", {
  h <- ui_html_32()
  expect_match(h, 'id="m"[^>]*value="18"')
  expect_match(h, 'id="d"[^>]*value="2.9"')
  expect_match(h, 'id="T"[^>]*value="24"')
  expect_true(grepl("Default: 18 × 2.9 Gy (52.2 Gy)", h, fixed = TRUE))
  expect_equal(schedule_time(18, 5, 1, integer()), 24)
  shiny::testServer(app_dir, {
    set_32(session)
    expect_equal(regimen()[c("m", "d", "T")], list(m = 18, d = 2.9, T = 24))
  })
})

# --- R04: every benchmark row shows its tumor BED difference --------------------------------

test_that("R04: all three benchmark rows carry their DI and DD differences, matching Compare regimens", {
  shiny::testServer(app_dir, {
    set_32(session)
    b <- output$benchmarks
    expect_false(grepl("not available", b, fixed = TRUE))
    for (v in c("-1.39", "-8.67", "+1.11", "-2.87", "-10.60", "-15.14")) expect_true(grepl(v, b, fixed = TRUE), info = v)
    expect_equal(tumor_delta_txt(c(70, 60, 50), 71.51), c("-1.51", "-11.51", "-21.51"))
    expect_equal(tumor_delta_txt(c(70, -1), 71.51)[2], "not evaluable (BED ≤ 0)")
    expect_error(tumor_delta_txt(c(70, 60, 50), c(71, 72)), "length 1 or length")
  })
})

# --- R05: an indeterminate DVH status stays indeterminate over the alpha/beta range ---------

test_that("R05: a straddling top bin is indeterminate at every alpha/beta, never 'within throughout'", {
  dvh <- list(mode = "dvh", value = data.frame(dose = c(0, 40, 50), volume = c(100, 1, 0)))
  r <- oar_row(lib, "cord", dvh, 20, 2.5, "QUANTEC", ab = 2, ab_range = c(0.87, 2))
  expect_equal(r$status, "indeterminate at this DVH resolution")
  expect_match(r$flip, "^indeterminate at this DVH resolution at every α/β")
  expect_false(grepl("within at every", r$flip))
  low <- list(mode = "dvh", value = data.frame(dose = c(0, 30, 35), volume = c(100, 1, 0)))
  expect_match(oar_row(lib, "cord", low, 20, 2.5, "QUANTEC", ab = 2, ab_range = c(0.87, 2))$flip, "^within at every")
  high <- list(mode = "dvh", value = data.frame(dose = c(0, 50, 60), volume = c(100, 1, 0)))
  expect_match(oar_row(lib, "cord", high, 20, 2.5, "QUANTEC", ab = 2, ab_range = c(0.87, 2))$flip, "^exceeds at every")
  g <- seq(0.87, 2, length.out = 6)
  expect_equal(flip_text(g, c(rep("indeterminate", 3), rep("within", 3))),
               "indeterminate for α/β 0.87-1.3 Gy (within otherwise)")
  expect_equal(flip_text(g, c(rep("exceeds", 3), rep("indeterminate", 3))),
               "exceeds for α/β 0.87-1.3 Gy; indeterminate for α/β 1.5-2 Gy")
  expect_equal(status_class(c("within", "within 5% of the source limit", "exceeds the source limit",
                              "indeterminate at this DVH resolution", "not assessed")),
               c("within", "within", "exceeds", "indeterminate", "not assessed"))
})

# --- R10: malformed rows inside an Eclipse data block are refused, not dropped -------------

eclipse_lines <- function(rows) c("Patient Name         : ", "Patient ID           : ", "Type                 : Cumulative Dose Volume Histogram",
                                  "", "Structure: Cord_PRV", "Volume [cm³]: 20.0", "",
                                  "  Dose [Gy]   Ratio of Total Structure Volume [%]", rows)

test_that("R10: NaN, Inf, missing and text values inside an Eclipse DVH block are refused with a line number", {
  ok <- read_dvh_upload(tmp_lines(eclipse_lines(c("0 100", "40 1", "50 0"))), 1000)
  expect_equal(nrow(ok[[1]]), 3)
  bad <- list(nan = c("0 100", "40 NaN", "50 0"), inf = c("0 100", "40 Inf", "50 0"),
              missing = c("0 100", "40", "50 0"), text = c("0 100", "40 abc", "50 0"))
  for (nm in names(bad)) {
    err <- tryCatch(read_dvh_upload(tmp_lines(eclipse_lines(bad[[nm]])), 1000), error = conditionMessage)
    expect_true(is.character(err) && grepl("line 10", err), info = nm)
  }
  # existing synthetic exports still parse
  fx <- function(f) testthat::test_path("fixtures", "dvh", f)
  expect_silent(read_dvh_upload(fx("eclipse_synthetic.txt"), 1000))
  expect_silent(read_dvh_upload(fx("eclipse_gy_cc.txt"), 1000))
  # metadata lines naming a dose in cGy before the table do not count as its header (audit 2026-10-26)
  meta <- eclipse_lines(c("0 100", "40 1", "50 0"))
  meta <- append(meta, c("Min Dose [Gy]: 0.0", "Max Dose [Gy]: 50.0", "Mean Dose [Gy]: 10.0"), after = 6)
  expect_equal(nrow(read_dvh_upload(tmp_lines(meta), 1000)[[1]]), 3)
})

# --- R03: one organ-mapping validator for both replay paths and the direct checkers -------

test_that("R03: unnamed, duplicate, blank, unknown, non-scalar and non-finite organ mappings are refused", {
  bad <- list(unnamed = list(100), duplicate = list(cord = 10, cord = 100), blank = setNames(list(100), " "),
              unknown = list(liver = 50), vector = list(cord = c(50, 60)), nonfinite = list(cord = Inf),
              text = list(cord = "50"), range = list(cord = 500), nan = list(cord = NaN))
  for (nm in names(bad)) {
    expect_error(validate_oar_percents(bad[[nm]]), info = nm)
    expect_error(build_oar_check(lib, bad[[nm]], "QUANTEC"), info = nm)
    expect_error(oar_draw_checker(lib, bad[[nm]], "QUANTEC"), info = nm)
  }
  expect_equal(validate_oar_percents(list(cord = 80, parotid = NA, brainstem = NULL)), list(cord = 80))
  expect_null(build_oar_check(lib, list(cord = NA), "QUANTEC"))     # all organs blank: no organ limit
  expect_true(is.function(build_oar_check(lib, list(cord = 80, brainstem = 85), "QUANTEC")))
})

test_that("R03: both replay paths refuse a malformed organ mapping and replay a well-formed one", {
  shiny::testServer(app_dir, {
    set_32(session, comparator = "hypno_test", min_gain = -3, `oar-v_cord` = 80, `oar-use_in_search` = TRUE)
    # uncertainty replay
    sc <- jsonlite::fromJSON(file.path(unzip_32(output$`unc-download`), "scenario.json"))
    expect_true(isTRUE(sc$oar$applied))
    expect_silent(in_app_32(rerun_uncertainty(sc)))
    for (bad in list(list(100), list(cord = 80, cord = 100), list(liver = 80))) {
      x <- sc; x$oar$percents <- bad
      expect_error(in_app_32(rerun_uncertainty(x)), "Cannot rerun")
    }
    # unnamed survives a JSON round trip and is still refused
    x <- sc; x$oar$percents <- list(100)
    y <- jsonlite::fromJSON(scenario_json(x))
    expect_error(in_app_32(rerun_uncertainty(y)), "Cannot rerun")
    # deterministic search replay
    sd <- jsonlite::fromJSON(file.path(unzip_32(output$download), "scenario.json"))
    expect_false(is.null(sd$search$oar))
    expect_silent(in_app_32(rerun_scenario(sd)))
    for (bad in list(list(100), list(cord = 80, cord = 100))) {
      x <- sd; x$search$oar$percents <- bad
      expect_error(in_app_32(rerun_scenario(x)), "Cannot rerun")
    }
  })
})

# --- replay grid bound (review phase 1) --------------------------------------------------

test_that("an uncertainty replay grid above the app's largest is refused before any calculation", {
  expect_equal(grid_size(2, 5:34, 0.05, 5), nrow(regimen_grid(list(m = 35, d = 2, T = 47), 5:34, 0.05, 5, 0)))
  expect_lte(grid_size(0.1, 1:99, 0.05, 30), UNC_MAX_GRID_ROWS)       # the largest the UI can ask for (audit 2026-10-26)
  expect_lte(grid_size(0.5, 1:34, 0.05, 30), UNC_MAX_GRID_ROWS)       # the audit's refused example A
  expect_lte(grid_size(2, 1:59, 0.05, 30), UNC_MAX_GRID_ROWS)         # the audit's refused example B
  expect_gt(grid_size(2.75, 1:100, 0.001, 30), UNC_MAX_GRID_ROWS)     # the review's 2.7-million probe
  shiny::testServer(app_dir, {
    set_32(session, comparator = "hypno_test", min_gain = -3, m_range = c(16, 18), d_max = 3.2)
    session$setInputs(`unc-run` = 1)
    sc <- jsonlite::fromJSON(file.path(unzip_32(output$`unc-download`), "scenario.json"))
    expect_false(is.null(sc$search))
    expect_silent(in_app_32(rerun_uncertainty(sc, "search")))
    x <- sc; x$search$m_range <- 1:100; x$search$d_step <- 0.001; x$search$d_max <- 30
    expect_error(in_app_32(rerun_uncertainty(x, "search")), "exceeds")
  })
})

# --- R09: the suite marker is written only for a resolvable, committed revision -----------

test_that("R09: release_revision fails closed", {
  h <- paste(rep("a", 40), collapse = "")
  fake <- function(head, status = 0L, dirty = character(), dstat = 0L)
    function(args) if (args[1] == "rev-parse") list(status = status, out = head) else list(status = dstat, out = dirty)
  expect_true(release_revision(fake(h))$ok)
  expect_equal(release_revision(fake(h))$rev, h)
  expect_false(release_revision(fake(h, status = 128L))$ok)                 # git failed
  expect_false(release_revision(fake(character()))$ok)                     # empty output
  expect_false(release_revision(fake("fatal: not a git repository"))$ok)   # malformed
  expect_false(release_revision(fake(h, dirty = " M app.R"))$ok)           # dirty tracked code
  expect_false(release_revision(fake(h, dstat = 1L))$ok)                   # git status failed
  expect_false(release_revision(function(args) list(status = 127L, out = "not found"))$ok)
  s <- readLines(file.path(app_abs, "tools", "run_tests.R"))
  expect_true(any(grepl("rr$ok && rr_end$ok && identical(rr$rev, rr_end$rev)", s, fixed = TRUE)))
})

# --- R01 on the app side: draws diagnostics are checked against the gate itself ------------

test_that("R01: the app checks recorded draws diagnostics against the gate, ESS included", {
  ok <- list(divergences = 0, treedepth_hits = 0, max_rhat = 1.003, min_ess_bulk_per_fit = 420,
             min_ess_tail_per_fit = 500, min_ebfmi = 0.8)
  expect_length(draws_gate_failures(ok), 0)
  expect_match(draws_gate_failures(modifyList(ok, list(min_ess_bulk_per_fit = 151))), "ess_bulk 151")
  expect_length(draws_gate_failures(modifyList(ok, list(min_ess_tail_per_fit = NULL))), 1)   # missing fails
  expect_length(draws_gate_failures(list(divergences = 0, treedepth_hits = 0, max_rhat = 1.001, min_ess_bulk = 900,
                                         min_ess_tail = 800, min_ebfmi = 0.9)), 0)
  # a 3.2 manifest whose DD diagnostics say gate_pass but record ESS 151 is refused
  td <- tempfile("draws32_"); dir.create(td)
  file.copy(file.path(app_abs, "draws", c("tumor_di_primary.csv", "tumor_dd_primary.csv")), td)
  man <- jsonlite::fromJSON(file.path(app_abs, "draws", "draws_manifest.json"), simplifyVector = FALSE)
  man$draws_version <- "3.2"
  man$files$tumor_dd_primary.csv$diagnostics <- modifyList(ok, list(gate_pass = TRUE, min_ess_bulk_per_fit = 151))
  man$files$tumor_di_primary.csv$diagnostics <- list(gate_pass = TRUE, divergences = 0, treedepth_hits = 0, max_rhat = 1.001,
                                                    min_ess_bulk = 900, min_ess_tail = 800, min_ebfmi = 0.9)
  jsonlite::write_json(man, file.path(td, "draws_manifest.json"), auto_unbox = TRUE, digits = NA)
  expect_error(load_draws("dd", "primary", td), "ess_bulk 151")
  expect_silent(load_draws("di", "primary", td))
  # draws older than 3.2, or with no version, are refused (audit 2026-10-26)
  for (v in list("3.1", NULL)) {
    man$draws_version <- v
    man$files$tumor_dd_primary.csv$diagnostics$min_ess_bulk_per_fit <- 420
    jsonlite::write_json(man, file.path(td, "draws_manifest.json"), auto_unbox = TRUE, digits = NA)
    expect_error(load_draws("di", "primary", td), "needs 3.2 or later")
  }
})

# --- R07: wording ---------------------------------------------------------------------

test_that("R07: the review's wording corrections are in place", {
  a <- paste(readLines(file.path(app_abs, "about.md"), encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl("all in 2-Gy fractions", a, fixed = TRUE))
  expect_false(grepl("other than nasopharyngeal cancer", a, fixed = TRUE))
  expect_false(grepl("The calculator has no dose-response relation", a, fixed = TRUE))
  expect_false(grepl("because it is a population-level slope", a, fixed = TRUE))
  expect_true(grepl("not a measurement of intrinsic cellular radiosensitivity", a, fixed = TRUE))
  expect_true(grepl("in-sample check", a, fixed = TRUE))
  expect_true(grepl("55 pairs among 11 regimens", a, fixed = TRUE))
  expect_true(grepl("may understate the toxicity of hypofractionation", a, fixed = TRUE))
  h <- ui_html_32()
  expect_false(grepl("2-Gy fractions, converted", h, fixed = TRUE))
  expect_false(grepl("to be re-checked", paste(deparse(toxicity_prior()), collapse = ""), fixed = TRUE))
  shiny::testServer(app_dir, {
    set_32(session)
    msg <- as.character(output$suggest_msg$html)
    expect_false(grepl("Model agreement tab", msg, fixed = TRUE))
    expect_true(grepl("The Both-model screen uses", msg, fixed = TRUE))
  })
  r <- oar_row(lib, "parotid", list(mode = "gy", value = 19), 5, 5, "QUANTEC", ab = 3)
  expect_match(r$flags, "physical source comparison only")
})
