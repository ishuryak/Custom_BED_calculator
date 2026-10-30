# Version 3.1, stage E (plan PLAN_codex_review_2026_10_25.md): the app uses the
# verified-trial refit. DI posterior draws; DD calibrated to DI at 35 x 2 Gy in
# 47 days and fitted conditional on each DI draw, so DD row i belongs to DI row i
# (review M1). Draws that fail the sampling-diagnostic gate are refused (M7).
# Variable-alpha/beta BEDs are not labeled Gy10 (M2), and the Monte Carlo error is
# not called an upper bound (M8).

app_abs <- normalizePath(testthat::test_path("..", ".."))
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
draws_dir <- file.path(app_abs, "draws")
man <- jsonlite::fromJSON(file.path(draws_dir, "draws_manifest.json"), simplifyVector = FALSE)

test_that("the shipped draws come from the corrected verified-trial refit and pass the gate", {
  expect_identical(man$draws_version, "3.2")   # 3.2: regenerated after the full-gate ESS repair (review R01)
  expect_identical(man$dataset, "refit_trials_2026_10_25c_main.csv")
  data_sha <- digest::digest(file.path(app_abs, "data", man$dataset), algo = "sha256", file = TRUE)
  expect_identical(man$data_sha256, data_sha)
  for (k in names(man$files)) expect_true(isTRUE(man$files[[k]]$diagnostics$gate_pass), info = k)
})

test_that("DD row i was fitted on DI row i, in both variants", {
  for (v in c("primary", "paper_restrictions")) {
    di <- utils::read.csv(file.path(draws_dir, sprintf("tumor_di_%s.csv", v)))
    dd <- utils::read.csv(file.path(draws_dir, sprintf("tumor_dd_%s.csv", v)))
    expect_equal(nrow(di), nrow(dd))
    expect_identical(di$di_draw, dd$di_draw)
    expect_equal(dd$alpha, di$alpha, tolerance = 0); expect_equal(dd$ab, di$ab, tolerance = 0)
    expect_equal(dd$g, di$g, tolerance = 0)
    # DD repopulation equals DI repopulation at the anchor in every row
    K0 <- 70 * (1 + 2 / di$ab)
    rep_di <- di$lambda * smooth_max0(47 - di$Tk)
    rep_dd <- dd$lambda * (1 - exp(-dd$alpha * K0 / 47)) * smooth_max0(47 - dd$c_over_alpha * 47 / K0)
    expect_lt(max(abs(rep_dd - rep_di)), 1e-6)
  }
  expect_true(all(utils::read.csv(file.path(draws_dir, "tumor_di_paper_restrictions.csv"))$ab == 10))
})

test_that("a draws file whose fit failed the diagnostic gate is refused", {
  tmp <- file.path(tempdir(), "draws_gate"); unlink(tmp, recursive = TRUE); dir.create(tmp)
  file.copy(list.files(draws_dir, full.names = TRUE), tmp)
  m2 <- man; m2$files[["tumor_dd_primary.csv"]]$diagnostics$gate_pass <- FALSE
  jsonlite::write_json(m2, file.path(tmp, "draws_manifest.json"), auto_unbox = TRUE, digits = NA)
  expect_error(load_draws("dd", "primary", tmp), "diagnostic gate")
  expect_s3_class(load_draws("di", "primary", tmp), "data.frame")
})

test_that("the tumor data range comes from the dataset named in the manifest", {
  old <- setwd(app_abs); on.exit(setwd(old))
  sup <- tumor_data_support()
  expect_identical(sup$file, man$dataset)
  expect_equal(c(sup$d_min, sup$d_max), c(1.15, 2.55))
  expect_equal(sup$n_arms, 45)
})

test_that("BED differences over draws are labeled Gy10 only when alpha/beta is 10 in every draw", {
  old <- setwd(app_abs); on.exit(setwd(old))
  expect_match(draws_bed_unit(load_draws("di", "primary")), "varies by draw")
  expect_identical(draws_bed_unit(load_draws("di", "paper_restrictions")), "Gy₁₀")
})

test_that("the refit point values are the calibrated-DD model at the DI medians", {
  pv <- jsonlite::fromJSON(file.path(draws_dir, "refit_point_values.json"))$values
  cal <- dd_calibration(list(m = 35, d = 2, T = 47), pv$di$alpha, pv$di$ab, pv$di$lambda, pv$dd$lambda, pv$di$Tk)
  expect_identical(cal$status, "ok")
  expect_equal(cal$c_over_alpha, pv$dd$c_over_alpha, tolerance = 1e-8)
  expect_identical(APP_VERSION, "4.2.1")   # 4.2.1: range explanation under Details; 4.2.0: main page 30 fx max, early/late columns, plot; 4.1.0: v4 review fixes; 4.0.0: layout only; 3.2.1: v3.2 review fixes; 3.2.0: 2026-10-26 review
})

test_that("superseded claims about the old fit are gone from the app text", {
  src <- paste(c(readLines(file.path(app_abs, "app.R"), encoding = "UTF-8"),
                 readLines(file.path(app_abs, "R", "mod_uncertainty.R"), encoding = "UTF-8")), collapse = "\n")
  for (pat in c("prior-dominated", "(at most)", "35 randomized", "Bayesian refit 2026-10-21", "only within-trial contrasts"))
    expect_false(grepl(pat, src, fixed = TRUE), info = pat)
})
