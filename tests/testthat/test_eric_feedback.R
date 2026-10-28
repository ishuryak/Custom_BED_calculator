# check report card and search settings
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
feedback_defaults <- list(
  m = 20, d = 2.75, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 26,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  cmp1_m = 17, cmp1_d = 3.05, cmp2_m = 15, cmp2_d = 3.1, cmp3_m = 30, cmp3_d = 2,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent",
  `oar-use_in_search` = FALSE, `unc-variant` = "primary", `unc-tumor_crit` = "both", `unc-tox_source` = "prior",
  `unc-tcp0` = 70, `unc-rule` = "min_q10", `unc-tox_share` = 90)
set_feedback <- function(session, ...) do.call(session$setInputs, modifyList(feedback_defaults, list(...)))

test_that("negative BEDs cannot meet the tumor threshold", {
  shiny::testServer(app_dir, {
    set_feedback(session, m = 18, d = 2.9, preset = "custom", alpha = 0.01,
                 g = 0.1, lambda_di = 0, lambda_dd = 0)
    b <- both()
    expect_equal(b$bed_di, c(-172.662, -386), tolerance = 1e-6)
    expect_equal(b$bed_dd, b$bed_di)
    r <- output$report
    expect_match(r, "not evaluable (BED ≤ 0)", fixed = TRUE)
    expect_equal(length(regmatches(r, gregexpr("threshold +1: not evaluable", r, fixed = TRUE))[[1]]), 2)
    expect_false(grepl("threshold +1: met", r, fixed = TRUE))
    expect_identical(tumor_threshold_txt(0, 1, -2), "threshold -2: not evaluable")
    expect_identical(tumor_threshold_txt(10, NA_real_, 1), "threshold +1: not evaluable")
    expect_identical(tumor_threshold_txt(NA_real_, 10, 1), "threshold +1: not evaluable")
    expect_identical(tumor_threshold_txt(12, 10, 1), "threshold +1: met")
    expect_identical(tumor_threshold_txt(10, 10, 1), "threshold +1: not met")
    expect_identical(tumor_threshold_txt(10, 10, 0), "threshold +0: met")
    session$setInputs(alpha = 0.35, g = 0, lambda_di = 0.23, lambda_dd = 0)
    expect_match(output$report, "threshold +1: not evaluable", fixed = TRUE)
    expect_false(grepl("threshold +1: met", output$report, fixed = TRUE))
    expect_match(output$report_notes$html, "It does not show uncertainty around the preset selected in the sidebar.", fixed = TRUE)
  })
})
test_that("uncertainty tab starts with both models", {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  ui <- env$unc_ui("unc")
  html <- htmltools::renderTags(ui)$html
  expect_match(html, '<option value="both" selected>', fixed = TRUE)
  expect_match(html, "Search within the toxicity limits", fixed = TRUE)
  note <- "These regimens meet the toxicity requirement, but they may still miss the tumor BED threshold."
  expect_match(html, note, fixed = TRUE)
  expect_lt(regexpr(note, html, fixed = TRUE)[1], regexpr('id="unc-search_table"', html, fixed = TRUE)[1])
})
test_that("source link uses the deployed commit", {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  f <- tempfile(); on.exit(unlink(f), add = TRUE)
  expect_match(env$deployed_revision_text(f), "Local copy", fixed = TRUE)
  rev <- "4876ae4000000000000000000000000000000000"
  writeLines(rev, f)
  expect_match(env$deployed_revision_text(f), paste0("/tree/", rev), fixed = TRUE)
  expect_match(env$deployed_revision_text(f), "[4876ae4]", fixed = TRUE)
  for (bad in c("", "4876ae4", "not-a-commit", "<script>")) {
    writeLines(bad, f)
    expect_match(env$deployed_revision_text(f), "Local copy", fixed = TRUE)
  }
})
