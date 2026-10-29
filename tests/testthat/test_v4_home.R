# Version 4 main page: the list of hypofractionated regimens comparable to
# 35 x 2 Gy. The BED formulas are checked against an oracle in
# test_bed_models.R; these tests check the v4 selection, summary and page.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)

p0 <- default_params()
std <- V4_STANDARD

# Brute force over the same grid, with the three conditions written out here.
brute_force <- function(p, margins = V4_MARGINS) {
  r0 <- compute_all(std$m, std$d, std$T, p, default_anchor())
  out <- list()
  for (m in V4_GRID$m_range) {
    T <- overall_time(m)
    if (T >= std$T) next
    for (k in 1:200) {
      d <- round(2 + k * V4_GRID$d_step, 6)
      if (d > V4_GRID$d_max + 1e-9) break
      a <- compute_all(m, d, T, p, default_anchor())
      di_ok <- is.finite(a$bed_di) && a$bed_di > 0 && a$bed_di >= r0$bed_di - margins$tumor - 1e-9
      dd_ok <- is.finite(a$bed_dd) && a$bed_dd > 0 && is.finite(r0$bed_dd) && a$bed_dd >= r0$bed_dd - margins$tumor - 1e-9
      if ((di_ok || dd_ok) && a$bed_early > 0 && a$bed_early <= r0$bed_early + margins$early + 1e-9 &&
          a$bed_late <= r0$bed_late + margins$late + 1e-9)
        out[[length(out) + 1]] <- c(m = m, d = d)
    }
  }
  do.call(rbind, out)
}

test_that("the standard is 35 x 2 Gy in 47 days on the Mon-Fri calendar", {
  expect_equal(overall_time(std$m), std$T)
  expect_equal(c(std$m, std$d), c(35, 2))
  expect_equal(unlist(V4_MARGINS), c(tumor = 5, early = 5, late = 5))
})

test_that("the list is exactly the regimens meeting all three conditions (defaults)", {
  h <- v4_regimens(p0)
  bf <- brute_force(p0)
  got <- h$rows[, c("m", "d")]
  expect_equal(nrow(got), nrow(bf))
  expect_equal(paste(got$m, sprintf("%.2f", got$d)), paste(bf[, "m"], sprintf("%.2f", bf[, "d"])))
  # Each listed row, checked against the standard's own BEDs.
  r0 <- h$standard
  x <- h$rows
  expect_true(all(x$bed_di >= r0$bed_di - 5 - 1e-9 | x$bed_dd >= r0$bed_dd - 5 - 1e-9))
  expect_true(all(x$bed_early <= r0$bed_early + 5 + 1e-9))
  expect_true(all(x$bed_late <= r0$bed_late + 5 + 1e-9))
  expect_true(all(x$T < std$T & x$d > std$d))
})

test_that("the list matches brute force under other tumor settings", {
  for (pr in c("refit", "hyhope", "calculator")) {
    p <- if (pr == "refit") {
      v <- jsonlite::fromJSON(file.path(app_abs, "draws", "refit_point_values.json"))$values
      modifyList(p0, list(alpha = v$di$alpha, ab = v$di$ab, lambda_di = v$di$lambda, Tk = v$di$Tk,
                          g = v$di$g, lambda_dd = v$dd$lambda))
    } else modifyList(p0, resolve_preset(pr))
    h <- v4_regimens(p)
    bf <- brute_force(p)
    expect_false(is.null(bf), info = pr)
    expect_equal(paste(h$rows$m, sprintf("%.2f", h$rows$d)), paste(bf[, "m"], sprintf("%.2f", bf[, "d"])), info = pr)
  }
})

test_that("pinned result at the defaults (scenario D18, calibrated DD)", {
  h <- v4_regimens(p0)
  hy <- compute_all(20, 2.75, 26, p0, default_anchor())$bed_early
  s <- v4_summary(h$rows, hypno_early = hy, data_d_max = 2.55)
  expect_equal(s$m, 16:34)
  r20 <- s[s$m == 20, ]
  expect_equal(r20$dose, "2.65–2.80"); expect_equal(r20$total, "53.0–56.0")
  expect_equal(c(r20$T, r20$saved), c(26, 21))
  expect_equal(r20$d_di, "−4.5 to +0.2")
  expect_match(r20$d_dd, "^−11.1 to −")
  expect_equal(r20$models, "DI only")
  expect_match(r20$notes, "outside trial data", fixed = TRUE)
  expect_match(r20$notes, "acute > HYPNO from 2.80 Gy", fixed = TRUE)
  expect_match(r20$notes, "HYPNO test arm (55 Gy in 20) inside", fixed = TRUE)
  expect_equal(s$models[s$m == 34], "both models")
  expect_equal(s$models[s$m == 22], "both at 2.70 Gy, else DI only")
  # Extrapolation: rows 16-20 lie wholly above 2.55 Gy, row 21 starts at it, rows from 26 fractions lie below it.
  expect_true(all(startsWith(s$notes[s$m <= 20], "outside trial data;")))
  expect_match(s$notes[s$m == 21], "outside trial data above 2.55 Gy", fixed = TRUE)
  expect_false(any(s$outside_data[s$m >= 26]))
  # The HYPNO arm itself is not "above HYPNO's" acute BED.
  x20 <- h$rows[h$rows$m == 20 & abs(h$rows$d - 2.75) < 1e-9, ]
  expect_equal(x20$bed_early, hy, tolerance = 1e-12)
  expect_equal(sum(nzchar(s$trial)), 1)
})

test_that("refit support: share of refit draws meeting the tumor condition", {
  di <- load_draws("di", dir = file.path(app_abs, "draws")); dd <- load_draws("dd", dir = file.path(app_abs, "draws"))
  regs <- data.frame(m = c(20, 21, 25), d = c(2.75, 2.55, 2.40)); regs$T <- overall_time(regs$m)
  sh <- v4_refit_support(regs, di, dd)
  expect_equal(round(sh, 2), c(0.95, 0.37, 1.00))
  # Independent count for one regimen.
  td <- tumor_deltas(regs[2, ], V4_STANDARD, di, dd)
  expect_equal(sh[2], mean((td$delta_di >= -5 & !is.na(td$delta_di)) | (td$delta_dd >= -5 & !is.na(td$delta_dd))))
  rows <- v4_regimens(p0)$rows
  rows$support <- v4_refit_support(rows, di, dd)
  s <- v4_summary(rows, data_d_max = 2.55)
  expect_match(s$notes[s$m == 21], "support < 80% at 2.55", fixed = TRUE)
  expect_true(all(grepl("^[0-9]+(–[0-9]+)?%$", s$support)))
})

test_that("summary wording follows which models qualify", {
  mk <- function(d, di, dd) data.frame(m = 20, d = d, T = 26, qual_di = di, qual_dd = dd,
                                       d_di = 0, d_dd = 0, bed_early = 50)
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(TRUE, TRUE), c(TRUE, TRUE)))$models, "both models")
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(TRUE, TRUE), c(FALSE, FALSE)))$models, "DI only")
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(FALSE, FALSE), c(TRUE, TRUE)))$models, "DD only")
  expect_equal(v4_summary(mk(c(2.7, 2.75, 2.8), c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE)))$models,
               "both at 2.75–2.80 Gy, else DI only")
  expect_equal(v4_summary(mk(c(2.7, 2.75, 2.8), c(FALSE, TRUE, TRUE), c(TRUE, TRUE, FALSE)))$models,
               "both at 2.75 Gy, else DI or DD only")
  expect_equal(v4_signed_range(c(-1.04, 0.33)), "−1.0 to +0.3")
  expect_equal(v4_signed_range(c(NA_real_, NaN)), "n/a")
  expect_equal(v4_summary(mk(2.7, TRUE, FALSE))$dose, "2.70")
  expect_equal(nrow(v4_summary(mk(2.7, TRUE, FALSE)[0, ])), 0)
})

test_that("wider margins list more regimens, zero margins fewer", {
  n <- function(mg) nrow(v4_regimens(p0, margins = list(tumor = mg, early = mg, late = mg))$rows)
  expect_lt(n(0), n(5))
  expect_lt(n(5), n(10))
})

test_that("an empty list carries the engine's reason", {
  p <- modifyList(p0, list(alpha = 0.01, lambda_di = 0.5, lambda_dd = 0.5))  # negative tumor BEDs
  h <- v4_regimens(p)
  expect_equal(nrow(h$rows), 0)
  expect_identical(h$reason, "negative_reference_both")
  expect_equal(nrow(v4_summary(h$rows)), 0)
})

# --- the page -----------------------------------------------------------------
v4_defaults <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE)

test_that("the main page lists the regimens and reacts to the model settings", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    expect_equal(home_rows()$m, 16:34)
    expect_match(as.character(output$home_settings$html), "not fitted to outcome data", fixed = TRUE)
    expect_match(as.character(output$home_msg$html), "19 schedules qualify, from 16 to 34 fractions", fixed = TRUE)
    expect_match(as.character(output$home_msg$html), "0.05 Gy", fixed = TRUE)
    expect_false(grepl("DD model is not available", as.character(output$home_msg$html), fixed = TRUE))
    tb <- output$home_table
    for (h in c("Tumor BED change, DI / DD (Gy₁₀)", "Refit support", "Notes")) expect_match(tb, h, fixed = TRUE)
    session$setInputs(home_table_rows_selected = 5)  # 20 fractions
    det <- as.character(output$home_detail$html)
    expect_match(det, "20 fractions in 26 days", fixed = TRUE)
    expect_match(det, "Standard 35 × 2 Gy in 47 days", fixed = TRUE)
    tb <- output$home_detail_table
    expect_match(tb, "2.65", fixed = TRUE); expect_match(tb, "2.80", fixed = TRUE)
    expect_false(grepl("2.85", tb, fixed = TRUE))
    expect_match(tb, "Refit support", fixed = TRUE)
    # Open in Report card: with the main-page rules the Report card agrees with the list.
    rc <- v4_report_inputs(20, 2.75)
    expect_equal(rc$min_gain, -5); expect_equal(c(rc$tol_early, rc$tol_late), c(5, 5))
    expect_identical(c(rc$early_mode, rc$late_mode, rc$comparator), c("reference", "reference", "custom"))
    do.call(session$setInputs, rc)
    r <- output$report
    expect_match(r, "threshold -5: met", fixed = TRUE)
    expect_false(grepl("exceeds", r, fixed = TRUE))
    # A dose just above the listed range breaks a toxicity limit there too.
    do.call(session$setInputs, v4_report_inputs(20, 2.85))
    expect_match(output$report, "exceeds", fixed = TRUE)
    do.call(session$setInputs, v4_defaults)
    # The main list ignores the Report card comparison.
    session$setInputs(comparator = "hypno_test")
    expect_equal(home_rows()$m, 16:34)
    # DD unavailable is announced.
    session$setInputs(lambda_dd = 0)
    expect_match(as.character(output$home_msg$html), "DD model is not available", fixed = TRUE)
    session$setInputs(lambda_dd = 0.4)
    # It follows the model settings.
    session$setInputs(ab_l = 2)
    expect_match(as.character(output$home_settings$html), "Model settings changed", fixed = TRUE)
    expect_false(identical(home_rows()$dose, v4_summary(v4_regimens(p0)$rows)$dose))
  })
})

test_that("the main page opens first and every menu page exists", {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  html <- htmltools::renderTags(env$ui)$html
  expect_match(html, "Hypofractionated alternatives to 70 Gy in 35 fractions", fixed = TRUE)
  expect_match(html, 'data-value="home"', fixed = TRUE)
  for (grp in env$V4_MENU) for (title in names(grp)) {
    expect_match(html, sprintf('id="go_%s"', grp[[title]]), fixed = TRUE, info = title)
    expect_match(html, sprintf('data-value="%s"', title), fixed = TRUE, info = title)
    expect_match(html, sprintf('id="back_%s"', grp[[title]]), fixed = TRUE, info = title)
  }
  # The home panel is the active one, and the menu starts closed.
  expect_match(html, 'class="tab-pane active" data-value="home"', fixed = TRUE)
  expect_match(html, 'data-open-desktop="closed"', fixed = TRUE)
  # No visible text still points to the old sidebar.
  vis <- gsub("<[^>]+>", " ", gsub("(?s)<script.*?</script>", " ", html, perl = TRUE))
  expect_false(grepl("sidebar", vis, ignore.case = TRUE))
})
