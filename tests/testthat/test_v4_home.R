# Version 4 main page: the list of hypofractionated regimens comparable to
# 35 x 2 Gy. The BED formulas are checked against an oracle in
# test_bed_models.R; these tests check the v4 selection, summary and page.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "v4_home.R"))
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
    n_bf <- if (is.null(bf)) 0 else nrow(bf)
    expect_equal(nrow(h$rows), n_bf, info = pr)
  }
})

test_that("pinned result at the defaults (scenario D18, calibrated DD)", {
  s <- v4_summary(v4_regimens(p0)$rows)
  expect_equal(s$m, 16:34)
  expect_equal(s$dose[s$m == 20], "2.65–2.80")
  expect_equal(s$total[s$m == 20], "53.0–56.0")
  expect_equal(s$T[s$m == 20], 26)
  expect_equal(s$saved[s$m == 20], 21)
  expect_equal(s$trial[s$m == 20], "HYPNO test arm, 55 Gy in 20")
  expect_equal(s$models[s$m == 20], "DI model only")
  expect_equal(s$models[s$m == 34], "both models")
  expect_equal(s$models[s$m == 22], "both models at 2.70 Gy, DI model only otherwise")
  expect_true(s$beyond[s$m == 16])
  expect_false(any(s$beyond[s$m >= 18]))
  expect_true(all(s$trial[s$m != 20] == ""))
})

test_that("summary wording follows which models qualify", {
  mk <- function(d, di, dd) data.frame(m = 20, d = d, T = 26, qual_di = di, qual_dd = dd)
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(TRUE, TRUE), c(TRUE, TRUE)))$models, "both models")
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(TRUE, TRUE), c(FALSE, FALSE)))$models, "DI model only")
  expect_equal(v4_summary(mk(c(2.7, 2.75), c(FALSE, FALSE), c(TRUE, TRUE)))$models, "DD model only")
  expect_equal(v4_summary(mk(c(2.7, 2.75, 2.8), c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE)))$models,
               "both models at 2.75–2.80 Gy, DI model only otherwise")
  expect_equal(v4_summary(mk(c(2.7, 2.75, 2.8), c(FALSE, TRUE, TRUE), c(TRUE, TRUE, FALSE)))$models,
               "both models at 2.75 Gy, DI or DD model only otherwise")
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
  expect_equal(nrow(v4_summary(h$rows)), 0)
})

# --- the page -----------------------------------------------------------------
v4_defaults <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3)

test_that("the main page lists the regimens and reacts to the model settings", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    expect_equal(home_rows()$m, 16:34)
    expect_match(as.character(output$home_settings$html), "defaults (tumor scenario D18)", fixed = TRUE)
    expect_match(as.character(output$home_msg$html), "19 schedules qualify, from 16 to 34 fractions", fixed = TRUE)
    session$setInputs(home_table_rows_selected = 5)  # 20 fractions
    det <- as.character(output$home_detail$html)
    expect_match(det, "20 fractions in 26 days", fixed = TRUE)
    expect_match(det, "Standard 35 × 2 Gy in 47 days", fixed = TRUE)
    tb <- output$home_detail_table
    expect_match(tb, "2.65", fixed = TRUE); expect_match(tb, "2.80", fixed = TRUE)
    expect_false(grepl("2.85", tb, fixed = TRUE))
    # The main list ignores the Report card comparison.
    session$setInputs(comparator = "hypno_test")
    expect_equal(home_rows()$m, 16:34)
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
