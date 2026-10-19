# Server-side tests of app.R, driven with shiny::testServer.

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)

defaults <- list(
  m = 20, d = 2.65, auto_T = TRUE, T = 26,
  ref_m = 35, ref_d = 2, ref_auto_T = TRUE, ref_T = 47,
  preset = "calculator", alpha = 0.2, ab = 10, lambda_di = 0.2, lambda_dd = 0.5, Tk = 28,
  ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "either", rank_by = "best", m_range = c(5, 34), min_gain = 1, min_days = 7,
  d_step = "0.05", d_max = 8, early_mode = "reference", tol_early = 0, tol_late = 0)

set_defaults <- function(session, ...) {
  args <- modifyList(defaults, list(...))
  do.call(session$setInputs, args)
}

err_msg <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

test_that("default inputs render results, comparisons and suggestions", {
  shiny::testServer(app_dir, {
    set_defaults(session)
    b <- both()
    expect_equal(round(b$bed_di, 2), c(67.04, 65.00))
    expect_true(grepl("Tumor BED, DI model", output$results))
    expect_true(grepl("HYPNO test arm", output$benchmarks))
    expect_true(grepl("55.87", output$benchmarks))   # HYPNO acute BED, elapsed clock
    s <- suggestions()
    expect_equal(nrow(s), 7)
    expect_true(all(s$d_di >= 1 - 1e-9))
    expect_true(all(s$T <= 40))
    v <- output$verdicts$html
    expect_true(grepl("higher under DI", v) && grepl("lower under DD", v))
    expect_true(grepl(">\\s*26\\s*<", output$results))
    expect_false(grepl("26\\.00", output$results))
    expect_true(grepl("+0.9 Gy under DI and -4.7 Gy under DD", output$hypno_note$html, fixed = TRUE))
    expect_true(grepl("55.9 Gy", output$hypno_cap_note$html) && grepl("above the current early limit of 53.1", output$hypno_cap_note$html))
    expect_false(grepl("bg-success|bg-danger", v))   # no green/red safety badges
    expect_true(grepl("59\u201361 Gy", output$acute_zone$html))
  })
})

test_that("DD rate 0: warning on both tabs, no NaN, no NA rows", {
  shiny::testServer(app_dir, {
    set_defaults(session, lambda_dd = 0, tol_early = 5)
    expect_true(grepl("DD model unavailable", output$dd_warning$html))
    expect_true(grepl("DD model unavailable", output$dd_warning2$html))
    expect_true(grepl("not available", output$results))
    expect_false(grepl("NaN", output$results))
    expect_true(grepl("not available under DD", output$verdicts$html))
    s <- suggestions()
    expect_false(any(is.na(s$m)))
    session$setInputs(criterion = "dd")
    expect_true(grepl("DD model is unavailable", output$suggest_msg$html))
  })
})

test_that("DI rate 0 is treated as no repopulation, without a warning", {
  shiny::testServer(app_dir, {
    set_defaults(session, lambda_di = 0)
    expect_null(output$dd_warning)
    b <- both()
    expect_equal(b$bed_dd, b$bed_di)
  })
})

test_that("inputs are validated with readable messages", {
  shiny::testServer(app_dir, {
    set_defaults(session, m = 20.5)
    expect_match(err_msg(regimen()), "Number of fractions must be a whole number")
    set_defaults(session, lambda_dd = 9)
    expect_match(err_msg(params()), "DD maximum repopulation rate must be between")
    set_defaults(session, Tp_e = NA)
    expect_match(err_msg(params()), "Mucosal doubling time")
  })
})

test_that("typed overall times cannot be shorter than the calendar allows", {
  shiny::testServer(app_dir, {
    set_defaults(session, auto_T = FALSE, T = 1)
    expect_match(err_msg(regimen()), "need at least 26 days")
    set_defaults(session, auto_T = FALSE, T = 30)
    expect_equal(regimen()$T, 30)
    set_defaults(session, auto_T = FALSE, T = 26)   # exactly the calendar minimum
    expect_equal(regimen()$T, 26)
    set_defaults(session, ref_auto_T = FALSE, ref_T = 40)
    expect_match(err_msg(ref()), "need at least 47 days")
  })
})

test_that("reference overall time follows the reference fractions", {
  shiny::testServer(app_dir, {
    set_defaults(session, ref_m = 30)
    expect_equal(ref()$T, 40)
    s <- suggestions()
    expect_true(all(s$m < 30))
    expect_true(all(s$T <= 40 - 7))
  })
})

test_that("Fowler limit option and chemotherapy change the early cap", {
  shiny::testServer(app_dir, {
    set_defaults(session, early_mode = "fowler")
    expect_equal(search()$cap, 59)
    s <- suggestions()
    expect_true(all(s$bed_early <= 59 + 1e-9))
    session$setInputs(chemo = TRUE)
    expect_equal(search()$cap, 54)
    expect_true(grepl("54.58 Gy", output$acute_zone$html))
  })
})

test_that("search inputs reach the search", {
  shiny::testServer(app_dir, {
    set_defaults(session, d_max = 2.5, min_gain = 0, min_days = 1, tol_early = 3, d_step = "0.1")
    s <- suggestions()
    expect_gt(nrow(s), 0)
    expect_true(all(s$d <= 2.5 + 1e-9))
    expect_true(all(abs(s$d * 10 - round(s$d * 10)) < 1e-9))
    set_defaults(session, min_gain = 3)
    expect_true(all(suggestions()$d_di >= 3 - 1e-9))
  })
})

read_zip <- function(f) {
  d <- tempfile(); utils::unzip(f, exdir = d)
  list(res = read.csv(file.path(d, "results.csv"), check.names = FALSE),
       all = read.csv(file.path(d, "all_qualifying.csv"), check.names = FALSE),
       sc = jsonlite::fromJSON(file.path(d, "scenario.json")),
       files = list.files(d))
}

test_that("download is a ZIP with full-precision results and a rerunnable scenario", {
  shiny::testServer(app_dir, {
    set_defaults(session)
    z <- read_zip(output$download)
    expect_setequal(z$files, c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt"))
    expect_equal(nrow(z$res), 7)
    expect_gte(nrow(z$all), nrow(z$res))
    s <- suggestions()
    expect_equal(z$res$bed_dd, s$bed_dd, tolerance = 1e-12)       # not rounded
    expect_equal(z$sc$anchor$m, 35)
    expect_equal(z$sc$clocks$acute, "elapsed")
    expect_true(nchar(z$sc$scenario_id) == 12)
    env <- new.env(); source(file.path(app_abs, "R", "bed_models.R"), local = env)
    source(file.path(app_abs, "R", "export.R"), local = env)
    r2 <- env$rerun_scenario(z$sc)
    expect_equal(r2$bed_di, z$res$bed_di, tolerance = 1e-8)
    expect_equal(r2$d, z$res$d, tolerance = 1e-12)
  })
})

test_that("every preset in the menu exists, and each one calibrates", {
  app_env <- new.env()
  old <- setwd(app_dir)
  on.exit(setwd(old))
  source("app.R", local = app_env)
  menu <- setdiff(unname(app_env$PRESET_LABELS), "custom")
  expect_setequal(menu, names(app_env$tumor_presets()))
  presets <- app_env$tumor_presets()
  setwd(old)
  shiny::testServer(app_dir, {
    for (nm in menu) {
      pr <- presets[[nm]]
      set_defaults(session, alpha = pr$alpha, lambda_di = pr$lambda_di, Tk = pr$Tk)
      expect_null(output$dd_warning)
      expect_true(is.finite(both()$bed_dd[1]))
    }
  })
})

test_that("allowed early excess reaches the search", {
  shiny::testServer(app_dir, {
    set_defaults(session, tol_early = 3)
    s <- suggestions()
    expect_true(any(s$d_early > 0))
    expect_true(all(s$d_early <= 3 + 1e-9))
  })
})

test_that("a regimen equal to the reference is reported as the same", {
  shiny::testServer(app_dir, {
    set_defaults(session, m = 35, d = 2)
    v <- output$verdicts$html
    expect_true(grepl("the same under DI", v) && grepl("the same under DD", v))
  })
})

test_that("changed acute parameters disable the Fowler zone and limit", {
  shiny::testServer(app_dir, {
    set_defaults(session, Tp_e = 10)
    expect_true(grepl("have been changed", output$acute_zone$html))
    set_defaults(session, Tp_e = 10, early_mode = "fowler")
    expect_true(search()$blocked)
    expect_true(grepl("applies only with Fowler", output$suggest_msg$html))
    expect_null(output$search_note)        # the message appears once, not on every output
    expect_null(output$hypno_cap_note)
  })
})

test_that("worst-case ranking with an unavailable DD model explains itself", {
  shiny::testServer(app_dir, {
    set_defaults(session, lambda_dd = 0, rank_by = "worst")
    expect_true(grepl("smaller of", output$suggest_msg$html))
  })
})

test_that("negative BEDs are flagged", {
  shiny::testServer(app_dir, {
    set_defaults(session)
    expect_null(output$neg_warning)
    set_defaults(session, Tp_e = 0.1)
    expect_true(grepl("negative", output$neg_warning$html))
  })
})

test_that("early BED unit follows the mucosal alpha/beta", {
  shiny::testServer(app_dir, {
    set_defaults(session, ab_e = 5)
    expect_true(grepl("Gy (\u03b1/\u03b2 = 5 Gy)", output$hypno_cap_note$html, fixed = TRUE))
    expect_false(grepl("This regimen: [0-9.]+ Gy\u2081\u2080", output$acute_zone$html))
  })
})

test_that("empty-search hint does not suggest the criterion already chosen", {
  shiny::testServer(app_dir, {
    set_defaults(session, alpha = 0.2, lambda_di = 0.144, Tk = 28)   # Withers preset: no rows
    expect_equal(nrow(suggestions()), 0)
    expect_false(grepl("require the change under", output$suggest_msg$html))
    set_defaults(session, alpha = 0.2, lambda_di = 0.144, Tk = 28, criterion = "di")
    expect_true(grepl("require the change under", output$suggest_msg$html))
  })
})

test_that("an unavailable DD column reads 'not available'", {
  shiny::testServer(app_dir, {
    set_defaults(session, lambda_dd = 0, tol_early = 5)
    x <- pretty_suggestions()
    expect_gt(nrow(x), 0)
    expect_true(all(x[["\u0394 BED DD"]] == "not available"))
    expect_true(all(is.na(read_zip(output$download)$res$d_dd)))   # NA, never zero
  })
})

test_that("comparison regimen does not move the DD model unless legacy coupling is chosen", {
  shiny::testServer(app_dir, {
    set_defaults(session, m = 20, d = 2.75)
    dd_std <- both()$bed_dd[1]
    expect_equal(dd_std, 58.02143210, tolerance = 1e-8)
    session$setInputs(comparator = "hypno_control")
    expect_equal(ref()$T, 38)
    expect_equal(both()$bed_dd[1], dd_std, tolerance = 1e-12)       # same DD model
    session$setInputs(coupled = TRUE)
    expect_false(isTRUE(all.equal(both()$bed_dd[1], dd_std)))        # legacy: model moves
    expect_true(grepl("Legacy mode", output$anchor_note$html))
  })
})

test_that("HYPNO acute and late limits reach the search", {
  shiny::testServer(app_dir, {
    set_defaults(session, early_mode = "hypno", late_mode = "hypno", min_gain = -5)
    expect_equal(search()$cap, 55.86597229, tolerance = 1e-8)
    expect_equal(search()$late_cap, 105.41666667, tolerance = 1e-8)
    s <- suggestions()
    expect_gt(nrow(s), 0)
    expect_true(all(s$bed_early <= 55.86597229 + 1e-9 & s$bed_late <= 105.41666667 + 1e-9))
    expect_true(all(s$d_di >= -5 - 1e-9 | s$d_dd >= -5 - 1e-9))
  })
})

test_that("a negative minimum change (non-inferiority margin) is accepted", {
  shiny::testServer(app_dir, {
    set_defaults(session, min_gain = -2)
    expect_equal(search_args()$min_gain, -2)
    set_defaults(session, min_gain = -25)
    expect_match(err_msg(search_args()), "between -20 and 30")
  })
})

test_that("results show late BED at alpha/beta 2 and 4, and the sensitivity tab renders", {
  shiny::testServer(app_dir, {
    set_defaults(session, m = 20, d = 2.9)
    expect_true(grepl("Late toxicity BED at \u03b1/\u03b2 2 Gy", output$results))
    expect_true(grepl("Late toxicity BED at \u03b1/\u03b2 4 Gy", output$results))
    st <- sens()
    expect_equal(nrow(st), 5 * 3)
    expect_setequal(unique(st$preset), c("calculator", "withers", "hartley", "hyhope", "d18"))
    expect_true(grepl("scenarios under DI", output$sens_summary$html))
    expect_true(grepl("+2.10", output$sens_tox, fixed = TRUE))      # 142.1 - 140 at alpha/beta 2
    expect_true(grepl("-2.60", output$sens_tox, fixed = TRUE))      # 114.07 - 116.67 at alpha/beta 3
  })
})

test_that("search table shows which model credits each row and flags large fractions", {
  shiny::testServer(app_dir, {
    set_defaults(session, tol_early = 10, tol_late = 10, d_max = 4)
    x <- pretty_suggestions()
    expect_true(all(c("Credited by", "Beyond precedent") %in% names(x)))
    s <- suggestions()
    expect_equal(x[["Beyond precedent"]] == "yes", s$d > 3.1 + 1e-9)
  })
})

test_that("HYHOPE benchmark row and HYPNO context are shown", {
  shiny::testServer(app_dir, {
    set_defaults(session)
    expect_true(grepl("HYHOPE", output$benchmarks))
    expect_true(grepl("50.7%", output$hypno_note$html) && grepl("10-percentage-point", output$hypno_note$html))
  })
})

test_that("Fenwick acute limit reaches the search and the Results tab", {
  shiny::testServer(app_dir, {
    set_defaults(session, early_mode = "fenwick", min_gain = -20, tol_late = 30, d_max = 4)
    expect_true(search_args()$fenwick)
    expect_null(search_args()$early_cap)
    s <- suggestions()
    expect_gt(nrow(s), 0)
    ck <- s$m * s$d * (1 + s$d / 10); lim <- 69.5 * (s$T / 32.2) / sin(s$T / 32.2) - 3.5
    expect_true(all(ck <= lim + 1e-9))
    expect_true(grepl("Fenwick", output$search_note$html))
    expect_null(output$hypno_cap_note)
    set_defaults(session, m = 20, d = 2.9)
    expect_true(grepl("this regimen 74.8 vs 74.2 Gy\u2081\u2080, above", output$acute_zone$html, fixed = TRUE))
    z <- read_zip(output$download)
    expect_true(isTRUE(z$sc$search$fenwick) || identical(z$sc$search$fenwick, FALSE))
  })
})

test_that("model agreement tab lists the regimens both models support against HYPNO", {
  shiny::testServer(app_dir, {
    set_defaults(session, alpha = 0.35, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28,
                 early_mode = "hypno", d_max = 5, agree_cmp = "hypno_test", agree_margin = -3, agree_days = 1)
    a <- agreement()
    expect_equal(a$m, c(17, 18, 19))
    expect_equal(a$d, c(3.05, 2.90, 2.80), tolerance = 1e-12)
    tab <- output$agree_table
    expect_true(grepl("15 of 15", tab, fixed = TRUE))
    expect_true(grepl("-2.79", tab, fixed = TRUE) && grepl("-1.96", tab, fixed = TRUE))
    session$setInputs(agree_margin = 0)
    expect_equal(nrow(agreement()), 0)
    expect_true(grepl("No regimen meets the margin", output$agree_msg$html))
    session$setInputs(agree_cmp = "sidebar", agree_margin = -3, agree_days = 7)
    expect_true(all(agreement()$T <= 47 - 7))
  })
})
