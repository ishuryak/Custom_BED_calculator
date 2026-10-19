# Version 2.5: fixes from the 2026-10-19 review of version 2.4
# (docs/development/2026_10_19_review_and_plan.md, findings F01-F08, N1-N6).
# Expected values are independent: they come from the probe script
# docs/development/2026_10_19_review_scripts/probe_v24.R and from the external
# review's oracle (docs/development/2026_10_19_external_review_codex.md).

app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
models_file <- Sys.getenv("BED_MODELS_FILE", file.path(app_abs, "R", "bed_models.R"))
source(models_file, local = TRUE)
source(file.path(app_abs, "R", "export.R"), local = TRUE)

p <- default_params()
an <- list(m = 35, d = 2, T = 47)
hypno <- list(m = 20, d = 2.75, T = 26)
ref150 <- list(m = 35, d = 2, T = 150)      # stress case of the external review (F01)
q0 <- modifyList(p, list(lambda_di = 0))
cap_hypno <- function(pp) compute_all(20, 2.75, 26, pp, an)$bed_early
read_text <- function(f) paste(readLines(file.path(app_abs, f), encoding = "UTF-8", warn = FALSE), collapse = "\n")

# --- F01: one validity rule ------------------------------------------------

test_that("F01: agreement counts only scenarios in which both comparisons are evaluable", {
  a <- agreement_regimens(q0, ref150, an, margin = -3, early_cap = cap_hypno(q0),
                          min_days_saved = 1, d_max = 5, m_range = 1:34)
  r <- a[a$m == 32, ]
  expect_equal(nrow(r), 1)
  expect_equal(r$n_scenarios, 15)
  expect_equal(r$n_valid_both, 3)
  expect_equal(r$n_invalid, 12)
  expect_equal(r$n_pass_both, 3)
  # The +103 Gy "gain" against a -38 Gy comparison is excluded. The largest
  # evaluable DI change is +67.33 Gy (HYHOPE at alpha/beta 7, comparison
  # 4.84 Gy): a positive BED is a validity floor, not a plausibility check.
  expect_equal(r$d_di_max, 67.33486, tolerance = 1e-5)
  expect_lt(r$d_di_max, 100)
})

test_that("F01: the default agreement rows keep 15 evaluable scenarios out of 15", {
  a <- agreement_regimens(p, hypno, an, margin = -3, early_cap = cap_hypno(p),
                          min_days_saved = 1, d_max = 5, m_range = 1:19)
  expect_equal(a$m, c(17, 18, 19))
  expect_equal(a$n_valid_both, rep(15, 3))
  expect_equal(a$n_pass_both, rep(15, 3))
  expect_equal(a$n_invalid, rep(0, 3))
})

test_that("F01: sensitivity table hides deltas against a negative comparison BED", {
  st <- sensitivity_table(list(m = 32, d = 2.10, T = 44), ref150, an, q0)
  expect_equal(sum(st$valid_di & st$valid_dd), 3)
  neg <- st$ref_di <= 0
  expect_equal(sum(neg), 12)
  expect_true(all(is.na(st$d_di[neg])))
  expect_true(all(st$reason_di[neg] == "negative_reference"))
  expect_equal(st$ref_di[st$preset == "calculator" & st$ab == 10], -38, tolerance = 1e-9)
})

test_that("F01: endpoint rule covers negative candidate, zero, unavailable and all-invalid", {
  e <- tumor_endpoint(c(-1, 0, 5, NA), 3)
  expect_equal(e$valid, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(e$reason, c("negative_candidate", "negative_candidate", "ok", "negative_candidate"))
  expect_equal(e$delta, c(NA, NA, 2, NA))
  e2 <- tumor_endpoint(c(5, 6), NA_real_, available = FALSE)
  expect_equal(e2$reason, c("dd_unavailable", "dd_unavailable"))
  expect_true(all(is.na(e2$delta)))
  e3 <- tumor_endpoint(5, -2)
  expect_equal(e3$reason, "negative_reference")
  st <- data.frame(d_di = c(NA, NA), d_dd = c(NA, NA), valid_di = c(FALSE, FALSE),
                   valid_dd = c(FALSE, FALSE))
  s <- scenario_summary(st, margin = -3)
  expect_equal(s$n_valid_both, 0)
  expect_equal(s$n_pass_both, 0)
  expect_true(is.na(s$d_di_min) && is.na(s$d_dd_max))
})

test_that("F01: Sensitivity tab and its summary count evaluable scenarios only", {
  shiny::testServer(app_dir, {
    session$setInputs(m = 32, d = 2.1, auto_T = TRUE, comparator = "custom", ref_m = 35, ref_d = 2,
                      ref_auto_T = FALSE, ref_T = 150, preset = "d18", alpha = 0.35, ab = 10,
                      lambda_di = 0, lambda_dd = 0.4, Tk = 28, g = 0, ab_e = 10, alpha_e = 0.35,
                      Tp_e = 2.5, Tk_e = 7, ab_l = 3, min_gain = 1)
    txt <- output$sens_summary$html
    expect_true(grepl("12 of 15 scenarios are not evaluable", txt, fixed = TRUE))
    expect_true(grepl("not evaluable", output$sens_table))
    expect_false(grepl("+103", output$sens_table, fixed = TRUE))
  })
})

# --- F02 and N1/N2: claims ----------------------------------------------------

withdrawn <- c("within the uncertainty of the model parameters", "within the parameter uncertainty",
               "pessimistic bound", "close to HYPNO in tumor effect", "as a non-inferiority margin",
               "as in a non-inferiority design", "can only lose dose", "cannot separate regimens 1",
               "does not support escalation", "was needed before any sparing")

test_that("F02: withdrawn claims are gone from the About page and the interface", {
  for (f in c("about.md", "app.R", "README.md")) {
    txt <- read_text(f)
    for (w in withdrawn) expect_false(grepl(w, txt, fixed = TRUE), info = paste(f, w))
  }
})

test_that("N1: the D18 HYPNO statement is conditional on a cited dose-response slope", {
  txt <- read_text("about.md")
  expect_true(grepl("3.2 (95% CI 1.6", txt, fixed = TRUE))
  expect_true(grepl("Stuschke", txt, fixed = TRUE))
})

test_that("F02: before the DI onset, background growth still rewards a shorter course", {
  x <- function(g) {
    a <- compute_all(c(17, 20), c(3.05, 2.75), c(23, 26), modifyList(p, list(g = g)), an)
    a$bed_di[1] - a$bed_di[2]
  }
  expect_equal(x(0), -2.46020777, tolerance = 1e-7)
  expect_equal(x(0.0154) - x(0), 0.0154 * 3 / 0.35, tolerance = 1e-9)
  expect_true(grepl("background growth", read_text("about.md"), fixed = TRUE))
  shiny::testServer(app_dir, {
    session$setInputs(m = 20, d = 2.65, auto_T = TRUE, comparator = "custom", ref_m = 35, ref_d = 2,
                      ref_auto_T = TRUE, preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23,
                      lambda_dd = 0.4, Tk = 28, g = 0, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5,
                      Tk_e = 7, ab_l = 3, agree_cmp = "hypno_test")
    expect_false(grepl("only lose", output$agree_note$html, fixed = TRUE))
    expect_false(grepl("Background growth g still gives", output$agree_note$html, fixed = TRUE))
    session$setInputs(g = 0.0154)
    expect_true(grepl("Background growth g still gives 0.044 Gy\u2081\u2080 per day saved", output$agree_note$html, fixed = TRUE))
  })
})

# --- F03 and N6: independent tabs -------------------------------------------

server_defaults <- list(
  m = 20, d = 2.65, auto_T = TRUE, T = 26, comparator = "custom",
  ref_m = 35, ref_d = 2, ref_auto_T = TRUE, ref_T = 47,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0,
  anchor_m = 35, anchor_d = 2, coupled = FALSE,
  ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "either", rank_by = "best", m_range = c(5, 34), min_gain = 1, min_days = 7,
  d_step = "0.05", d_max = 5, early_mode = "hypno", late_mode = "reference",
  tol_early = 0, tol_late = 0, agree_cmp = "hypno_test", agree_margin = -3, agree_days = 1)
set_srv <- function(session, ...) do.call(session$setInputs, modifyList(server_defaults, list(...)))

test_that("F03: agreement does not depend on candidate-only inputs", {
  shiny::testServer(app_dir, {
    set_srv(session)
    base <- agreement()
    expect_equal(base$m, c(17, 18, 19))
    set_srv(session, min_gain = NA); expect_equal(agreement()$m, base$m)
    set_srv(session, min_days = NA); expect_equal(agreement()$m, base$m)
    set_srv(session, m_range = c(36, 40)); expect_equal(agreement()$m, base$m)
    set_srv(session, criterion = "dd", rank_by = "dd"); expect_equal(agreement()$d_dd, base$d_dd)
  })
})

test_that("F03: the candidate export does not depend on the Results regimen", {
  shiny::testServer(app_dir, {
    set_srv(session, m = NA)
    f <- output$download
    expect_true(file.exists(f))
    d <- tempfile(); utils::unzip(f, exdir = d)
    expect_equal(nrow(read.csv(file.path(d, "results.csv"))), 3)
  })
})

test_that("F03: in fixed-anchor mode the sidebar comparison does not move the agreement", {
  shiny::testServer(app_dir, {
    set_srv(session)
    a1 <- agreement()
    set_srv(session, comparator = "hypno_control")
    a2 <- agreement()
    expect_equal(a2$d_di, a1$d_di, tolerance = 1e-12)
    expect_equal(a2$d_dd, a1$d_dd, tolerance = 1e-12)
    expect_equal(a2$bed_dd, a1$bed_dd, tolerance = 1e-12)
  })
})

test_that("N6: legacy coupling calibrates the agreement at its own comparison and says so", {
  shiny::testServer(app_dir, {
    set_srv(session, coupled = TRUE)
    expect_equal(agree_anchor()[c("m", "d", "T")], list(m = 20, d = 2.75, T = 26))
    expect_true(grepl("20 × 2.75 Gy in 26 days", output$agree_note$html, fixed = TRUE))
    a <- agreement()
    if (nrow(a) > 0) {
      # DD equals DI at the anchor, which is the agreement comparison
      x <- compute_all(20, 2.75, 26, p, list(m = 20, d = 2.75, T = 26))
      expect_equal(x$bed_dd, x$bed_di, tolerance = 1e-9)
    }
  })
})

test_that("F03: an empty agreement explains its own reason", {
  shiny::testServer(app_dir, {
    set_srv(session, lambda_dd = 0)          # DD unavailable at these values
    expect_equal(nrow(agreement()), 0)
    msg <- output$agree_msg$html
    expect_true(grepl("DD model is unavailable", msg, fixed = TRUE))
    expect_false(grepl("See the candidate tab", msg, fixed = TRUE))
  })
})

# --- F04 and F05: export and replay ---------------------------------------

unzip_all <- function(f) {
  d <- tempfile(); utils::unzip(f, exdir = d)
  list(dir = d, files = list.files(d), sc = jsonlite::fromJSON(file.path(d, "scenario.json")))
}

test_that("F04: the agreement analysis downloads and replays in a fresh process", {
  shiny::testServer(app_dir, {
    set_srv(session)
    z <- unzip_all(output$download_agree)
    expect_true(all(c("results.csv", "all_qualifying.csv", "scenario_detail.csv", "scenario.json",
                      "README.txt") %in% z$files))
    expect_equal(z$sc$analysis, "agreement")
    expect_equal(z$sc$comparison$m, 20)
    expect_equal(z$sc$agreement$margin, -3)
    expect_equal(nrow(z$sc$scenario_set), 15)
    det <- read.csv(file.path(z$dir, "scenario_detail.csv"))
    expect_equal(nrow(det), 3 * 15)
    res <- read.csv(file.path(z$dir, "results.csv"))
    rr <- callr::r(function(root, js) {
      source(file.path(root, "R", "bed_models.R")); source(file.path(root, "R", "export.R"))
      rerun_scenario(jsonlite::fromJSON(js))
    }, args = list(root = app_abs, js = file.path(z$dir, "scenario.json")))
    expect_equal(rr$m, res$m)
    expect_equal(rr$d_di, res$d_di, tolerance = 1e-8)
    expect_equal(rr$d_dd, res$d_dd, tolerance = 1e-8)
    expect_equal(rr$n_pass_both, res$n_pass_both)
    expect_equal(attr(rr, "scenario_detail")$d_dd, det$d_dd, tolerance = 1e-8)
  })
})

test_that("F04: download names carry the app version and the analysis", {
  # testServer exposes the content file, not the filename, so the filename
  # builder is a plain function that is tested directly.
  app_env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = app_env)
  f <- app_env$export_filename("agreement", "abc123def456", as.Date("2026-10-19"))
  expect_equal(f, sprintf("bed_v%s_agreement_abc123def456_2026-10-19.zip", APP_VERSION))
  expect_false(grepl("bed_v2.1", read_text("app.R"), fixed = TRUE))
})

cand_scenario <- function() {
  build_scenario(p, default_ref(), an,
                 search = list(m_range = 5:34, d_step = 0.05, d_max = 5, criterion = "either",
                               rank_by = "best", tol_early = 0, tol_late = 0, early_cap = cap_hypno(p),
                               late_cap = NULL, min_gain = 1, min_days_saved = 7, fenwick = FALSE),
                 context = list())
}

test_that("F05: replay rejects unsupported model versions, hinge and clocks", {
  sc <- jsonlite::fromJSON(scenario_json(cand_scenario()))
  expect_equal(nrow(rerun_scenario(sc)), 3)
  bad <- function(f) { x <- sc; f(x) }
  expect_error(rerun_scenario(bad(function(x) { x$equation_version <- "9.9"; x })), "equation version")
  expect_error(rerun_scenario(bad(function(x) { x$schema_version <- "99"; x })), "schema")
  expect_error(rerun_scenario(bad(function(x) { x$selection_version <- "0.1"; x })), "selection version")
  expect_error(rerun_scenario(bad(function(x) { x$hinge$type <- "hard"; x })), "hinge")
  expect_error(rerun_scenario(bad(function(x) { x$hinge$k_per_day <- 10; x })), "hinge")
  expect_error(rerun_scenario(bad(function(x) { x$clocks$tumor <- "elapsed"; x })), "clock")
  expect_error(rerun_scenario(bad(function(x) { x$clocks$acute <- "weeks"; x })), "clock")
  expect_error(rerun_scenario(bad(function(x) { x$params$alpha <- NULL; x })), "alpha")
  expect_error(rerun_scenario(bad(function(x) { x$params$alpha <- -1; x })), "alpha")
  expect_error(rerun_scenario(bad(function(x) { x$params$ab <- "ten"; x })), "ab")
  expect_error(rerun_scenario(bad(function(x) { x$search$criterion <- "any"; x })), "criterion")
  expect_error(rerun_scenario(bad(function(x) { x$calibration$c_over_alpha <- 1; x })), "calibration")
})

test_that("F05: schema 1 candidate files are migrated explicitly", {
  sc <- jsonlite::fromJSON(scenario_json(cand_scenario()))
  sc$schema_version <- "1"; sc$analysis <- NULL; sc$selection_version <- NULL
  sc$scenario_set_version <- NULL
  r <- rerun_scenario(sc)
  expect_equal(nrow(r), 3)
  expect_equal(attr(r, "migrated_from"), "1")
})

test_that("F05: the scenario ID ignores key order but not scientific inputs", {
  sc <- cand_scenario()
  rev_sc <- sc; rev_sc$params <- rev(rev_sc$params); rev_sc$search <- rev(rev_sc$search)
  rev_sc <- rev_sc[rev(names(rev_sc))]
  expect_identical(scenario_id(rev_sc), scenario_id(sc))
  sc2 <- sc; sc2$params$alpha <- 0.36
  expect_false(identical(scenario_id(sc2), scenario_id(sc)))
  sc3 <- sc; sc3$search$m_range <- rev(sc3$search$m_range)   # arrays keep their order
  expect_false(identical(scenario_id(sc3), scenario_id(sc)))
  rt <- jsonlite::fromJSON(scenario_json(sc))
  expect_identical(scenario_id(rt), scenario_id(sc))
  sc4 <- sc; sc4$params$g <- 0.01
  expect_false(identical(scenario_id(sc4), scenario_id(sc)))
})

# --- F06, N3-N5: units, labels, provenance --------------------------------

test_that("F06: the HYPNO acute cap is labeled with the mucosal alpha/beta in use", {
  shiny::testServer(app_dir, {
    set_srv(session, ab_e = 5)
    txt <- output$search_note$html
    expect_true(grepl("71.0 Gy (α/β = 5 Gy)", txt, fixed = TRUE))
    expect_false(grepl("71.0 Gy₁₀", txt, fixed = TRUE))
    set_srv(session)
    expect_true(grepl("55.9 Gy₁₀", output$search_note$html, fixed = TRUE))
  })
})

test_that("F06: Fenwick is not described as an elapsed-day model", {
  src <- read_text("R/bed_models.R")
  expect_false(grepl("Fowler, Fenwick, Hartley and Meade count elapsed", src, fixed = TRUE))
})

test_that("F06: the fraction-size flag names the benchmark set literally", {
  shiny::testServer(app_dir, {
    set_srv(session, early_mode = "reference", tol_early = 10, tol_late = 10, min_gain = -20, d_max = 4)
    x <- pretty_suggestions()
    expect_true("Above benchmark fraction sizes (> 3.1 Gy)" %in% names(x))
    expect_false("Beyond precedent" %in% names(x))
  })
})

test_that("N3: an infinite DD onset is shown as 'none'", {
  shiny::testServer(app_dir, {
    set_srv(session, lambda_di = 0)
    expect_gt(nrow(agreement()), 0)
    tab <- output$agree_table
    expect_false(grepl("Inf", tab, fixed = TRUE))
    expect_true(grepl("none", tab, fixed = TRUE))
  })
})

test_that("Wide tables scroll horizontally instead of being cut off (walkthrough finding)", {
  src <- read_text("app.R")
  for (id in c("agree_table", "sens_table", "results", "benchmarks"))
    expect_true(grepl(sprintf('overflow-x: auto;", tableOutput("%s")', id), src, fixed = TRUE), info = id)
})

test_that("F01: the Results and benchmark tables apply the validity rule too", {
  shiny::testServer(app_dir, {
    set_srv(session, preset = "custom", alpha = 0.2, lambda_di = 2)   # comparison DI BED -106 Gy
    b <- both()
    expect_lt(b$bed_di[2], 0)
    r <- output$results
    expect_true(grepl("not evaluable (BED \u2264 0)", r, fixed = TRUE))
    expect_false(grepl(sprintf("%+.2f", b$bed_di[1] - b$bed_di[2]), r, fixed = TRUE))
    expect_true(grepl("not evaluable", output$benchmarks, fixed = TRUE))
    expect_true(grepl("not evaluable", output$verdicts$html, fixed = TRUE))
    set_srv(session)
    expect_false(grepl("not evaluable", output$results, fixed = TRUE))
  })
})

test_that("N5: Results rows carry BED units", {
  shiny::testServer(app_dir, {
    set_srv(session)
    r <- output$results
    expect_true(grepl("Early toxicity BED (Gy₁₀)", r, fixed = TRUE))
    expect_true(grepl("Late toxicity BED (Gy₃)", r, fixed = TRUE))
  })
})

test_that("F06: the sensitivity record states inherited values", {
  st <- sensitivity_table(list(m = 20, d = 2.75, T = 26), default_ref(), an, p)
  expect_true(all(c("lambda_dd_inherited", "g", "anchor_m") %in% names(st)))
  expect_true(all(st$lambda_dd_inherited[st$preset %in% c("withers", "hartley", "hyhope")]))
  expect_false(any(st$lambda_dd_inherited[st$preset %in% c("calculator", "d18")]))
})

# --- F07: records and versions ---------------------------------------------

test_that("F07: versions agree and source statements are current", {
  desc <- read.dcf(file.path(app_abs, "DESCRIPTION"))
  expect_equal(unname(desc[1, "Version"]), APP_VERSION)
  for (f in c("R/bed_models.R", "R/export.R")) expect_false(grepl("version 2.1", read_text(f), fixed = TRUE), info = f)
  ab <- read_text("about.md")
  expect_false(grepl("33 trial arms", ab, fixed = TRUE))
  expect_true(grepl("35 regimens", ab, fixed = TRUE))
  expect_true(grepl("version 4.2", ab, fixed = TRUE))
  rr <- read_text("References/README.md")
  expect_false(grepl("which the calculator omits", rr, fixed = TRUE))
  expect_true(grepl("0.7 Gy EQD2 per day", rr, fixed = TRUE))
})
