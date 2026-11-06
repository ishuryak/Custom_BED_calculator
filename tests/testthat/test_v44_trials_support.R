# Version 4.4.0 (2026-11-06): the Refit support column is always shown on the main page, and
# Details opens a list of the trials the Bayesian refit of the DI and DD tumor models used.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "v4_home.R"))
  source(file.path(app_abs, "R", f), local = TRUE)
v4_defaults <- list(
  m = 18, d = 2.9, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE, T = 24,
  comparator = "custom", ref_m = 35, ref_d = 2, ref_fx_week = "5", ref_auto_T = TRUE, ref_T = 47, min_gain = 1,
  preset = "d18", alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.4, Tk = 28, g = 0, dd_mode = "calibrated",
  anchor_m = 35, anchor_d = 2, ab_e = 10, alpha_e = 0.35, Tp_e = 2.5, Tk_e = 7, chemo = FALSE, ab_l = 3,
  criterion = "both", rank_by = "worst", m_range = c(5, 34), min_days = 7, d_step = "0.05", d_max = 5,
  early_mode = "hypno", late_mode = "reference", tol_early = 0, tol_late = 0,
  agree_cmp = "sidebar", agree_threshold_mode = "global", agree_margin = -3, agree_days = 1,
  `oar-library` = "QUANTEC", `oar-ab_cord` = 2, `oar-ab_oar` = 3, `oar-mode` = "percent", `oar-use_in_search` = FALSE)

test_that("the Refit support column and its footnote are shown without Details", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    session$setInputs(home_details = FALSE)
    tb <- output$home_table
    expect_match(tb, "Refit support", fixed = TRUE)
    # The model column and Notes stay under Details.
    for (h in c("Tumor condition holds in", "\"Notes\"")) expect_false(grepl(h, tb, fixed = TRUE), info = h)
    expect_match(as.character(output$home_footnotes$html), "Refit support.", fixed = TRUE)
    session$setInputs(home_details = TRUE)
    tb <- output$home_table
    for (h in c("Refit support", "Tumor condition holds in", "Notes")) expect_match(tb, h, fixed = TRUE)
    # Opening the trial list does not fail.
    expect_silent(session$setInputs(home_trials = 1))
  })
})

test_that("Details offers the fitted-trial list", {
  html <- as.character(htmltools::renderTags(eval(parse(text = "v4_home"), envir = local({
    e <- new.env(); sys.source(file.path(app_abs, "app.R"), envir = e, chdir = TRUE); e }))())$html)
  expect_match(html, "id=\"home_trials\"", fixed = TRUE)
  expect_match(html, "data-display-if=\"input.home_details\"", fixed = TRUE)
  expect_match(html, "Trials the tumor models were refitted to", fixed = TRUE)
})

test_that("the fitted-trial list follows the fitting data named in the draws manifest", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    x <- utils::read.csv(file.path("data", draws_manifest()$dataset))
    expect_equal(nrow(t), nrow(x))
    expect_equal(nrow(t), 45)
    expect_equal(length(unique(t$pmid)), 19)
    expect_setequal(t$arm, x$arm)
    # Values come from the fitting file, formatted for reading.
    ch <- t[t$trial == "CHART (Saunders 2010)" & t$arm == "CHART", ]
    expect_identical(c(ch$dose, ch$fractions, ch$time, ch$total), c("1.5", "36", "12", "54.0"))
    gh <- t[t$trial == "Ghoshal 2008" & grepl("boost", t$arm), ]
    expect_identical(gh$dose, "25 × 1.8 + 15 × 1.5")
    kr <- t[t$trial == "KROG 0201 (Moon 2014)" & t$arm == "CONV", ]
    expect_identical(kr$fractions, "≈ 33.2")
    # Drugs: no arm had chemotherapy; DAHANCA 6&7 gave nimorazole in both arms.
    expect_true(all(startsWith(t$chemotherapy, "none")))
    expect_true(all(grepl("nimorazole", t$other_drugs[t$trial == "DAHANCA 6&7 (Overgaard 2003)"])))
    expect_equal(sum(nzchar(t$other_drugs)), 2)
    expect_match(t$citation[t$pmid == "14511925"][1], "^Overgaard et al\\. Five compared with six fractions")
  })
})

test_that("a trial list that does not match the fitting data is refused", {
  withr::with_dir(app_abs, {
    ref <- utils::read.csv(V4_FIT_TRIAL_LIST, stringsAsFactors = FALSE, colClasses = c(pmid = "character"))
    tmp <- tempfile(fileext = ".csv")
    utils::write.csv(ref[-1, ], tmp, row.names = FALSE)
    expect_error(v4_fit_trials(list_path = tmp), "does not match the fitting data")
    utils::write.csv(rbind(ref, ref[1, ]), tmp, row.names = FALSE)
    expect_error(v4_fit_trials(list_path = tmp), "Repeated arm")
    r2 <- ref; r2$arm[1] <- "renamed arm"
    utils::write.csv(r2, tmp, row.names = FALSE)
    expect_error(v4_fit_trials(list_path = tmp), "does not match the fitting data")
  })
})

test_that("the trial table gives each trial's reference once, with footnotes for drugs and new trials", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    h <- as.character(v4_fit_trials_table(t))
    for (k in c("Trial and reference", "Overall time (days)", "Time basis")) expect_match(h, k, fixed = TRUE)
    # No drug columns: the few trials with a drug note get a footnote instead.
    for (k in c("<th>Chemotherapy</th>", "<th>Other drugs</th>")) expect_false(grepl(k, h, fixed = TRUE), info = k)
    expect_equal(lengths(regmatches(h, gregexpr("<th>", h, fixed = TRUE))), 7)
    expect_equal(lengths(regmatches(h, gregexpr("v4-trial-ref", h, fixed = TRUE))), 19)
    expect_match(h, "rowspan=\"4\"", fixed = TRUE)    # RTOG 9003 and Yamazaki have 4 arms
    expect_match(h, "https://pubmed.ncbi.nlm.nih.gov/14511925/", fixed = TRUE)
    expect_match(h, "https://doi.org/10.1016/s0140-6736(03)14361-9", fixed = TRUE)
    expect_equal(lengths(regmatches(h, gregexpr("<tr>", h, fixed = TRUE))), 1 + 45)
    fn <- v4_fit_trial_notes(t)
    # Drug notes: DAHANCA 6&7 (nimorazole) and RTOG 83-13 (chemotherapy only long before randomization).
    expect_identical(fn$marks[["DAHANCA 6&7 (Overgaard 2003)"]], "a")
    expect_identical(fn$marks[["RTOG 83-13 (Cox 1990)"]], "b,*")
    expect_match(h, "DAHANCA 6&amp;7 (Overgaard 2003): nimorazole", fixed = TRUE)
    expect_match(h, "RTOG 83-13 (Cox 1990): Chemotherapy: no chemotherapy within 6 weeks of randomization", fixed = TRUE)
    # The eight trials added in the refit, from the 2018 paper's references, are starred; no others.
    new <- names(fn$marks)[grepl("\\*", fn$marks)]
    expect_setequal(new, c("EORTC 22811 (Van den Bogaert 1995)", "Ghoshal 2008", "IAEA-ACC (Overgaard 2010)",
                           "KROG 0201 (Moon 2014)", "Princess Margaret (O'Sullivan 2020)", "RTOG 83-13 (Cox 1990)",
                           "RTOG 88-09 (Fu 1995)", "RTOG 9512 (Trotti 2014)"))
    expect_match(h, "Added in the refit: not among the trials of the earlier fits, Shuryak I, Hall EJ, Brenner DJ. Dose dependence", fixed = TRUE)
    expect_match(h, "Int J Radiat Oncol Biol Phys 2019;104:272", fixed = TRUE)
    expect_equal(length(fn$notes), 3)
  })
})

# The three main-page limits can be changed (version 4.4.0).
test_that("margins: checked, and out-of-range values are refused with a message", {
  expect_identical(v4_check_margins(3, 4, 6.5), list(tumor = 3, early = 4, late = 6.5))
  expect_error(v4_check_margins(-1, 5, 5), "tumor BED loss must be a number from 0 to 20")
  expect_error(v4_check_margins(5, NA, 5), "early BED increase")
  expect_error(v4_check_margins(5, 5, 21), "late BED increase")
})

test_that("the main page follows the limits the user sets, and passes them on", {
  shiny::testServer(app_dir, {
    do.call(session$setInputs, v4_defaults)
    n5 <- nrow(home()$rows)
    session$setInputs(home_m_tumor = 3, home_m_early = 4, home_m_late = 6)
    mg <- list(tumor = 3, early = 4, late = 6)
    expect_identical(home()$margins, mg)
    # The same list as the engine gives directly for these limits, and different from 5/5/5.
    expect_identical(home()$rows, v4_regimens(params(), anchor(), margins = mg)$rows)
    expect_false(isTRUE(all.equal(nrow(home()$rows), n5)))
    intro <- as.character(output$home_intro$html)
    expect_match(intro, "at most 3 Gy₁₀ below the standard", fixed = TRUE)
    expect_match(intro, "BED at most 4 Gy₁₀ above the standard", fixed = TRUE)
    expect_match(intro, "BED at most 6 Gy₃ above the standard", fixed = TRUE)
    # Support uses the tumor limit chosen.
    dr <- v4_draws()
    expect_false(is.null(dr)); expect_gt(nrow(home_points()), 0)
    expect_equal(home_points()$support, v4_refit_support(home()$rows, dr$di, dr$dd, margin = 3))
    # One limit tighter and another looser at once.
    session$setInputs(home_m_tumor = 2, home_m_early = 8, home_m_late = 5)
    expect_identical(home()$rows, v4_regimens(params(), anchor(), margins = list(tumor = 2, early = 8, late = 5))$rows)
    # An empty or invalid box (for example while a number is retyped) keeps the last valid
    # limits, so the list and plot stay, and a message says so beside the boxes.
    good <- home()$rows
    for (bad in list(NA, -1, 25)) {
      session$setInputs(home_m_tumor = bad)
      expect_identical(home()$margins, list(tumor = 2, early = 8, late = 5))
      expect_identical(home()$rows, good)
      expect_match(as.character(output$home_plot_ui$html), "home_plot", fixed = TRUE)
      expect_match(as.character(output$home_m_msg$html), "tumor BED loss must be a number from 0 to 20", fixed = TRUE)
      expect_match(as.character(output$home_m_msg$html), "Still using 2, 8 and 5.", fixed = TRUE)
    }
    session$setInputs(home_m_tumor = 4)
    expect_identical(home()$margins, list(tumor = 4, early = 8, late = 5))
    expect_null(output$home_m_msg$html)
  })
})

test_that("the Report card receives the limits in use", {
  rc <- v4_report_inputs(20, 2.75, list(tumor = 3, early = 4, late = 6))
  expect_equal(c(rc$min_gain, rc$tol_early, rc$tol_late), c(-3, 4, 6))
  expect_match(v4_transfer_scope(list(tumor = 3, early = 4, late = 6)),
               "threshold \u22123 and the early and late limits (standard + 4 and standard + 6)", fixed = TRUE)
  expect_equal(v4_report_inputs(20, 2.75)$min_gain, -5)
})

test_that("a main-list export with changed limits replays exactly", {
  for (f in c("oar.R", "export.R", "v4_export.R")) source(file.path(app_abs, "R", f), local = TRUE)
  withr::with_dir(app_abs, {
    p <- default_params(); anchor <- default_anchor(); dom <- tumor_data_support()
    draws <- list(di = load_draws("di"), dd = load_draws("dd"))
    mg <- list(tumor = 3, early = 4, late = 6)
    sc <- main_scenario(p, anchor, draws, dom, margins = mg)
    expect_identical(sc$margins, mg); expect_equal(sc$support$margin, 3)
    r <- main_compute(p, anchor, draws, dom$d_max, margins = mg)
    expect_identical(r$points[, c("m", "d")], v4_regimens(p, anchor, margins = mg)$rows[, c("m", "d")])
    zf <- tempfile(fileext = ".zip"); write_main_zip(zf, sc, r)
    ex <- tempfile(); utils::unzip(zf, exdir = ex)
    back <- jsonlite::fromJSON(file.path(ex, "scenario.json"), simplifyVector = TRUE)
    rr <- main_rerun(back)
    expect_equal(rr$points, r$points, tolerance = 1e-12)
    expect_equal(rr$summary, r$summary)
    # The id differs from the 5/5/5 list, and a tampered support margin is refused.
    expect_false(identical(main_id(sc), main_id(main_scenario(p, anchor, draws, dom))))
    bad <- back; bad$support$margin <- 5
    expect_error(main_rerun(bad), "support margin differs from its tumor margin")
    bad <- back; bad$margins$late <- 30
    expect_error(main_rerun(bad), "margins are not valid")
  })
})
