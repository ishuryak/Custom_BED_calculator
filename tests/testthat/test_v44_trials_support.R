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

test_that("the trial table gives each trial's reference once, spanning its arms, with links", {
  withr::with_dir(app_abs, {
    t <- v4_fit_trials()
    h <- as.character(v4_fit_trials_table(t))
    for (k in c("Trial and reference", "Overall time (days)", "Time basis", "Chemotherapy", "Other drugs"))
      expect_match(h, k, fixed = TRUE)
    expect_equal(lengths(regmatches(h, gregexpr("v4-trial-ref", h, fixed = TRUE))), 19)
    expect_match(h, "rowspan=\"4\"", fixed = TRUE)    # RTOG 9003 and Yamazaki have 4 arms
    expect_match(h, "https://pubmed.ncbi.nlm.nih.gov/14511925/", fixed = TRUE)
    expect_match(h, "https://doi.org/10.1016/s0140-6736(03)14361-9", fixed = TRUE)
    expect_equal(lengths(regmatches(h, gregexpr("<tr>", h, fixed = TRUE))), 1 + 45)
  })
})
