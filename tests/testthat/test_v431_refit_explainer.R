# Version 4.3.1 (2026-11-04): the refit explainer is served by the site and linked from the
# main page and the About section. Each test was run against 4.3.0 (20175de) and failed there.
app_dir <- testthat::test_path("..", "..")
app_abs <- normalizePath(app_dir)
md_path <- file.path(app_abs, "docs", "refit_explainer", "HNSCC_BED_v4_refit_explainer_2026_11_04.md")
html_path <- file.path(app_abs, "www", "refit_explainer.html")
docx_path <- file.path(app_abs, "www", "refit_explainer.docx")
read_utf8 <- function(f) paste(readLines(f, encoding = "UTF-8", warn = FALSE), collapse = "\n")
# Whitespace collapsed, since pandoc wraps long lines of the page.
read_flat <- function(f) gsub("\\s+", " ", read_utf8(f))

test_that("the explainer is served from www/ as a self-contained page and a Word file", {
  expect_true(file.exists(html_path))
  expect_true(file.exists(docx_path))
  html <- read_flat(html_path)
  expect_match(html, "how the tumor model was refitted to randomized trials", fixed = TRUE)
  expect_match(html, "Back to the calculator", fixed = TRUE)
  expect_gte(lengths(regmatches(html, gregexpr("<math", html, fixed = TRUE))), 4)   # equations as MathML
  expect_false(grepl("$$", html, fixed = TRUE))                                    # no unconverted LaTeX
  # Nothing is loaded from elsewhere: no external scripts or stylesheets.
  expect_false(grepl("<script[^>]+src=", html, perl = TRUE))
  expect_false(grepl("<link[^>]+href=\"http", html, perl = TRUE))
  # The Word file is a real document with its equations.
  files <- utils::unzip(docx_path, list = TRUE)$Name
  expect_true("word/document.xml" %in% files)
  xml_dir <- tempfile(); on.exit(unlink(xml_dir, recursive = TRUE))
  utils::unzip(docx_path, files = "word/document.xml", exdir = xml_dir)
  xml <- read_utf8(file.path(xml_dir, "word", "document.xml"))
  expect_gte(lengths(regmatches(xml, gregexpr("<m:oMath", xml, fixed = TRUE))), 4)
})

test_that("the served page is built from the committed source", {
  md <- read_flat(md_path)
  html <- read_flat(html_path)
  # A sentence that changed in the 2026-11-04 audit appears in both.
  key <- "The calibrated DD fits for the sensitivity sets used 200 DI draws each"
  expect_match(md, key, fixed = TRUE)
  expect_match(html, key, fixed = TRUE)
})

test_that("the main page links to the explainer", {
  env <- new.env()
  old <- setwd(app_abs); on.exit(setwd(old))
  source("app.R", local = env)
  html <- htmltools::renderTags(env$ui)$html
  expect_match(html, 'href="refit_explainer.html"', fixed = TRUE)
  expect_match(html, 'href="refit_explainer.docx"', fixed = TRUE)
  # Every linked file exists where Shiny serves it.
  for (f in c(env$REFIT_EXPLAINER_HTML, env$REFIT_EXPLAINER_DOCX))
    expect_true(file.exists(file.path(app_abs, "www", f)), info = f)
})

test_that("the About section links to the explainer", {
  shiny::testServer(app_dir, {
    about <- as.character(output$about$html)
    expect_match(about, 'href="refit_explainer.html"', fixed = TRUE)
    expect_match(about, 'href="refit_explainer.docx"', fixed = TRUE)
  })
})

test_that("the explainer's sensitivity table matches the fit summaries", {
  md <- readLines(md_path, encoding = "UTF-8", warn = FALSE)
  sets <- c("Main" = "main/dd_calibrated_onset7_n1000",
            "Effective N halved" = "effn_half/dd_calibrated_onset7_n200",
            "Planned instead of delivered overall time" = "nominal_clock/dd_calibrated_onset7_n200",
            "Kaplan–Meier units only" = "km_only/dd_calibrated_onset7_n200",
            "RTOG 9003 rates from the 2014 update" = "rtog9003_beitler/dd_calibrated_onset7_n200",
            "Plus ARTSCAN" = "artscan/dd_calibrated_onset7_n200",
            "Plus three later trials" = "new_trials/dd_calibrated_onset7_n200",
            "Plus same-chemotherapy trials" = "chemo/dd_calibrated_onset7_n200",
            "Intercept prior sd 10 instead of 3" = "prior_c10/dd_calibrated_onset7_n200")
  # The table writes every value with its sign, using the Unicode minus.
  fmt <- function(x) sub("^-", "−", sprintf("%+.1f", x))
  for (lab in names(sets)) {
    row <- grep(sprintf("^\\| %s \\|", lab), md, value = TRUE)
    expect_length(row, 1)
    cells <- trimws(strsplit(row, "|", fixed = TRUE)[[1]])
    s <- jsonlite::fromJSON(file.path(app_abs, "fit", paste0("output_refit_2026_10_25c_", sets[[lab]]), "summary.json"),
                            simplifyVector = FALSE)
    hy <- Filter(function(r) startsWith(r$regimen, "HYPNO"), s$regimens)[[1]]$dd
    want <- sprintf("%s (%s to %s)", fmt(hy$median), fmt(hy$q10), fmt(hy$q90))
    # Columns: "", label, arms, DI, calibrated DD, alpha/beta.
    expect_identical(cells[5], want, info = lab)
  }
})
