# Headless-browser walkthrough of version 2.7 (a colleague's comments, and the
# independent review of the plan, section 5F). Run from the app folder after
# tools/restore_library.R:
#   Rscript docs/development/2026_10_20_v27/walkthrough_v27.R              # local app
#   Rscript docs/development/2026_10_20_v27/walkthrough_v27.R --live URL   # deployed app
#   Rscript docs/development/2026_10_20_v27/walkthrough_v27.R --self-test-fail
# Local mode starts the app in the release library; live mode only opens the
# URL. The opening page is checked before any control is touched. Each run
# writes a log and screenshots to its own folder under walkthrough/. Any failed
# check, or a timeout, gives exit status 1. --self-test-fail adds one check
# that must fail, to confirm that a failure reaches the exit status.
library(chromote)
`%||%` <- function(a, b) if (is.null(a)) b else a

main <- function(args) {
  root <- normalizePath(".")
  source(file.path(root, "tools", "env.R"))
  rel_lib <- bed_library()
  stopifnot(dir.exists(file.path(rel_lib, "shiny")))
  live <- "--live" %in% args
  base <- if (live) sub("/$", "", args[which(args == "--live") + 1]) else NULL
  self_fail <- "--self-test-fail" %in% args
  expected <- local({ e <- new.env(); sys.source(file.path(root, "R", "bed_models.R"), e)
                      list(app = e$APP_VERSION, sel = e$SELECTION_VERSION, eq = e$EQUATION_VERSION) })
  short_ver <- sub("\\.0$", "", expected$app)
  out <- file.path(root, "docs", "development", "2026_10_20_v27", "walkthrough",
                   paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), if (live) "_live" else "_local"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  log <- character()
  check <- function(name, cond, detail = "") {
    line <- sprintf("%s  %s%s", if (isTRUE(cond)) "PASS" else "FAIL", name,
                    if (nzchar(detail)) paste0("  [", detail, "]") else "")
    log <<- c(log, line); cat(line, "\n")
  }
  app <- NULL; b <- NULL
  on.exit({
    if (!is.null(b)) try(b$close(), silent = TRUE)
    if (!is.null(app)) try(app$kill(), silent = TRUE)
  }, add = TRUE)

  if (!live) {
    port <- 8795
    app <- callr::r_bg(function(root, port) { setwd(root); shiny::runApp(root, port = port, launch.browser = FALSE) },
                       args = list(root = root, port = port), libpath = c(rel_lib, .Library))
    base <- sprintf("http://127.0.0.1:%d", port)
    ok <- FALSE
    for (i in 1:90) { ok <- tryCatch({ suppressWarnings(readLines(base, warn = FALSE)); TRUE }, error = function(e) FALSE); if (ok) break; Sys.sleep(1) }
    if (!ok) stop("local app did not start within 90 s")
  }

  b <- ChromoteSession$new(width = 1400, height = 1000)
  console <- character()
  invisible(b$Runtime$enable())
  b$Runtime$consoleAPICalled(callback_ = function(x) if (x$type %in% c("error", "assert"))
    console <<- c(console, paste(vapply(x$args, function(a) as.character(a$value %||% a$description %||% ""), ""), collapse = " ")))
  b$Runtime$exceptionThrown(callback_ = function(x) console <<- c(console, paste("exception:", x$exceptionDetails$text)))
  js <- function(expr) b$Runtime$evaluate(expr, returnByValue = TRUE)$result$value
  idle <- function(min_wait = 1.2) {
    Sys.sleep(min_wait)
    for (i in 1:120) { if (isTRUE(js("!$('html').hasClass('shiny-busy')"))) return(invisible(TRUE)); Sys.sleep(0.5) }
    stop("page stayed busy for 60 s")
  }
  errors <- function() js("Array.from(document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')).map(e => e.id + ': ' + e.innerText).join(' | ')")
  tab <- function(label) { js(sprintf("Array.from(document.querySelectorAll('.nav-link')).find(a => a.innerText.trim().startsWith('%s')).click(); true", label)); idle() }
  text_of <- function(id) js(sprintf("(document.getElementById('%s') || {innerText: ''}).innerText", id))
  value_of <- function(id) js(sprintf("$('#%s').val()", id))
  set_val <- function(id, v) {
    vv <- if (is.character(v)) shQuote(v) else v
    js(sprintf("var e=$('#%s'); if (e[0] && e[0].selectize) e[0].selectize.setValue(%s); else e.val(%s).trigger('change'); true", id, vv, vv))
    idle()
    if (isTRUE(js(sprintf("!!($('#%s')[0] || {}).selectize", id))) && !identical(value_of(id), as.character(v)))
      stop("menu ", id, " did not take value ", v)
  }
  set_check <- function(id, v) { js(sprintf("var e=$('#%s'); if (e.prop('checked') !== %s) e.click(); true", id, tolower(v))); idle() }
  set_radio <- function(id, v) { js(sprintf("$('input[name=%s][value=%s]').click(); true", id, v)); idle() }
  shot <- function(name) b$screenshot(file.path(out, paste0(name, ".png")), show = FALSE)
  overflow_px <- function() { w <- js("[document.documentElement.scrollWidth, document.documentElement.clientWidth]"); as.numeric(w[[1]]) - as.numeric(w[[2]]) }
  # Rows of the candidate table once its output has rendered ("" when empty).
  rows <- function() {
    for (i in 1:40) { if (isTRUE(js("!!document.querySelector('#suggest_table table') || document.getElementById('suggest_msg').innerText.length > 0"))) break; Sys.sleep(0.5) }
    v <- rows_now(); if (is.null(v)) "<not readable>" else v
  }
  # Poll until the table shows the expected rows (a live server can lag behind
  # the busy flag); returns what it last read, so a mismatch is still reported.
  rows_until <- function(want, timeout = 20) {
    v <- rows()
    for (i in seq_len(timeout * 2)) { if (identical(v, want)) break; Sys.sleep(0.5); v <- rows_now() %||% "<not readable>" }
    v
  }
  rows_now <- function() js("Array.from(document.querySelectorAll('#suggest_table tbody tr')).filter(r => !r.querySelector('.dataTables_empty')).map(r => r.cells[0].innerText + 'x' + r.cells[1].innerText).join(',')")
  visible <- function(id) isTRUE(js(sprintf("(function(){var e=document.getElementById('%s'); if(!e) return false; var r=e.getBoundingClientRect(); return e.offsetParent !== null && r.width > 0 && r.height > 0;})()", id)))
  fetch <- function(id) {
    href <- js(sprintf("document.getElementById('%s').getAttribute('href')", id))
    f <- tempfile(fileext = ".zip"); download.file(paste0(base, "/", href), f, mode = "wb", quiet = TRUE)
    d <- tempfile(); utils::unzip(f, exdir = d); list(files = list.files(d), dir = d)
  }
  replay <- function(dir) callr::r(function(root, js) {
    source(file.path(root, "R", "bed_models.R")); source(file.path(root, "R", "export.R"))
    rerun_scenario(jsonlite::fromJSON(js))
  }, args = list(root = root, js = file.path(dir, "scenario.json")), libpath = c(rel_lib, .Library))
  fresh <- function() { invisible(b$Page$navigate(base)); invisible(b$Page$loadEventFired()); idle(4) }

  # 1. Fresh page, before any control is touched -------------------------------------
  fresh()
  check(sprintf("1 title shows version %s", short_ver), grepl(paste("version", short_ver), js("document.title + document.body.innerText"), fixed = TRUE))
  check("1 opening criterion is both models", identical(value_of("criterion"), "both"), value_of("criterion"))
  check("1 opening ranking is the smaller change", identical(value_of("rank_by"), "worst"), value_of("rank_by"))
  check("1 opening minimum change is +1 Gy", isTRUE(as.numeric(value_of("min_gain")) == 1), value_of("min_gain"))
  check("1 visible labels agree", grepl("both models", js("$('#criterion').parent().text()")) &&
          grepl("smaller of", js("$('#rank_by').parent().text()")))
  check("1 Results caveat visible", visible("caveat-results"))
  check("1 Results: no output errors", !nzchar(errors()), errors())
  shot("01_fresh_results")
  tab("Candidate regimens")
  check("1 Candidate caveat visible", visible("caveat-candidates"))
  check("1 opening candidate table is empty", identical(rows(), ""), rows())
  msg <- text_of("suggest_msg")
  check("1 empty message: both models, valid outcome, separate agreement comparison",
        grepl("under both models", msg) && grepl("valid outcome", msg) && grepl("explores a separate comparison", msg), substr(msg, 1, 160))
  check("1 candidate context names the sidebar comparison", grepl("Comparison: 35 × 2 Gy in 47 days (sidebar)", text_of("search_context"), fixed = TRUE), text_of("search_context"))
  shot("02_fresh_candidates")
  tab("Model agreement")
  check("1 Agreement caveat visible, above the table", visible("caveat-agreement") &&
          isTRUE(js("!!(document.getElementById('caveat-agreement').compareDocumentPosition(document.getElementById('agree_table')) & Node.DOCUMENT_POSITION_FOLLOWING)")))
  check("1 agreement context names HYPNO, -3 Gy and 1 day", grepl("Comparison: 20 × 2.75 Gy in 26 days (HYPNO test arm)", text_of("agree_context"), fixed = TRUE) &&
          grepl("at least -3 Gy under both models", text_of("agree_context"), fixed = TRUE))
  at <- text_of("agree_table")
  check("1 agreement rows 17, 18, 19 with 15 of 15 settings", grepl("15 of 15 settings", at) && grepl("3.05", at) && grepl("2.90", at) && grepl("2.80", at))
  check("1 agreement headers say settings", grepl("Selected parameter settings with both", at) && grepl("over evaluable settings", at) && !grepl("Scenarios with both", at))
  shot("03_fresh_agreement")
  if (self_fail) check("self-test: this check must fail", FALSE, "deliberate")

  # 2. All tabs and sidebar groups ------------------------------------------------------
  js("document.querySelectorAll('.accordion-button.collapsed').forEach(x => x.click()); true"); idle()
  check("2 all six sidebar panels open", js("document.querySelectorAll('.accordion-collapse.show').length") == 6)
  for (t in c("Results", "Candidate regimens", "Model agreement", "Sensitivity", "Model and sources")) {
    tab(t); check(sprintf("2 %s opens without output errors", t), !nzchar(errors()), errors())
  }
  idle(3)
  check("2 equations rendered by MathJax", js("document.querySelectorAll('mjx-container, .MathJax, .MathJax_Display').length") > 0)
  href <- js("(Array.from(document.querySelectorAll('#about a')).find(a => a.innerText.includes('v2-toxicity-hypofractionation')) || {}).href")
  check("2 source link targets the version 2 branch", identical(href, "https://github.com/ishuryak/Custom_BED_calculator/tree/v2-toxicity-hypofractionation"), href %||% "missing")
  check("2 version footer from the code", grepl(sprintf("Version %s of the Custom BED calculator (equations %s, selection rules %s", expected$app, expected$eq, expected$sel),
                                               text_of("about"), fixed = TRUE))
  shot("04_model")

  # 3. Candidate fixtures, chosen explicitly ---------------------------------------------
  tab("Candidate regimens")
  set_val("min_gain", -3)
  check("3 both/smaller/-3 Gy lists 30, 29, 28 fractions in that order", identical(r3 <- rows_until("30x2.20,29x2.25,28x2.30"), "30x2.20,29x2.25,28x2.30"), r3)
  check("3 plot renders", isTRUE(js("!!document.querySelector('#suggest_plot img')")))
  z <- fetch("download")
  res <- read.csv(file.path(z$dir, "results.csv")); rr <- replay(z$dir)
  sc <- jsonlite::fromJSON(file.path(z$dir, "scenario.json"))
  check("3 download has its four files", setequal(z$files, c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt")))
  check("3 download replays to full precision", identical(rr$m, res$m) && isTRUE(all.equal(rr$d_dd, res$d_dd, tolerance = 1e-8)))
  check("3 export records version, settings and selection rules",
        identical(sc$app_version, expected$app) && identical(sc$selection_version, expected$sel) &&
          identical(sc$search$criterion, "both") && identical(sc$search$rank_by, "worst") && isTRUE(as.numeric(sc$search$min_gain) == -3))
  if (!live) check("3 export code unchanged since load", identical(sc$code$changed_since_load, FALSE))
  shot("05_candidates_m3")
  set_val("min_gain", 1); set_val("criterion", "either"); set_val("rank_by", "best")
  check("3 either/larger/+1 Gy lists 29, 28, 30", identical(r4 <- rows_until("29x2.25,28x2.30,30x2.20"), "29x2.25,28x2.30,30x2.20"), r4)
  set_val("criterion", "both"); set_val("rank_by", "worst")
  check("3 back to the opening settings: empty again, no stale rows", identical(r5 <- rows_until(""), ""), r5)
  ze <- fetch("download")
  check("3 empty-state download still has its files", setequal(ze$files, c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt")))

  # 4. Agreement downloads and its own comparison --------------------------------------
  tab("Model agreement")
  za <- fetch("download_agree")
  det <- read.csv(file.path(za$dir, "scenario_detail.csv")); ra <- read.csv(file.path(za$dir, "results.csv"))
  rra <- replay(za$dir)
  check("4 agreement download replays, 45 setting rows", nrow(det) == 45 && identical(rra$m, ra$m))
  set_val("agree_cmp", "sidebar")
  check("4 sidebar comparison switches the agreement context", grepl("Comparison: 35 × 2 Gy in 47 days (the sidebar comparison)", text_of("agree_context"), fixed = TRUE))
  set_val("agree_cmp", "hypno_test")

  # 5. Validity states ----------------------------------------------------------------
  tab("Candidate regimens")
  set_val("lambda_dd", 0)
  check("5 DD unavailable: DD-specific explanation, no generic text", grepl("DD model is unavailable", text_of("suggest_msg")) && !grepl("valid outcome", text_of("suggest_msg")))
  check("5 DD unavailable: no stale rows", identical(rows(), ""))
  set_val("lambda_dd", 0.4)
  set_val("m", 20.5); tab("Results")
  check("5 invalid fraction count gives a validation message", grepl("whole number", text_of("results")))
  set_val("m", 20)

  # 6. Chemotherapy checkbox ----------------------------------------------------------
  r0 <- text_of("results"); z0 <- text_of("acute_zone")
  set_check("chemo", TRUE)
  check("6 chemotherapy lowers the Fowler zone", grepl("54–58", text_of("acute_zone")) && !identical(z0, text_of("acute_zone")))
  check("6 chemotherapy leaves every BED unchanged", identical(r0, text_of("results")))
  tab("Candidate regimens")
  check("6 chemotherapy leaves the HYPNO limit unchanged", grepl("55.87 Gy₁₀ (HYPNO 55 Gy/20, fixed)", text_of("search_context"), fixed = TRUE))
  set_check("chemo", FALSE)

  # 7. Sensitivity wording ----------------------------------------------------------------
  tab("Sensitivity")
  s <- text_of("sens_summary")
  check("7 summary counts parameter settings and says what they are not",
        grepl("evaluable parameter settings under DI", s) && grepl("not independent confirmations or a probability sample", s) &&
          grepl("Uncertainty intervals are not calculated", s), substr(s, 1, 200))
  check("7 sensitivity: 15 rows", js("document.querySelectorAll('#sens_table tbody tr').length") == 15)
  set_val("min_gain", 0); check("7 summary follows the candidate threshold", grepl("at least 0 Gy", text_of("sens_summary"))); set_val("min_gain", 1)
  shot("06_sensitivity")

  # 8. Narrow screen ------------------------------------------------------------------
  b$Emulation$setDeviceMetricsOverride(width = 420, height = 900, deviceScaleFactor = 1, mobile = TRUE); idle(2)
  for (t in c("Results", "Candidate regimens", "Model agreement", "Sensitivity", "Model and sources")) {
    tab(t); check(sprintf("8 narrow 420 px: %s has no page-level horizontal scroll", t), overflow_px() <= 1, paste(overflow_px(), "px"))
  }
  tab("Candidate regimens")
  check("8 narrow: caveat and criterion menu visible", visible("caveat-candidates") && (visible("criterion-selectized") || visible("criterion")))
  shot("07_narrow_candidates")
  b$Emulation$clearDeviceMetricsOverride()

  # 9. Reload: the opening state returns ------------------------------------------------
  fresh()
  tab("Candidate regimens")
  check("9 after reload: both, smaller change, +1 Gy, empty table",
        identical(value_of("criterion"), "both") && identical(value_of("rank_by"), "worst") && isTRUE(as.numeric(value_of("min_gain")) == 1) && identical(rows(), ""),
        paste(vapply(list(value_of("criterion"), value_of("rank_by"), value_of("min_gain"), rows()),
                     function(x) paste(deparse(x), collapse = ""), ""), collapse = " | "))
  check("browser console: no errors", length(console) == 0, paste(head(console, 5), collapse = " | "))

  rev <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "?")
  dirty <- tryCatch(length(system("git status --porcelain -- R app.R about.md DESCRIPTION", intern = TRUE)) > 0, error = function(e) NA)
  writeLines(c(sprintf("Walkthrough %s; target %s (%s); local revision %s%s; expected version %s",
                       format(Sys.time()), base, if (live) "live" else "local app started by this script", rev,
                       if (isTRUE(dirty)) " with uncommitted app changes" else "", expected$app),
               log), file.path(out, "walkthrough_log.txt"))
  sum(startsWith(log, "FAIL"))
}

status <- tryCatch({
  n_fail <- main(commandArgs(trailingOnly = TRUE))
  cat(n_fail, "failed\n")
  if (n_fail > 0) 1L else 0L
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); 1L })
quit(status = status)
