# Headless-browser walkthrough of version 3 (PLAN_v3 phase 5 checkpoint). Run from the app folder after
# tools/restore_library.R:
#   Rscript docs/development/2026_10_21_v3/walkthrough_v3.R              # local app
#   Rscript docs/development/2026_10_21_v3/walkthrough_v3.R --live URL   # deployed app
#   Rscript docs/development/2026_10_21_v3/walkthrough_v3.R --self-test-fail
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
  out <- file.path(root, "docs", "development", "2026_10_21_v3", "walkthrough",
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
    port <- 8796
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

  num_input <- function(id, v) { js(sprintf("var e=$('#%s'); e.val(%s); e.trigger('change'); true", id, v)); idle() }
  table_text <- function(id) for (i in 1:40) { t <- text_of(id); if (nzchar(t)) return(t); Sys.sleep(0.5) }

  # 1. Fresh page --------------------------------------------------------------------
  fresh()
  check(sprintf("1 title shows version %s", short_ver), grepl(paste("version", short_ver), js("document.title + document.body.innerText"), fixed = TRUE))
  bn <- text_of("banner")
  check("1 banner: HYPNO 55/20 evaluated against 35 x 2 Gy, +1 Gy, DD calibrated",
        grepl("20 × 2.75 Gy in 26 days", bn, fixed = TRUE) && grepl("35 × 2 Gy in 47 days", bn, fixed = TRUE) &&
          grepl("+1 Gy", bn, fixed = TRUE) && grepl("DD calibrated", bn, fixed = TRUE), bn)
  # 3.1 (2026-10-26): column name, DD status wording and fitted-DD banner text updated to the 3.1 app.
  rc <- table_text("report")
  check("1 report card opens first, with a named difference and the refit column",
        grepl("Δ vs 35 × 2 Gy in 47 days", rc, fixed = TRUE) && grepl("Verified-trial refit 2026-10-25", rc, fixed = TRUE), substr(rc, 1, 200))
  check("1 report card flags the extrapolation (2.75 Gy > 2.5 Gy)", grepl("outside the tumor data", text_of("report_notes")))
  check("1 how-to-read card visible", grepl("How to read this", js("document.body.innerText"), fixed = TRUE))
  check("1 Report card: no output errors", !nzchar(errors()), errors())
  shot("01_report_card")
  if (self_fail) check("self-test: this check must fail", FALSE, "deliberate")

  # 2. Every tab opens without errors --------------------------------------------------
  for (t in c("Results", "Compare regimens", "Search", "Both-model screen", "Sensitivity", "Uncertainty", "Organs at risk", "Model and sources")) {
    tab(t); idle(2); check(sprintf("2 %s opens without output errors", t), !nzchar(errors()), errors())
  }
  tab("Search")
  check("2 early-limit value shown beside its control", grepl("Value now: 55.87", text_of("early_limit_value"), fixed = TRUE), text_of("early_limit_value"))
  check("2 search opens empty at the defaults", identical(rows(), ""), rows())
  tab("Both-model screen")
  check("2 both-model screen follows the global comparison (no override banner)", !nzchar(text_of("agree_override")))
  shot("02_both_model_screen")

  # 3. Global comparison: HYPNO and -3 Gy restore 17-19 fractions -----------------------
  set_val("comparator", "hypno_test"); num_input("min_gain", -3)
  tab("Both-model screen"); idle(3)
  at <- table_text("agree_table")
  check("3 HYPNO, -3 Gy: 17, 18, 19 fractions, 15 of 15 settings", grepl("15 of 15 settings", at) && grepl("3.05", at) && grepl("2.90", at) && grepl("2.80", at), substr(at, 1, 160))
  set_radio("agree_threshold_mode", "own"); num_input("agree_margin", 0)
  check("3 an override is highlighted", grepl("Override on this tab", text_of("agree_override"), fixed = TRUE), text_of("agree_override"))
  set_radio("agree_threshold_mode", "global")
  shot("03_override")

  # 4. Uncertainty tab --------------------------------------------------------------------
  tab("Uncertainty"); idle(3)
  us <- table_text("unc-summary")
  check("4 uncertainty summary against HYPNO", grepl("vs 20 × 2.75 Gy in 26 days", us, fixed = TRUE) && grepl("80% range", us, fixed = TRUE), substr(us, 1, 160))
  check("4 refit status names the DD limitation", grepl("DD is not an independent check of DI", text_of("unc-refit_status"), fixed = TRUE))
  js("document.getElementById('unc-run').click(); true"); idle(6)
  st <- table_text("unc-search_table")
  check("4 search under uncertainty lists regimens and neighbors", grepl("listed", st) && grepl("neighbor", st), substr(st, 1, 120))
  shot("04_uncertainty")
  z <- fetch("unc-download")
  check("4 uncertainty ZIP has scenario, summary, search and README",
        all(c("scenario.json", "evaluated_regimen_summary.csv", "search_with_neighbors.csv", "README.txt") %in% z$files), paste(z$files, collapse = ","))
  rr <- callr::r(function(root, js) {
    setwd(root); for (f in c("bed_models.R", "schedule.R", "uncertainty.R", "oar.R", "export.R")) source(file.path("R", f))
    rerun_uncertainty(jsonlite::fromJSON(js))
  }, args = list(root = root, js = file.path(z$dir, "scenario.json")), libpath = c(rel_lib, .Library))
  saved <- read.csv(file.path(z$dir, "evaluated_regimen_summary.csv"))
  check("4 uncertainty download replays exactly", isTRUE(all.equal(rr$di_median, saved$di_median, tolerance = 1e-10)) && isTRUE(all.equal(as.numeric(rr$share_tox), as.numeric(saved$share_tox))))

  # 5. Organs at risk ---------------------------------------------------------------------
  tab("Organs at risk")
  num_input("oar-v_cord", 75); num_input("oar-v_parotid", 40); idle(2)
  ot <- table_text("oar-table")
  check("5 organ table shows converted limits and statuses", grepl("Spinal cord", ot) && grepl("within", ot) && grepl("not assessed", ot), substr(ot, 1, 160))
  check("5 the binding constraint is named", grepl("Binds first", text_of("oar-binding"), fixed = TRUE), text_of("oar-binding"))
  tab("Report card")
  check("5 organ rows reach the report card", grepl("Spinal cord, maximum EQD2", table_text("report"), fixed = TRUE))
  tab("Organs at risk")
  set_radio("oar-mode", "dvh"); idle()
  js(sprintf("var s=document.querySelector('#oar-dvh_file'); true"))
  node <- b$DOM$getDocument()$root$nodeId
  q <- b$DOM$querySelector(node, "#oar-dvh_file")$nodeId
  fx <- normalizePath(file.path(root, "tests", "testthat", "fixtures", "dvh"))
  invisible(b$DOM$setFileInputFiles(files = list(file.path(fx, "eclipse_with_identifiers.txt")), nodeId = q)); idle(3)
  check("5 a DVH file with patient identifiers is refused", grepl("patient identifier fields", text_of("oar-dvh_msg"), fixed = TRUE), text_of("oar-dvh_msg"))
  invisible(b$DOM$setFileInputFiles(files = list(file.path(fx, "eclipse_synthetic.txt")), nodeId = q)); idle(3)
  check("5 a synthetic DVH without identifiers is read", grepl("Read 2 structures", text_of("oar-dvh_msg"), fixed = TRUE), text_of("oar-dvh_msg"))
  shot("05_organs")

  # 6. Calendar and fitted DD --------------------------------------------------------------
  set_val("comparator", "custom"); num_input("min_gain", 1)
  num_input("m", 33); num_input("d", 2); set_radio("fx_week", 6); idle(2)
  check("6 33 fractions at 6 per week take 38 days", grepl("33 fractions, 38 days", text_of("T_shown"), fixed = TRUE), text_of("T_shown"))
  set_radio("fx_week", 5); num_input("m", 20); num_input("d", 2.75)
  set_radio("dd_mode", "fitted"); idle(2)
  check("6 fitted DD named in the banner", grepl("DD values entered directly (not calibrated)", text_of("banner"), fixed = TRUE), text_of("banner"))
  set_radio("dd_mode", "calibrated")
  check("6 no output errors after the calendar and DD changes", !nzchar(errors()), errors())

  # 7. Narrow screen ------------------------------------------------------------------------
  b$Emulation$setDeviceMetricsOverride(width = 390, height = 844, deviceScaleFactor = 1, mobile = TRUE); idle(2)
  tab("Report card")
  check("7 narrow: no horizontal page scroll", overflow_px() <= 1, overflow_px())
  shot("07_narrow_report")
  b$Emulation$clearDeviceMetricsOverride()
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
