# Headless-browser walkthrough of version 2.6, following the browser acceptance
# procedure (section 10) of the external review of 2.5.1. Run from the app
# folder after tools/restore_library.R:
#   Rscript docs/development/2026_10_19_v26/walkthrough_v26.R
# The app runs in the release library; the driver (chromote) runs in the
# packages installed in R. Writes a log and screenshots to
# docs/development/2026_10_19_v26/walkthrough/.
library(chromote)
root <- normalizePath(".")
source(file.path(root, "tools", "env.R"))
rel_lib <- bed_library()
stopifnot(dir.exists(file.path(rel_lib, "shiny")))
out <- file.path(root, "docs", "development", "2026_10_19_v26", "walkthrough")
dir.create(out, showWarnings = FALSE)
port <- 8793
app <- callr::r_bg(function(root, port) { setwd(root); shiny::runApp(root, port = port, launch.browser = FALSE) },
                   args = list(root = root, port = port), libpath = c(rel_lib, .Library))
on.exit(app$kill(), add = TRUE)
base <- sprintf("http://127.0.0.1:%d", port)
for (i in 1:90) { ok <- tryCatch({ suppressWarnings(readLines(base, warn = FALSE)); TRUE }, error = function(e) FALSE); if (ok) break; Sys.sleep(1) }
stopifnot(ok)

b <- ChromoteSession$new(width = 1400, height = 1000)
console <- character()
invisible(b$Runtime$enable())
b$Runtime$consoleAPICalled(callback_ = function(x) if (x$type %in% c("error", "assert"))
  console <<- c(console, paste(vapply(x$args, function(a) as.character(a$value %||% a$description %||% ""), ""), collapse = " ")))
b$Runtime$exceptionThrown(callback_ = function(x) console <<- c(console, paste("exception:", x$exceptionDetails$text)))
`%||%` <- function(a, b) if (is.null(a)) b else a
js <- function(expr) b$Runtime$evaluate(expr, returnByValue = TRUE)$result$value
idle <- function(min_wait = 1.2) {
  Sys.sleep(min_wait)
  for (i in 1:80) { if (isTRUE(js("!$('html').hasClass('shiny-busy')"))) break; Sys.sleep(0.5) }
  Sys.sleep(0.4)
}
log <- character()
check <- function(name, cond, detail = "") {
  line <- sprintf("%s  %s%s", if (isTRUE(cond)) "PASS" else "FAIL", name, if (nzchar(detail)) paste0("  [", detail, "]") else "")
  log <<- c(log, line); cat(line, "\n")
}
errors <- function() js("Array.from(document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')).map(e => e.id + ': ' + e.innerText).join(' | ')")
tab <- function(label) { js(sprintf("Array.from(document.querySelectorAll('.nav-link')).find(a => a.innerText.trim().startsWith('%s')).click(); true", label)); idle() }
# Select menus are selectize widgets: set them through selectize, and confirm
# the value took, so a check can never pass on an unchanged menu.
set_val <- function(id, v) {
  vv <- if (is.character(v)) shQuote(v) else v
  js(sprintf("var e=$('#%s'); if (e[0] && e[0].selectize) e[0].selectize.setValue(%s); else e.val(%s).trigger('change'); true", id, vv, vv))
  idle()
  if (isTRUE(js(sprintf("!!($('#%s')[0] || {}).selectize", id))) && !identical(value_of(id), as.character(v)))
    stop("menu ", id, " did not take value ", v)
}
set_check <- function(id, v) { js(sprintf("var e=$('#%s'); if (e.prop('checked') !== %s) e.click(); true", id, tolower(v))); idle() }
set_radio <- function(id, v) { js(sprintf("$('input[name=%s][value=%s]').click(); true", id, v)); idle() }
set_slider <- function(id, from, to) { js(sprintf("$('#%s').data('ionRangeSlider').update({from: %d, to: %d}); true", id, from, to)); idle() }
text_of <- function(id) js(sprintf("(document.getElementById('%s') || {innerText: ''}).innerText", id))
value_of <- function(id) js(sprintf("$('#%s').val()", id))
shot <- function(name) b$screenshot(file.path(out, paste0(name, ".png")), show = FALSE)
overflow_px <- function() { w <- js("[document.documentElement.scrollWidth, document.documentElement.clientWidth]"); as.numeric(w[[1]]) - as.numeric(w[[2]]) }
no_overflow <- function() overflow_px() <= 1
widest <- function() js("(function(){var cw=document.documentElement.clientWidth, r=[]; document.querySelectorAll('body *').forEach(function(e){var b=e.getBoundingClientRect(); if (b.right > cw + 1 && getComputedStyle(e).position !== 'fixed') r.push((e.id||e.className||e.tagName).toString().slice(0,40) + ':' + Math.round(b.right));}); return r.slice(0,6).join(' ; ');})()")
fetch <- function(id) {
  href <- js(sprintf("document.getElementById('%s').getAttribute('href')", id))
  f <- tempfile(fileext = ".zip"); download.file(paste0(base, "/", href), f, mode = "wb", quiet = TRUE)
  d <- tempfile(); utils::unzip(f, exdir = d); list(files = list.files(d), dir = d)
}
replay <- function(dir) callr::r(function(root, js) {
  source(file.path(root, "R", "bed_models.R")); source(file.path(root, "R", "export.R"))
  rerun_scenario(jsonlite::fromJSON(js))
}, args = list(root = root, js = file.path(dir, "scenario.json")), libpath = c(rel_lib, .Library))

invisible(b$Page$navigate(base)); invisible(b$Page$loadEventFired()); idle(4)
check("title shows version 2.6", grepl("version 2.6", js("document.title + document.body.innerText")))

# 1. Accordions, invalid fraction counts ---------------------------------------
js("document.querySelectorAll('.accordion-button.collapsed').forEach(x => x.click()); true"); idle()
check("1 all six sidebar panels open", js("document.querySelectorAll('.accordion-collapse.show').length") == 6)
check("1 Results: no output errors at defaults", !nzchar(errors()), errors())
set_val("m", 20.5)
check("1 fractional fraction count gives a validation message", grepl("whole number", text_of("results")))
set_val("m", "''")
check("1 blank fraction count gives a validation message, no stale table", grepl("Enter a number", text_of("results")) &&
        js("document.querySelectorAll('#results table').length") == 0)
set_val("m", 500)
check("1 out-of-range fraction count gives a validation message", grepl("between 1 and 100", text_of("results")))
set_val("m", 20)

# 2. Calendar and comparators ------------------------------------------------------
for (mm in c(35, 20, 15)) { set_val("m", mm)
  want <- c(`35` = 47, `20` = 26, `15` = 19)[[as.character(mm)]]
  check(sprintf("2 %d fractions give %d days", mm, want), grepl(sprintf("%d fractions, %d days", mm, want), text_of("T_shown"), fixed = TRUE)) }
set_val("m", 20)
set_check("auto_T", FALSE); set_val("T", 20)
check("2 too-short manual time is rejected", grepl("need at least 26 days", text_of("results")))
set_val("T", 26); set_check("auto_T", TRUE)
set_val("comparator", "hypno_test")
check("2 HYPNO test comparison is 20 fractions, 26 days", grepl("20 fractions, 26 days, total dose 55.00", text_of("ref_T_shown"), fixed = TRUE))
set_val("comparator", "hypno_control")
check("2 HYPNO control comparison is 33 fractions, 38 days", grepl("33 fractions, 38 days, total dose 66.00", text_of("ref_T_shown"), fixed = TRUE))
check("2 fraction slider upper bound follows the comparison (32)", js("$('#m_range').data('ionRangeSlider').options.max") == 32)
check("2 custom comparison fields hidden for a preset comparison", isTRUE(js("$('#ref_m').is(':hidden')")))
dd_ctl <- text_of("results")
set_val("comparator", "custom")
check("2 fixed anchor: changing the comparison leaves this regimen's DD BED unchanged",
      identical(sub(".*Tumor BED, DD model[^\n]*?\t([0-9.]+).*", "\\1", dd_ctl), sub(".*Tumor BED, DD model[^\n]*?\t([0-9.]+).*", "\\1", text_of("results"))))
shot("01_results")

# 3. Presets ---------------------------------------------------------------------
for (pr in c("calculator", "withers", "hartley", "hyhope", "d18")) {
  set_val("preset", pr)
  check(sprintf("3 preset %s keeps itself selected and has a note", pr), identical(value_of("preset"), pr) && nzchar(text_of("preset_note")))
}
set_val("preset", "withers"); set_val("ab", 5)
check("3 Withers at α/β 5 gives 0.84 Gy/day", grepl("costs 0.84 Gy", text_of("tumor_note"), fixed = TRUE), text_of("tumor_note"))
set_val("ab", 10)
check("3 Withers at α/β 10 gives 0.72 Gy/day", grepl("costs 0.72 Gy", text_of("tumor_note"), fixed = TRUE), text_of("tumor_note"))
set_val("preset", "hartley"); set_val("lambda_dd", 0.6)
check("3 editing an inherited DD value keeps the Hartley preset", identical(value_of("preset"), "hartley"))
set_val("g", 0.01)
check("3 editing g keeps the Hartley preset", identical(value_of("preset"), "hartley"))
set_val("alpha", 0.31)
check("3 editing a preset-controlled value switches to Custom", identical(value_of("preset"), "custom"))
set_val("g", 0); set_val("preset", "d18")
check("3 D18 restored", identical(value_of("alpha"), "0.35") && identical(value_of("lambda_dd"), "0.4"), paste(value_of("alpha"), value_of("lambda_dd")))

# 4. Edge states, and R1 through the controls -------------------------------------------
set_val("lambda_di", 0)
check("4 zero DI: no output errors", !nzchar(errors()), errors())
set_val("lambda_di", 0.23); set_val("lambda_dd", 0)
check("4 zero DD with positive DI: DD unavailable warning", grepl("DD model unavailable", text_of("dd_warning")))
tab("Model agreement")
check("4 zero DD: agreement explains", grepl("DD model is unavailable", text_of("agree_msg")))
tab("Sensitivity")
check("4 zero DD: sensitivity labels unavailable values", grepl("not available (no DD calibration)", text_of("sens_table"), fixed = TRUE))
set_val("lambda_dd", 0.4)
tab("Candidate regimens")
set_val("preset", "custom"); set_val("alpha", 0.2); set_val("lambda_di", 0.4); set_val("lambda_dd", 5)
set_radio("early_mode", "reference"); set_val("tol_early", 30); set_val("tol_late", 30)
set_val("criterion", "either"); set_val("rank_by", "worst"); set_val("min_gain", 0); set_val("min_days", 1)
set_val("d_max", 2.5); set_slider("m_range", 32, 32)
ct <- js("document.getElementById('suggest_table').innerText")
check("4 R1 through the controls: conservative ranking lists 32 × 2.10 Gy", grepl("32\t2.10", ct, fixed = TRUE), substr(gsub("\n", " / ", ct), 1, 200))
check("4 R1: 32 × 2.50 Gy is not listed", !grepl("32\t2.50", ct, fixed = TRUE))
set_val("rank_by", "dd"); set_val("d_step", "0.5")
check("4 R1: DD ranking with no evaluable DD explains, lists nothing", grepl("none can be ranked", text_of("suggest_msg")))
shot("02_r1_no_rankable")
set_val("d_step", "0.05"); set_val("rank_by", "best"); set_val("min_gain", 1); set_val("min_days", 7); set_val("d_max", 5)
set_slider("m_range", 5, 34); set_val("tol_early", 0); set_val("tol_late", 0); set_radio("early_mode", "hypno")
set_val("preset", "d18")
idle(2)
check("4 defaults restored", identical(value_of("preset"), "d18") &&
        isTRUE(all.equal(as.numeric(c(value_of("alpha"), value_of("lambda_di"), value_of("lambda_dd"))), c(0.35, 0.23, 0.4))), paste(value_of("preset"), value_of("alpha"), value_of("lambda_di"), value_of("lambda_dd")))

# 5. Results with nonstandard alpha/beta ---------------------------------------------
tab("Results")
set_val("ab_e", 7); set_val("ab_l", 2.5)
r <- text_of("results")
check("5 units follow nonstandard α/β", grepl("Early toxicity BED (Gy (α/β = 7 Gy))", r, fixed = TRUE) &&
        grepl("Late toxicity BED (Gy (α/β = 2.5 Gy))", r, fixed = TRUE))
check("5 Fowler zone flagged as not applicable", grepl("have been changed", text_of("acute_zone")))
check("5 benchmark heading says recomputed", grepl("recomputed with the current parameters", js("document.body.innerText")))
set_val("ab_e", 10); set_val("ab_l", 3)

# 6. Candidate tab, every choice ------------------------------------------------------
tab("Candidate regimens")
n_rows <- js("document.querySelectorAll('#suggest_table tbody tr').length")
check("6 default candidate rows: 3", identical(as.integer(n_rows), 3L), as.character(n_rows))
ct <- js("document.getElementById('suggest_table').innerText")
check("6 default rows 29 × 2.25, 28 × 2.30, 30 × 2.20", grepl("29\t2.25", ct) && grepl("28\t2.30", ct) && grepl("30\t2.20", ct))
check("6 slack and fractions-saved columns", grepl("Early limit minus regimen", ct) && grepl("Fractions saved", ct))
check("6 context strip names comparison and limits", grepl("Comparison: 35 × 2 Gy in 47 days (sidebar)", text_of("search_context"), fixed = TRUE) &&
        grepl("55.87 Gy₁₀ (HYPNO 55 Gy/20, fixed)", text_of("search_context"), fixed = TRUE), text_of("search_context"))
bad <- character()
for (cr in c("either", "both", "di", "dd")) for (rk in c("best", "worst", "di", "dd")) {
  set_val("criterion", cr); set_val("rank_by", rk)
  if (nzchar(errors())) bad <- c(bad, paste(cr, rk, errors()))
}
check("6 all 16 criterion/ranking pairs render without errors", length(bad) == 0, paste(bad, collapse = "; "))
set_val("criterion", "either"); set_val("rank_by", "best")
bad <- character()
for (st in c("0.1", "0.2", "0.25", "0.5", "0.05")) { set_val("d_step", st); if (nzchar(errors())) bad <- c(bad, st) }
check("6 every dose step renders", length(bad) == 0, paste(bad, collapse = ","))
for (em in c("reference", "fowler", "fenwick", "hypno")) {
  set_radio("early_mode", em)
  check(sprintf("6 early limit %s renders", em), !nzchar(errors()), errors())
  if (em == "reference") check("6 reference mode shows the excess input", isTRUE(js("$('#tol_early').is(':visible')")))
  if (em == "fowler") {
    check("6 Fowler limit shown as 59", grepl("59.00 Gy₁₀ (Fowler zone lower edge, fixed)", text_of("search_context"), fixed = TRUE), text_of("search_context"))
    set_check("chemo", TRUE)
    check("6 chemotherapy lowers the Fowler limit to 54", grepl("54.00 Gy₁₀ (Fowler", text_of("search_context"), fixed = TRUE))
    set_check("chemo", FALSE)
    set_val("Tp_e", 3)
    check("6 changed acute parameters block the Fowler limit", grepl("applies only with Fowler", text_of("suggest_msg")))
    set_val("Tp_e", 2.5)
  }
}
set_radio("late_mode", "hypno")
check("6 HYPNO late limit shown", grepl("105.42 Gy₃ (HYPNO 55 Gy/20, fixed)", text_of("search_context"), fixed = TRUE))
set_radio("late_mode", "reference")
set_val("min_gain", 30)
check("6 empty search explains itself", grepl("No regimen in the searched range", text_of("suggest_msg")))
ze <- fetch("download")
check("6 empty-state download has its files", setequal(ze$files, c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt")))
set_val("min_gain", 1)
js("Array.from(document.querySelectorAll('#suggest_table th')).find(t => t.innerText.startsWith('Days')).click(); true"); idle()
check("6 table sorts by a column click", grepl("sorting_asc|sorting_desc", js("Array.from(document.querySelectorAll('#suggest_table th')).find(t => t.innerText.startsWith('Days')).className")))
check("6 plot renders", isTRUE(js("!!document.querySelector('#suggest_plot img')")))
z <- fetch("download")
res <- read.csv(file.path(z$dir, "results.csv"))
rr <- replay(z$dir)
check("9 candidate download replays in a fresh process to full precision",
      identical(rr$m, res$m) && isTRUE(all.equal(rr$d_di, res$d_di, tolerance = 1e-8)) && isTRUE(all.equal(rr$d_dd, res$d_dd, tolerance = 1e-8)))
sc <- jsonlite::fromJSON(file.path(z$dir, "scenario.json"))
check("9 export records 2.6 versions and unchanged code", identical(sc$app_version, "2.6.0") && identical(sc$selection_version, "2.6") &&
        identical(sc$code$changed_since_load, FALSE))
shot("03_candidates")

# 7. Agreement -----------------------------------------------------------------------
tab("Model agreement")
at <- text_of("agree_table")
check("7 default agreement rows with 15 of 15", grepl("15 of 15", at) && grepl("3.05", at) && grepl("2.90", at) && grepl("2.80", at))
check("7 agreement context names its own comparison and limits", grepl("Comparison: 20 × 2.75 Gy in 26 days (HYPNO test arm)", text_of("agree_context"), fixed = TRUE))
check("7 range columns say evaluable scenarios", grepl("over evaluable scenarios", at))
set_val("min_gain", "''"); set_val("min_days", "''"); tab("Model agreement")
check("7 agreement survives cleared candidate-only inputs", grepl("15 of 15", text_of("agree_table")))
set_val("min_gain", 1); set_val("min_days", 7)
set_val("agree_cmp", "sidebar")
check("7 sidebar comparison switches the agreement comparison", grepl("Comparison: 35 × 2 Gy in 47 days (the sidebar comparison)", text_of("agree_context"), fixed = TRUE))
set_val("agree_cmp", "hypno_test")
wide <- js("(function(){var w=document.querySelector('#agree_table').closest('.table-scroll'); return getComputedStyle(w).overflowX;})()")
check("7 agreement table scrolls inside its box", identical(wide, "auto") && no_overflow())
za <- fetch("download_agree")
det <- read.csv(file.path(za$dir, "scenario_detail.csv")); ra <- read.csv(file.path(za$dir, "results.csv"))
rra <- replay(za$dir)
check("9 agreement download replays to full precision, details included",
      nrow(det) == 45 && identical(rra$m, ra$m) && isTRUE(all.equal(attr(rra, "scenario_detail")$d_dd, det$d_dd, tolerance = 1e-8)))
ok_rng <- all(vapply(seq_len(nrow(ra)), function(i) { d <- det[det$m == ra$m[i] & det$valid_di & det$valid_dd, ]
  isTRUE(all.equal(c(ra$d_dd_min[i], ra$d_dd_max[i]), range(d$d_dd), tolerance = 1e-10)) }, TRUE))
check("9 agreement ranges reconstruct from scenario_detail", ok_rng)
shot("04_agreement")

# 8. Sensitivity ---------------------------------------------------------------------
tab("Sensitivity")
check("8 sensitivity: 15 rows, no errors", js("document.querySelectorAll('#sens_table tbody tr').length") == 15 && !nzchar(errors()))
check("8 sensitivity context: the sidebar regimen and the threshold", grepl("20 × 2.65 Gy in 26 days (the sidebar regimen", text_of("sens_context"), fixed = TRUE) &&
        grepl("at least 1 Gy", text_of("sens_context"), fixed = TRUE))
check("8 inherited DD rates marked", grepl("(sidebar)", text_of("sens_table"), fixed = TRUE))
shot("05_sensitivity")

# 10. Model and sources ----------------------------------------------------------------
tab("Model and sources")
idle(3)
mj <- js("document.querySelectorAll('mjx-container, .MathJax, .MathJax_Display').length")
check("10 equations rendered by MathJax", mj > 0, as.character(mj))
check("10 no raw TeX left in the page", !grepl("\\\\frac|\\$\\$", js("document.getElementById('about').innerText")))
check("10 version footer from the code", grepl("Version 2.6.0 of the Custom BED calculator (equations 2.2, selection rules 2.6", text_of("about"), fixed = TRUE))
check("10 conditional illustration present", grepl("Conditional illustration, not a prediction", text_of("about"), fixed = TRUE))
nl <- js("document.querySelectorAll('#about a[href^=\"http\"]').length")
check("10 source links present", nl >= 16, as.character(nl))
shot("06_model")

# Narrow width ---------------------------------------------------------------------
b$Emulation$setDeviceMetricsOverride(width = 420, height = 900, deviceScaleFactor = 1, mobile = TRUE); idle(2)
for (t in c("Results", "Candidate regimens", "Model agreement", "Sensitivity", "Model and sources")) {
  tab(t)
  check(sprintf("narrow 420 px: %s has no page-level horizontal scroll", t), no_overflow(), paste(overflow_px(), "px;", widest()))
  shot(paste0("07_narrow_", gsub("[^a-z]+", "_", tolower(t))))
}
b$Emulation$clearDeviceMetricsOverride()

check("browser console: no errors", length(console) == 0, paste(head(console, 5), collapse = " | "))
b$close()
ver <- callr::r(function() c(R.version.string, as.character(packageVersion("shiny"))), libpath = c(rel_lib, .Library))
rev <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "?")
dirty <- tryCatch(length(system("git status --porcelain -- R app.R about.md", intern = TRUE)) > 0, error = function(e) NA)
writeLines(c(sprintf("Walkthrough %s; revision %s%s; app in release library %s (%s, shiny %s); chromote %s, Chrome %s",
                     format(Sys.time()), rev, if (isTRUE(dirty)) " with uncommitted app changes" else "", rel_lib, ver[1], ver[2],
                     packageVersion("chromote"),
                     tryCatch(chromote::default_chromote_object()$get_browser()$get_path(), error = function(e) "?")),
             log), file.path(out, "walkthrough_log.txt"))
cat(sum(startsWith(log, "PASS")), "passed,", sum(startsWith(log, "FAIL")), "failed\n")
