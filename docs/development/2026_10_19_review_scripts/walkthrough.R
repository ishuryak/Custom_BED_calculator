# Headless-browser walkthrough of the app (chromote). Run from the app folder:
#   Rscript docs/development/2026_10_19_review_scripts/walkthrough.R
# Starts the app on a local port, drives every tab, changes and clears inputs,
# fetches both downloads and writes screenshots and a log to
# docs/development/2026_10_19_walkthrough/.
library(chromote)
root <- normalizePath(".")
out <- file.path(root, "docs", "development", "2026_10_19_walkthrough")
dir.create(out, showWarnings = FALSE)
port <- 8791
app <- callr::r_bg(function(root, port) shiny::runApp(root, port = port, launch.browser = FALSE),
                   args = list(root = root, port = port))
on.exit(app$kill(), add = TRUE)
base <- sprintf("http://127.0.0.1:%d", port)
for (i in 1:60) { ok <- tryCatch({ readLines(base, warn = FALSE); TRUE }, error = function(e) FALSE); if (ok) break; Sys.sleep(1) }
stopifnot(ok)

b <- ChromoteSession$new(width = 1400, height = 1000)
js <- function(expr) b$Runtime$evaluate(expr, returnByValue = TRUE)$result$value
idle <- function(min_wait = 1.5) {
  Sys.sleep(min_wait)
  for (i in 1:60) { if (isTRUE(js("!$('html').hasClass('shiny-busy')"))) break; Sys.sleep(0.5) }
  Sys.sleep(0.5)
}
log <- character()
check <- function(name, cond, detail = "") {
  line <- sprintf("%s  %s%s", if (isTRUE(cond)) "PASS" else "FAIL", name, if (nzchar(detail)) paste0("  [", detail, "]") else "")
  log <<- c(log, line); cat(line, "\n")
}
errors <- function() js("Array.from(document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')).map(e => e.id + ': ' + e.innerText).join(' | ')")
tab <- function(label) { js(sprintf("Array.from(document.querySelectorAll('.nav-link')).find(a => a.innerText.trim().startsWith('%s')).click()", label)); idle() }
set_val <- function(id, v) { js(sprintf("$('#%s').val(%s).trigger('change')", id, if (is.character(v)) shQuote(v) else v)); idle() }
text_of <- function(id) js(sprintf("(document.getElementById('%s') || {innerText: ''}).innerText", id))
shot <- function(name) b$screenshot(file.path(out, paste0(name, ".png")), show = FALSE)
fetch <- function(id) {
  href <- js(sprintf("document.getElementById('%s').getAttribute('href')", id))
  f <- tempfile(fileext = ".zip"); download.file(paste0(base, "/", href), f, mode = "wb", quiet = TRUE)
  d <- tempfile(); utils::unzip(f, exdir = d); list(files = list.files(d), dir = d)
}

b$Page$navigate(base); b$Page$loadEventFired(); idle(4)
check("app title shows version 2.5", grepl("version 2.5", js("document.body.innerText")))
check("Results tab: no output errors", !nzchar(errors()), errors())
check("Results tab: units on BED rows", grepl("Early toxicity BED (Gy₁₀)", text_of("results"), fixed = TRUE))
shot("01_results")

tab("Candidate regimens")
check("candidate tab: no output errors", !nzchar(errors()), errors())
n_rows <- js("document.querySelectorAll('#suggest_table tbody tr').length")
check("candidate tab: 3 default rows", identical(as.integer(n_rows), 3L), as.character(n_rows))
check("candidate tab: renamed fraction-size column",
      grepl("Above benchmark fraction sizes", js("document.getElementById('suggest_table').innerText")))
check("candidate tab: no uncertainty claim", !grepl("within the uncertainty", js("document.body.innerText")))
z <- fetch("download")
check("candidate download: 4 files", setequal(z$files, c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt")),
      paste(z$files, collapse = ","))
check("candidate download: 3 rows", nrow(read.csv(file.path(z$dir, "results.csv"))) == 3)
shot("02_candidates")

tab("Model agreement")
check("agreement tab: no output errors", !nzchar(errors()), errors())
at <- text_of("agree_table")
check("agreement tab: 17, 18, 19 fractions with 15 of 15", grepl("15 of 15", at) && grepl("3.05", at) && grepl("2.90", at) && grepl("2.80", at))
check("agreement note states the anchor", grepl("DD anchor: 35 × 2 Gy in 47 days (fixed)", text_of("agree_note"), fixed = TRUE))
wide <- js("(function(){var w=document.querySelector('#agree_table').closest('.table-scroll'); var t=w.querySelector('table'); return [t.scrollWidth, w.clientWidth, getComputedStyle(w).overflowX, document.documentElement.scrollWidth, document.documentElement.clientWidth];})()")
check("agreement table: wider than its box only inside a horizontal scroller, page does not overflow",
      identical(wide[[3]], "auto") && as.numeric(wide[[4]]) <= as.numeric(wide[[5]]) + 1, paste(unlist(wide), collapse = " "))
js("document.querySelector('#agree_table').closest('.table-scroll').scrollLeft = 10000"); Sys.sleep(0.5)
shot("03b_agreement_scrolled")
za <- fetch("download_agree")
check("agreement download: includes scenario_detail.csv", "scenario_detail.csv" %in% za$files, paste(za$files, collapse = ","))
check("agreement download: 45 scenario rows", nrow(read.csv(file.path(za$dir, "scenario_detail.csv"))) == 45)
shot("03_agreement")
set_val("min_gain", "''")                     # clear a candidate-only input
tab("Model agreement")
check("agreement survives a cleared candidate minimum", grepl("15 of 15", text_of("agree_table")))
set_val("min_gain", 1)
js("document.querySelectorAll('.accordion-button.collapsed').forEach(x => x.click())"); idle()
check("no legacy coupling switch", isTRUE(js("document.getElementById('coupled') === null")))
check("anchor note explains the anchor", grepl("same tumor BED for this anchor course", text_of("anchor_note"), fixed = TRUE))
check("preset note explains the default", grepl("The default. Scenario D18", text_of("preset_note"), fixed = TRUE))
set_val("lambda_dd", 0)
check("DD unavailable is explained on the agreement tab", grepl("DD model is unavailable", text_of("agree_msg")))
set_val("lambda_dd", 0.4)

tab("Candidate regimens")
set_val("ab_e", 5)
check("custom mucosal alpha/beta: cap labeled in BED Gy at alpha/beta 5",
      grepl("71.0 Gy (α/β = 5 Gy)", text_of("search_note"), fixed = TRUE))
set_val("ab_e", 10)

tab("Sensitivity")
check("sensitivity tab: no output errors", !nzchar(errors()), errors())
check("sensitivity tab: 15 scenario rows", js("document.querySelectorAll('#sens_table tbody tr').length") == 15)
shot("04_sensitivity")

tab("Model and sources")
mj <- js("document.querySelectorAll('mjx-container, .MathJax, .MathJax_Display').length")
check("Model page: equations rendered by MathJax", mj > 0, as.character(mj))
check("Model page: provenance table", grepl("Origin of preset values", js("document.body.innerText")))
shot("05_model")

b$close()
writeLines(c(sprintf("Walkthrough %s, app at %s, chromote %s, Chrome %s", format(Sys.time()), root,
                     packageVersion("chromote"), tryCatch(chromote::default_chromote_object()$get_browser()$get_path(), error = function(e) "?")),
             log), file.path(out, "walkthrough_log.txt"))
cat(sum(startsWith(log, "PASS")), "passed,", sum(startsWith(log, "FAIL")), "failed\n")
