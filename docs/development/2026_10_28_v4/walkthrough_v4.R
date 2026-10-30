# Headless-browser walkthrough of version 4 (main page + menu). Run from the app
# folder after tools/restore_library.R:
#   Rscript docs/development/2026_10_28_v4/walkthrough_v4.R              # local app
#   Rscript docs/development/2026_10_28_v4/walkthrough_v4.R --live URL   # deployed app
#   Rscript docs/development/2026_10_28_v4/walkthrough_v4.R --self-test-fail
# Writes a log and screenshots to walkthrough/<time>; any failed check gives exit 1.
library(chromote)
`%||%` <- function(a, b) if (is.null(a)) b else a

main <- function(args) {
  root <- normalizePath(".")
  source(file.path(root, "tools", "env.R"))
  rel_lib <- bed_library()
  stopifnot(dir.exists(file.path(rel_lib, "shiny")))
  live <- "--live" %in% args
  base <- if (live) sub("/$", "", args[which(args == "--live") + 1]) else NULL
  out <- file.path(root, "docs", "development", "2026_10_28_v4", "walkthrough",
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
    port <- 8797
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
    invisible(FALSE)
  }
  shot <- function(name) b$screenshot(file.path(out, paste0(name, ".png")), selector = "body", show = FALSE)
  visible_errors <- function() js(paste0(
    "Array.from(document.querySelectorAll('.shiny-output-error')).filter(e => e.offsetParent !== null)",
    ".map(e => e.textContent.trim()).filter(t => t.length > 0).join(' | ')"))

  b$Page$navigate(base); Sys.sleep(3); idle(2)
  # The main page, before anything is touched.
  check("main page is shown first", isTRUE(js("$('.tab-pane.active').attr('data-value') === 'home'")))
  check("main heading", isTRUE(js("document.body.innerText.includes('Hypofractionated alternatives to 70 Gy in 35 fractions')")))
  n <- js("$('#home_table table.dataTable tbody tr').length")
  check("main table has 15 rows at the defaults", identical(as.integer(n), 15L), paste("rows:", n))
  first <- js("$('#home_table table.dataTable tbody tr:first td').map((i,e)=>e.innerText).get().join(' | ')")
  check("first row is 30 fractions", grepl("^30 \\|", first), first)
  check("early and late toxicity columns", isTRUE(js("$('#home_table').text().includes('Early toxicity BED change') && $('#home_table').text().includes('Late toxicity BED change')")))
  check("details hidden by default", isTRUE(js("!$('#home_table').text().includes('Refit support') && !$('#home_table').text().includes('Notes')")))
  check("BED change plot below the table", isTRUE(js("$('#home_plot img').length === 1 && $('#home_plot').offset().top > $('#home_table').offset().top")))
  shot("01a_main_no_details")
  js("$('#home_details').prop('checked', true).trigger('change')"); idle(2)
  check("HYPNO arm marked", isTRUE(js("$('#home_table').text().includes('HYPNO test arm (55 Gy in 20) inside')")))
  check("extrapolation note on the main page", isTRUE(js("$('#home_table').text().includes('outside trial data')")))
  check("refit support column", isTRUE(js("$('#home_table').text().includes('Refit support')")))
  js("$('#home_details').prop('checked', false).trigger('change')"); idle(2)
  check("caution above the table", isTRUE(js("$('.v4-caution').offset().top < $('#home_table').offset().top")))
  check("menu starts closed", isTRUE(js("$('.bslib-sidebar-layout').hasClass('sidebar-collapsed')")))
  check("no banner on the main page", isTRUE(js("$('.bed-banner').is(':hidden')")))
  check("no output errors on the main page", !nzchar(visible_errors() %||% ""), visible_errors() %||% "")
  shot("01_main")

  # A real mouse click on the 20-fraction row (synthetic jQuery clicks do not reach DT's selection).
  js("document.querySelectorAll('#home_table table.dataTable tbody tr')[10].scrollIntoView({block: 'center'})"); Sys.sleep(0.5)
  xy <- js("(() => { const r = document.querySelectorAll('#home_table table.dataTable tbody tr')[10].cells[1].getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; })()")
  for (ty in c("mousePressed", "mouseReleased"))
    b$Input$dispatchMouseEvent(type = ty, x = xy[[1]], y = xy[[2]], button = "left", clickCount = 1)
  idle(2)
  det <- js("$('#home_detail').text()")
  check("row click shows the dose-by-dose card", grepl("20 fractions in 26 days", det %||% ""), substr(det %||% "", 1, 80))
  top <- js("document.getElementById('home_detail').getBoundingClientRect().top")
  check("dose card is scrolled into view", is.numeric(top) && top >= -5 && top < 900, paste("top:", round(top %||% NA)))
  shot("02_main_row_selected")

  # Open the selected regimen (20 x 2.75 Gy) on the Report card with the main-page rules.
  js("$('#home_open_dose').val('2.75'); Shiny.setInputValue('home_open_dose', '2.75')"); idle()
  bx <- js("(() => { const el = document.getElementById('home_open_report'); el.scrollIntoView({block: 'center'}); const r = el.getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; })()")
  for (ty in c("mousePressed", "mouseReleased"))
    b$Input$dispatchMouseEvent(type = ty, x = bx[[1]], y = bx[[2]], button = "left", clickCount = 1)
  idle(3)
  check("Open on the Report card goes to the Report card", identical(js("$('.tab-pane.active').attr('data-value')"), "Report card"))
  idle(2)
  md <- c(js("String($('#m').val())") %||% "NULL", js("String($('#d').val())") %||% "NULL")
  check("Report card got the regimen", identical(md, c("20", "2.75")), paste(md, collapse = " x "))
  check("Report card opens at the top", isTRUE(js("window.scrollY < 50")), paste("scrollY:", js("window.scrollY")))
  rep_txt <- js("$('#report').text()") %||% ""
  check("Report card agrees with the main list", grepl("threshold -5: met", rep_txt, fixed = TRUE) && !grepl("exceeds", rep_txt),
        substr(rep_txt, 1, 100))
  shot("02b_report_from_main")
  js("$('#back_report').click()"); idle()

  # Open the menu and visit every page.
  js("$('#open_menu').click()"); idle()
  check("menu opens", isTRUE(js("!$('.bslib-sidebar-layout').hasClass('sidebar-collapsed')")))
  shot("03_menu_open")
  pages <- c(report = "Report card", results = "Results", compare = "Compare regimens", search = "Search",
             agree = "Both-model screen", sens = "Sensitivity", unc = "Uncertainty", oar = "Organs at risk",
             settings = "Model settings", about = "Model and sources")
  k <- 4
  for (slug in names(pages)) {
    js(sprintf("$('#go_%s').click()", slug)); idle(2)
    act <- js("$('.tab-pane.active').attr('data-value')")
    check(sprintf("menu opens %s", pages[[slug]]), identical(act, pages[[slug]]), act %||% "")
    err <- visible_errors() %||% ""
    check(sprintf("%s has no output errors", pages[[slug]]), !nzchar(err), substr(err, 1, 120))
    shot(sprintf("%02d_%s", k, slug)); k <- k + 1
  }
  check("banner shown on a detail page", isTRUE(js("$('.bed-banner').is(':visible')")))
  # A setting changes the main list.
  js("$('#go_settings').click()"); idle()
  js("$('#ab_l').val(2).trigger('change')"); idle(2)
  js("$('#back_settings').click()"); idle(2)
  check("back link returns to the main page", isTRUE(js("$('.tab-pane.active').attr('data-value') === 'home'")))
  check("main page reports changed settings", isTRUE(js("$('#home_settings').text().includes('Model settings changed')")))
  n2 <- js("$('#home_table table.dataTable tbody tr').length")
  check("list follows the late α/β setting", !identical(as.integer(n2), 15L) ||
          !identical(js("$('#home_table').text()"), first), paste("rows:", n2))

  # Phone width.
  b2 <- ChromoteSession$new(width = 390, height = 844)
  b2$Page$navigate(base); Sys.sleep(4)
  js2 <- function(e) b2$Runtime$evaluate(e, returnByValue = TRUE)$result$value
  wide <- js2("document.documentElement.scrollWidth > window.innerWidth + 1")
  check("no horizontal page scroll at 390 px", !isTRUE(wide))
  check("first column frozen on the phone", isTRUE(js2("document.querySelectorAll('#home_table .dtfc-fixed-left, #home_table .dtfc-fixed-start').length > 0")))
  js2("$('#open_menu').click()"); Sys.sleep(1.5)
  js2("$('#go_report').click()"); Sys.sleep(2)
  check("phone: menu closes after choosing a page", isTRUE(js2("$('.bslib-sidebar-layout').first().hasClass('sidebar-collapsed')")))
  js2("$('#back_report').click()"); Sys.sleep(1.5)
  b2$screenshot(file.path(out, "20_phone_main.png"), selector = "body", show = FALSE)
  b2$close()

  check("no browser console errors", length(console) == 0, paste(head(console, 3), collapse = " | "))
  if ("--self-test-fail" %in% args) check("self-test: this check must fail", FALSE)
  writeLines(log, file.path(out, "walkthrough.log"))
  n_fail <- sum(startsWith(log, "FAIL"))
  cat(sprintf("\n%d checks, %d failed. Log and screenshots: %s\n", length(log), n_fail, out))
  if (n_fail > 0) quit(status = 1)
}
main(commandArgs(trailingOnly = TRUE))
