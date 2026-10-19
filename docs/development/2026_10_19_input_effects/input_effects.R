# Input-effect audit: change every menu box of the running app, one at a time,
# in a headless browser, and record which outputs change. Run from the app folder:
#   Rscript docs/development/2026_10_19_input_effects/input_effects.R
# Each case reloads the page (all defaults), applies optional base actions,
# takes a snapshot of every output on every tab and of both downloads, applies
# the change, and takes a second snapshot. A case passes when every output in
# `must` changed, no output in `keep` changed, no output shows an error, and
# something changed at all. Writes input_effects_log.txt and input_effects.csv.
library(chromote)
root <- normalizePath(".")
out <- file.path(root, "docs", "development", "2026_10_19_input_effects")
port <- 8793
app <- callr::r_bg(function(root, port) shiny::runApp(root, port = port, launch.browser = FALSE),
                   args = list(root = root, port = port))
on.exit(app$kill(), add = TRUE)  # inside a function only; top level relies on the finalizer below
reg.finalizer(environment(), function(e) try(app$kill(), silent = TRUE), onexit = TRUE)
base <- sprintf("http://127.0.0.1:%d", port)
for (i in 1:60) { ok <- tryCatch({ readLines(base, warn = FALSE); TRUE }, error = function(e) FALSE); if (ok) break; Sys.sleep(1) }
stopifnot(ok)

b <- ChromoteSession$new(width = 1400, height = 1000)
js <- function(expr) b$Runtime$evaluate(expr, returnByValue = TRUE)$result$value
idle <- function(min_wait = 1) {
  Sys.sleep(min_wait)
  for (i in 1:120) { if (isTRUE(js("!$('html').hasClass('shiny-busy')"))) break; Sys.sleep(0.5) }
  Sys.sleep(0.4)
}
tab <- function(label) { js(sprintf("Array.from(document.querySelectorAll('.nav-link')).find(a => a.innerText.trim().startsWith('%s')).click()", label)); idle() }
load <- function() {
  # a page load occasionally times out in headless Chrome; retry it
  for (k in 1:3) {
    ok <- tryCatch({ b$Page$navigate(base); b$Page$loadEventFired(timeout_ = 60); TRUE }, error = function(e) FALSE)
    if (ok) break
  }
  stopifnot(ok); idle(3)
  # open every sidebar panel so that its outputs are not suspended
  js("document.querySelectorAll('.accordion-button.collapsed').forEach(x => x.click())"); idle()
}

# Actions on inputs, through the same widgets a user operates.
num <- function(id, v) sprintf("$('#%s').val(%s).trigger('change')", id, v)
sel <- function(id, v) sprintf("$('#%s')[0].selectize.setValue('%s')", id, v)
chk <- function(id, v) sprintf("$('#%s').prop('checked', %s).trigger('change')", id, tolower(v))
rad <- function(id, v) sprintf("$('input[name=%s][value=%s]').click()", id, v)
sld <- function(id, a, z) sprintf("$('#%s').data('ionRangeSlider').update({from: %d, to: %d}); $('#%s').trigger('change')", id, a, z, id)
apply_actions <- function(a) for (x in a) { js(sprintf("void(function(){ %s })()", x)); idle() }

text_of <- function(id) js(sprintf("(document.getElementById('%s') || {innerText: '<absent>'}).innerText", id))
hash <- function(x) digest::digest(x, algo = "md5")
# A download gives two keys: its result rows (CSV files) and its scenario file.
# The scenario also records interface context (for example the Results-tab
# regimen) that does not change any result, so the two are kept apart.
fetch <- function(id, csvs) {
  href <- js(sprintf("document.getElementById('%s').getAttribute('href')", id))
  f <- tempfile(fileext = ".zip")
  ok <- tryCatch({ download.file(paste0(base, "/", href), f, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
  if (!ok) return(c(rows = "<download failed>", scenario = "<download failed>"))
  d <- tempfile(); utils::unzip(f, exdir = d)
  rows <- vapply(csvs, function(x) {
    p <- file.path(d, x); if (file.exists(p)) paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n") else "<missing>"
  }, "")
  sc <- jsonlite::fromJSON(file.path(d, "scenario.json"), simplifyVector = FALSE)
  sc$exported_at <- NULL
  c(rows = hash(rows), scenario = hash(jsonlite::toJSON(sc, auto_unbox = TRUE, digits = NA)))
}
errors <- function() js("Array.from(document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')).map(e => e.id).join(',')")

snapshot <- function() {
  s <- list()
  for (id in c("T_shown", "ref_T_shown", "preset_note", "tumor_note", "anchor_note", "lambda_e_shown"))
    s[[id]] <- hash(text_of(id))
  tab("Results")
  rows <- js("Array.from(document.querySelectorAll('#results tbody tr')).map(r => Array.from(r.cells).slice(1).map(c => c.innerText).join('|'))")
  rn <- c("r_T", "r_total", "r_DI", "r_DD", "r_early", "r_eqd2e", "r_late", "r_eqd2l", "r_late_s1", "r_late_s2")
  for (i in seq_along(rn)) s[[rn[i]]] <- hash(if (i <= length(rows)) rows[[i]] else "<none>")
  for (id in c("dd_warning", "neg_warning", "verdicts", "acute_zone", "benchmarks", "hypno_note")) s[[id]] <- hash(text_of(id))
  err <- errors()
  tab("Candidate regimens")
  for (id in c("search_note", "hypno_cap_note", "suggest_msg", "suggest_table")) s[[id]] <- hash(text_of(id))
  s$suggest_plot <- hash(js("(document.querySelector('#suggest_plot img') || {src: ''}).src"))
  z <- fetch("download", c("results.csv", "all_qualifying.csv")); s$download <- z[["rows"]]; s$download_scenario <- z[["scenario"]]
  err <- paste(err, errors())
  tab("Model agreement")
  for (id in c("agree_note", "agree_msg", "agree_table")) s[[id]] <- hash(text_of(id))
  z <- fetch("download_agree", c("results.csv", "scenario_detail.csv")); s$download_agree <- z[["rows"]]; s$download_agree_scenario <- z[["scenario"]]
  err <- paste(err, errors())
  tab("Sensitivity")
  for (id in c("sens_summary", "sens_table", "sens_tox")) s[[id]] <- hash(text_of(id))
  err <- trimws(paste(err, errors()))
  list(h = unlist(s), err = err)
}

C <- function(name, change, must = character(), keep = character(), base = character(), note = "")
  list(name = name, change = change, must = must, keep = keep, base = base, note = note)
TAB_C <- c("suggest_table", "download")
TAB_A <- c("agree_table", "download_agree")
cases <- list(
  # Regimen evaluated (Results tab only)
  C("m (fractions)", num("m", 18), c("r_T", "r_DI", "T_shown", "sens_table"), c(TAB_C, TAB_A, "benchmarks", "anchor_note")),
  C("d (dose per fraction)", num("d", 3), c("r_DI", "r_DD", "r_early", "r_late", "sens_table"), c(TAB_C, TAB_A, "benchmarks")),
  C("auto_T off + T = 33", c(chk("auto_T", FALSE), num("T", 33)), c("r_T", "r_DI", "T_shown"), c(TAB_C, TAB_A, "r_late")),
  # Comparison regimen
  C("comparator = HYPNO test", sel("comparator", "hypno_test"), c("r_DI", "benchmarks", "ref_T_shown", TAB_C), c(TAB_A, "anchor_note")),
  C("comparator = HYPNO control", sel("comparator", "hypno_control"), c("r_DI", "benchmarks", "ref_T_shown", TAB_C), c(TAB_A, "anchor_note")),
  C("ref_m", num("ref_m", 30), c("r_DI", "benchmarks", "ref_T_shown", TAB_C), c(TAB_A, "anchor_note")),
  C("ref_d", num("ref_d", 2.2), c("r_DI", "benchmarks", "ref_T_shown", TAB_C), c(TAB_A, "anchor_note")),
  C("ref_auto_T off + ref_T = 50", c(chk("ref_auto_T", FALSE), num("ref_T", 50)), c("r_DI", "ref_T_shown", TAB_C), c(TAB_A, "anchor_note", "r_late")),
  # Tumor parameters
  C("preset = version 1 values", sel("preset", "calculator"), c("preset_note", "tumor_note", "anchor_note", "r_DI", "r_DD", TAB_C, TAB_A), c("r_early", "r_late")),
  C("preset = Withers", sel("preset", "withers"), c("preset_note", "tumor_note", "r_DI", TAB_C), c("r_early", "r_late")),
  C("preset = Hartley", sel("preset", "hartley"), c("preset_note", "tumor_note", "r_DI", TAB_C), c("r_early", "r_late")),
  C("preset = HYHOPE", sel("preset", "hyhope"), c("preset_note", "tumor_note", "r_DI", TAB_C), c("r_early", "r_late")),
  C("alpha", num("alpha", 0.3), c("preset_note", "tumor_note", "anchor_note", "r_DI", "r_DD", TAB_C, TAB_A), c("r_early", "r_late")),
  C("ab (tumor)", num("ab", 7), c("preset_note", "tumor_note", "r_DI", "r_DD", TAB_C), c("r_early", "r_late")),
  C("lambda_di", num("lambda_di", 0.3), c("preset_note", "tumor_note", "anchor_note", "r_DI", TAB_C, TAB_A), c("r_early", "r_late")),
  C("lambda_dd", num("lambda_dd", 0.5), c("preset_note", "anchor_note", "r_DD", TAB_C, TAB_A, "sens_table"), c("r_DI", "tumor_note", "r_early")),
  C("Tk", num("Tk", 21), c("preset_note", "tumor_note", "anchor_note", "r_DI", TAB_C), c("r_early", "r_late")),
  C("g", num("g", 0.01), c("tumor_note", "r_DI", "r_DD", "sens_table", "agree_note", TAB_C), c("r_early", "r_late", "anchor_note")),
  # DD calibration anchor
  C("anchor_m", num("anchor_m", 33), c("anchor_note", "r_DD", TAB_A), c("r_DI", "r_early", "r_late", "tumor_note")),
  C("anchor_d", num("anchor_d", 2.2), c("anchor_note", "r_DD", TAB_A), c("r_DI", "r_early", "r_late", "tumor_note")),
  # Early toxicity
  C("ab_e", num("ab_e", 8), c("r_early", "acute_zone", "search_note", TAB_C), c("r_DI", "r_DD", "r_late")),
  C("alpha_e", num("alpha_e", 0.3), c("r_early", "acute_zone", TAB_C), c("r_DI", "r_DD", "r_late")),
  C("Tp_e", num("Tp_e", 3.5), c("lambda_e_shown", "r_early", "acute_zone", TAB_C), c("r_DI", "r_late")),
  C("Tk_e", num("Tk_e", 10), c("r_early", "acute_zone", TAB_C), c("r_DI", "r_late")),
  C("chemo (HYPNO early limit)", chk("chemo", TRUE), c("acute_zone", "download_scenario"), c(TAB_C, "r_early", "r_DI"),
    note = "with the default HYPNO limit chemotherapy changes only the Fowler zone text and the saved scenario"),
  C("chemo (Fowler early limit)", chk("chemo", TRUE), c("acute_zone", "search_note", TAB_C), c("r_early", "r_DI"), base = rad("early_mode", "fowler")),
  # Late toxicity
  C("ab_l", num("ab_l", 4), c("r_late", "r_late_s2"), c("r_DI", "r_DD", "r_early")),
  # Candidate tab
  C("criterion = both", sel("criterion", "both"), TAB_C, c(TAB_A, "r_DI")),
  C("criterion = DD", sel("criterion", "dd"), TAB_C, c(TAB_A, "r_DI")),
  C("criterion = DI", sel("criterion", "di"), "download_scenario", c(TAB_A, "r_DI"),
    note = "with the defaults no regimen qualifies only under DD, so the rows may stay the same"),
  C("rank_by = worst", sel("rank_by", "worst"), "download_scenario", c(TAB_A, "r_DI")),
  C("rank_by = DD", sel("rank_by", "dd"), "download_scenario", c(TAB_A, "r_DI")),
  # The default list is 28-30 fractions at 2.2-2.3 Gy, so 15-30 fractions or a 2.75 Gy maximum
  # would not change it: these cases use values that cut into it.
  C("m_range 29-30", sld("m_range", 29, 30), c(TAB_C, "download_scenario"), c(TAB_A, "r_DI")),
  C("min_gain = -3", num("min_gain", -3), c(TAB_C, "sens_summary"), c(TAB_A, "r_DI")),
  C("min_days = 7 (no effect expected)", num("min_days", 7), character(), c(TAB_C, TAB_A),
    note = "re-entering the default value must change nothing"),
  C("min_days = 3 (late excess 2 Gy)", num("min_days", 3), TAB_C, c(TAB_A, "r_DI"), base = num("tol_late", 2),
    note = "31 x 2.2 and 32 x 2.15 Gy save 4 and 3 days but exceed the default late limit by 1.4-1.5 Gy"),
  C("d_step = 0.25", sel("d_step", "0.25"), c(TAB_C, TAB_A), "r_DI"),
  C("d_max = 2.75", num("d_max", 2.75), TAB_A, c(TAB_C, "r_DI"),
    note = "removes the agreement regimens above 2.75 Gy; candidates are 2.2-2.3 Gy"),
  C("d_max = 2.25", num("d_max", 2.25), TAB_C, "r_DI"),
  C("early_mode = reference", rad("early_mode", "reference"), c("search_note", TAB_C), c(TAB_A, "r_DI"),
    note = "the agreement comparison is HYPNO itself, so 'comparison plus 0' equals the HYPNO limit there"),
  C("early_mode = Fowler", rad("early_mode", "fowler"), c("search_note", TAB_C, TAB_A), "r_DI"),
  C("early_mode = Fenwick", rad("early_mode", "fenwick"), c("search_note", TAB_C, TAB_A), "r_DI"),
  C("tol_early = 3", num("tol_early", 3), c(TAB_C), "r_DI", base = rad("early_mode", "reference")),
  C("late_mode = HYPNO", rad("late_mode", "hypno"), c("download_scenario", "search_note"), "r_DI"),
  C("tol_late = 5", num("tol_late", 5), c("download_scenario"), "r_DI"),
  # Agreement tab
  C("agree_cmp = sidebar", sel("agree_cmp", "sidebar"), c("agree_note", TAB_A), c(TAB_C, "r_DI")),
  C("agree_margin = -1", num("agree_margin", -1), TAB_A, c(TAB_C, "r_DI")),
  C("agree_days = 3", num("agree_days", 3), TAB_A, c(TAB_C, "r_DI"))
)

only <- Sys.getenv("INPUT_EFFECTS_ONLY")
if (nzchar(only)) cases <- Filter(function(cs) grepl(only, cs$name), cases)
load(); default_snap <- snapshot()
log <- character(); rows <- list()
for (cs in cases) {
  load()
  s0 <- if (length(cs$base)) { apply_actions(cs$base); snapshot() } else default_snap
  apply_actions(cs$change)
  s1 <- snapshot()
  changed <- names(s1$h)[s1$h != s0$h[names(s1$h)]]
  miss <- setdiff(cs$must, changed); moved <- intersect(cs$keep, changed)
  pass <- (length(changed) > 0 || !length(cs$must)) && !length(miss) && !length(moved) && !nzchar(s1$err)
  line <- sprintf("%s  %-30s changed: %s%s%s%s", if (pass) "PASS" else "FAIL", cs$name,
                  if (length(changed)) paste(changed, collapse = ",") else "<nothing>",
                  if (length(miss)) paste0("  | did NOT change: ", paste(miss, collapse = ",")) else "",
                  if (length(moved)) paste0("  | changed but should not: ", paste(moved, collapse = ",")) else "",
                  if (nzchar(s1$err)) paste0("  | output errors: ", s1$err) else "")
  cat(line, "\n"); log <- c(log, line)
  rows[[length(rows) + 1]] <- data.frame(case = cs$name, pass = pass, changed = paste(changed, collapse = " "),
                                         not_changed_but_required = paste(miss, collapse = " "),
                                         changed_but_should_not = paste(moved, collapse = " "),
                                         errors = s1$err, note = cs$note)
}
b$close()
rev <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "?")
dirty <- tryCatch(length(system2("git", c("status", "--porcelain", "--", "app.R", "R"), stdout = TRUE)) > 0, error = function(e) NA)
writeLines(c(sprintf("Input-effect audit %s, revision %s%s, chromote %s", format(Sys.time()), rev,
                     if (isTRUE(dirty)) " plus uncommitted app changes" else "", packageVersion("chromote")),
             sprintf("Default snapshot output errors: %s", if (nzchar(default_snap$err)) default_snap$err else "none"),
             log, sprintf("%d passed, %d failed", sum(startsWith(log, "PASS")), sum(startsWith(log, "FAIL")))),
           file.path(out, "input_effects_log.txt"))
write.csv(do.call(rbind, rows), file.path(out, "input_effects.csv"), row.names = FALSE)
cat(sum(startsWith(log, "PASS")), "passed,", sum(startsWith(log, "FAIL")), "failed\n")
