# Export and exact replay of the version 4 main list (review F08 of 2026-10-31).
# A Candidate ZIP records the Search page's settings (one ranked dose per number
# of fractions); it is not a record of the main list. This file writes the main
# analysis itself: every qualifying grid point at full precision, the
# per-fraction summary rows, the refit support with the draws it came from, and
# every input needed to recompute them. The biological equations are those of
# R/bed_models.R; only the display analysis gains a schema.
# Requires R/bed_models.R, R/uncertainty.R, R/export.R and R/v4_home.R.

MAIN_SCHEMA_VERSION <- "main-1"
MAIN_ANALYSIS_VERSION <- "1.0"   # v4_regimens + v4_refit_support + v4_summary, as in app 4.3.0
MAIN_SUPPORT_VARIANT <- "primary"

# Fields that determine the result (hashed for the analysis ID). App version,
# code identity, export time and interface context (for example the Details
# toggle) are excluded.
MAIN_SCIENTIFIC_FIELDS <- c("schema_version", "main_version", "equation_version", "selection_version",
                            "analysis", "hinge", "clocks", "params", "dd", "anchor", "standard", "margins",
                            "grid", "calendar", "rule", "support", "data_domain", "trial_arms")

# Draws used for refit support: list(di, dd) from load_draws(), or NULL.
main_support_block <- function(draws, margin = V4_MARGINS$tumor, threshold = V4_SUPPORT_MIN) {
  if (is.null(draws))
    return(list(status = "unavailable", reason = "the refit draws could not be read when the list was computed"))
  list(status = "available", variant = MAIN_SUPPORT_VARIANT,
       files = list(di = list(file = attr(draws$di, "file"), sha256 = attr(draws$di, "sha256")),
                    dd = list(file = attr(draws$dd, "file"), sha256 = attr(draws$dd, "sha256"))),
       n_rows = nrow(draws$di),
       pairing = "row k of the DI draws paired with row k of the conditional calibrated-DD draws; every row weighted equally",
       event = sprintf("tumor BED change >= -%g under DI or DD, each draw at its own alpha/beta", margin),
       margin = margin, flag_below = threshold, engine = "v4_refit_support",
       meaning = "share of refit draws; not a probability of tumor control, safety or non-inferiority")
}

# The scenario of the main list for the current settings. draws: as above.
main_scenario <- function(p, anchor, draws, data_domain = tumor_data_support(), context = list(),
                          code = NULL, code_root = ".") {
  pp <- lapply(setNames(PARAM_NAMES, PARAM_NAMES), function(n) if (is.null(p[[n]]) && n == "g") 0 else p[[n]])
  cal <- dd_setup(p, anchor)
  sc <- list(
    schema_version = MAIN_SCHEMA_VERSION, main_version = MAIN_ANALYSIS_VERSION,
    app_version = APP_VERSION, equation_version = EQUATION_VERSION, selection_version = SELECTION_VERSION,
    analysis = "main",
    hinge = list(type = "softplus", k_per_day = SMOOTH_K),
    clocks = list(input = "inclusive calendar days, first to last treatment day",
                  tumor = "inclusive", acute = if (is.null(p$clock_e)) "elapsed" else p$clock_e),
    params = pp, dd = dd_block(p),
    anchor = anchor[c("m", "d", "T")],
    calibration = list(status = cal$status,
                       c_over_alpha = if (is.finite(cal$c_over_alpha)) cal$c_over_alpha else NULL,
                       onset_anchor = if (is.finite(cal$onset_ref)) cal$onset_ref else NULL),
    standard = V4_STANDARD, margins = V4_MARGINS,
    grid = list(m_range = V4_GRID$m_range, d_step = V4_GRID$d_step, d_max = V4_GRID$d_max, min_days_saved = 1),
    calendar = list(fx_per_week = 5, start_weekday = 1),
    rule = "listed when tumor (DI or DD), early and late conditions all hold against the standard",
    support = main_support_block(draws),
    data_domain = list(file = data_domain$file, d_max = data_domain$d_max,
                       n_arms = data_domain$n_arms, n_studies = data_domain$n_studies),
    trial_arms = V4_TRIAL_ARMS)
  now <- code_identity(code_root)
  if (is.null(code)) code <- now
  code$changed_since_load <- !identical(code$files, now$files)
  sc$code <- code
  sc$context <- context
  sc
}

main_id <- function(sc) {
  x <- sc[intersect(MAIN_SCIENTIFIC_FIELDS, names(sc))]
  js <- enc2utf8(as.character(scenario_json(canonicalize(x), pretty = FALSE)))
  substr(digest::digest(js, algo = "sha256", serialize = FALSE), 1, 12)
}

# The main list for given inputs: every qualifying point (with support when
# draws are given) and the summary rows.
main_compute <- function(p, anchor, draws, data_d_max) {
  h <- v4_regimens(p, anchor)
  pts <- h$rows
  if (!is.null(draws) && nrow(pts) > 0) pts$support <- v4_refit_support(pts, draws$di, draws$dd)
  hy <- comparator_presets()$hypno_test
  hypno_early <- compute_all(hy$m, hy$d, hy$T, p, anchor)$bed_early
  list(points = pts, summary = v4_summary(pts, hypno_early = hypno_early, data_d_max = data_d_max),
       reason = h$reason)
}

main_readme <- function(sc, id) {
  c(sprintf("HNSCC BED calculator %s, main list (analysis %s %s), ID %s", sc$app_version,
            sc$schema_version, sc$main_version, id),
    sprintf("Standard: %g x %g Gy in %g days. DD: %s, anchor %g x %g Gy in %g days.", sc$standard$m, sc$standard$d,
            sc$standard$T, sc$dd$mode, sc$anchor$m, sc$anchor$d, sc$anchor$T),
    sprintf("Listed when tumor BED >= standard - %g under DI or DD, early BED <= standard + %g, late BED <= standard + %g;",
            sc$margins$tumor, sc$margins$early, sc$margins$late),
    sprintf("  %g-%g fractions, 5 a week from a Monday, doses above 2 Gy in %g Gy steps up to %g Gy.",
            min(sc$grid$m_range), max(sc$grid$m_range), sc$grid$d_step, sc$grid$d_max),
    "",
    "points.csv: every qualifying (fractions, dose) point at full precision: absolute BEDs (bed_*), changes",
    "  from the standard (d_*; NA where a model is not usable), qual_di / qual_dd, and support (share of the",
    "  refit draws meeting the tumor condition under either model; absent when the draws were unavailable).",
    "summary.csv: the rows of the main page (one per number of fractions), with the ranges as displayed.",
    "  Ranges are the lowest and highest values over the doses kept, not uncertainty intervals.",
    "scenario.json: every input, the DD calibration, the draws files used for support (SHA-256), versions.",
    if (identical(sc$support$status, "available")) "Refit support: computed from the draws named in scenario.json." else
      paste("Refit support: NOT AVAILABLE in this export:", sc$support$reason),
    "",
    "To rerun in R from the app folder (refused if a version, input or draws file differs):",
    "  for (f in c('bed_models','schedule','uncertainty','oar','export','v4_home','v4_export')) source(sprintf('R/%s.R', f))",
    "  sc <- jsonlite::fromJSON('scenario.json', simplifyVector = TRUE); r <- main_rerun(sc)",
    "  r$points and r$summary match points.csv and summary.csv.",
    "",
    "Model outputs for trial design and hypothesis generation, not treatment recommendations.")
}

# Write the ZIP. r: main_compute() of exactly the inputs in sc.
write_main_zip <- function(file, sc, r) {
  id <- main_id(sc)
  sc$analysis_id <- id
  sc$result <- list(points = nrow(r$points), rows = nrow(r$summary), empty_reason = r$reason)
  sc$exported_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  dir <- tempfile("bed_main_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  utils::write.csv(r$points, file.path(dir, "points.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(r$summary, file.path(dir, "summary.csv"), row.names = FALSE, na = "NA")
  writeLines(enc2utf8(as.character(scenario_json(sc))), file.path(dir, "scenario.json"), useBytes = TRUE)
  writeLines(main_readme(sc, id), file.path(dir, "README.txt"))
  zip::zip(file, files = c("points.csv", "summary.csv", "scenario.json", "README.txt"), root = dir)
  invisible(id)
}

# Strict validation of a main scenario read back with jsonlite::fromJSON().
# Returns the inputs to recompute. Refuses another version, a changed fixed
# rule, a different trial-data domain, or support draws that are not the ones
# recorded.
validate_main_scenario <- function(sc, draws_loader = function() list(di = load_draws("di"), dd = load_draws("dd")),
                                   data_domain = tumor_data_support()) {
  same <- function(a, b) isTRUE(all.equal(unlist(a), unlist(b), check.attributes = FALSE, tolerance = 1e-12))
  for (k in c("schema_version", "main_version", "equation_version", "selection_version")) {
    cur <- list(schema_version = MAIN_SCHEMA_VERSION, main_version = MAIN_ANALYSIS_VERSION,
                equation_version = EQUATION_VERSION, selection_version = SELECTION_VERSION)[[k]]
    if (!identical(as.character(sc[[k]]), cur))
      replay_fail(k, " '", paste(sc[[k]], collapse = ","), "' is not this calculator's (", cur, ").")
  }
  if (!identical(sc$analysis, "main")) replay_fail("this is not a main-list export.")
  if (!same(sc$standard, V4_STANDARD) || !same(sc$margins, V4_MARGINS) ||
      !same(sc$grid[c("m_range", "d_step", "d_max")], V4_GRID[c("m_range", "d_step", "d_max")]) ||
      !same(sc$grid$min_days_saved, 1) || !same(sc$calendar, list(fx_per_week = 5, start_weekday = 1)))
    replay_fail("its fixed main-list rules (standard, margins, grid or calendar) differ from this calculator's.")
  if (!same(sc$data_domain$d_max, data_domain$d_max) || !identical(sc$data_domain$file, data_domain$file))
    replay_fail("its trial-data domain (", sc$data_domain$file, ", largest dose ", sc$data_domain$d_max,
                " Gy) differs from this calculator's, so the notes could differ.")
  # Shared input checks (hinge, clocks, parameters, DD block, anchor, calibration).
  sc2 <- sc; sc2$comparison <- sc$standard
  sc2$search <- list(criterion = "either", rank_by = "best", m_range = V4_GRID$m_range, d_step = V4_GRID$d_step,
                     d_max = V4_GRID$d_max, tol_early = V4_MARGINS$early, tol_late = V4_MARGINS$late,
                     min_gain = -V4_MARGINS$tumor, min_days_saved = 1, fenwick = FALSE)
  v <- check_inputs(sc2)
  st <- check_choice(sc$support$status, "support status", c("available", "unavailable"))
  draws <- NULL
  if (st == "available") {
    draws <- tryCatch(draws_loader(), error = function(e) replay_fail("the refit draws cannot be read: ", conditionMessage(e)))
    for (m in c("di", "dd"))
      if (!identical(as.character(sc$support$files[[m]]$sha256), attr(draws[[m]], "sha256")))
        replay_fail("its ", m, " support draws (SHA-256 ", sc$support$files[[m]]$sha256,
                    ") are not the draws of this calculator, so the support could differ.")
    if (!same(sc$support$margin, V4_MARGINS$tumor) || !same(sc$support$flag_below, V4_SUPPORT_MIN))
      replay_fail("its support margin or flag threshold differs from this calculator's.")
  }
  list(p = v$p, anchor = v$anchor, draws = draws)
}

# Exact rerun of a main export: the same points, summary and support.
main_rerun <- function(sc, ...) {
  v <- validate_main_scenario(sc, ...)
  main_compute(v$p, v$anchor, v$draws, sc$data_domain$d_max)
}
