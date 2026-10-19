# Export of a complete, rerunnable scenario (version 2.1).
# One ZIP download holds:
#   results.csv       the listed regimens (one per number of fractions), full precision
#   all_qualifying.csv every grid regimen that met all criteria
#   scenario.json     every input, the calibration, versions and a scenario ID
#   README.txt        what the files are and how to rerun them
# Requires R/bed_models.R to be sourced first.

# Scientific content of a scenario. `search` holds the suggest_regimens()
# arguments; `context` holds the interface choices that produced them.
build_scenario <- function(p, ref, anchor, search, context) {
  cal <- dd_calibration(anchor, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  c_over_alpha <- if (is.finite(cal$c_over_alpha)) cal$c_over_alpha else NULL
  list(
    schema_version = "1",
    app_version = APP_VERSION,
    equation_version = EQUATION_VERSION,
    hinge = list(type = "softplus", k_per_day = SMOOTH_K),
    clocks = list(input = "inclusive calendar days, first to last treatment day",
                  tumor = "inclusive", acute = if (is.null(p$clock_e)) "elapsed" else p$clock_e),
    params = p[c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g",
                 "ab_e", "alpha_e", "lambda_e", "Tk_e", "ab_l")],
    comparison = ref[c("m", "d", "T")],
    anchor = anchor[c("m", "d", "T")],
    calibration = list(status = cal$status, c_over_alpha = c_over_alpha,
                       onset_anchor = if (is.finite(cal$onset_ref)) cal$onset_ref else NULL),
    search = search,
    context = context)
}

# Canonical JSON (full precision, NULL kept as null).
scenario_json <- function(sc) {
  jsonlite::toJSON(sc, auto_unbox = TRUE, digits = NA, null = "null", pretty = TRUE)
}

# Deterministic ID from the scientific content (no dates, no software paths).
scenario_id <- function(sc) {
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f))
  writeLines(scenario_json(sc), f, useBytes = TRUE)
  substr(unname(tools::md5sum(f)), 1, 12)
}

# Rerun the search of a scenario read back with jsonlite::fromJSON().
rerun_scenario <- function(sc) {
  p <- sc$params
  p$clock_e <- sc$clocks$acute
  s <- sc$search
  nul <- function(x) if (is.null(x) || length(x) == 0) NULL else x
  suggest_regimens(p, sc$comparison, anchor = sc$anchor,
                   m_range = as.integer(s$m_range), d_step = s$d_step, d_max = s$d_max,
                   criterion = s$criterion, rank_by = s$rank_by,
                   tol_early = s$tol_early, tol_late = s$tol_late,
                   early_cap = nul(s$early_cap), late_cap = nul(s$late_cap),
                   min_gain = s$min_gain, min_days_saved = s$min_days_saved)
}

export_readme <- function(sc, id) {
  c(sprintf("HNSCC BED calculator %s, scenario %s", sc$app_version, id),
    "",
    "results.csv: the listed regimens, one per number of fractions (the highest-ranked dose per fraction).",
    "all_qualifying.csv: every searched regimen that met all criteria.",
    "scenario.json: every input, the DD calibration, the clocks and versions.",
    "",
    "Columns: m fractions, d Gy per fraction, T overall days (inclusive), bed_* absolute BEDs (Gy),",
    "d_* differences from the comparison regimen, d_min the smaller of d_di and d_dd,",
    "qual_di / qual_dd whether that model credits the regimen, beyond_precedent d > 3.1 Gy.",
    "Empty cells (NA) mean 'not available', never zero.",
    "",
    "To rerun in R from the app folder:",
    "  source('R/bed_models.R'); source('R/export.R')",
    "  sc <- jsonlite::fromJSON('scenario.json'); rerun_scenario(sc)",
    "",
    "Model outputs for hypothesis generation, not treatment recommendations.")
}

# Write the ZIP to `file`. `best` and `feasible` come from suggest_regimens().
write_export_zip <- function(file, sc, best, feasible) {
  id <- scenario_id(sc)
  sc$scenario_id <- id
  sc$exported_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  dir <- tempfile("bed_export_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  write.csv(best, file.path(dir, "results.csv"), row.names = FALSE, na = "NA")
  write.csv(if (is.null(feasible)) best[0, ] else feasible,
            file.path(dir, "all_qualifying.csv"), row.names = FALSE, na = "NA")
  writeLines(scenario_json(sc), file.path(dir, "scenario.json"), useBytes = TRUE)
  writeLines(export_readme(sc, id), file.path(dir, "README.txt"))
  zip::zip(file, files = c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt"),
           root = dir)
  invisible(id)
}
