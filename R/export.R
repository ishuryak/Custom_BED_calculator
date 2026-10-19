# Export of a complete, rerunnable analysis (scenario schema 2, app 2.5).
# One ZIP download holds:
#   results.csv          the listed regimens (one per number of fractions), full precision
#   all_qualifying.csv   every grid regimen that met all criteria
#   scenario_detail.csv  agreement only: every listed regimen under every sensitivity scenario
#   scenario.json        every input, the calibration, versions, code identity and an ID
#   README.txt           what the files are and how to rerun them
# Requires R/bed_models.R to be sourced first.

SCHEMA_VERSION <- "2"
SUPPORTED_SCHEMAS <- c("1", "2")          # "1": candidate files of versions 2.1-2.4
SUPPORTED_SELECTION <- c(SELECTION_VERSION)
PARAM_NAMES <- c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g",
                 "ab_e", "alpha_e", "lambda_e", "Tk_e", "ab_l")

# SHA-256 of the source files that compute the results, plus the R and
# package versions. Provenance only: it is not part of the scenario ID.
code_identity <- function(root = ".") {
  files <- c("R/bed_models.R", "R/export.R", "app.R")
  sha <- vapply(files, function(f) {
    path <- file.path(root, f)
    if (file.exists(path)) digest::digest(file = path, algo = "sha256") else NA_character_
  }, character(1))
  pk <- c("shiny", "jsonlite", "zip", "digest")
  list(files = as.list(sha), r_version = R.version.string,
       packages = as.list(vapply(pk, function(x) tryCatch(as.character(utils::packageVersion(x)),
                                                          error = function(e) NA_character_), character(1))))
}

# Scientific content of an analysis. analysis: "candidate" (the candidate
# search) or "agreement" (the two-model agreement). `search` holds the
# suggest_regimens() arguments; `agreement` holds margin and ab_values;
# `context` holds interface choices that do not change the result.
build_scenario <- function(p, ref, anchor, search, context, analysis = "candidate",
                           agreement = NULL, anchor_mode = "fixed", code_root = ".") {
  cal <- dd_calibration(anchor, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  pp <- lapply(setNames(PARAM_NAMES, PARAM_NAMES), function(n) {
    v <- p[[n]]
    if (is.null(v) && n == "g") 0 else v
  })
  sc <- list(
    schema_version = SCHEMA_VERSION,
    app_version = APP_VERSION,
    equation_version = EQUATION_VERSION,
    selection_version = SELECTION_VERSION,
    analysis = analysis,
    hinge = list(type = "softplus", k_per_day = SMOOTH_K),
    clocks = list(input = "inclusive calendar days, first to last treatment day",
                  tumor = "inclusive", acute = if (is.null(p$clock_e)) "elapsed" else p$clock_e),
    params = pp,
    comparison = ref[c("m", "d", "T")],
    anchor = c(anchor[c("m", "d", "T")], list(mode = anchor_mode)),
    calibration = list(status = cal$status,
                       c_over_alpha = if (is.finite(cal$c_over_alpha)) cal$c_over_alpha else NULL,
                       onset_anchor = if (is.finite(cal$onset_ref)) cal$onset_ref else NULL),
    search = search)
  if (analysis == "agreement") {
    abv <- if (is.null(agreement$ab_values)) SCENARIO_AB_VALUES else agreement$ab_values
    sc$agreement <- list(
      margin = agreement$margin, ab_values = abv,
      margin_meaning = "tumor BED change (Gy) at each scenario's own tumor alpha/beta",
      margin_eqd2 = as.list(setNames(agreement$margin / (1 + 2 / abv), paste0("ab_", abv))),
      selection = "one dose per fraction per number of fractions, chosen with the current parameters by the smaller of the two changes; scenarios evaluated afterwards")
    sc$scenario_set_version <- SCENARIO_SET_VERSION
    sc$scenario_set <- scenario_set(p, abv)
  }
  sc$code <- code_identity(code_root)
  sc$context <- context
  sc
}

# Canonical form for the ID: names of every named list sorted recursively,
# unnamed lists and vectors kept in order, data frames as lists of rows.
canonicalize <- function(x) {
  if (is.data.frame(x)) x <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
  if (is.list(x)) {
    if (!is.null(names(x)) && all(names(x) != "")) {
      x <- x[order(names(x), method = "radix")]
    }
    return(lapply(x, canonicalize))
  }
  x
}

# JSON (full precision, NULL kept as null).
scenario_json <- function(sc, pretty = TRUE) {
  jsonlite::toJSON(sc, auto_unbox = TRUE, digits = NA, null = "null", pretty = pretty)
}

# Deterministic ID from the fields that determine the result: versions,
# hinge, clocks, parameters, comparison, anchor, search, agreement and
# scenario set. App version, calibration (derived), code identity and
# interface context are excluded. SHA-256 of canonical UTF-8 JSON.
SCIENTIFIC_FIELDS <- c("schema_version", "equation_version", "selection_version", "analysis",
                       "hinge", "clocks", "params", "comparison", "anchor", "search",
                       "agreement", "scenario_set_version", "scenario_set")
scenario_id <- function(sc) {
  x <- sc[intersect(SCIENTIFIC_FIELDS, names(sc))]
  js <- enc2utf8(as.character(scenario_json(canonicalize(x), pretty = FALSE)))
  substr(digest::digest(js, algo = "sha256", serialize = FALSE), 1, 12)
}

# --- strict replay -----------------------------------------------------------

replay_fail <- function(...) stop("Cannot rerun this scenario: ", ..., call. = FALSE)

check_num <- function(x, name, lo = -Inf, hi = Inf, integer = FALSE) {
  if (is.null(x) || length(x) != 1 || !is.numeric(x) || !is.finite(x))
    replay_fail(name, " must be one finite number.")
  if (x < lo || x > hi) replay_fail(name, " must be between ", lo, " and ", hi, ".")
  if (integer && x != round(x)) replay_fail(name, " must be a whole number.")
  x
}

check_regimen <- function(r, what) {
  if (!is.list(r)) replay_fail(what, " is missing.")
  list(m = check_num(r$m, paste(what, "m"), 1, 100, TRUE), d = check_num(r$d, paste(what, "d"), 0.01, 30),
       T = check_num(r$T, paste(what, "T"), 1, 1000))
}

# Validate a scenario read back with jsonlite::fromJSON() and return the
# inputs to recompute it. Schema 1 (candidate files of versions 2.1-2.4) is
# migrated explicitly; anything else unsupported is refused.
validate_scenario <- function(sc) {
  schema <- as.character(sc$schema_version)
  if (length(schema) != 1 || !schema %in% SUPPORTED_SCHEMAS)
    replay_fail("unsupported schema version '", paste(schema, collapse = ","), "'.")
  migrated <- NULL
  if (schema == "1") {
    # Versions 2.1-2.4 wrote candidate searches only. Equation 2.1 equals
    # equation 2.2 without the optional Fenwick boundary, which such files cannot use.
    migrated <- "1"
    sc$analysis <- "candidate"
    sc$selection_version <- SELECTION_VERSION
    if (identical(as.character(sc$equation_version), "2.1") && !isTRUE(sc$search$fenwick))
      sc$equation_version <- EQUATION_VERSION
  }
  if (!identical(as.character(sc$equation_version), EQUATION_VERSION))
    replay_fail("equation version '", sc$equation_version, "' is not this calculator's (", EQUATION_VERSION, ").")
  if (!identical(as.character(sc$selection_version), SUPPORTED_SELECTION))
    replay_fail("selection version '", sc$selection_version, "' is not supported.")
  if (!sc$analysis %in% c("candidate", "agreement")) replay_fail("unknown analysis '", sc$analysis, "'.")
  if (!identical(sc$hinge$type, "softplus") || !isTRUE(all.equal(as.numeric(sc$hinge$k_per_day), SMOOTH_K)))
    replay_fail("hinge (", sc$hinge$type, ", k ", sc$hinge$k_per_day, ") is not the softplus hinge with k = ",
                SMOOTH_K, " that this calculator computes.")
  if (!identical(sc$clocks$tumor, "inclusive")) replay_fail("tumor clock must be 'inclusive'.")
  if (!isTRUE(sc$clocks$acute %in% c("elapsed", "inclusive"))) replay_fail("acute clock must be 'elapsed' or 'inclusive'.")

  pr <- sc$params
  if (!is.list(pr)) replay_fail("params are missing.")
  lim <- list(alpha = c(1e-6, 10), ab = c(0.01, 1000), lambda_di = c(0, 10), lambda_dd = c(0, 10),
              Tk = c(0, 1000), g = c(0, 1), ab_e = c(0.01, 1000), alpha_e = c(1e-6, 10),
              lambda_e = c(0, 100), Tk_e = c(0, 1000), ab_l = c(0.01, 1000))
  p <- lapply(PARAM_NAMES, function(n) check_num(pr[[n]], n, lim[[n]][1], lim[[n]][2]))
  names(p) <- PARAM_NAMES
  p$clock_e <- sc$clocks$acute

  ref <- check_regimen(sc$comparison, "comparison")
  anchor <- check_regimen(sc$anchor, "anchor")

  cal <- dd_calibration(anchor, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  saved <- sc$calibration
  if (!is.null(saved)) {
    if (!identical(saved$status, cal$status))
      replay_fail("saved calibration status '", saved$status, "' differs from the recomputed '", cal$status, "'.")
    if (is.finite(cal$c_over_alpha) &&
        !(length(saved$c_over_alpha) == 1 && isTRUE(abs(saved$c_over_alpha - cal$c_over_alpha) <=
                                                     1e-8 * max(1, abs(cal$c_over_alpha)))))
      replay_fail("saved calibration C/alpha does not match the one recomputed from the inputs.")
  }

  s <- sc$search
  if (!is.list(s)) replay_fail("search settings are missing.")
  if (!isTRUE(s$criterion %in% c("either", "both", "di", "dd"))) replay_fail("unknown criterion '", s$criterion, "'.")
  if (!isTRUE(s$rank_by %in% c("best", "worst", "di", "dd"))) replay_fail("unknown rank_by '", s$rank_by, "'.")
  mr <- s$m_range
  if (!is.numeric(mr) || length(mr) < 1 || any(!is.finite(mr)) || any(mr < 1) || any(mr != round(mr)))
    replay_fail("m_range must be whole numbers of at least 1.")
  opt <- function(x, name) if (is.null(x) || length(x) == 0) NULL else check_num(x, name)
  args <- list(m_range = as.integer(mr), d_step = check_num(s$d_step, "d_step", 1e-3, 30),
               d_max = check_num(s$d_max, "d_max", 0.01, 30), criterion = s$criterion,
               rank_by = s$rank_by, tol_early = check_num(s$tol_early, "tol_early", 0, 1000),
               tol_late = check_num(s$tol_late, "tol_late", 0, 1000),
               early_cap = opt(s$early_cap, "early_cap"), late_cap = opt(s$late_cap, "late_cap"),
               min_gain = check_num(s$min_gain, "min_gain", -1000, 1000),
               min_days_saved = check_num(s$min_days_saved, "min_days_saved", 0, 1000, TRUE),
               fenwick = isTRUE(s$fenwick))

  agree <- NULL
  if (sc$analysis == "agreement") {
    ag <- sc$agreement
    if (!is.list(ag)) replay_fail("agreement settings are missing.")
    if (!identical(as.character(sc$scenario_set_version), SCENARIO_SET_VERSION))
      replay_fail("scenario set version '", sc$scenario_set_version, "' is not supported.")
    abv <- as.numeric(ag$ab_values)
    if (length(abv) < 1 || any(!is.finite(abv)) || any(abv <= 0)) replay_fail("agreement ab_values are invalid.")
    agree <- list(margin = check_num(ag$margin, "agreement margin", -1000, 1000), ab_values = abv)
    now <- scenario_set(p, abv)
    old <- as.data.frame(sc$scenario_set)
    cols <- c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g")
    if (nrow(old) != nrow(now) || !identical(as.character(old$preset), now$preset) ||
        !isTRUE(all.equal(as.matrix(old[cols]), as.matrix(now[cols]), tolerance = 1e-10, check.attributes = FALSE)))
      replay_fail("the saved scenario set does not match the one this calculator builds from the inputs.")
  }
  list(analysis = sc$analysis, p = p, ref = ref, anchor = anchor, args = args, agreement = agree,
       migrated_from = migrated)
}

# Rerun the analysis of a scenario read back with jsonlite::fromJSON().
rerun_scenario <- function(sc) {
  v <- validate_scenario(sc)
  out <- if (v$analysis == "agreement") {
    args <- v$args
    args$min_gain <- NULL; args$criterion <- NULL; args$rank_by <- NULL
    do.call(agreement_regimens, c(list(p = v$p, ref = v$ref, anchor = v$anchor,
                                       margin = v$agreement$margin, ab_values = v$agreement$ab_values), args))
  } else {
    do.call(suggest_regimens, c(list(p = v$p, ref = v$ref, anchor = v$anchor), v$args))
  }
  if (!is.null(v$migrated_from)) attr(out, "migrated_from") <- v$migrated_from
  out
}

export_readme <- function(sc, id) {
  agree <- identical(sc$analysis, "agreement")
  c(sprintf("HNSCC BED calculator %s, %s analysis, scenario %s", sc$app_version,
            if (agree) "model agreement" else "candidate regimen", id),
    sprintf("Comparison: %g x %g Gy in %g days. DD anchor: %g x %g Gy in %g days (%s).",
            sc$comparison$m, sc$comparison$d, sc$comparison$T, sc$anchor$m, sc$anchor$d, sc$anchor$T,
            if (identical(sc$anchor$mode, "coupled")) "legacy: the comparison regimen" else "fixed"),
    "",
    "results.csv: the listed regimens, one per number of fractions.",
    "all_qualifying.csv: every searched regimen that met all criteria.",
    if (agree) c("scenario_detail.csv: each listed regimen under each sensitivity scenario: raw BEDs,",
                 "  validity (valid_di, valid_dd), the reason when not evaluable, differences, pass_both.",
                 sprintf("Agreement margin: %g Gy of tumor BED at each scenario's own alpha/beta (not an outcome margin).",
                         sc$agreement$margin),
                 "Scenario counts: n_valid_both evaluable scenarios, n_pass_both of them meet the margin,",
                 "  n_invalid are not evaluable (a comparison or candidate BED not positive, or DD unavailable).",
                 "  Ranges are deterministic minimum and maximum over evaluable scenarios, not confidence intervals."),
    "scenario.json: every input, the DD calibration, clocks, versions and code identity.",
    "",
    "Columns: m fractions, d Gy per fraction, T overall days (inclusive), bed_* absolute BEDs (Gy),",
    "d_* differences from the comparison regimen, d_min the smaller of the available d_di and d_dd,",
    "qual_di / qual_dd whether that model credits the regimen,",
    "beyond_precedent d above 3.1 Gy, the largest fraction size among the app's benchmark regimens.",
    "Empty cells (NA) mean 'not available', never zero. A difference against a comparison or",
    "candidate BED that is not positive is NA. If no regimen qualified, scenario.json gives the",
    "reason under result$empty_reason (no value means the limits excluded every regimen).",
    "",
    "To rerun in R from the app folder (refused if the versions, hinge or clocks differ):",
    "  source('R/bed_models.R'); source('R/export.R')",
    "  sc <- jsonlite::fromJSON('scenario.json'); rerun_scenario(sc)",
    "",
    "Model outputs for hypothesis generation, not treatment recommendations.")
}

# Write the ZIP to `file`. `best` comes from suggest_regimens() or
# agreement_regimens(); `feasible` is its attribute "feasible".
write_export_zip <- function(file, sc, best, feasible) {
  id <- scenario_id(sc)
  sc$scenario_id <- id
  sc$result <- list(rows_listed = nrow(best),
                    rows_qualifying = if (is.null(feasible)) nrow(best) else nrow(feasible),
                    empty_reason = attr(best, "reason"))
  sc$exported_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  dir <- tempfile("bed_export_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- c("results.csv", "all_qualifying.csv", "scenario.json", "README.txt")
  write.csv(best, file.path(dir, "results.csv"), row.names = FALSE, na = "NA")
  write.csv(if (is.null(feasible)) best[0, ] else feasible,
            file.path(dir, "all_qualifying.csv"), row.names = FALSE, na = "NA")
  if (identical(sc$analysis, "agreement")) {
    det <- attr(best, "scenario_detail")
    write.csv(if (is.null(det)) data.frame() else det, file.path(dir, "scenario_detail.csv"),
              row.names = FALSE, na = "NA")
    files <- c(files, "scenario_detail.csv")
  }
  writeLines(enc2utf8(as.character(scenario_json(sc))), file.path(dir, "scenario.json"), useBytes = TRUE)
  writeLines(export_readme(sc, id), file.path(dir, "README.txt"))
  zip::zip(file, files = files, root = dir)
  invisible(id)
}
