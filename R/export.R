# Export of a complete, rerunnable analysis (scenario schema 3, app 3.0).
# One ZIP download holds:
#   results.csv          the listed regimens (one per number of fractions), full precision
#   all_qualifying.csv   every grid regimen that met all criteria
#   scenario_detail.csv  agreement only: every listed regimen under every sensitivity scenario
#   scenario.json        every input, the calibration, versions, code identity and an ID
#   README.txt           what the files are and how to rerun them
# Requires R/bed_models.R to be sourced first.

SCHEMA_VERSION <- "3"
# Files that this engine replays exactly: schema 3 of this version, and
# schema 2 of versions 2.6-2.7 (selection 2.6, equations 2.2), whose
# calibrated-DD equations version 3.0 computes unchanged.
REPLAYABLE <- list(c(schema = "3", selection = "2.6", equation = "3.0"),
                   c(schema = "2", selection = "2.6", equation = "2.2"))
# Exact replay needs the engine that wrote the file. This calculator keeps
# only its own engine, so rerun_scenario() replays files of the current
# schema, equation and selection versions and refuses others.
# migrate_scenario() turns an older file (schema 1 of versions 2.1-2.4, or
# schema 2 of version 2.5) into a new scenario under the current rules,
# recording where it came from; its results may differ from the original.
MIGRATABLE_SCHEMAS <- c("1", "2", "3")
MIGRATABLE_SELECTION <- c("2.5", "2.6")   # schema 1 files carry none
# Practical bound on the search grid of a replayed file (rows evaluated).
MAX_GRID_ROWS <- 200000
PARAM_NAMES <- c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g",
                 "ab_e", "alpha_e", "lambda_e", "Tk_e", "ab_l")

# SHA-256 of the source files that compute the results, plus the R and
# package versions. Provenance only: it is not part of the scenario ID.
# The bytes are read with readBin(): digest(file = ...) first calls
# file.access(), which on some Windows and OneDrive setups reports a readable
# file as unreadable. The SHA-256 is the same either way.
# Version 3 hashes every R file, the app, the About page and the data and
# draws files, so a new module cannot be invisible to the record (plan E).
code_identity_files <- function(root = ".") {
  f <- c(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = FALSE),
         list.files(file.path(root, "data"), pattern = "\\.(csv|json)$", full.names = FALSE),
         list.files(file.path(root, "draws"), pattern = "\\.(csv|json)$", full.names = FALSE))
  dirs <- c(rep("R", length(list.files(file.path(root, "R"), pattern = "\\.R$"))),
            rep("data", length(list.files(file.path(root, "data"), pattern = "\\.(csv|json)$"))),
            rep("draws", length(list.files(file.path(root, "draws"), pattern = "\\.(csv|json)$"))))
  c(file.path(dirs, f), "app.R", "about.md")
}
code_identity <- function(root = ".") {
  files <- code_identity_files(root)
  sha <- vapply(files, function(f) {
    path <- file.path(root, f)
    if (!file.exists(path)) return(NA_character_)
    bytes <- readBin(path, "raw", n = file.size(path))
    digest::digest(bytes, algo = "sha256", serialize = FALSE)
  }, character(1))
  pk <- c("shiny", "jsonlite", "zip", "digest")
  list(files = as.list(sha), r_version = R.version.string,
       packages = as.list(vapply(pk, function(x) tryCatch(as.character(utils::packageVersion(x)),
                                                          error = function(e) NA_character_), character(1))))
}

# Scientific content of an analysis. analysis: "candidate" (the candidate
# search) or "agreement" (the two-model agreement). `search` holds the
# suggest_regimens() arguments; `agreement` holds margin and ab_values;
# `context` holds interface choices that do not change the result. `code` is
# the code_identity() taken when the app loaded its R code; if the files on
# disk have changed since, code$changed_since_load is TRUE.
build_scenario <- function(p, ref, anchor, search, context, analysis = "candidate",
                           agreement = NULL, anchor_mode = "fixed", code_root = ".", code = NULL) {
  cal <- dd_setup(p, anchor)
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
    dd = dd_block(p),
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
  now <- code_identity(code_root)
  if (is.null(code)) code <- now
  code$changed_since_load <- !identical(code$files, now$files)
  sc$code <- code
  sc$context <- context
  sc
}

# The DD model's settings: calibrated (default) or fitted with its own values.
dd_block <- function(p) {
  if (!identical(p$dd_mode, "fitted")) return(list(mode = "calibrated"))
  list(mode = "fitted", alpha_dd = p$alpha_dd, c_over_alpha_dd = p$c_over_alpha_dd,
       g_dd = if (is.null(p$g_dd)) p$g else p$g_dd)
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
                       "hinge", "clocks", "params", "dd", "comparison", "anchor", "search",
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

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) replay_fail(name, " must be true or false.")
  x
}
check_choice <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1 || !x %in% choices)
    replay_fail(name, " must be one of: ", paste(choices, collapse = ", "), ".")
  x
}
same_num <- function(a, b, tol = 1e-8) length(a) == 1 && is.numeric(a) && isTRUE(abs(a - b) <= tol * max(1, abs(b)))

# Checks shared by replay and migration: hinge, clocks, parameters, regimens,
# calibration and search settings. Returns the inputs to recompute.
check_inputs <- function(sc) {
  if (!identical(sc$hinge$type, "softplus") || !isTRUE(all.equal(as.numeric(sc$hinge$k_per_day), SMOOTH_K)))
    replay_fail("hinge (", sc$hinge$type, ", k ", sc$hinge$k_per_day, ") is not the softplus hinge with k = ",
                SMOOTH_K, " that this calculator computes.")
  if (!identical(sc$clocks$tumor, "inclusive")) replay_fail("tumor clock must be 'inclusive'.")
  if (!isTRUE(sc$clocks$acute %in% c("elapsed", "inclusive")) || length(sc$clocks$acute) != 1)
    replay_fail("acute clock must be 'elapsed' or 'inclusive'.")

  pr <- sc$params
  if (!is.list(pr)) replay_fail("params are missing.")
  lim <- list(alpha = c(1e-6, 10), ab = c(0.01, 1000), lambda_di = c(0, 10), lambda_dd = c(0, 10),
              Tk = c(0, 1000), g = c(0, 1), ab_e = c(0.01, 1000), alpha_e = c(1e-6, 10),
              lambda_e = c(0, 100), Tk_e = c(0, 1000), ab_l = c(0.01, 1000))
  p <- lapply(PARAM_NAMES, function(n) check_num(pr[[n]], n, lim[[n]][1], lim[[n]][2]))
  names(p) <- PARAM_NAMES
  p$clock_e <- sc$clocks$acute
  ddb <- sc$dd
  if (!is.null(ddb)) {
    mode <- check_choice(ddb$mode, "dd mode", c("calibrated", "fitted"))
    if (mode == "fitted") {
      p$dd_mode <- "fitted"
      p$alpha_dd <- check_num(ddb$alpha_dd, "alpha_dd", 1e-6, 10)
      p$c_over_alpha_dd <- check_num(ddb$c_over_alpha_dd, "c_over_alpha_dd", 1e-6, 1e6)
      p$g_dd <- check_num(ddb$g_dd, "g_dd", 0, 1)
    }
  }

  ref <- check_regimen(sc$comparison, "comparison")
  anchor <- check_regimen(sc$anchor, "anchor")

  # Every saved calibration field must describe the inputs.
  cal <- dd_setup(p, anchor)
  saved <- sc$calibration
  if (!is.null(saved)) {
    if (!identical(saved$status, cal$status))
      replay_fail("saved calibration status '", saved$status, "' differs from the recomputed '", cal$status, "'.")
    if (is.finite(cal$c_over_alpha) && !same_num(saved$c_over_alpha, cal$c_over_alpha))
      replay_fail("saved calibration C/alpha does not match the one recomputed from the inputs.")
    if (is.finite(cal$onset_ref) && !is.null(saved$onset_anchor) && !same_num(saved$onset_anchor, cal$onset_ref))
      replay_fail("saved calibration onset_anchor does not match the one recomputed from the inputs.")
  }

  s <- sc$search
  if (!is.list(s)) replay_fail("search settings are missing.")
  criterion <- check_choice(s$criterion, "criterion", c("either", "both", "di", "dd"))
  rank_by <- check_choice(s$rank_by, "rank_by", c("best", "worst", "di", "dd"))
  mr <- s$m_range
  if (!is.numeric(mr) || length(mr) < 1 || length(mr) > 100 || any(!is.finite(mr)) || any(mr < 1) ||
      any(mr > 100) || any(mr != round(mr)) || anyDuplicated(mr))
    replay_fail("m_range must be 1 to 100 distinct whole numbers between 1 and 100.")
  opt <- function(x, name) if (is.null(x) || length(x) == 0) NULL else check_num(x, name, 0, 1000)
  args <- list(m_range = as.integer(mr), d_step = check_num(s$d_step, "d_step", 1e-3, 30),
               d_max = check_num(s$d_max, "d_max", 0.01, 30), criterion = criterion,
               rank_by = rank_by, tol_early = check_num(s$tol_early, "tol_early", 0, 1000),
               tol_late = check_num(s$tol_late, "tol_late", 0, 1000),
               early_cap = opt(s$early_cap, "early_cap"), late_cap = opt(s$late_cap, "late_cap"),
               min_gain = check_num(s$min_gain, "min_gain", -1000, 1000),
               min_days_saved = check_num(s$min_days_saved, "min_days_saved", 0, 1000, TRUE),
               fenwick = check_flag(s$fenwick, "fenwick"))
  if (!is.null(s$calendar)) {
    fpw <- check_num(s$calendar$fx_per_week, "calendar fx_per_week", 5, 6, TRUE)
    args$calendar <- list(fx_per_week = fpw,
                          start_weekday = check_num(s$calendar$start_weekday, "calendar start_weekday", 1, fpw, TRUE))
  }
  n_d <- max(0, floor((args$d_max - ref$d) / args$d_step) + 1)
  if (n_d * length(mr) > MAX_GRID_ROWS)
    replay_fail("the search grid (", n_d * length(mr), " regimens) exceeds ", MAX_GRID_ROWS, ".")
  list(p = p, ref = ref, anchor = anchor, args = args)
}

# Agreement settings: the margin, the scenario set, and the search fields the
# agreement fixes (criterion "both", rank_by "worst", min_gain = margin).
check_agreement <- function(sc, p, args) {
  ag <- sc$agreement
  if (!is.list(ag)) replay_fail("agreement settings are missing.")
  if (!identical(as.character(sc$scenario_set_version), SCENARIO_SET_VERSION))
    replay_fail("scenario set version '", sc$scenario_set_version, "' is not supported.")
  abv <- ag$ab_values
  if (!is.numeric(abv) || length(abv) < 1 || any(!is.finite(abv)) || any(abv <= 0))
    replay_fail("agreement ab_values are invalid.")
  margin <- check_num(ag$margin, "agreement margin", -1000, 1000)
  if (!identical(args$criterion, "both") || !identical(args$rank_by, "worst") ||
      !isTRUE(abs(args$min_gain - margin) <= 1e-12))
    replay_fail("an agreement file must have search criterion 'both', rank_by 'worst' and min_gain ",
                "equal to the agreement margin; this one has '", args$criterion, "', '", args$rank_by,
                "' and ", args$min_gain, " (margin ", margin, ").")
  if (!is.null(ag$margin_eqd2)) {
    want <- margin / (1 + 2 / abv)
    got <- unlist(ag$margin_eqd2)
    if (length(got) != length(want) || !isTRUE(all.equal(unname(got), want, tolerance = 1e-10)))
      replay_fail("agreement margin_eqd2 does not match the margin and ab_values.")
  }
  now <- scenario_set(p, abv)
  old <- as.data.frame(sc$scenario_set)
  cols <- c("alpha", "ab", "lambda_di", "lambda_dd", "Tk", "g")
  if (nrow(old) != nrow(now) || !all(cols %in% names(old)) || !identical(as.character(old$preset), now$preset) ||
      !isTRUE(all.equal(as.matrix(old[cols]), as.matrix(now[cols]), tolerance = 1e-10, check.attributes = FALSE)))
    replay_fail("the saved scenario set does not match the one this calculator builds from the inputs.")
  list(margin = margin, ab_values = abv)
}

# Strict validation of a scenario read back with jsonlite::fromJSON(). Only
# files of this calculator's schema, equation and selection versions are
# accepted; an older file must go through migrate_scenario() first.
validate_scenario <- function(sc) {
  schema <- as.character(sc$schema_version)
  if (length(schema) != 1 || !schema %in% MIGRATABLE_SCHEMAS)
    replay_fail("unsupported schema version '", paste(schema, collapse = ","), "'.")
  old_sel <- if (is.null(sc$selection_version)) "none (schema 1)" else as.character(sc$selection_version)
  key <- c(schema = schema, selection = old_sel, equation = as.character(sc$equation_version)[1])
  replayable <- any(vapply(REPLAYABLE, function(r) identical(unname(r), unname(key)), logical(1)))
  if (!replayable && schema == SCHEMA_VERSION && identical(old_sel, SELECTION_VERSION))
    replay_fail("equation version '", sc$equation_version, "' is not this calculator's (", EQUATION_VERSION, ").")
  if (!replayable) {
    if (schema == SCHEMA_VERSION && !old_sel %in% MIGRATABLE_SELECTION)
      replay_fail("selection version '", old_sel, "' is not supported.")
    replay_fail("this file was written by version ", sc$app_version, " (schema ", schema, ", selection ",
                old_sel, "). Exact replay needs that version's engine, which this calculator (selection ",
                SELECTION_VERSION, ") does not include. To reanalyze it with the current rules, run ",
                "rerun_scenario(migrate_scenario(sc)); its results may differ from the original export.")
  }
  analysis <- check_choice(sc$analysis, "analysis", c("candidate", "agreement"))
  v <- check_inputs(sc)
  agree <- if (analysis == "agreement") check_agreement(sc, v$p, v$args) else NULL
  list(analysis = analysis, p = v$p, ref = v$ref, anchor = v$anchor, args = v$args, agreement = agree,
       migration = sc$migration)
}

# Turn an older scenario file into a new scenario under the current rules.
# The original versions and ID are recorded under `migration`; the new file
# gets its own ID. Equation 2.1 (versions 2.1) equals 2.2 without the optional
# Fenwick boundary, so only non-Fenwick 2.1 files can be migrated. Schema 1
# files (versions 2.1-2.4) hold candidate searches; a missing Fenwick flag in
# them means "not used" and is recorded in migration$defaults_filled.
migrate_scenario <- function(sc) {
  schema <- as.character(sc$schema_version)
  if (length(schema) != 1 || !schema %in% MIGRATABLE_SCHEMAS)
    replay_fail("unsupported schema version '", paste(schema, collapse = ","), "'.")
  sel <- if (is.null(sc$selection_version)) NA_character_ else as.character(sc$selection_version)
  if (schema == SCHEMA_VERSION && identical(sel, SELECTION_VERSION))
    replay_fail("this file is already of the current versions; rerun it with rerun_scenario().")
  if (schema == "3") replay_fail("schema 3 files of other selection versions cannot be migrated.")
  if (schema == "2" && !sel %in% MIGRATABLE_SELECTION)
    replay_fail("selection version '", sel, "' cannot be migrated.")
  eq <- as.character(sc$equation_version)
  filled <- character()
  if (schema == "1") {
    if (is.null(sc$search$fenwick)) { sc$search$fenwick <- FALSE; filled <- c(filled, "search$fenwick = false") }
    sc$analysis <- "candidate"
  }
  if (identical(eq, "2.1")) {
    if (!identical(sc$search$fenwick, FALSE))
      replay_fail("equation 2.1 has no Fenwick boundary; this file cannot be migrated.")
  } else if (!eq %in% c("2.2", EQUATION_VERSION)) {
    replay_fail("equation version '", eq, "' cannot be migrated.")
  }
  old_id <- if (!is.null(sc$scenario_id)) as.character(sc$scenario_id) else NA_character_
  from <- list(schema_version = schema, app_version = sc$app_version, equation_version = eq,
               selection_version = if (is.na(sel)) NULL else sel, scenario_id = old_id,
               result = sc$result)
  sc$schema_version <- SCHEMA_VERSION
  sc$equation_version <- EQUATION_VERSION
  sc$selection_version <- SELECTION_VERSION
  sc$app_version <- APP_VERSION
  sc$scenario_id <- NULL; sc$result <- NULL; sc$exported_at <- NULL
  if (identical(sc$analysis, "agreement") && is.null(sc$scenario_set_version)) sc$scenario_set_version <- SCENARIO_SET_VERSION
  sc$migration <- list(from = from,
                       to = list(schema_version = SCHEMA_VERSION, app_version = APP_VERSION,
                                 equation_version = EQUATION_VERSION, selection_version = SELECTION_VERSION),
                       defaults_filled = filled,
                       note = paste("Reanalysis of an older scenario with the current selection rules.",
                                    "This is not an exact replay; results may differ from the original export."))
  validate_scenario(sc)  # the migrated file must pass the strict checks
  sc$scenario_id <- scenario_id(sc)
  sc
}

# Rerun the analysis of a scenario read back with jsonlite::fromJSON(). Strict:
# files of other versions are refused (see migrate_scenario()).
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
  if (!is.null(v$migration)) attr(out, "migration") <- v$migration
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
                 "  n_valid_di and n_valid_dd give each model's own count. Both ranges are the minimum and maximum",
                 "  over the n_valid_both evaluable scenarios, not confidence intervals.",
                 "scenario_detail.csv covers the listed regimens only, not every qualifying regimen."),
    "scenario.json: every input, the DD calibration, clocks, versions and code identity.",
    "",
    "Columns: m fractions, d Gy per fraction, T overall days (inclusive), bed_* absolute BEDs (Gy),",
    "d_* differences from the comparison regimen, d_min the smaller of d_di and d_dd (NA unless both),",
    "score the ranking value; all_qualifying.csv has rank_eligible (FALSE: the score needs a model",
    "  whose difference is not evaluable, so the row qualified but could not be ranked),",
    "fractions_saved, early_slack and late_slack (limit minus the regimen's BED; for the Fenwick",
    "  option, the boundary minus the cell-kill BED),",
    "qual_di / qual_dd whether that model credits the regimen,",
    "beyond_precedent d above 3.1 Gy, the largest fraction size among the app's benchmark regimens.",
    "Empty cells (NA) mean 'not available', never zero. A difference against a comparison or",
    "candidate BED that is not positive is NA. If no regimen qualified, scenario.json gives the",
    "reason under result$empty_reason (no value means the limits excluded every regimen).",
    "",
    "To rerun in R from the app folder (refused if the versions, hinge or clocks differ):",
    "  source('R/bed_models.R'); source('R/export.R')",
    "  sc <- jsonlite::fromJSON('scenario.json', simplifyVector = TRUE); rerun_scenario(sc)",
    "A file from an older version is refused by rerun_scenario(). rerun_scenario(migrate_scenario(sc))",
    "reanalyzes it with the current rules and records the original versions; results may differ.",
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
