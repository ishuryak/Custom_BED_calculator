# Organ-at-risk (OAR) module (version 3).
#
# Constraints stated for one fractionation are converted to another with the
# linear-quadratic isoeffect and an organ alpha/beta, as the HYPNO trial
# describes (Bentzen et al. 2026: alpha/beta 3 Gy, 2 Gy for spinal cord and
# brainstem; for mean doses, ideally bin by bin over the dose-volume
# histogram). Every constraint in data/oar_constraints.csv carries its source,
# page and a verbatim quote. Nothing here is a statement that a plan is safe:
# the output is "within" or "exceeds the converted constraint under the
# stated alpha/beta".
#
# Version 3.1 (2026-10-25 review, OAR-01 to OAR-05):
# - A mean-dose limit is a physical mean from its source. It is compared with
#   the planned physical mean, with the source fractionation stated. The mean
#   EQD2 is shown for illustration only and never decides a status.
# - Maximum-dose limits keep the LQ conversion (the HYPNO method); the cord's
#   2-Gy-per-fraction anchor stays an EQD2 comparison.
# - Each limit has its source operator: HYHOPE's limits are strict ("<").
# - HYHOPE's mandatory PRVs (cord + 5 mm, brainstem + 5 mm) are organs of
#   their own, and unassessed mandatory metrics are listed.
# - A DVH maximum is the interval of the highest occupied bin, never a bin
#   midpoint; a bin that straddles the limit is "indeterminate".
# - An entered organ without a verified constraint is never a pass.

OAR_VERSION <- "3.1"
NEAR_FRACTION <- 0.05          # "within 5% of" the converted limit
LQ_POINT_DOSE_FLAG_GY <- 5     # dose per fraction above which a point-dose conversion is flagged

load_constraints <- function(path = file.path("data", "oar_constraints.csv")) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8")
  stopifnot(all(c("id", "organ", "metric", "limit_gy", "limit_lo_gy", "n_fractions", "level",
                  "ab_class", "library", "source_short", "page", "quote", "operator", "scope_note") %in% names(x)))
  x$scope_note[is.na(x$scope_note)] <- ""
  stopifnot(all(x$operator %in% c("<", "\u2264")))
  x
}

# Version and clinical scope of each library, shown with the results.
OAR_LIBRARY_META <- list(
  QUANTEC = list(version = "QUANTEC 2010 (Int J Radiat Oncol Biol Phys 76(3) Suppl)",
                 scope = paste("Conventional fractionation (about 1.8-2 Gy per fraction to the target).",
                               "Mean-dose limits are physical means from those courses.")),
  HYHOPE = list(version = "HYHOPE protocol v4.2 (30 Jan 2024), sec. 4.1.4.2; body headers say v4.1",
                scope = paste("Tables for 12, 15 and 20 fractions. The 12-fraction level was removed by",
                              "amendment and is kept as historical modeling context.")))

# Organs offered in the app, each tied to the constraint used by default in
# each library. QUANTEC values are stated in 2-Gy fractions (n_fractions 0).
# HYHOPE values are stated for 12, 15 and 20 fractions. Organs whose
# conventional constraint was not found in a verified source have no QUANTEC
# entry (NA) and can be assessed with the HYHOPE library only.
OAR_ORGANS <- data.frame(
  key = c("cord", "cord_prv", "brainstem", "brainstem_prv", "mandible", "parotid", "submandibular",
          "constrictors", "larynx", "oral_cavity", "esophagus", "cochlea", "plexus"),
  label = c("Spinal cord", "Spinal cord + 5 mm", "Brainstem", "Brainstem + 5 mm", "Mandible",
            "Parotid (spared side)", "Submandibular gland (spared side)",
            "Pharyngeal constrictors", "Larynx", "Oral cavity (outside target)", "Cervical esophagus",
            "Cochlea", "Brachial plexus"),
  metric = c("Dmax", "Dmax", "Dmax", "Dmax", "Dmax", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmax"),
  ab_class = c("cord", "cord", "cord", "cord", "oar", "oar", "oar", "oar", "oar", "oar", "oar", "oar", "oar"),
  quantec_id = c("OAR079", NA, "OAR089", NA, NA, "OAR086", "OAR088", NA, "OAR094", NA, NA, "OAR099", NA),
  hyhope_organ = c("Spinal cord", "Spinal cord + 5 mm", "Brainstem", "Brainstem + 5 mm", "Mandible",
                   "Contralateral parotid", "Contralateral submandibular gland", "Pharyngeal constrictors", "Larynx",
                   "Oral cavity (excluding PTV)", "Cervical esophagus", "Cochlea", "Brachial plexus"),
  stringsAsFactors = FALSE)

# Default organ alpha/beta (Gy) by class, as used by HYPNO.
default_oar_ab <- function() c(cord = 2, oar = 3)

# --- LQ conversion -----------------------------------------------------------

# EQD2 of a uniform physical dose D given in n fractions.
eqd2_physical <- function(D, n, ab) D * (D / n + ab) / (2 + ab)

# Physical dose in n fractions with the given EQD2 (inverse of eqd2_physical).
physical_from_eqd2 <- function(E, n, ab) n / 2 * (-ab + sqrt(ab^2 + 4 * E * (2 + ab) / n))

# EQD2 of a constraint row: limits stated in 2-Gy fractions (n_fractions 0)
# are EQD2 already; others are converted from their own fractionation.
constraint_eqd2 <- function(limit_gy, n_fractions, ab) {
  ifelse(n_fractions == 0, limit_gy, eqd2_physical(limit_gy, n_fractions, ab))
}

# Physical limit at m fractions of a constraint (limit_gy at n_fractions).
convert_limit <- function(limit_gy, n_fractions, m, ab) {
  physical_from_eqd2(constraint_eqd2(limit_gy, n_fractions, ab), m, ab)
}

# Does a planned value meet a limit under the source operator? Strict "<"
# fails at equality; a 1e-9 Gy tolerance absorbs rounding only.
limit_met <- function(planned, limit, op = "\u2264") {
  if (identical(op, "<")) planned < limit - 1e-9 else planned <= limit + 1e-9
}

# The one organ check shared by the table and both searches. row: a constraint
# row; metric: "Dmax" or "Dmean"; phys: planned physical dose in m fractions
# (vector over draws or regimens); ab: organ alpha/beta (same length or 1).
# A mean-dose limit is a physical mean: compared on physical dose, whatever
# alpha/beta. A maximum is compared as EQD2 (LQ conversion).
organ_meets <- function(row, metric, phys, m, ab) {
  if (metric == "Dmean") return(limit_met(phys, row$limit_gy, row$operator))
  limit_met(eqd2_physical(phys, m, ab), constraint_eqd2(row$limit_gy, row$n_fractions, ab), row$operator)
}

no_constraint_error <- function(organ_key, library) {
  stop(sprintf("%s has no verified constraint in the %s library. Clear its dose or choose another library.",
               OAR_ORGANS$label[OAR_ORGANS$key == organ_key], library), call. = FALSE)
}

# The constraint row an organ uses for a library at m fractions. HYHOPE:
# the per-protocol (mandatory) row or the recommended row, from the table
# with the nearest number of fractions (12, 15 or 20; ties go to the larger).
pick_constraint <- function(lib, organ_key, library = c("QUANTEC", "HYHOPE"), m = 20, metric = NULL) {
  library <- match.arg(library)
  o <- OAR_ORGANS[OAR_ORGANS$key == organ_key, ]
  if (nrow(o) != 1) stop("unknown organ: ", organ_key)
  met <- if (is.null(metric)) o$metric else metric
  if (library == "QUANTEC") {
    if (is.na(o$quantec_id)) return(NULL)
    return(lib[lib$id == o$quantec_id, , drop = FALSE])
  }
  rows <- lib[lib$library == "HYHOPE" & lib$organ == o$hyhope_organ & lib$metric == met &
                lib$level %in% c("per_protocol", "recommended"), , drop = FALSE]
  if (nrow(rows) == 0) return(NULL)
  ns <- sort(unique(rows$n_fractions))
  n_use <- ns[which.min(abs(ns - m) - 1e-6 * ns)]
  rows[rows$n_fractions == n_use, , drop = FALSE][1, , drop = FALSE]
}

# --- dose-volume histograms --------------------------------------------------

DVH_MAX_DOSE_GY <- 300   # a larger dose means a unit error (cGy read as Gy)

# Shared validator of a cumulative DVH, before any arithmetic (files and the
# public functions). dose: bin edges in Gy, strictly increasing; volume:
# cumulative volume at or above each dose, non-increasing, ending at zero (a
# closed tail). unit "%" also requires the curve to start at 100%.
validate_dvh <- function(dose, volume, unit = NA, name = "The DVH") {
  fail <- function(...) stop(name, ": ", ..., call. = FALSE)
  if (!is.numeric(dose) || !is.numeric(volume) || length(dose) != length(volume))
    fail("dose and volume must be numbers of the same length.")
  if (length(dose) < 2) fail("at least two rows are needed.")
  if (any(!is.finite(dose)) || any(!is.finite(volume))) fail("every dose and volume must be a finite number.")
  if (any(dose < 0) || any(volume < 0)) fail("doses and volumes must not be negative.")
  if (any(diff(dose) <= 0)) fail("doses must increase from row to row (no repeated or unordered doses).")
  if (any(diff(volume) > 1e-9 * max(1, max(volume)))) stop("A cumulative DVH must not increase with dose.", call. = FALSE)
  if (volume[1] <= 0) fail("the total volume must be positive.")
  if (volume[length(volume)] > 1e-6 * volume[1])
    fail(sprintf("the curve ends at %g (not zero volume), so its maximum dose is unknown. Export the full curve.",
                 volume[length(volume)]))
  if (identical(unit, "%") && abs(volume[1] - 100) > 0.5)
    fail(sprintf("a %% DVH must start at 100%% of the structure; this one starts at %g%%.", volume[1]))
  if (max(dose) > DVH_MAX_DOSE_GY) fail(sprintf("doses above %g Gy: check the dose unit (Gy or cGy).", DVH_MAX_DOSE_GY))
  invisible(TRUE)
}

# Differential DVH from a cumulative one (validated first). Returns, per
# occupied bin, its lower and upper dose edges, its midpoint (for mean-dose
# quadrature only) and its volume.
cumulative_to_differential <- function(dose, volume) {
  validate_dvh(dose, volume)
  k <- length(dose)
  v <- volume[-k] - volume[-1]
  keep <- v > 0
  data.frame(dose = ((dose[-k] + dose[-1]) / 2)[keep], volume = v[keep],
             lo = dose[-k][keep], hi = dose[-1][keep])
}

# Mean physical dose and mean EQD2 (bin by bin, midpoint quadrature) of a
# differential DVH whose doses were delivered in n fractions, and the
# maximum as the interval [max_lo, max_hi] of the highest occupied bin: the
# histogram does not locate the maximum inside that bin.
dvh_means <- function(diff_dvh, n, ab) {
  w <- diff_dvh$volume / sum(diff_dvh$volume)
  top <- which.max(diff_dvh$hi)
  list(mean_dose = sum(w * diff_dvh$dose),
       mean_eqd2 = sum(w * eqd2_physical(diff_dvh$dose, n, ab)),
       max_lo = diff_dvh$lo[top], max_hi = diff_dvh$hi[top])
}

# --- DVH files ---------------------------------------------------------------

DVH_MAX_BYTES <- 2 * 1024^2
ID_PATTERN <- "(patient\\s*name|patient\\s*id|\\bmrn\\b|medical\\s*record|date\\s*of\\s*birth|\\bdob\\b|birth\\s*date)\\s*[:=,;\\t]\\s*\\S"

# Refuse text that carries a patient identifier with a value. Returns the
# offending line numbers (never the content) or integer(0).
identifier_lines <- function(lines) which(grepl(ID_PATTERN, lines, ignore.case = TRUE, perl = TRUE))

# CSV columns whose header names an identifier, whatever its quoting, byte
# order mark, case, spacing or separators ("Patient ID", " patient_id ", MRN).
identifier_columns <- function(header) {
  n <- tolower(gsub("[^A-Za-z0-9]", "", gsub("\ufeff", "", header, fixed = TRUE)))
  which(grepl("patientname|patientid|medicalrecord|dateofbirth|birthdate", n) |
          n %in% c("mrn", "dob", "id", "name", "patient", "birthday"))
}

dvh_unit_scale <- function(header) {
  if (grepl("cGy", header, fixed = TRUE)) 0.01 else if (grepl("\\bGy\\b", header, perl = TRUE)) 1 else NA_real_
}

# Parse an Eclipse-style text DVH export (one or more structures, each with a
# "Structure:" line and a table whose header contains "Dose [Gy]" or
# "Dose [cGy]" and a volume column in [%] or [cm³]/[cc]). Written to the
# layout of Eclipse text exports and tested on synthetic files; check the
# first real export against the planning system's own mean and maximum.
parse_dvh_eclipse <- function(lines) {
  starts <- grep("^\\s*Structure\\s*:", lines)
  if (length(starts) == 0) stop("No 'Structure:' line found: this does not look like an Eclipse DVH text export.")
  ends <- c(starts[-1] - 1, length(lines))
  out <- list(); first_line <- integer()
  for (k in seq_along(starts)) {
    blk <- lines[starts[k]:ends[k]]
    name <- gsub("\\s+", " ", trimws(sub("^\\s*Structure\\s*:", "", blk[1])))
    # The table header names both columns; metadata lines such as "Max Dose [cGy]: ..." do not.
    h <- grep("Dose \\[c?Gy\\].*Volume \\[", blk)
    if (length(h) == 0) next
    # A blank or repeated structure name is refused before anything is stored (review F01 of
    # 2026-10-26): a second block of the same name used to replace the first silently.
    if (!nzchar(name)) stop("Line ", starts[k], ": the structure has no name.")
    if (name %in% names(out))
      stop("Structure '", name, "' appears more than once (lines ", first_line[[name]], " and ", starts[k],
           "). Export each structure once, or rename one of them.")
    first_line[[name]] <- starts[k]
    hdr <- blk[h[1]]
    # Columns are the header's bracketed units, in order: "Dose [Gy]", optional "Relative dose [%]",
    # "... Volume [%]". Every data row must have exactly that many fields (review F01).
    cols <- regmatches(hdr, gregexpr("[^\\[\\]]*\\[[^]]*\\]", hdr, perl = TRUE))[[1]]
    n_fields <- length(cols)
    dcol <- grep("(^|\\s)Dose \\[c?Gy\\]", cols)[1]
    vcol <- grep("Volume \\[", cols)[1]
    scale <- dvh_unit_scale(regmatches(hdr, regexpr("Dose \\[c?Gy\\]", hdr)))
    vol_unit <- if (grepl("Volume \\[%\\]", hdr)) "%" else if (grepl("Volume \\[(cm\u00b3|cm3|cc)\\]", hdr)) "cc" else NA
    if (is.na(vol_unit)) stop("Structure '", name, "': volume column must be in % or cm\u00b3.")
    # Every non-blank line after the column header is a data row (review R10): a malformed or
    # non-finite row is refused with its line number, never dropped.
    idx <- if (h[1] < length(blk)) (h[1] + 1):length(blk) else integer()
    rows <- blk[idx]; keep <- nzchar(trimws(rows)); rows <- rows[keep]; line_no <- starts[k] - 1 + idx[keep]
    for (i in seq_along(rows)) {
      tok <- strsplit(trimws(rows[i]), "\\s+")[[1]]
      v <- suppressWarnings(as.numeric(tok))
      if (length(tok) < 2 || any(is.na(v)) || any(!is.finite(v)))
        stop("Structure '", name, "', line ", line_no[i], ": not a row of finite numbers ('", trimws(rows[i]), "').")
      if (length(tok) != n_fields)
        stop("Structure '", name, "', line ", line_no[i], ": ", length(tok), " fields where the column header has ",
             n_fields, " fields ('", trimws(rows[i]), "').")
    }
    if (length(rows) < 2) stop("Structure '", name, "': fewer than two DVH rows.")
    num <- matrix(as.numeric(unlist(strsplit(trimws(rows), "\\s+"))), ncol = n_fields, byrow = TRUE)
    out[[name]] <- data.frame(dose = num[, dcol] * scale, volume = num[, vcol])
    validate_dvh(out[[name]]$dose, out[[name]]$volume, vol_unit, sprintf("Structure '%s'", name))
    attr(out[[name]], "volume_unit") <- vol_unit
  }
  if (length(out) == 0) stop("No DVH table found.")
  out
}

# Parse a two-column CSV (dose, cumulative volume), optionally with a
# structure column. Units come from the header: "cGy" or "Gy"; "%" or
# "cc"/"cm3". A generic format for exports from any planning system.
parse_dvh_csv <- function(lines) {
  # Every line must have the header's field count: read.csv would otherwise pad a short row with NA or
  # wrap a long one into a new row (review F01 of 2026-10-26).
  body <- lines[nzchar(trimws(lines))]
  nf <- utils::count.fields(textConnection(body), sep = ",", quote = "\"", blank.lines.skip = FALSE)
  bad <- which(is.na(nf) | nf != nf[1])
  if (length(bad))
    stop(sprintf("Line %d has %s fields where the header has %d.", which(nzchar(trimws(lines)))[bad[1]],
                 if (is.na(nf[bad[1]])) "an unreadable number of" else nf[bad[1]], nf[1]))
  x <- utils::read.csv(text = paste(lines, collapse = "\n"), check.names = FALSE, stringsAsFactors = FALSE,
                       colClasses = "character")
  nm <- names(x)
  idc <- identifier_columns(nm)
  idc <- idc[vapply(idc, function(j) any(nzchar(trimws(x[[j]]))), logical(1))]
  if (length(idc))
    stop(sprintf(paste("The file has patient identifier columns with values (column%s %s). It was not read.",
                       "Remove the column, or use the local version of the calculator."),
                 if (length(idc) > 1) "s" else "", paste(idc, collapse = ", ")))
  dcol <- grep("dose", nm, ignore.case = TRUE)
  vcol <- grep("vol", nm, ignore.case = TRUE)
  scol <- grep("struct|roi|organ", nm, ignore.case = TRUE)
  if (length(dcol) != 1 || length(vcol) != 1)
    stop("The CSV needs exactly one dose column and one volume column (headers containing 'dose' and 'vol').")
  scale <- dvh_unit_scale(nm[dcol])
  if (is.na(scale)) stop("The dose column header must give the unit, Gy or cGy (e.g. 'Dose (Gy)').")
  vol_unit <- if (grepl("%", nm[vcol])) "%" else if (grepl("cc|cm3|cm\u00b3", nm[vcol], ignore.case = TRUE)) "cc" else NA
  if (is.na(vol_unit)) stop("The volume column header must give the unit, % or cc (e.g. 'Volume (%)').")
  if (length(scol) == 1) {
    s <- gsub("\\s+", " ", trimws(x[[scol]]))
    if (any(!nzchar(s))) stop("Data row ", which(!nzchar(s))[1], ": the structure has no name.")
    runs <- rle(s)$values
    dup <- unique(runs[duplicated(runs)])
    if (length(dup)) stop("Structure '", dup[1], "' appears more than once, in separate blocks of rows. ",
                          "Export each structure once, or rename one of them.")
    groups <- split(seq_len(nrow(x)), factor(s, levels = unique(s)))
  } else groups <- list(structure = seq_len(nrow(x)))
  lapply(stats::setNames(names(groups), names(groups)), function(g) {
    ix <- groups[[g]]
    d <- suppressWarnings(data.frame(dose = as.numeric(x[[dcol]][ix]) * scale, volume = as.numeric(x[[vcol]][ix])))
    if (anyNA(d)) stop("The DVH contains non-numeric values.")
    validate_dvh(d$dose, d$volume, vol_unit, if (length(groups) > 1) sprintf("Structure '%s'", g) else "The DVH")
    attr(d, "volume_unit") <- vol_unit
    d
  })
}

# Read an uploaded DVH file: size limit, identifier refusal, format choice,
# curve validation. Shiny stores an upload as a temporary file on the server;
# this function writes nothing and logs nothing, and the caller deletes the
# temporary file at once and keeps only the parsed curves in memory.
read_dvh_upload <- function(path, size) {
  if (size > DVH_MAX_BYTES) stop("The file is larger than 2 MB.")
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  bad <- identifier_lines(lines)
  if (length(bad))
    stop(sprintf(paste("The file has patient identifier fields with values (line%s %s). It was not read.",
                       "Remove or blank the patient name and ID fields, or use the local version of the calculator."),
                 if (length(bad) > 1) "s" else "", paste(utils::head(bad, 5), collapse = ", ")))
  if (any(grepl("^\\s*Structure\\s*:", lines))) parse_dvh_eclipse(lines) else parse_dvh_csv(lines)
}

# --- feasibility ---------------------------------------------------------------

# Status of a planned value against a limit, in words. basis "EQD2": a
# maximum compared with its LQ-converted constraint; "physical": a mean
# compared with the source's physical mean limit.
oar_status <- function(planned, limit, op = "≤", basis = "EQD2") {
  what <- if (identical(basis, "physical")) "the source limit" else "the converted constraint"
  ifelse(is.na(planned), "not assessed",
    ifelse(!limit_met(planned, limit, op), paste("exceeds", what),
      ifelse(planned > limit * (1 - NEAR_FRACTION), paste("within 5% of", what), "within")))
}
status_exceeds <- function(status) grepl("^exceeds", status)
# Four-way class of a status (review R05): an indeterminate or unassessed status is never "within".
status_class <- function(status) {
  ifelse(is.na(status) | status == "not assessed", "not assessed",
    ifelse(grepl("^exceeds", status), "exceeds",
      ifelse(grepl("^indeterminate", status), "indeterminate", "within")))
}

# Description of a source's fractionation and metric.
source_fractionation_text <- function(row) {
  if (row$n_fractions == 0) "conventional fractionation (about 1.8-2 Gy per fraction)" else
    sprintf("%d fractions", as.integer(row$n_fractions))
}
source_metric_text <- function(row, metric) {
  if (metric == "Dmean") return("Mean dose (physical)")
  if (identical(row$library, "HYHOPE")) "Max (volume definition not specified by the protocol)" else "Max"
}

# One organ at one regimen (m fractions, prescription dose per fraction d).
# input: list(mode = "percent" | "gy" | "dvh", value = number or DVH data frame)
#   percent: planned dose as % of the prescription total dose (plan scales
#            linearly with the prescription)
#   gy:      planned physical max or mean dose (Gy) in m fractions
#   dvh:     cumulative DVH (Gy) from a plan of m fractions; mean bin by bin,
#            maximum as the interval of the highest occupied bin
# A maximum is compared as EQD2 with the LQ-converted constraint. A mean is
# compared as a physical dose with the source's physical mean limit; its EQD2
# is shown for illustration only. ab: organ alpha/beta; ab_range: c(lo, hi)
# for the status-flip check.
oar_row <- function(lib, organ_key, input, m, d, library = "QUANTEC", ab = NULL,
                    ab_range = NULL) {
  o <- OAR_ORGANS[OAR_ORGANS$key == organ_key, ]
  if (nrow(o) != 1) stop("unknown organ: ", organ_key)
  if (is.null(ab)) ab <- default_oar_ab()[[o$ab_class]]
  row <- pick_constraint(lib, organ_key, library, m)
  entered <- !(is.null(input) || is.null(input$value) || (is.numeric(input$value) && length(input$value) == 1 && is.na(input$value)))
  base <- data.frame(organ = o$label, metric = o$metric, library = library, ab = ab,
                     source = NA_character_, source_metric = NA_character_, source_fractionation = NA_character_,
                     basis = NA_character_, operator = NA_character_, source_n = NA_real_,
                     limit_eqd2 = NA_real_, limit_physical = NA_real_, planned_physical = NA_real_, planned_eqd2 = NA_real_,
                     planned_max_lo = NA_real_, planned_max_hi = NA_real_,
                     margin_gy = NA_real_, status = "not assessed", flags = "", flip = "", entered = entered,
                     stringsAsFactors = FALSE)
  if (is.null(row)) {
    base$status <- if (entered) "not assessed: no verified constraint in this library" else "no verified constraint in this library"
    return(base)
  }
  mean_metric <- o$metric == "Dmean"
  basis <- if (mean_metric) "physical" else "EQD2"
  base$source <- sprintf("%s (%s, %s)", row$source_short, row$id, row$page)
  base$source_metric <- source_metric_text(row, o$metric)
  base$source_fractionation <- source_fractionation_text(row)
  base$basis <- if (mean_metric) "physical mean (source metric)" else sprintf("EQD2 (LQ-converted, α/β %g Gy)", ab)
  base$operator <- row$operator
  base$source_n <- row$n_fractions
  if (mean_metric) {
    base$limit_physical <- row$limit_gy
  } else {
    base$limit_eqd2 <- constraint_eqd2(row$limit_gy, row$n_fractions, ab)
    base$limit_physical <- convert_limit(row$limit_gy, row$n_fractions, m, ab)
  }
  flags <- character()
  if (row$limit_lo_gy < row$limit_gy)
    flags <- c(flags, sprintf("source gives a range %g-%g Gy; the upper end is this calculator's choice, not a single source limit",
                              row$limit_lo_gy, row$limit_gy))
  if (mean_metric) {
    flags <- c(flags, sprintf(paste("source physical mean stated for %s; physical source comparison only: the effect of",
                                    "changed fractionation is not validated and may understate hypofractionation toxicity"),
                              base$source_fractionation))
    if (row$n_fractions > 0 && row$n_fractions != m)
      flags <- c(flags, sprintf("the selected HYHOPE table is for %d fractions, not %g", row$n_fractions, m))
  } else if (row$n_fractions > 0 && row$n_fractions != m) {
    flags <- c(flags, sprintf("converted from HYHOPE's %d-fraction table", row$n_fractions))
  }
  if (identical(row$library, "HYHOPE") && row$n_fractions == 12)
    flags <- c(flags, "the 12-fraction level was removed by amendment (historical modeling context)")
  if (nzchar(row$scope_note)) flags <- c(flags, row$scope_note)
  if (!entered) {
    base$flags <- paste(flags, collapse = "; ")
    return(base)
  }
  lim <- if (mean_metric) base$limit_physical else base$limit_eqd2
  pl <- switch(input$mode,
    percent = list(phys = input$value / 100 * m * d, eqd2 = eqd2_physical(input$value / 100 * m * d, m, ab), uniform = TRUE),
    gy = list(phys = input$value, eqd2 = eqd2_physical(input$value, m, ab), uniform = TRUE),
    dvh = {
      dd <- cumulative_to_differential(input$value$dose, input$value$volume)
      mm <- dvh_means(dd, m, ab)
      if (o$metric == "Dmax") list(phys = mm$max_hi, eqd2 = eqd2_physical(mm$max_hi, m, ab), uniform = FALSE,
                                   lo = mm$max_lo, hi = mm$max_hi)
      else list(phys = mm$mean_dose, eqd2 = mm$mean_eqd2, uniform = FALSE)
    },
    stop("unknown input mode: ", input$mode))
  if (mean_metric) {
    flags <- c(flags, if (pl$uniform) "mean EQD2 shown for illustration only, converted as if uniform (it understates a non-uniform dose)" else
      "mean EQD2 (bin by bin) shown for illustration only")
  }
  if (o$metric == "Dmax" && pl$phys / m > LQ_POINT_DOSE_FLAG_GY)
    flags <- c(flags, sprintf("dose per fraction above %g Gy: LQ conversion is an extrapolation", LQ_POINT_DOSE_FLAG_GY))
  base$planned_physical <- pl$phys
  base$planned_eqd2 <- pl$eqd2
  planned <- if (mean_metric) pl$phys else pl$eqd2
  base$margin_gy <- lim - planned
  base$status <- oar_status(planned, lim, row$operator, basis)
  if (!is.null(pl$lo)) {
    # A DVH maximum lies somewhere in its top bin [lo, hi]: definite only when
    # both bounds give the same answer.
    base$planned_max_lo <- pl$lo; base$planned_max_hi <- pl$hi
    flags <- c(flags, sprintf("maximum lies in the DVH bin %g-%g Gy; the status uses both bounds", pl$lo, pl$hi))
    if (limit_met(eqd2_physical(pl$hi, m, ab), lim, row$operator)) {
      # the upper bound meets the limit: the status above stands
    } else if (!limit_met(eqd2_physical(pl$lo, m, ab), lim, row$operator)) {
      base$status <- "exceeds the converted constraint"
    } else {
      base$status <- "indeterminate at this DVH resolution"
    }
  }
  if (!is.null(ab_range)) {
    if (mean_metric) {
      base$flip <- "physical mean: α/β does not change the status"
    } else {
      grid <- seq(ab_range[1], ab_range[2], length.out = 41)
      st <- vapply(grid, function(a) {
        r <- oar_row(lib, organ_key, input, m, d, library, a, NULL)
        status_class(r$status)
      }, character(1))
      base$flip <- flip_text(grid, st)
    }
  }
  base$flags <- paste(flags, collapse = "; ")
  base
}

# The status over an alpha/beta grid, in words. st: status classes per grid value ("within",
# "exceeds", "indeterminate", "not assessed"); a logical st is read as exceeds (TRUE) / within.
# "within at every" is said only when every value is definitively within (review R05).
flip_text <- function(grid, st) {
  if (is.logical(st)) st <- ifelse(st, "exceeds", "within")
  rng <- sprintf("%.2g to %.2g Gy", min(grid), max(grid))
  if (length(unique(st)) == 1) {
    return(switch(st[1],
      within = paste("within at every α/β from", rng),
      exceeds = paste("exceeds at every α/β from", rng),
      indeterminate = paste("indeterminate at this DVH resolution at every α/β from", rng),
      paste("not assessed at any α/β from", rng)))
  }
  r <- rle(st); hi <- cumsum(r$lengths); lo <- hi - r$lengths + 1
  seg <- function(cls) {
    k <- which(r$values == cls)
    if (!length(k)) return(NULL)
    paste(cls, "for α/β", paste(sprintf("%.2g-%.2g Gy", grid[lo[k]], grid[hi[k]]), collapse = " and "))
  }
  parts <- unlist(lapply(c("exceeds", "indeterminate", "not assessed"), seg))
  paste0(paste(parts, collapse = "; "), if (any(st == "within")) " (within otherwise)" else "")
}

# Feasibility table over organs. inputs: named list organ_key -> input.
oar_table <- function(lib, inputs, m, d, library = "QUANTEC", ab = default_oar_ab(), ab_ranges = NULL) {
  rows <- lapply(OAR_ORGANS$key, function(k) {
    a <- ab[[OAR_ORGANS$ab_class[OAR_ORGANS$key == k]]]
    rng <- if (is.null(ab_ranges)) NULL else ab_ranges[[OAR_ORGANS$ab_class[OAR_ORGANS$key == k]]]
    oar_row(lib, k, inputs[[k]], m, d, library, a, rng)
  })
  out <- do.call(rbind, rows)
  # Which assessed constraint binds first: the smallest margin relative to the
  # limit, each in its own basis (EQD2 for maxima, physical for means).
  rel <- out$margin_gy / ifelse(out$metric == "Dmean", out$limit_physical, out$limit_eqd2)
  out$binds_first <- FALSE
  if (any(is.finite(rel))) out$binds_first[which.min(ifelse(is.finite(rel), rel, Inf))] <- TRUE
  out
}

# What a selection of entered organs leaves unassessed in a library: the
# source's mandatory metrics (HYHOPE per-protocol limits) not entered, and
# source metrics this calculator does not offer. complete is FALSE while a
# mandatory metric is unassessed, so no result may be read as "all met".
oar_coverage <- function(lib, entered, library = "QUANTEC") {
  keys <- names(entered)[!vapply(entered, function(v) is.null(v) || (is.atomic(v) && length(v) == 1 && is.na(v)), logical(1))]
  if (library == "HYHOPE") {
    h <- lib[lib$library == "HYHOPE", ]
    mand <- unique(h$organ[h$level == "per_protocol"])
    mand_labels <- OAR_ORGANS$label[match(mand, OAR_ORGANS$hyhope_organ)]
    mandatory_unassessed <- mand_labels[!OAR_ORGANS$key[match(mand, OAR_ORGANS$hyhope_organ)] %in% keys]
    offered <- paste(OAR_ORGANS$hyhope_organ, OAR_ORGANS$metric)
    extra <- unique(h[!paste(h$organ, h$metric) %in% offered & h$level %in% c("per_protocol", "recommended"), c("organ", "metric")])
    other <- sprintf("%s (%s)", extra$organ, ifelse(extra$metric == "Dmax", "maximum", "mean"))
  } else {
    mandatory_unassessed <- character()
    other <- c("Larynx V50Gy ≤ 27% (only the mean is checked)",
               "Parotid: both glands below about 25 Gy (the one-gland limit is used)",
               "Cochlea: the more conservative 35 Gy mean")
  }
  list(entered = OAR_ORGANS$label[OAR_ORGANS$key %in% keys], mandatory_unassessed = mandatory_unassessed,
       other_unassessed = other, complete = length(mandatory_unassessed) == 0)
}

# One sentence on coverage, for every place that reports organ limits.
oar_coverage_text <- function(cov) {
  s <- "These are selected constraint metrics, not full protocol compliance."
  if (length(cov$mandatory_unassessed))
    s <- paste(s, sprintf("Mandatory metrics not assessed: %s, so the limits cannot all be called met.",
                          paste(cov$mandatory_unassessed, collapse = ", ")))
  if (length(cov$other_unassessed))
    s <- paste(s, sprintf("Other source metrics not assessed here: %s.", paste(cov$other_unassessed, collapse = "; ")))
  s
}

# For the search under uncertainty: organs entered as % of the prescription.
# Returns a function(regs, tx) giving a logical regimens x draws matrix that
# is TRUE where every entered organ meets its limit: a maximum under that
# draw's organ alpha/beta (tx$ab_cord or tx$ab_oar), a mean on physical dose.
# An entered organ without a constraint is an error, never a pass.
# One organ-mapping validator (review R03, 2026-10-26), used by both replay paths and the
# direct checkers: a named mapping of unique, nonblank, supported organ keys to one finite
# percent of the prescription (0-200) each. An unnamed list, a duplicate or unknown key, or a
# non-scalar value is refused, never read as "no organ limit". NULL or NA values are organs left
# blank (not assessed) and are dropped after the names are checked.
validate_oar_percents <- function(percents) {
  if (is.null(percents) || length(percents) == 0) return(list())
  nm <- names(percents)
  if (is.null(nm) || any(is.na(nm) | !nzchar(trimws(nm))))
    stop("organ limits must name each organ (an unnamed or blank organ key was given).", call. = FALSE)
  if (anyDuplicated(nm)) stop("organ '", nm[anyDuplicated(nm)], "' is given more than once.", call. = FALSE)
  bad <- setdiff(nm, OAR_ORGANS$key)
  if (length(bad)) stop("unknown organ key(s): ", paste(bad, collapse = ", "), ".", call. = FALSE)
  out <- lapply(nm, function(k) {
    v <- percents[[k]]
    if (is.null(v) || (length(v) == 1 && is.na(v) && !is.nan(v))) return(NULL)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 0 || v > 200)
      stop("organ '", k, "': the dose must be one number from 0 to 200% of the prescription.", call. = FALSE)
    v
  })
  names(out) <- nm
  out[!vapply(out, is.null, logical(1))]
}

oar_draw_checker <- function(lib, percents, library = "QUANTEC") {
  percents <- validate_oar_percents(percents)
  if (length(percents) == 0) return(NULL)
  function(regs, tx) {
    ok <- matrix(TRUE, nrow(regs), nrow(tx))
    for (k in names(percents)) {
      o <- OAR_ORGANS[OAR_ORGANS$key == k, ]
      ab <- if (o$ab_class == "cord") tx$ab_cord else tx$ab_oar
      phys <- percents[[k]] / 100 * regs$m * regs$d
      for (i in seq_len(nrow(regs))) {
        row <- pick_constraint(lib, k, library, regs$m[i])
        if (is.null(row)) no_constraint_error(k, library)
        ok[i, ] <- ok[i, ] & organ_meets(row, o$metric, phys[i], regs$m[i], ab)
      }
    }
    ok
  }
}

# Deterministic organ check for suggest_regimens(): organs entered as % of the
# prescription, fixed alpha/beta per class. Returns function(m, d) -> logical
# (TRUE where every entered organ meets its limit), with its settings attached
# for the scenario record, or NULL when nothing is entered.
build_oar_check <- function(lib, percents, library = "QUANTEC", ab = default_oar_ab()) {
  percents <- validate_oar_percents(percents)
  if (length(percents) == 0) return(NULL)
  f <- function(m, d) {
    ok <- rep(TRUE, length(m))
    for (k in names(percents)) {
      o <- OAR_ORGANS[OAR_ORGANS$key == k, ]
      a <- ab[[o$ab_class]]
      for (i in seq_along(m)) {
        row <- pick_constraint(lib, k, library, m[i])
        if (is.null(row)) no_constraint_error(k, library)
        ok[i] <- ok[i] && organ_meets(row, o$metric, percents[[k]] / 100 * m[i] * d[i], m[i], a)
      }
    }
    ok
  }
  attr(f, "settings") <- list(library = library, ab = as.list(ab), percents = percents, oar_version = OAR_VERSION)
  f
}
