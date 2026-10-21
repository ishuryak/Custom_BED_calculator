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

OAR_VERSION <- "3.0"
NEAR_FRACTION <- 0.05          # "within 5% of" the converted limit
LQ_POINT_DOSE_FLAG_GY <- 5     # dose per fraction above which a point-dose conversion is flagged

load_constraints <- function(path = file.path("data", "oar_constraints.csv")) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8")
  stopifnot(all(c("id", "organ", "metric", "limit_gy", "limit_lo_gy", "n_fractions", "level",
                  "ab_class", "library", "source_short", "page", "quote") %in% names(x)))
  x
}

# Organs offered in the app, each tied to the constraint used by default in
# each library. QUANTEC values are stated in 2-Gy fractions (n_fractions 0).
# HYHOPE values are stated for 12, 15 and 20 fractions. Organs whose
# conventional constraint was not found in a verified source have no QUANTEC
# entry (NA) and can be assessed with the HYHOPE library only.
OAR_ORGANS <- data.frame(
  key = c("cord", "brainstem", "mandible", "parotid", "submandibular", "constrictors",
          "larynx", "oral_cavity", "esophagus", "cochlea", "plexus"),
  label = c("Spinal cord", "Brainstem", "Mandible", "Parotid (spared side)", "Submandibular gland (spared side)",
            "Pharyngeal constrictors", "Larynx", "Oral cavity (outside target)", "Cervical esophagus",
            "Cochlea", "Brachial plexus"),
  metric = c("Dmax", "Dmax", "Dmax", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmean", "Dmax"),
  ab_class = c("cord", "cord", "oar", "oar", "oar", "oar", "oar", "oar", "oar", "oar", "oar"),
  quantec_id = c("OAR079", "OAR089", NA, "OAR086", "OAR088", NA, "OAR094", NA, NA, "OAR099", NA),
  hyhope_organ = c("Spinal cord", "Brainstem", "Mandible", "Contralateral parotid",
                   "Contralateral submandibular gland", "Pharyngeal constrictors", "Larynx",
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

# Differential DVH from a cumulative one. dose: increasing bin edges (Gy);
# volume: cumulative volume at or above each dose (any unit, non-increasing).
# Returns bin mid doses and bin volumes (the last edge's volume is assigned
# to that dose, so a DVH ending above zero volume keeps it).
cumulative_to_differential <- function(dose, volume) {
  o <- order(dose); dose <- dose[o]; volume <- volume[o]
  if (any(diff(volume) > 1e-9 * max(1, max(volume))))
    stop("A cumulative DVH must not increase with dose.")
  k <- length(dose)
  v <- c(volume[-k] - volume[-1], volume[k])
  mid <- c((dose[-k] + dose[-1]) / 2, dose[k])
  keep <- v > 0
  data.frame(dose = mid[keep], volume = v[keep])
}

# Mean physical dose and mean EQD2 (bin by bin) of a differential DVH whose
# doses were delivered in n fractions.
dvh_means <- function(diff_dvh, n, ab) {
  w <- diff_dvh$volume / sum(diff_dvh$volume)
  list(mean_dose = sum(w * diff_dvh$dose),
       mean_eqd2 = sum(w * eqd2_physical(diff_dvh$dose, n, ab)),
       max_dose = max(diff_dvh$dose))
}

# --- DVH files ---------------------------------------------------------------

DVH_MAX_BYTES <- 2 * 1024^2
ID_PATTERN <- "(patient\\s*name|patient\\s*id|\\bmrn\\b|medical\\s*record|date\\s*of\\s*birth|\\bdob\\b|birth\\s*date)\\s*[:=,;\\t]\\s*\\S"

# Refuse text that carries a patient identifier with a value. Returns the
# offending line numbers (never the content) or integer(0).
identifier_lines <- function(lines) which(grepl(ID_PATTERN, lines, ignore.case = TRUE, perl = TRUE))

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
  out <- list()
  for (k in seq_along(starts)) {
    blk <- lines[starts[k]:ends[k]]
    name <- trimws(sub("^\\s*Structure\\s*:", "", blk[1]))
    h <- grep("Dose \\[c?Gy\\]", blk)
    if (length(h) == 0) next
    hdr <- blk[h[1]]
    scale <- dvh_unit_scale(regmatches(hdr, regexpr("Dose \\[c?Gy\\]", hdr)))
    vol_unit <- if (grepl("Volume \\[%\\]", hdr)) "%" else if (grepl("Volume \\[(cm\u00b3|cm3|cc)\\]", hdr)) "cc" else NA
    if (is.na(vol_unit)) stop("Structure '", name, "': volume column must be in % or cm\u00b3.")
    rows <- blk[(h[1] + 1):length(blk)]
    rows <- rows[grepl("^\\s*[0-9.eE+-]+(\\s+[0-9.eE+-]+)+\\s*$", rows)]
    if (length(rows) < 2) stop("Structure '", name, "': fewer than two DVH rows.")
    num <- do.call(rbind, lapply(strsplit(trimws(rows), "\\s+"), as.numeric))
    out[[name]] <- data.frame(dose = num[, 1] * scale, volume = num[, ncol(num)])
    attr(out[[name]], "volume_unit") <- vol_unit
  }
  if (length(out) == 0) stop("No DVH table found.")
  out
}

# Parse a two-column CSV (dose, cumulative volume), optionally with a
# structure column. Units come from the header: "cGy" or "Gy"; "%" or
# "cc"/"cm3". A generic format for exports from any planning system.
parse_dvh_csv <- function(lines) {
  x <- utils::read.csv(text = paste(lines, collapse = "\n"), check.names = FALSE, stringsAsFactors = FALSE)
  nm <- names(x)
  dcol <- grep("dose", nm, ignore.case = TRUE)
  vcol <- grep("vol", nm, ignore.case = TRUE)
  scol <- grep("struct|roi|organ", nm, ignore.case = TRUE)
  if (length(dcol) != 1 || length(vcol) != 1)
    stop("The CSV needs exactly one dose column and one volume column (headers containing 'dose' and 'vol').")
  scale <- dvh_unit_scale(nm[dcol])
  if (is.na(scale)) stop("The dose column header must give the unit, Gy or cGy (e.g. 'Dose (Gy)').")
  vol_unit <- if (grepl("%", nm[vcol])) "%" else if (grepl("cc|cm3|cm\u00b3", nm[vcol], ignore.case = TRUE)) "cc" else NA
  if (is.na(vol_unit)) stop("The volume column header must give the unit, % or cc (e.g. 'Volume (%)').")
  groups <- if (length(scol) == 1) split(seq_len(nrow(x)), x[[scol]]) else list(structure = seq_len(nrow(x)))
  lapply(groups, function(ix) {
    d <- data.frame(dose = as.numeric(x[[dcol]][ix]) * scale, volume = as.numeric(x[[vcol]][ix]))
    if (anyNA(d)) stop("The DVH contains non-numeric values.")
    attr(d, "volume_unit") <- vol_unit
    d
  })
}

# Read an uploaded DVH file: size limit, identifier refusal, format choice.
# Nothing is written to disk or logged; the caller keeps the result in memory.
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

# Status of a planned value against a converted limit, in words.
oar_status <- function(planned, limit) {
  ifelse(is.na(planned), "not assessed",
    ifelse(planned > limit + 1e-9, "exceeds the converted constraint",
      ifelse(planned > limit * (1 - NEAR_FRACTION), "within 5% of the converted constraint", "within")))
}

# One organ at one regimen (m fractions, prescription dose per fraction d).
# input: list(mode = "percent" | "gy" | "dvh", value = number or DVH data frame)
#   percent: planned dose as % of the prescription total dose (plan scales
#            linearly with the prescription; mean-dose organs treated as
#            uniform, flagged)
#   gy:      planned physical max or mean dose (Gy) in m fractions (mean
#            treated as uniform, flagged)
#   dvh:     cumulative DVH (Gy) from a plan of m fractions; max and mean
#            EQD2 computed bin by bin
# ab: organ alpha/beta; ab_range: c(lo, hi) for the status-flip check.
oar_row <- function(lib, organ_key, input, m, d, library = "QUANTEC", ab = NULL,
                    ab_range = NULL) {
  o <- OAR_ORGANS[OAR_ORGANS$key == organ_key, ]
  if (is.null(ab)) ab <- default_oar_ab()[[o$ab_class]]
  row <- pick_constraint(lib, organ_key, library, m)
  base <- data.frame(organ = o$label, metric = o$metric, library = library, ab = ab,
                     source = NA_character_, source_n = NA_real_, limit_eqd2 = NA_real_,
                     limit_physical = NA_real_, planned_physical = NA_real_, planned_eqd2 = NA_real_,
                     margin_gy = NA_real_, status = "not assessed", flags = "", flip = "",
                     stringsAsFactors = FALSE)
  if (is.null(row)) {
    base$status <- "no verified constraint in this library"
    return(base)
  }
  base$source <- sprintf("%s (%s, %s)", row$source_short, row$id, row$page)
  base$source_n <- row$n_fractions
  base$limit_eqd2 <- constraint_eqd2(row$limit_gy, row$n_fractions, ab)
  base$limit_physical <- convert_limit(row$limit_gy, row$n_fractions, m, ab)
  flags <- character()
  if (row$limit_lo_gy < row$limit_gy) flags <- c(flags, sprintf("source gives a range %g-%g Gy; the upper value is used", row$limit_lo_gy, row$limit_gy))
  if (row$n_fractions > 0 && row$n_fractions != m) flags <- c(flags, sprintf("converted from HYHOPE's %d-fraction table", row$n_fractions))
  if (is.null(input) || is.null(input$value) || (is.numeric(input$value) && is.na(input$value))) {
    base$flags <- paste(flags, collapse = "; ")
    return(base)
  }
  pl <- switch(input$mode,
    percent = list(phys = input$value / 100 * m * d, eqd2 = eqd2_physical(input$value / 100 * m * d, m, ab), uniform = TRUE),
    gy = list(phys = input$value, eqd2 = eqd2_physical(input$value, m, ab), uniform = TRUE),
    dvh = {
      dd <- cumulative_to_differential(input$value$dose, input$value$volume)
      mm <- dvh_means(dd, m, ab)
      if (o$metric == "Dmax") list(phys = mm$max_dose, eqd2 = eqd2_physical(mm$max_dose, m, ab), uniform = FALSE)
      else list(phys = mm$mean_dose, eqd2 = mm$mean_eqd2, uniform = FALSE)
    },
    stop("unknown input mode: ", input$mode))
  if (o$metric == "Dmean" && pl$uniform)
    flags <- c(flags, "mean dose converted as if uniform: this understates the EQD2 of a non-uniform dose (use a DVH)")
  if (o$metric == "Dmax" && pl$phys / m > LQ_POINT_DOSE_FLAG_GY)
    flags <- c(flags, sprintf("dose per fraction above %g Gy: LQ conversion is an extrapolation", LQ_POINT_DOSE_FLAG_GY))
  base$planned_physical <- pl$phys
  base$planned_eqd2 <- pl$eqd2
  base$margin_gy <- base$limit_eqd2 - pl$eqd2
  base$status <- oar_status(pl$eqd2, base$limit_eqd2)
  if (!is.null(ab_range)) {
    grid <- seq(ab_range[1], ab_range[2], length.out = 41)
    st <- vapply(grid, function(a) {
      r <- oar_row(lib, organ_key, input, m, d, library, a, NULL)
      r$status == "exceeds the converted constraint"
    }, logical(1))
    base$flip <- flip_text(grid, st)
  }
  base$flags <- paste(flags, collapse = "; ")
  base
}

# Where over an alpha/beta grid a constraint is exceeded (st TRUE), in words.
flip_text <- function(grid, st) {
  rng <- sprintf("%.2g to %.2g Gy", min(grid), max(grid))
  if (all(st)) return(paste("exceeds at every α/β from", rng))
  if (!any(st)) return(paste("within at every α/β from", rng))
  r <- rle(st); hi <- cumsum(r$lengths); lo <- hi - r$lengths + 1
  parts <- sprintf("%.2g-%.2g Gy", grid[lo[r$values]], grid[hi[r$values]])
  paste("exceeds for α/β", paste(parts, collapse = " and "), "(within otherwise)")
}

# Feasibility table over organs. inputs: named list organ_key -> input.
oar_table <- function(lib, inputs, m, d, library = "QUANTEC", ab = default_oar_ab(), ab_ranges = NULL) {
  rows <- lapply(OAR_ORGANS$key, function(k) {
    a <- ab[[OAR_ORGANS$ab_class[OAR_ORGANS$key == k]]]
    rng <- if (is.null(ab_ranges)) NULL else ab_ranges[[OAR_ORGANS$ab_class[OAR_ORGANS$key == k]]]
    oar_row(lib, k, inputs[[k]], m, d, library, a, rng)
  })
  out <- do.call(rbind, rows)
  # Which assessed constraint binds first: the smallest margin relative to the limit.
  rel <- out$margin_gy / out$limit_eqd2
  out$binds_first <- FALSE
  if (any(is.finite(rel))) out$binds_first[which.min(ifelse(is.finite(rel), rel, Inf))] <- TRUE
  out
}

# For the search under uncertainty: organs entered as % of the prescription.
# Returns a function(regs, tx) giving a logical regimens x draws matrix that
# is TRUE where every entered organ is within its converted constraint under
# that draw's organ alpha/beta (tx$ab_cord or tx$ab_oar).
oar_draw_checker <- function(lib, percents, library = "QUANTEC") {
  percents <- percents[!vapply(percents, function(v) is.null(v) || is.na(v), logical(1))]
  if (length(percents) == 0) return(NULL)
  function(regs, tx) {
    ok <- matrix(TRUE, nrow(regs), nrow(tx))
    for (k in names(percents)) {
      o <- OAR_ORGANS[OAR_ORGANS$key == k, ]
      ab <- if (o$ab_class == "cord") tx$ab_cord else tx$ab_oar
      phys <- percents[[k]] / 100 * regs$m * regs$d
      for (i in seq_len(nrow(regs))) {
        row <- pick_constraint(lib, k, library, regs$m[i])
        if (is.null(row)) next
        lim <- constraint_eqd2(row$limit_gy, row$n_fractions, ab)
        ok[i, ] <- ok[i, ] & eqd2_physical(phys[i], regs$m[i], ab) <= lim + 1e-9
      }
    }
    ok
  }
}

# Deterministic organ check for suggest_regimens(): organs entered as % of the
# prescription, fixed alpha/beta per class. Returns function(m, d) -> logical
# (TRUE where every entered organ is within its converted constraint), with
# its settings attached for the scenario record, or NULL when nothing is entered.
build_oar_check <- function(lib, percents, library = "QUANTEC", ab = default_oar_ab()) {
  percents <- percents[!vapply(percents, function(v) is.null(v) || is.na(v), logical(1))]
  if (length(percents) == 0) return(NULL)
  f <- function(m, d) {
    ok <- rep(TRUE, length(m))
    for (k in names(percents)) {
      a <- ab[[OAR_ORGANS$ab_class[OAR_ORGANS$key == k]]]
      for (i in seq_along(m)) {
        row <- pick_constraint(lib, k, library, m[i])
        if (is.null(row)) next
        ok[i] <- ok[i] && eqd2_physical(percents[[k]] / 100 * m[i] * d[i], m[i], a) <=
          constraint_eqd2(row$limit_gy, row$n_fractions, a) + 1e-9
      }
    }
    ok
  }
  attr(f, "settings") <- list(library = library, ab = as.list(ab), percents = percents)
  f
}
