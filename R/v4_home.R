# Version 4 main page: hypofractionated regimens comparable to the standard
# 35 x 2 Gy in 47 days. Uses the version 3 engine (suggest_regimens) unchanged.
#
# A regimen is listed when, compared with the standard:
#   1. tumor BED (DI or DD model) >= the standard's BED under that model - 5 Gy,
#   2. early (acute mucosal) BED <= the standard's early BED + 5 Gy,
#   3. late BED <= the standard's late BED + 5 Gy,
# and it has at least 5 fractions fewer than the standard (5 fractions per week
# from a Monday, dose per fraction above 2 Gy).

V4_STANDARD <- list(m = 35, d = 2, T = 47)
V4_MARGINS <- list(tumor = 5, early = 5, late = 5)
V4_MIN_FEWER_FRACTIONS <- 5  # listed regimens have at most 35 - 5 = 30 fractions
V4_GRID <- list(m_range = 5:(V4_STANDARD$m - V4_MIN_FEWER_FRACTIONS), d_step = 0.05, d_max = 8)

# Every qualifying (fractions, dose per fraction) pair, the standard's BEDs and,
# when nothing can be listed, the engine's reason.
v4_regimens <- function(p, anchor = default_anchor(), standard = V4_STANDARD,
                        margins = V4_MARGINS, grid = V4_GRID) {
  s <- suggest_regimens(p, standard, anchor = anchor, m_range = grid$m_range,
                        d_step = grid$d_step, d_max = grid$d_max,
                        criterion = "either", rank_by = "best",
                        tol_early = margins$early, tol_late = margins$late,
                        min_gain = -margins$tumor, min_days_saved = 1)
  rows <- attr(s, "feasible")
  if (is.null(rows)) rows <- empty_suggestions()[, setdiff(names(empty_suggestions()), "next_step")]
  rows <- rows[order(rows$m, rows$d), , drop = FALSE]
  rownames(rows) <- NULL
  list(rows = rows, standard = compute_all(standard$m, standard$d, standard$T, p, anchor),
       reason = attr(s, "reason"), margins = margins)
}

# "2.65–2.80", or one value when the range has a single dose.
v4_range <- function(x, digits = 2) {
  lo <- min(x); hi <- max(x)
  f <- paste0("%.", digits, "f")
  if (abs(hi - lo) < 1e-9) sprintf(f, lo) else sprintf(paste0(f, "–", f), lo, hi)
}

# Signed range with a true minus sign, "−4.5 to −0.3"; "n/a" when not evaluable.
v4_signed_range <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("n/a")
  f <- function(v) gsub("-", "−", sprintf(paste0("%+.", digits, "f"), v), fixed = TRUE)
  lo <- round(min(x), digits) + 0; hi <- round(max(x), digits) + 0  # + 0 turns -0 into 0
  if (lo == hi) f(lo) else paste(f(lo), "to", f(hi))
}

# Doses with refit support below this share are marked (decision 2, 2026-10-28).
V4_SUPPORT_MIN <- 0.8

# Share of the refit's primary draws (DI posterior draws paired with calibrated-DD
# draws) in which a regimen meets the tumor condition under the DI or the DD model.
# It uses the refit's own parameters, not the Model settings values.
v4_refit_support <- function(rows, di, dd, standard = V4_STANDARD, margin = V4_MARGINS$tumor) {
  if (nrow(rows) == 0) return(numeric())
  td <- tumor_deltas(rows[, c("m", "d", "T")], standard, di, dd)
  ok <- (td$valid_di & td$delta_di >= -margin - 1e-9) | (td$valid_dd & td$delta_dd >= -margin - 1e-9)
  ok[is.na(ok)] <- FALSE
  rowMeans(ok)
}

# Settings that make the Report card and the other menu pages use the main-page rules
# for one regimen: comparison 35 x 2 Gy, tumor threshold -5 Gy, early and late limits
# equal to the comparison plus 5 Gy, 5 fractions a week from a Monday.
v4_report_inputs <- function(m, d) {
  list(m = m, d = d, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE,
       comparator = "custom", ref_m = V4_STANDARD$m, ref_d = V4_STANDARD$d, ref_fx_week = "5", ref_auto_T = TRUE,
       min_gain = -V4_MARGINS$tumor, early_mode = "reference", tol_early = V4_MARGINS$early,
       late_mode = "reference", tol_late = V4_MARGINS$late)
}

# Published trial arm that the main list can contain (5 fractions a week).
V4_TRIAL_ARMS <- data.frame(m = 20, d = 2.75, T = 26, label = "HYPNO test arm (55 Gy in 20) inside")

# One row per number of fractions, most fractions first: dose and total-dose ranges,
# time, the tumor BED change range under each model, the early and late BED change
# ranges, where both models meet the tumor condition, the
# refit support range, and notes (outside the trial data, acute BED above HYPNO's,
# low refit support, a published trial arm). rows: v4_regimens()$rows, with column
# `support` when refit support is shown. hypno_early: HYPNO 55 Gy/20 early BED
# under the same parameters. data_d_max: largest dose per fraction in the trial data.
v4_summary <- function(rows, standard_T = V4_STANDARD$T, hypno_early = Inf, data_d_max = Inf,
                       trials = V4_TRIAL_ARMS, support_min = V4_SUPPORT_MIN) {
  if (nrow(rows) == 0)
    return(data.frame(m = integer(), dose = character(), total = character(), T = numeric(),
                      saved = numeric(), d_di = character(), d_dd = character(),
                      d_early = character(), d_late = character(), models = character(),
                      support = character(), notes = character(), outside_data = logical(),
                      above_hypno = logical(), low_support = logical(), trial = character()))
  has_support <- "support" %in% names(rows)
  out <- do.call(rbind, lapply(split(rows, rows$m), function(x) {
    both <- x$qual_di & x$qual_dd
    one <- if (all(x$qual_di | both)) "DI" else if (all(x$qual_dd | both)) "DD" else "DI or DD"
    models <- if (all(both)) "both models" else if (any(both))
      sprintf("both at %s Gy, else %s only", v4_range(x$d[both]), one) else sprintf("%s only", one)
    outside <- x$d > data_d_max + 1e-9
    hot <- x$bed_early > hypno_early + 1e-9
    low <- if (has_support) x$support < support_min - 1e-12 else rep(FALSE, nrow(x))
    trial <- paste(trials$label[trials$m == x$m[1] & trials$T == x$T[1] &
                                  vapply(trials$d, function(d) any(abs(x$d - d) < 1e-9), logical(1))],
                   collapse = "; ")
    notes <- c(
      if (all(outside)) "outside trial data" else if (any(outside))
        sprintf("outside trial data above %.2f Gy", data_d_max),
      if (any(hot)) sprintf("acute > HYPNO from %.2f Gy", min(x$d[hot])),
      if (any(low)) sprintf("support < %d%% at %s Gy", round(100 * support_min), v4_range(x$d[low])),
      if (nzchar(trial)) trial)
    data.frame(m = x$m[1], dose = v4_range(x$d), total = v4_range(x$m * x$d, 1),
               T = x$T[1], saved = standard_T - x$T[1],
               d_di = v4_signed_range(x$d_di), d_dd = v4_signed_range(x$d_dd),
               d_early = v4_signed_range(x$d_early), d_late = v4_signed_range(x$d_late), models = models,
               support = if (has_support) {
                 lo <- round(100 * min(x$support)); hi <- round(100 * max(x$support))
                 if (lo == hi) sprintf("%d%%", lo) else sprintf("%d–%d%%", lo, hi)
               } else "",
               notes = paste(notes, collapse = "; "),
               outside_data = any(outside), above_hypno = any(hot), low_support = any(low), trial = trial)
  }))
  out <- out[order(-out$m), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Long table for the main-page plot: for each fraction number and endpoint, the
# lowest and highest BED change (regimen minus standard) over the listed doses.
v4_plot_data <- function(rows) {
  ends <- c(d_di = "Tumor (DI model)", d_dd = "Tumor (DD model)",
            d_early = "Early toxicity", d_late = "Late toxicity")
  if (nrow(rows) == 0)
    return(data.frame(m = integer(), endpoint = factor(levels = ends), lo = numeric(), hi = numeric()))
  out <- do.call(rbind, lapply(split(rows, rows$m), function(x) do.call(rbind, lapply(names(ends), function(k) {
    v <- x[[k]][is.finite(x[[k]])]
    if (length(v) == 0) return(NULL)
    data.frame(m = x$m[1], endpoint = ends[[k]], lo = min(v), hi = max(v))
  }))))
  out$endpoint <- factor(out$endpoint, levels = ends)
  rownames(out) <- NULL
  out
}

# Refit draws for the support column, loaded once per R process. NULL when the
# draws cannot be read (the column then says "not available").
V4_DRAWS_CACHE <- new.env()
v4_draws <- function() {
  if (!exists("x", envir = V4_DRAWS_CACHE, inherits = FALSE))
    assign("x", tryCatch(list(di = load_draws("di"), dd = load_draws("dd")), error = function(e) NULL),
           envir = V4_DRAWS_CACHE)
  get("x", envir = V4_DRAWS_CACHE)
}

# Why the main list is empty, in the main page's terms (the standard, not a chosen comparison).
v4_empty_reason <- function(reason) {
  if (is.null(reason)) return(NULL)
  switch(reason,
    negative_reference_di = , negative_reference_dd = , negative_reference_both =
      "The standard's tumor BED is not positive with these settings, so no schedule can meet the tumor condition. Check the repopulation settings.",
    negative_reference_early = , negative_reference_late =
      "The standard's side-effect BED is not positive with these settings. Check the toxicity settings.",
    dd_unavailable = "The DD model is not available with these settings.",
    paste("No search was possible:", reason))
}
