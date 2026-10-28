# Version 4 main page: hypofractionated regimens comparable to the standard
# 35 x 2 Gy in 47 days. Uses the version 3 engine (suggest_regimens) unchanged.
#
# A regimen is listed when, compared with the standard:
#   1. tumor BED (DI or DD model) >= the standard's BED under that model - 5 Gy,
#   2. early (acute mucosal) BED <= the standard's early BED + 5 Gy,
#   3. late BED <= the standard's late BED + 5 Gy,
# and it is shorter than the standard (5 fractions per week from a Monday,
# dose per fraction above 2 Gy).

V4_STANDARD <- list(m = 35, d = 2, T = 47)
V4_MARGINS <- list(tumor = 5, early = 5, late = 5)
V4_GRID <- list(m_range = 5:34, d_step = 0.05, d_max = 8)

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

# Published trial arms that the main list can contain (5 fractions a week).
V4_TRIAL_ARMS <- data.frame(m = c(20, 15), d = c(2.75, 3.1), T = c(26, 19),
                            label = c("HYPNO test arm, 55 Gy in 20", "HYHOPE, 46.5 Gy in 15"))

# One row per number of fractions: the range of doses per fraction that meet all
# three conditions, whether both tumor models agree (and where), and any
# published trial arm inside the range.
v4_summary <- function(rows, standard_T = V4_STANDARD$T, precedent_max_d = PRECEDENT_MAX_D,
                       trials = V4_TRIAL_ARMS) {
  if (nrow(rows) == 0)
    return(data.frame(m = integer(), dose = character(), total = character(), T = numeric(),
                      saved = numeric(), models = character(), beyond = logical(), trial = character()))
  out <- do.call(rbind, lapply(split(rows, rows$m), function(x) {
    both <- x$qual_di & x$qual_dd
    one <- if (all(x$qual_di | both)) "DI" else if (all(x$qual_dd | both)) "DD" else "DI or DD"
    models <- if (all(both)) "both models" else if (any(both))
      sprintf("both models at %s Gy, %s model only otherwise", v4_range(x$d[both]), one) else
      sprintf("%s model only", one)
    data.frame(m = x$m[1], dose = v4_range(x$d), total = v4_range(x$m * x$d, 1),
               T = x$T[1], saved = standard_T - x$T[1], models = models,
               beyond = any(x$d > precedent_max_d + 1e-9),
               trial = paste(trials$label[trials$m == x$m[1] & trials$T == x$T[1] &
                                            vapply(trials$d, function(d) any(abs(x$d - d) < 1e-9), logical(1))],
                             collapse = "; "))
  }))
  rownames(out) <- NULL
  out
}
