# BED models for the HNSCC BED calculator, version 2.1.
#
# Notation: m = number of fractions, d = dose per fraction (Gy),
# T = overall treatment time (days, first to last treatment day inclusive),
# ab = alpha/beta ratio (Gy), alpha = linear LQ parameter (1/Gy),
# lambda = repopulation rate (1/day), g = background tumor growth (1/day).
#
# All max(0, x) terms use the same smooth version as the v1 website,
# log(1 + exp(3x)) / 3, so v1 and v2 agree at the v1 defaults.

APP_VERSION <- "2.3.0"
EQUATION_VERSION <- "2.2"  # 2.1 + optional Fenwick 2008 acute boundary

SMOOTH_K <- 3  # sharpness (1/day) of the smooth max(0, x), as in v1

# Fowler et al. 2003 acute mucosal tolerance zone (Gy10, radiotherapy alone)
# and the suggested reduction with concurrent chemotherapy (Gy10). Valid only
# with Fowler's acute parameters (alpha 0.35, alpha/beta 10, Tk 7, Tp 2.5) and
# his elapsed-day clock.
FOWLER_ZONE <- c(59, 61)
FOWLER_CHEMO_REDUCTION <- c(3, 5)

# Fenwick et al. 2008 acute tolerance boundary on cell-kill BED (Gy10, no
# repopulation term): B(T) = 69.5 (T/32.2)/sin(T/32.2) - 3.5, T inclusive
# days, valid for T <= 50 days. Cell-kill BED uses alpha/beta 10 Gy and, for
# once-daily fractions, complete repair.
FENWICK <- list(b0 = 69.5, tc = 32.2, shift = 3.5, t_max = 50, ab = 10)

fenwick_limit <- function(T) {
  x <- T / FENWICK$tc
  out <- FENWICK$b0 * x / sin(x) - FENWICK$shift
  out[!(T > 0 & T <= FENWICK$t_max)] <- NA_real_
  out
}

fenwick_cell_kill <- function(m, d) bed_lq(m, d, FENWICK$ab)

# Largest dose per fraction with curative precedent in the benchmark trials
# (HYHOPE 46.5 Gy in 15 fractions). Larger fractions are flagged.
PRECEDENT_MAX_D <- 3.1

# Smooth max(0, x); switches to x where exp() would overflow.
smooth_max0 <- function(x, k = SMOOTH_K) {
  z <- k * x
  ifelse(z > 30, x, log1p(exp(pmin(z, 30))) / k)
}

# Inverse of smooth_max0 for y > 0: the x with smooth_max0(x) = y.
inv_smooth_max0 <- function(y, k = SMOOTH_K) {
  z <- k * y
  ifelse(z > 30, y, log(expm1(z)) / k)
}

# Plain linear-quadratic BED, no repopulation.
bed_lq <- function(m, d, ab) m * d * (1 + d / ab)

# Equivalent dose in 2-Gy fractions.
eqd2 <- function(bed, ab) bed / (1 + 2 / ab)

# Overall time (days, first to last treatment day inclusive) for m fractions
# given once daily Monday to Friday, starting on a Monday.
# 35 fractions gives 47 days.
overall_time <- function(m) {
  w <- (m - 1) %/% 5
  7 * w + (m - 1) %% 5 + 1
}

# Convert the app's inclusive overall time to the clock of a tissue model.
# Fowler, Fenwick, Hartley and Meade count elapsed days (last day minus first
# day), which is the inclusive count minus 1.
acute_model_time <- function(T, clock = "elapsed") {
  switch(clock, elapsed = T - 1, inclusive = T, stop("unknown clock: ", clock))
}

# Mean fraction of tumor cells killed per day (2018 paper, Eq. 3).
p_kill <- function(K, T, alpha) -expm1(-alpha * K / T)

# Tumor BED, dose-independent (DI) repopulation model (2018 paper, Eq. 1):
# repopulation at fixed rate lambda_di starting at fixed time Tk, plus
# optional background growth g over the whole course.
bed_di <- function(m, d, T, alpha, ab, lambda_di, Tk, g = 0) {
  bed_lq(m, d, ab) - lambda_di * smooth_max0(T - Tk) / alpha - g * T / alpha
}

# Calibration of the dose-dependent (DD) model at a fixed anchor regimen.
# Repopulation starts when -ln(surviving fraction) reaches C, at time
# C*T/(alpha*K) with K = bed_lq(m, d, ab). C is solved so that
# BED_DD = BED_DI for the anchor. Background growth g is shared by both
# models and cancels. Returns c_over_alpha = C/alpha (Gy), the onset day at
# the anchor, and a status:
#   "ok"               calibration found, onset after treatment start
#   "no_repopulation"  DI repopulation rate is zero, so DD has none either
#   "infeasible"       no positive onset reproduces DI at the anchor, so
#                      DD values are not defined and are returned as NA.
dd_calibration <- function(anchor, alpha, ab, lambda_di, lambda_dd, Tk) {
  K0 <- bed_lq(anchor$m, anchor$d, ab)
  T0 <- anchor$T
  if (lambda_di <= 0)
    return(list(c_over_alpha = Inf, onset_ref = Inf, status = "no_repopulation"))
  if (lambda_dd <= 0)
    return(list(c_over_alpha = NA_real_, onset_ref = NA_real_, status = "infeasible"))
  y <- lambda_di * smooth_max0(T0 - Tk) / (lambda_dd * p_kill(K0, T0, alpha))
  x0 <- inv_smooth_max0(y)  # = T0 - onset time at the anchor
  onset <- T0 - x0
  if (!is.finite(onset) || onset <= 0)
    return(list(c_over_alpha = NA_real_, onset_ref = onset, status = "infeasible"))
  list(c_over_alpha = K0 * onset / T0, onset_ref = onset, status = "ok")
}

# Tumor BED, DD repopulation model (2018 paper, Eq. 2-5): the repopulation
# rate is lambda_dd times the mean fraction of cells killed per day, and the
# onset is reached at a fixed level of cumulative cell killing.
bed_dd <- function(m, d, T, alpha, ab, lambda_dd, c_over_alpha, g = 0) {
  K <- bed_lq(m, d, ab)
  if (is.na(c_over_alpha)) return(rep(NA_real_, length(K)))
  growth <- g * T / alpha
  if (is.infinite(c_over_alpha)) return(K - growth)  # no accelerated repopulation
  onset <- c_over_alpha * T / K
  K - lambda_dd * p_kill(K, T, alpha) * smooth_max0(T - onset) / alpha - growth
}

# DD onset day and repopulation rate for given regimens (diagnostics).
dd_onset_rate <- function(m, d, T, alpha, ab, lambda_dd, c_over_alpha) {
  K <- bed_lq(m, d, ab)
  list(onset = c_over_alpha * T / K, rate = lambda_dd * p_kill(K, T, alpha))
}

# Early (acute mucosal) toxicity BED: LQ with mucosal repopulation at rate
# lambda_e (= ln 2 / doubling time) starting at Tk_e. T must already be on
# the model's clock (see acute_model_time); compute_all() converts once.
bed_early <- function(m, d, T, ab_e, alpha_e, lambda_e, Tk_e) {
  bed_lq(m, d, ab_e) - lambda_e * smooth_max0(T - Tk_e) / alpha_e
}

# Late toxicity BED: LQ only, no repopulation. alpha does not enter.
bed_late <- function(m, d, ab_l) bed_lq(m, d, ab_l)

# Every quantity for one or more regimens (vectorized over m, d, T; T is the
# inclusive calendar count). p: list with alpha, ab, lambda_di, lambda_dd, Tk,
# g, ab_e, alpha_e, lambda_e, Tk_e, clock_e, ab_l. anchor: list(m, d, T), the
# regimen at which the DD model is calibrated.
compute_all <- function(m, d, T, p, anchor) {
  g <- if (is.null(p$g)) 0 else p$g
  clock_e <- if (is.null(p$clock_e)) "elapsed" else p$clock_e
  cal <- dd_calibration(anchor, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  early <- bed_early(m, d, acute_model_time(T, clock_e), p$ab_e, p$alpha_e, p$lambda_e, p$Tk_e)
  late <- bed_late(m, d, p$ab_l)
  data.frame(
    m = m, d = d, T = T, total_dose = m * d,
    bed_di = bed_di(m, d, T, p$alpha, p$ab, p$lambda_di, p$Tk, g),
    bed_dd = bed_dd(m, d, T, p$alpha, p$ab, p$lambda_dd, cal$c_over_alpha, g),
    bed_early = early, eqd2_early = eqd2(early, p$ab_e),
    bed_late = late, eqd2_late = eqd2(late, p$ab_l)
  )
}

# Fowler tolerance zone (Gy10), lowered with concurrent chemotherapy by the
# suggested 3-5 Gy10 (lower edge minus 5, upper edge minus 3).
fowler_zone <- function(chemo = FALSE) {
  if (chemo) FOWLER_ZONE - rev(FOWLER_CHEMO_REDUCTION) else FOWLER_ZONE
}

# Search for hypofractionated regimens, given once daily Monday to Friday.
# Dose per fraction runs over multiples of d_step above the comparison dose
# per fraction, up to d_max. Overall time must be at least min_days_saved
# shorter than the comparison regimen `ref`. The DD model is calibrated at
# `anchor`, never at `ref` unless the caller passes anchor = ref (legacy).
# A regimen is kept when
#   early BED <= early_cap (default: comparison early BED + tol_early),
#   late BED  <= late_cap  (default: comparison late BED + tol_late),
#   or, with fenwick = TRUE, cell-kill BED <= fenwick_limit(T) instead,
#   tumor BED change >= min_gain (negative = allowed loss) under the
#   criterion "either" (DI or DD), "both", "di" or "dd".
# A model whose comparison BED is negative or not finite cannot credit any
# regimen, and a candidate with a negative BED cannot qualify on it.
# For each number of fractions, the dose per fraction with the highest score
# is listed (rank_by: "best" = larger of dDI and dDD, "worst" = smaller,
# "di", "dd"). Column next_step says what the next dose step would break.
# All qualifying grid rows are attached as attribute "feasible".
# Empty results carry attribute "reason":
#   "dd_unavailable" / "negative_reference_di" / "negative_reference_dd" /
#   "negative_reference_both" /
#   "negative_reference_early" / "negative_reference_late": the criterion,
#     ranking or toxicity limit needs a model that is unavailable;
#   "empty_dose_range": no dose per fraction between the comparison dose and d_max;
#   "no_shorter_regimen": no fraction number saves enough days.
suggest_regimens <- function(p, ref, anchor = default_anchor(), m_range = 5:34,
                             d_step = 0.05, d_max = 8,
                             criterion = "either", rank_by = "best",
                             tol_early = 0, tol_late = 0, early_cap = NULL,
                             late_cap = NULL, min_gain = 0, min_days_saved = 1,
                             fenwick = FALSE) {
  eps <- 1e-9
  out <- empty_suggestions()
  blocked <- function(reason) { attr(out, "reason") <- reason; out }
  cal <- dd_calibration(anchor, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  r0 <- compute_all(ref$m, ref$d, ref$T, p, anchor)
  di_usable <- is.finite(r0$bed_di) && r0$bed_di > 0
  dd_usable <- cal$status != "infeasible" && is.finite(r0$bed_dd) && r0$bed_dd > 0
  needs_di <- criterion %in% c("both", "di") || rank_by %in% c("worst", "di")
  needs_dd <- criterion %in% c("both", "dd") || rank_by %in% c("worst", "dd")
  if (cal$status == "infeasible" && needs_dd) return(blocked("dd_unavailable"))
  if (!dd_usable && needs_dd) return(blocked("negative_reference_dd"))
  if (!di_usable && needs_di) return(blocked("negative_reference_di"))
  if (!di_usable && !dd_usable) return(blocked("negative_reference_both"))
  if (!fenwick && is.null(early_cap) && !(is.finite(r0$bed_early) && r0$bed_early > 0))
    return(blocked("negative_reference_early"))
  if (is.null(late_cap) && !(is.finite(r0$bed_late) && r0$bed_late > 0))
    return(blocked("negative_reference_late"))
  cap_early <- if (is.null(early_cap)) r0$bed_early + tol_early else early_cap
  cap_late <- if (is.null(late_cap)) r0$bed_late + tol_late else late_cap

  d_start <- ceiling((ref$d + 1e-9) / d_step) * d_step
  if (d_start > d_max + eps) return(blocked("empty_dose_range"))
  d_grid <- round(seq(d_start, d_max + eps, by = d_step), 6)
  g <- expand.grid(d = d_grid, m = m_range)
  g$T <- overall_time(g$m)
  g <- g[g$T <= ref$T - min_days_saved, , drop = FALSE]
  if (nrow(g) == 0) return(blocked("no_shorter_regimen"))

  evaluate <- function(m, d, T) {
    a <- compute_all(m, d, T, p, anchor)
    a$d_di <- a$bed_di - r0$bed_di
    a$d_dd <- a$bed_dd - r0$bed_dd
    a$d_early <- a$bed_early - r0$bed_early
    a$d_late <- a$bed_late - r0$bed_late
    a$qual_di <- di_usable & a$bed_di > 0 & a$d_di >= min_gain - eps
    a$qual_dd <- dd_usable & !is.na(a$d_dd) & a$bed_dd > 0 & a$d_dd >= min_gain - eps
    a$tumor_ok <- switch(criterion,
      either = a$qual_di | a$qual_dd, both = a$qual_di & a$qual_dd,
      di = a$qual_di, dd = a$qual_dd,
      stop("unknown criterion: ", criterion))
    a$early_ok <- if (fenwick) {
      lim <- fenwick_limit(T)
      !is.na(lim) & fenwick_cell_kill(m, d) <= lim + eps
    } else a$bed_early > 0 & a$bed_early <= cap_early + eps
    a$late_ok <- a$bed_late <= cap_late + eps
    a
  }
  a <- evaluate(g$m, g$d, g$T)
  a <- a[a$tumor_ok & a$early_ok & a$late_ok, , drop = FALSE]
  if (nrow(a) == 0) return(out)
  # A model that is not usable shows no difference at all (NA, not a number).
  if (!di_usable) a$d_di <- NA_real_
  if (!dd_usable) a$d_dd <- NA_real_

  # A model that is not usable does not enter the score.
  a$score <- switch(rank_by,
    best = pmax(a$d_di, a$d_dd, na.rm = TRUE),
    worst = pmin(a$d_di, a$d_dd, na.rm = TRUE),
    di = a$d_di,
    dd = a$d_dd,
    stop("unknown rank_by: ", rank_by))
  a$d_min <- pmin(a$d_di, a$d_dd, na.rm = TRUE)
  a$credited_by <- ifelse(a$qual_di & a$qual_dd, "DI and DD",
                          ifelse(a$qual_di, "DI only", "DD only"))
  a$beyond_precedent <- a$d > PRECEDENT_MAX_D + eps
  a$weeks_saved <- (ref$T - a$T) / 7
  feasible <- a[order(a$m, -a$score), setdiff(names(a), c("tumor_ok", "early_ok", "late_ok")),
                drop = FALSE]
  rownames(feasible) <- NULL

  a <- a[order(a$m, -a$score), , drop = FALSE]
  a <- a[!duplicated(a$m), , drop = FALSE]
  a <- a[order(-a$score), , drop = FALSE]

  # What would the next dose step break?
  nxt_d <- round(a$d + d_step, 6)
  nx <- evaluate(a$m, nxt_d, a$T)
  broken <- mapply(function(e, l, t) {
    parts <- c(if (!e) "early toxicity", if (!l) "late toxicity", if (!t) "tumor criterion")
    if (length(parts) == 0) "none (a higher dose ranks lower)" else paste(parts, collapse = " and ")
  }, nx$early_ok, nx$late_ok, nx$tumor_ok)
  a$next_step <- ifelse(nxt_d > d_max + eps, "end of dose search range", broken)
  a$tumor_ok <- NULL
  a$early_ok <- NULL
  a$late_ok <- NULL
  rownames(a) <- NULL
  a <- a[, names(empty_suggestions()), drop = FALSE]
  attr(a, "feasible") <- feasible
  a
}

empty_suggestions <- function() {
  data.frame(m = integer(), d = numeric(), T = numeric(), total_dose = numeric(),
             bed_di = numeric(), bed_dd = numeric(), bed_early = numeric(),
             eqd2_early = numeric(), bed_late = numeric(), eqd2_late = numeric(),
             d_di = numeric(), d_dd = numeric(), d_early = numeric(),
             d_late = numeric(), d_min = numeric(), qual_di = logical(),
             qual_dd = logical(), credited_by = character(),
             beyond_precedent = logical(), score = numeric(),
             next_step = character(), weeks_saved = numeric())
}

# Regimen at which the DD model is calibrated, unless the user changes it.
default_anchor <- function() list(m = 35, d = 2, T = 47)

# Comparison (control) regimens offered as presets. HYPNO (Bentzen et al.,
# IJROBP 2026): 55 Gy/20 fx at 5 fx/week and 66 Gy/33 fx at 6 fx/week, with
# overall times for a Monday start.
comparator_presets <- function() {
  list(
    hypno_test = list(m = 20, d = 2.75, T = 26,
                      label = "HYPNO test arm: 55 Gy in 20 fractions, 5 fractions/week (26 days)"),
    hypno_control = list(m = 33, d = 2, T = 38,
                         label = "HYPNO control arm: 66 Gy in 33 fractions, 6 fractions/week (38 days)"))
}

# Published regimens shown for comparison, each with the regimen it is
# compared against ("own_control" = the trial's own control arm, else the
# user's comparison regimen). Overall times are inclusive, Monday start.
benchmark_regimens <- function() {
  data.frame(
    name = c("HYPNO test arm: 55 Gy in 20 fractions, 5 fractions/week",
             "HYPNO control arm: 66 Gy in 33 fractions, 6 fractions/week",
             "HYHOPE: 46.5 Gy in 15 fractions (favorable-risk HPV-positive oropharynx)"),
    m = c(20, 33, 15), d = c(2.75, 2, 3.1), T = c(26, 38, 19),
    stringsAsFactors = FALSE)
}

# Withers et al. 1988: about 0.6 Gy per day, in 2-Gy fractions, compensates
# repopulation after about 4 weeks. In BED at tumor alpha/beta ab this is
# 0.6 * (1 + 2/ab) Gy per day, so lambda_di = alpha * 0.6 * (1 + 2/ab).
withers_lambda <- function(alpha, ab) alpha * 0.6 * (1 + 2 / ab)

# Tumor parameter presets (alpha in 1/Gy, ab in Gy, lambda in 1/day, Tk in
# inclusive days). Each preset sets only the fields it lists; a preset
# without lambda_dd leaves the DD maximum rate unchanged.
#   calculator: the version 1 website defaults (all five values)
#   withers:    a rule, lambda_di = withers_lambda(alpha, ab), with alpha 0.2
#               and Tk 28 days; lambda_di below is the value at ab = 10
#   hartley:    Hartley et al. 2011, alpha 0.3, ab 10, Tp 3 d, Tk 22 elapsed
#               days (23 inclusive)
#   hyhope:     HYHOPE protocol, alpha 0.35, ab 10, Tp 3 d, Tk 21 d
#   d18:        analyst scenario proposed 2026-10-18 (not a fitted set)
tumor_presets <- function() {
  list(
    calculator = list(alpha = 0.2, ab = 10, lambda_di = 0.2, lambda_dd = 0.5, Tk = 28),
    withers = list(alpha = 0.2, lambda_di = withers_lambda(0.2, 10), Tk = 28, rule = "withers"),
    hartley = list(alpha = 0.3, ab = 10, lambda_di = log(2) / 3, Tk = 23),
    hyhope = list(alpha = 0.35, ab = 10, lambda_di = log(2) / 3, Tk = 21),
    d18 = list(alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.40, Tk = 28))
}

# Numeric values a preset sets, given the current tumor alpha/beta (needed
# by the Withers rule).
resolve_preset <- function(name, ab = 10) {
  pr <- tumor_presets()[[name]]
  if (identical(pr$rule, "withers")) pr$lambda_di <- withers_lambda(pr$alpha, ab)
  pr$rule <- NULL
  pr
}

# Defaults (version 2.3, chosen by Igor Shuryak on 2026-10-19): tumor values
# are scenario D18; toxicity values are literature defaults (see
# References/README.md). The default early-toxicity limit is the acute BED of
# HYPNO's 55 Gy in 20 fractions.
default_params <- function() {
  list(alpha = 0.35, ab = 10, lambda_di = 0.23, lambda_dd = 0.40, Tk = 28, g = 0,
       ab_e = 10, alpha_e = 0.35, lambda_e = log(2) / 2.5, Tk_e = 7,
       clock_e = "elapsed", ab_l = 3)
}

# Version 1 tumor values with the same toxicity values (the "calculator" preset).
legacy_params <- function() modifyList(default_params(), resolve_preset("calculator"))

DEFAULT_EARLY_MODE <- "hypno"
DEFAULT_PRESET <- "d18"

default_ref <- function() list(m = 35, d = 2, T = 47)

# Fraction search range shown after the reference number of fractions
# changes: from 5 (or fewer, if the reference is very short) up to one
# fraction below the reference.
search_range_for_ref <- function(ref_m) {
  hi <- max(1, ref_m - 1)
  c(min(5, hi), hi)
}

# Name of the preset whose values equal the given ones, or "custom". A preset
# matches when every field it sets agrees (the Withers rule is evaluated at
# the given ab); fields it does not set are ignored.
preset_matching <- function(alpha, lambda_di, Tk, ab = 10, lambda_dd = NULL) {
  vals <- list(alpha = alpha, lambda_di = lambda_di, Tk = Tk, ab = ab, lambda_dd = lambda_dd)
  given <- unlist(vals)
  if (!is.numeric(given) || length(given) < 4 || any(!is.finite(given))) return("custom")
  for (nm in names(tumor_presets())) {
    pr <- resolve_preset(nm, ab)
    keys <- intersect(names(pr), names(vals))
    keys <- keys[!vapply(vals[keys], is.null, logical(1))]
    if (!all(c("alpha", "lambda_di", "Tk") %in% keys)) next
    if (isTRUE(all.equal(unlist(vals[keys]), unlist(pr[keys]), tolerance = 1e-6))) return(nm)
  }
  "custom"
}

# Deterministic sensitivity: tumor BED differences of one regimen from the
# comparison regimen, over the tumor presets and a few tumor alpha/beta
# values. DI-only presets keep the current lambda_dd. These are selected
# scenarios, not a probability distribution.
sensitivity_table <- function(reg, ref, anchor, p, ab_values = c(7, 10, 15)) {
  rows <- list()
  for (nm in names(tumor_presets())) for (ab in ab_values) {
    q <- modifyList(p, resolve_preset(nm, ab))
    q$ab <- ab
    cal <- dd_calibration(anchor, q$alpha, q$ab, q$lambda_di, q$lambda_dd, q$Tk)
    a <- compute_all(c(reg$m, ref$m), c(reg$d, ref$d), c(reg$T, ref$T), q, anchor)
    rows[[length(rows) + 1]] <- data.frame(
      preset = nm, ab = ab, alpha = q$alpha, lambda_di = q$lambda_di,
      lambda_dd = q$lambda_dd, Tk = q$Tk, dd_status = cal$status,
      d_di = a$bed_di[1] - a$bed_di[2], d_dd = a$bed_dd[1] - a$bed_dd[2],
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
