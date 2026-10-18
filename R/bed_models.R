# BED models for the HNSCC BED calculator v2.
#
# Notation: m = number of fractions, d = dose per fraction (Gy),
# T = overall treatment time (days), ab = alpha/beta ratio (Gy),
# alpha = linear LQ parameter (1/Gy), lambda = repopulation rate (1/day).
#
# All max(0, x) terms use the same smooth version as the v1 website,
# log(1 + exp(3x)) / 3, so v1 and v2 agree at the v1 defaults.

SMOOTH_K <- 3  # sharpness (1/day) of the smooth max(0, x), as in v1

# Fowler et al. 2003 acute mucosal tolerance zone (Gy10, radiotherapy alone)
# and the suggested reduction with concurrent chemotherapy (Gy10). Valid only
# with Fowler's acute parameters (alpha 0.35, alpha/beta 10, Tk 7, Tp 2.5).
FOWLER_ZONE <- c(59, 61)
FOWLER_CHEMO_REDUCTION <- c(3, 5)

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

# Tumor BED, dose-independent (DI) repopulation model:
# repopulation at fixed rate lambda_di starting at fixed time Tk.
bed_di <- function(m, d, T, alpha, ab, lambda_di, Tk) {
  bed_lq(m, d, ab) - lambda_di * smooth_max0(T - Tk) / alpha
}

# Calibration of the dose-dependent (DD) model.
# Repopulation starts when -ln(surviving fraction) reaches C, at time
# C*T/(alpha*K) with K = bed_lq(m, d, ab). C is solved so that
# BED_DD = BED_DI for the reference regimen. Returns c_over_alpha = C/alpha
# (Gy), the onset day at the reference, and a status:
#   "ok"               calibration found, onset after treatment start
#   "no_repopulation"  DI repopulation rate is zero, so DD has none either
#   "infeasible"       no positive onset reproduces DI at the reference, so
#                      DD values are not defined and are returned as NA.
dd_calibration <- function(ref, alpha, ab, lambda_di, lambda_dd, Tk) {
  K0 <- bed_lq(ref$m, ref$d, ab)
  T0 <- ref$T
  if (lambda_di <= 0)
    return(list(c_over_alpha = Inf, onset_ref = Inf, status = "no_repopulation"))
  if (lambda_dd <= 0)
    return(list(c_over_alpha = NA_real_, onset_ref = NA_real_, status = "infeasible"))
  p_kill0 <- 1 - exp(-alpha * K0 / T0)
  y <- lambda_di * smooth_max0(T0 - Tk) / (lambda_dd * p_kill0)
  x0 <- inv_smooth_max0(y)  # = T0 - onset time at the reference
  onset <- T0 - x0
  if (!is.finite(onset) || onset <= 0)
    return(list(c_over_alpha = NA_real_, onset_ref = onset, status = "infeasible"))
  list(c_over_alpha = K0 * onset / T0, onset_ref = onset, status = "ok")
}

# Tumor BED, DD repopulation model: the repopulation rate is lambda_dd times
# the mean fraction of cells killed per day, and the onset is reached at a
# fixed level of cumulative cell killing (see dd_calibration).
bed_dd <- function(m, d, T, alpha, ab, lambda_dd, c_over_alpha) {
  K <- bed_lq(m, d, ab)
  if (is.na(c_over_alpha)) return(rep(NA_real_, length(K)))
  if (is.infinite(c_over_alpha)) return(K)  # no repopulation
  p_kill <- 1 - exp(-alpha * K / T)
  onset <- c_over_alpha * T / K
  K - lambda_dd * p_kill * smooth_max0(T - onset) / alpha
}

# Early (acute mucosal) toxicity BED: LQ with mucosal repopulation at rate
# lambda_e (= ln 2 / doubling time) starting at Tk_e.
bed_early <- function(m, d, T, ab_e, alpha_e, lambda_e, Tk_e) {
  bed_lq(m, d, ab_e) - lambda_e * smooth_max0(T - Tk_e) / alpha_e
}

# Late toxicity BED: LQ only, no repopulation. alpha does not enter.
bed_late <- function(m, d, ab_l) bed_lq(m, d, ab_l)

# Every quantity for one or more regimens (vectorized over m, d, T).
# p: list with alpha, ab, lambda_di, lambda_dd, Tk, ab_e, alpha_e, lambda_e,
#    Tk_e, ab_l. ref: list(m, d, T) used to calibrate the DD model.
compute_all <- function(m, d, T, p, ref) {
  cal <- dd_calibration(ref, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  early <- bed_early(m, d, T, p$ab_e, p$alpha_e, p$lambda_e, p$Tk_e)
  late <- bed_late(m, d, p$ab_l)
  data.frame(
    m = m, d = d, T = T, total_dose = m * d,
    bed_di = bed_di(m, d, T, p$alpha, p$ab, p$lambda_di, p$Tk),
    bed_dd = bed_dd(m, d, T, p$alpha, p$ab, p$lambda_dd, cal$c_over_alpha),
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
# Dose per fraction runs over multiples of d_step above the reference dose per
# fraction, up to d_max. Overall time must be at least min_days_saved shorter
# than the reference. A regimen is kept when
#   early BED <= early_cap (default: reference early BED + tol_early),
#   late BED  <= reference late BED + tol_late,
#   tumor BED gain >= min_gain under the chosen criterion:
#   "either" (DI or DD), "both", "di" or "dd".
# For each number of fractions, the dose per fraction with the highest score
# is listed (rank_by: "best" = larger of dDI and dDD, "worst" = smaller,
# "di", "dd"). Column next_step says what the next dose step would break.
# Empty results carry attribute "reason": "dd_unavailable" when a DD-based
# criterion or ranking ("both", "dd", rank "worst" or "dd") is requested but
# the DD calibration is infeasible; "empty_dose_range" when no dose per
# fraction lies between the reference dose and d_max; "no_shorter_regimen"
# when no fraction number saves enough days. With an infeasible DD model,
# criterion "either" and rank "best" use DI alone.
suggest_regimens <- function(p, ref, m_range = 5:34, d_step = 0.05, d_max = 8,
                             criterion = "either", rank_by = "best",
                             tol_early = 0, tol_late = 0, early_cap = NULL,
                             min_gain = 0, min_days_saved = 1) {
  eps <- 1e-9
  out <- empty_suggestions()
  cal <- dd_calibration(ref, p$alpha, p$ab, p$lambda_di, p$lambda_dd, p$Tk)
  if (cal$status == "infeasible" && (criterion %in% c("both", "dd") || rank_by %in% c("worst", "dd"))) {
    attr(out, "reason") <- "dd_unavailable"
    return(out)
  }
  r0 <- compute_all(ref$m, ref$d, ref$T, p, ref)
  cap_early <- if (is.null(early_cap)) r0$bed_early + tol_early else early_cap
  cap_late <- r0$bed_late + tol_late

  d_start <- ceiling((ref$d + 1e-9) / d_step) * d_step
  if (d_start > d_max + eps) {
    attr(out, "reason") <- "empty_dose_range"
    return(out)
  }
  d_grid <- round(seq(d_start, d_max + eps, by = d_step), 6)
  g <- expand.grid(d = d_grid, m = m_range)
  g$T <- overall_time(g$m)
  g <- g[g$T <= ref$T - min_days_saved, , drop = FALSE]
  if (nrow(g) == 0) {
    attr(out, "reason") <- "no_shorter_regimen"
    return(out)
  }

  evaluate <- function(m, d, T) {
    a <- compute_all(m, d, T, p, ref)
    a$d_di <- a$bed_di - r0$bed_di
    a$d_dd <- a$bed_dd - r0$bed_dd
    a$d_early <- a$bed_early - r0$bed_early
    a$d_late <- a$bed_late - r0$bed_late
    di_ok <- a$d_di >= min_gain - eps
    dd_ok <- !is.na(a$d_dd) & a$d_dd >= min_gain - eps
    a$tumor_ok <- switch(criterion,
      either = di_ok | dd_ok, both = di_ok & dd_ok, di = di_ok, dd = dd_ok,
      stop("unknown criterion: ", criterion))
    a$early_ok <- a$bed_early <= cap_early + eps
    a$late_ok <- a$bed_late <= cap_late + eps
    a
  }
  a <- evaluate(g$m, g$d, g$T)
  a <- a[a$tumor_ok & a$early_ok & a$late_ok, , drop = FALSE]
  if (nrow(a) == 0) return(out)

  a$score <- switch(rank_by,
    best = pmax(a$d_di, a$d_dd, na.rm = TRUE),
    worst = pmin(a$d_di, a$d_dd, na.rm = TRUE),
    di = a$d_di,
    dd = a$d_dd,
    stop("unknown rank_by: ", rank_by))
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
  a$weeks_saved <- (ref$T - a$T) / 7
  a$tumor_ok <- NULL
  a$early_ok <- NULL
  a$late_ok <- NULL
  rownames(a) <- NULL
  a
}

empty_suggestions <- function() {
  data.frame(m = integer(), d = numeric(), T = numeric(), total_dose = numeric(),
             bed_di = numeric(), bed_dd = numeric(), bed_early = numeric(),
             eqd2_early = numeric(), bed_late = numeric(), eqd2_late = numeric(),
             d_di = numeric(), d_dd = numeric(), d_early = numeric(),
             d_late = numeric(), score = numeric(), next_step = character(),
             weeks_saved = numeric())
}

# Published regimens shown for comparison. Overall times assume a Monday
# start: HYPNO (Bentzen et al., IJROBP 2026) 55 Gy/20 fx at 5 fx/week and
# 66 Gy/33 fx at 6 fx/week (Monday to Saturday).
benchmark_regimens <- function() {
  data.frame(
    name = c("HYPNO test arm: 55 Gy in 20 fractions, 5 fractions/week",
             "HYPNO control arm: 66 Gy in 33 fractions, 6 fractions/week"),
    m = c(20, 33), d = c(2.75, 2), T = c(26, 38))
}

# Tumor parameter presets (alpha in 1/Gy, lambda_di in 1/day, Tk in days).
#   calculator: the version 1 website defaults
#   withers:    0.6 Gy/day in 2-Gy fractions lost to repopulation after
#               about 4 weeks (Withers et al. 1988). In BED10 this is
#               0.6 * 1.2 = 0.72 Gy/day, so lambda_di/alpha = 0.72, with the
#               calculator's alpha
#   hyhope:     alpha 0.35, Tp 3 d, Tk 21 d (HYHOPE protocol)
tumor_presets <- function() {
  list(
    calculator = list(alpha = 0.2, lambda_di = 0.2, Tk = 28),
    withers = list(alpha = 0.2, lambda_di = 0.72 * 0.2, Tk = 28),
    hyhope = list(alpha = 0.35, lambda_di = log(2) / 3, Tk = 21))
}

# Defaults: tumor values are the v1 website defaults; toxicity values are
# literature defaults (see References/README.md).
default_params <- function() {
  list(alpha = 0.2, ab = 10, lambda_di = 0.2, lambda_dd = 0.5, Tk = 28,
       ab_e = 10, alpha_e = 0.35, lambda_e = log(2) / 2.5, Tk_e = 7,
       ab_l = 3)
}

default_ref <- function() list(m = 35, d = 2, T = 47)

# Fraction search range shown after the reference number of fractions
# changes: from 5 (or fewer, if the reference is very short) up to one
# fraction below the reference.
search_range_for_ref <- function(ref_m) {
  hi <- max(1, ref_m - 1)
  c(min(5, hi), hi)
}

# Name of the preset whose tumor values equal the given ones, or "custom".
preset_matching <- function(alpha, lambda_di, Tk) {
  if (!is.numeric(alpha) || !is.numeric(lambda_di) || !is.numeric(Tk) ||
      any(!is.finite(c(alpha, lambda_di, Tk)))) return("custom")
  pr <- tumor_presets()
  for (nm in names(pr)) {
    if (isTRUE(all.equal(c(alpha, lambda_di, Tk), c(pr[[nm]]$alpha, pr[[nm]]$lambda_di, pr[[nm]]$Tk),
                         tolerance = 1e-6))) return(nm)
  }
  "custom"
}
