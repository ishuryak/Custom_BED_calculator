# Parameter uncertainty (version 3). Requires R/bed_models.R.
#
# Tumor parameters come as joint posterior draws from the Bayesian refit of
# the 2018 Table 1 data (fit/fit_refit.py, files in draws/). DI and DD are
# separate models with separate draws; they are never averaged. Toxicity
# parameters (acute mucosa, late tissue, organ alpha/beta) have no joint data
# and are drawn from stated, assumption-based priors (toxicity_prior()).
#
# Every regimen and the comparison are evaluated under the SAME draws, so a
# difference carries the parameter correlation. Matrices are regimens (rows)
# x draws (columns). A draw in which a BED is not positive is not evaluable
# for that model: it is counted, reported, and never counted as meeting a
# criterion.

UNCERTAINTY_VERSION <- "3.0"
DRAWS_DIR <- "draws"

# Range of dose per fraction in the tumor data the draws were fitted to,
# derived from the data file (not typed), so a change of data moves it.
tumor_data_support <- function(path = file.path("data", "shuryak2018_table1.csv")) {
  x <- utils::read.csv(path)
  list(d_min = min(x$dose_per_fraction_gy), d_max = max(x$dose_per_fraction_gy),
       T_min = min(x$total_time_days), T_max = max(x$total_time_days))
}

# Read a draws file and check its columns and SHA-256 against the manifest.
load_draws <- function(model = c("di", "dd"), variant = "primary", dir = DRAWS_DIR) {
  model <- match.arg(model)
  f <- file.path(dir, sprintf("tumor_%s_%s.csv", model, variant))
  man <- jsonlite::fromJSON(file.path(dir, "manifest.json"), simplifyVector = FALSE)
  key <- basename(f)
  bytes <- readBin(f, "raw", n = file.size(f))
  sha <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  if (is.null(man$files[[key]]) || !identical(man$files[[key]]$sha256, sha))
    stop("Draws file ", key, " does not match draws/manifest.json (SHA-256).")
  x <- utils::read.csv(f)
  need <- if (model == "di") c("alpha", "ab", "lambda", "Tk", "g") else
    c("alpha", "ab", "lambda", "c_over_alpha", "g")
  if (!all(need %in% names(x))) stop("Draws file ", key, " lacks columns: ", paste(setdiff(need, names(x)), collapse = ", "))
  if (any(!is.finite(as.matrix(x[need])))) stop("Draws file ", key, " has non-finite values.")
  x <- x[need]   # subset first: subsetting a data frame drops attributes
  attr(x, "sha256") <- sha
  attr(x, "file") <- key
  attr(x, "ess") <- man$files[[key]]$min_ess_bulk
  x
}

# --- vectorized BEDs: regimens (m, d, T vectors) x draws -------------------

bed_lq_mat <- function(m, d, ab) outer(m * d, rep(1, length(ab))) + outer(m * d^2, 1 / ab)

# DI tumor BED for every regimen and draw. dr: data frame with alpha, ab,
# lambda, Tk and g (one row per draw).
bed_di_draws <- function(m, d, T, dr) {
  K <- bed_lq_mat(m, d, dr$ab)
  hinge <- smooth_max0(outer(T, dr$Tk, "-"))
  K - sweep(hinge, 2, dr$lambda / dr$alpha, "*") - outer(T, dr$g / dr$alpha)
}

# DD tumor BED with the fitted onset constant c_over_alpha (= C/alpha) of each draw.
bed_dd_draws <- function(m, d, T, dr) {
  K <- bed_lq_mat(m, d, dr$ab)
  Tm <- matrix(T, nrow(K), ncol(K))
  pk <- -expm1(-sweep(K, 2, dr$alpha, "*") / Tm)
  onset <- sweep(Tm / K, 2, dr$c_over_alpha, "*")
  K - sweep(pk * smooth_max0(Tm - onset), 2, dr$lambda / dr$alpha, "*") -
    outer(T, dr$g / dr$alpha)
}

# Acute (Fowler-type) BED for every regimen and draw; T is inclusive and is
# converted to the elapsed clock as in compute_all().
bed_early_draws <- function(m, d, T, tx) {
  K <- bed_lq_mat(m, d, tx$ab_e)
  K - sweep(smooth_max0(outer(acute_model_time(T), tx$Tk_e, "-")), 2, tx$lambda_e / tx$alpha_e, "*")
}

bed_late_draws <- function(m, d, ab_l) bed_lq_mat(m, d, ab_l)

# --- toxicity priors (assumption-based) -------------------------------------

# Stated priors for parameters without joint data. Each is lognormal or
# uniform with a source note. Organ alpha/beta values feed the OAR module.
toxicity_prior <- function() {
  list(
    ab_e = list(dist = "lognormal", median = 10, sd_log = 0.354,
                note = "Acute mucosa α/β, 95% 5-20 Gy around Fowler's 10 Gy (assumption)."),
    alpha_e = list(dist = "uniform", lo = 0.30, hi = 0.40,
                   note = "Mucosal α: Hartley 2011 used 0.3/Gy, Fowler 2003 0.35/Gy (assumed range)."),
    Tp_e = list(dist = "uniform", lo = 2.0, hi = 4.0,
                note = "Mucosal doubling time: Fowler 2003 2.5 d; range assumed."),
    Tk_e = list(dist = "uniform", lo = 7, hi = 19,
                note = "Mucosal onset (elapsed days): Fowler 2003 day 7; Pardo-Montero 2021 refit near day 19 (as reported in the v3 plan review, to be re-checked)."),
    ab_l = list(dist = "lognormal", median = 3, sd_log = 0.207,
                note = "Late head and neck tissue α/β, 95% 2-4.5 Gy around 3 Gy (Stuschke 1999 reported 3.3-4.0 for grade 2+ late effects)."),
    ab_cord = list(dist = "lognormal", median = 1.32, sd_log = 0.212,
                   note = "Spinal cord and brainstem α/β, 95% 0.87-2 Gy: Kirkpatrick 2010 (QUANTEC) 0.87 Gy; HYPNO used 2 Gy."),
    ab_oar = list(dist = "lognormal", median = 3, sd_log = 0.207,
                  note = "Other organs at risk α/β, 95% 2-4.5 Gy around HYPNO's 3 Gy (assumption)."))
}

draw_from <- function(spec, n) {
  switch(spec$dist,
    lognormal = exp(log(spec$median) + spec$sd_log * stats::rnorm(n)),
    uniform = stats::runif(n, spec$lo, spec$hi),
    stop("unknown prior: ", spec$dist))
}

# n toxicity draws with a fixed seed and generator (recorded in exports).
toxicity_draws <- function(n, seed = 20261021, prior = toxicity_prior()) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  out <- data.frame(ab_e = draw_from(prior$ab_e, n), alpha_e = draw_from(prior$alpha_e, n),
                    Tp_e = draw_from(prior$Tp_e, n), Tk_e = draw_from(prior$Tk_e, n),
                    ab_l = draw_from(prior$ab_l, n), ab_cord = draw_from(prior$ab_cord, n),
                    ab_oar = draw_from(prior$ab_oar, n))
  out$lambda_e <- log(2) / out$Tp_e
  out
}

# Toxicity draws that reproduce one deterministic parameter set exactly
# (used for the "fixed toxicity" option and for tests).
fixed_toxicity_draws <- function(p, n) {
  data.frame(ab_e = rep(p$ab_e, n), alpha_e = p$alpha_e, Tp_e = log(2) / p$lambda_e,
             Tk_e = p$Tk_e, ab_l = p$ab_l, ab_cord = 2, ab_oar = 3, lambda_e = p$lambda_e)
}

# --- propagation -------------------------------------------------------------

# Tumor BED differences of regimens from a comparison under each model's
# draws. Returns matrices delta_di and delta_dd (NA where not evaluable) and
# the raw BEDs. ref: list(m, d, T).
tumor_deltas <- function(regs, ref, di, dd) {
  b_di <- bed_di_draws(regs$m, regs$d, regs$T, di)
  r_di <- bed_di_draws(ref$m, ref$d, ref$T, di)[1, ]
  b_dd <- bed_dd_draws(regs$m, regs$d, regs$T, dd)
  r_dd <- bed_dd_draws(ref$m, ref$d, ref$T, dd)[1, ]
  ok <- function(b, r) is.finite(b) & b > 0 & matrix(is.finite(r) & r > 0, nrow(b), ncol(b), byrow = TRUE)
  v_di <- ok(b_di, r_di); v_dd <- ok(b_dd, r_dd)
  d_di <- sweep(b_di, 2, r_di); d_di[!v_di] <- NA
  d_dd <- sweep(b_dd, 2, r_dd); d_dd[!v_dd] <- NA
  list(delta_di = d_di, delta_dd = d_dd, valid_di = v_di, valid_dd = v_dd,
       ref_di = r_di, ref_dd = r_dd)
}

# Toxicity check per regimen and draw. limits: list(early = list(mode, tol),
# late = list(mode, tol)); modes "hypno" (55 Gy/20 under the same draw) or
# "reference" (the comparison under the same draw, plus tol). oar: optional
# matrix (regimens x draws) of logical OAR status from oar_status_draws().
toxicity_ok_draws <- function(regs, ref, tx, limits, oar = NULL) {
  hy <- comparator_presets()$hypno_test
  base_early <- switch(limits$early$mode,
    hypno = bed_early_draws(hy$m, hy$d, hy$T, tx)[1, ],
    reference = bed_early_draws(ref$m, ref$d, ref$T, tx)[1, ] + limits$early$tol,
    stop("Under uncertainty the early limit must be 'hypno' or 'reference' (the Fowler and Fenwick limits are tied to fixed parameter values)."))
  base_late <- switch(limits$late$mode,
    hypno = bed_late_draws(hy$m, hy$d, tx$ab_l)[1, ],
    reference = bed_late_draws(ref$m, ref$d, tx$ab_l)[1, ] + limits$late$tol,
    stop("unknown late limit: ", limits$late$mode))
  e <- bed_early_draws(regs$m, regs$d, regs$T, tx)
  l <- bed_late_draws(regs$m, regs$d, tx$ab_l)
  ok_e <- e > 0 & sweep(e, 2, base_early + 1e-9, "<=")
  ok_l <- sweep(l, 2, base_late + 1e-9, "<=")
  out <- ok_e & ok_l
  if (!is.null(oar)) out <- out & oar
  list(ok = out, early_ok = ok_e, late_ok = ok_l)
}

q_or_na <- function(x, p) if (all(is.na(x))) NA_real_ else unname(stats::quantile(x, p, na.rm = TRUE, type = 7))

# Summary per regimen. threshold: minimum tumor BED change (Gy) under each
# model. Shares are over ALL draws: a draw that is not evaluable does not
# meet any criterion and is counted in n_invalid_*.
summarize_draws <- function(regs, td, tox = NULL, threshold = 0) {
  n <- ncol(td$delta_di)
  meet_di <- td$valid_di & td$delta_di >= threshold - 1e-9; meet_di[is.na(meet_di)] <- FALSE
  meet_dd <- td$valid_dd & td$delta_dd >= threshold - 1e-9; meet_dd[is.na(meet_dd)] <- FALSE
  both <- meet_di & meet_dd
  probs <- c(.025, .10, .5, .90, .975)
  qs <- function(x) t(apply(x, 1, function(r) if (all(is.na(r))) rep(NA_real_, length(probs)) else
    stats::quantile(r, probs, na.rm = TRUE, names = FALSE, type = 7)))
  qi <- qs(td$delta_di); qd <- qs(td$delta_dd)
  dmin <- pmin(td$delta_di, td$delta_dd)
  out <- data.frame(m = regs$m, d = regs$d, T = regs$T, n_draws = n,
    di_q025 = qi[, 1], di_q10 = qi[, 2], di_median = qi[, 3], di_q90 = qi[, 4], di_q975 = qi[, 5],
    dd_q025 = qd[, 1], dd_q10 = qd[, 2], dd_median = qd[, 3], dd_q90 = qd[, 4], dd_q975 = qd[, 5],
    min_q10 = apply(dmin, 1, function(r) if (all(is.na(r))) NA_real_ else
      stats::quantile(r, .10, na.rm = TRUE, names = FALSE, type = 7)),
    share_di = rowMeans(meet_di), share_dd = rowMeans(meet_dd), share_both = rowMeans(both),
    n_invalid_di = rowSums(!td$valid_di), n_invalid_dd = rowSums(!td$valid_dd))
  if (!is.null(tox)) {
    out$share_tox <- rowMeans(tox$ok)
    out$share_joint <- rowMeans(both & tox$ok)
  }
  out
}

# Tumor control under each draw, CONDITIONAL on a stated control rate tcp0 for
# the comparison: with TCP = exp(-K S), a BED difference scales S by
# exp(-alpha * delta), so TCP = tcp0 ^ exp(-alpha * delta). alpha is the refit's
# (population-level) value of the same draw.
tcp_given_control <- function(delta, alpha, tcp0) {
  stopifnot(tcp0 > 0, tcp0 < 1)
  tcp0^exp(-sweep(delta, 2, alpha, "*"))
}

# Display rules: shares to 5 percentage points, "<5%" and ">95%" at the ends.
fmt_share <- function(x) {
  ifelse(is.na(x), "not available",
    ifelse(x < 0.05, "<5%", ifelse(x > 0.95, ">95%", sprintf("%d%%", as.integer(round(x * 20) * 5)))))
}
# BEDs in summaries to 0.5 Gy.
round_half <- function(x) round(x * 2) / 2
fmt_half <- function(x, plus = TRUE) ifelse(is.na(x), "n/a", sprintf(if (plus) "%+.1f" else "%.1f", round_half(x)))

# Monte Carlo standard error of a share with effective sample size ess.
share_mcse <- function(p, ess) sqrt(p * (1 - p) / ess)

# --- search under uncertainty ------------------------------------------------

# Candidate grid like suggest_regimens(), on any calendar.
regimen_grid <- function(ref, m_range, d_step, d_max, min_days_saved, calendar = default_calendar()) {
  d_start <- ceiling((ref$d + 1e-9) / d_step) * d_step
  if (d_start > d_max + 1e-9) return(data.frame(m = integer(), d = numeric(), T = numeric()))
  d_grid <- round(seq(d_start, d_max + 1e-9, by = d_step), 6)
  g <- expand.grid(d = d_grid, m = m_range)
  Tm <- schedule_time(m_range, calendar$fx_per_week, calendar$start_weekday, calendar$skip_days)
  g$T <- Tm[match(g$m, m_range)]
  g <- g[g$T <= ref$T - min_days_saved, c("m", "d", "T"), drop = FALSE]
  rownames(g) <- NULL
  g
}

# Selection rules, stated before the search is run:
#   "min_q10"     highest 10th percentile of the smaller of dDI and dDD (per draw)
#   "share_joint" highest share of draws meeting the tumor threshold under both
#                 models AND every toxicity limit
SELECTION_RULES <- c("min_q10", "share_joint")
# The interactive search uses a fixed subset of the draws (every 4th, so 1000
# of 4000; Monte Carlo standard error of a share at most 1.6 points); the
# listed regimens are then summarized with every draw.
SEARCH_DRAWS <- 1000
search_subset <- function(x, n = SEARCH_DRAWS) {
  if (nrow(x) <= n) return(x)
  x[round(seq(1, nrow(x), length.out = n)), , drop = FALSE]
}

# Uncertainty-aware search. For each number of fractions it lists the dose per
# fraction chosen by `rule` among regimens whose share of draws meeting all
# toxicity limits is at least tox_share_min, and returns the full grid
# summary as attribute "grid" so neighbors can be shown. The per-draw maximum
# is never reported (it is optimistic by construction).
uncertain_search <- function(ref, di, dd, tx, limits, threshold = 0, rule = "min_q10",
                             m_range = 5:34, d_step = 0.05, d_max = 5, min_days_saved = 1,
                             calendar = default_calendar(), tox_share_min = 0.9, chunk = 250,
                             oar_fun = NULL) {
  stopifnot(rule %in% SELECTION_RULES)
  g <- regimen_grid(ref, m_range, d_step, d_max, min_days_saved, calendar)
  if (nrow(g) == 0) return(structure(data.frame(), reason = "empty_grid"))
  parts <- lapply(split(seq_len(nrow(g)), ceiling(seq_len(nrow(g)) / chunk)), function(ix) {
    r <- g[ix, , drop = FALSE]
    td <- tumor_deltas(r, ref, di, dd)
    oar <- if (is.null(oar_fun)) NULL else oar_fun(r, tx)
    tox <- toxicity_ok_draws(r, ref, tx, limits, oar)
    summarize_draws(r, td, tox, threshold)
  })
  s <- do.call(rbind, parts)
  rownames(s) <- NULL
  s$score <- if (rule == "min_q10") s$min_q10 else s$share_joint
  s$eligible <- s$share_tox >= tox_share_min - 1e-12 & is.finite(s$score)
  sup <- tumor_data_support()
  s$outside_data <- s$d > sup$d_max + 1e-9
  best <- s[s$eligible, , drop = FALSE]
  best <- best[order(best$m, -best$score, best$d), , drop = FALSE]
  best <- best[!duplicated(best$m), , drop = FALSE]
  best <- best[order(-best$score), , drop = FALSE]
  rownames(best) <- NULL
  attr(best, "grid") <- s
  attr(best, "rule") <- rule
  if (nrow(best) == 0) attr(best, "reason") <- "no_eligible_regimen"
  best
}

# The listed regimen's neighbors on the dose grid (one step either side).
grid_neighbors <- function(grid, m, d, d_step) {
  grid[grid$m == m & abs(grid$d - d) <= d_step + 1e-9, , drop = FALSE]
}
