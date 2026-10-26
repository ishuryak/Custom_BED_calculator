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

# 3.1 (2026-10-25 review): conservative score over the full ensemble with
# invalid draws as failures (M3); the screen's choice re-checked on all draws
# (M4); regular candidate calendar (S5); fixed toxicity uses the organ
# alpha/beta of the Organs at risk tab (S8); the organ policy is exported (S2).
# 3.2 (2026-10-26 review): draws regenerated after every conditional DD fit met the full gate
# (ESS included, R01) with an enforced identity chain (R02); organ mappings validated (R03).
UNCERTAINTY_VERSION <- "3.2"
DRAWS_MIN_VERSION <- "3.2"
DRAWS_DIR <- "draws"
UNC_VARIANT_KEYS <- c("primary", "paper_restrictions")
TOXICITY_SEED <- 20261021
TOXICITY_RNG <- "Mersenne-Twister/Inversion/Rejection"

draws_manifest <- function(dir = DRAWS_DIR)
  jsonlite::fromJSON(file.path(dir, "draws_manifest.json"), simplifyVector = FALSE)

# Range of dose per fraction in the tumor data the draws were fitted to,
# derived from the data file named in the draws manifest (not typed), so a
# change of data moves it.
tumor_data_support <- function(path = NULL, dir = DRAWS_DIR) {
  if (is.null(path)) {
    ds <- draws_manifest(dir)$dataset
    path <- file.path("data", if (is.null(ds)) "shuryak2018_table1.csv" else ds)
  }
  x <- utils::read.csv(path)
  list(d_min = min(x$dose_per_fraction_gy), d_max = max(x$dose_per_fraction_gy),
       T_min = min(x$total_time_days), T_max = max(x$total_time_days),
       n_arms = nrow(x), n_studies = length(unique(x$first_author)), file = basename(path))
}

# Read a draws file and check its columns and SHA-256 against the manifest.
load_draws <- function(model = c("di", "dd"), variant = "primary", dir = DRAWS_DIR) {
  model <- match.arg(model)
  f <- file.path(dir, sprintf("tumor_%s_%s.csv", model, variant))
  man <- draws_manifest(dir)
  key <- basename(f)
  bytes <- readBin(f, "raw", n = file.size(f))
  sha <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  if (is.null(man$files[[key]]) || !identical(man$files[[key]]$sha256, sha))
    stop("Draws file ", key, " does not match draws/draws_manifest.json (SHA-256).")
  # Draws from version 3.1 on carry their sampling diagnostics; a file whose fit
  # failed the release gate (divergences, R-hat, ESS, E-BFMI) is never used.
  if (!is.null(man$draws_version) && !isTRUE(man$files[[key]]$diagnostics$gate_pass))
    stop("Draws file ", key, " failed the sampling-diagnostic gate and is not used.")
  # From draws version 3.2 the recorded diagnostics are checked here against the gate itself,
  # ESS included, rather than trusting the gate_pass flag (review R01, 2026-10-26).
  # Draws older than 3.2 were accepted under a gate without ESS (review R01) and are refused.
  if (is.null(man$draws_version) || numeric_version(man$draws_version) < DRAWS_MIN_VERSION)
    stop("Draws file ", key, " is from draws version ", if (is.null(man$draws_version)) "unknown" else man$draws_version,
         "; this calculator needs ", DRAWS_MIN_VERSION, " or later (full diagnostic gate, ESS included).")
  {
    why <- draws_gate_failures(man$files[[key]]$diagnostics)
    if (length(why)) stop("Draws file ", key, " does not meet the sampling-diagnostic gate: ", paste(why, collapse = "; "), ".")
  }
  x <- utils::read.csv(f)
  need <- if (model == "di") c("alpha", "ab", "lambda", "Tk", "g") else
    c("alpha", "ab", "lambda", "c_over_alpha", "g")
  if (!all(need %in% names(x))) stop("Draws file ", key, " lacks columns: ", paste(setdiff(need, names(x)), collapse = ", "))
  if (any(!is.finite(as.matrix(x[need])))) stop("Draws file ", key, " has non-finite values.")
  x <- x[need]   # subset first: subsetting a data frame drops attributes
  attr(x, "sha256") <- sha
  attr(x, "file") <- key
  # Monte Carlo benchmark size: the number of distinct DI draws (3.1 files repeat
  # each DI draw across its paired DD draws), else the fit's minimum bulk ESS.
  attr(x, "ess") <- if (!is.null(man$files[[key]]$n_distinct_di_draws)) man$files[[key]]$n_distinct_di_draws else man$files[[key]]$min_ess_bulk
  x
}

# The diagnostic gate of fit/gate.py, applied to a manifest's recorded diagnostics (DI keys,
# or the per-fit minima of the conditional DD fits). Missing values fail.
DRAWS_GATE <- list(max_rhat = 1.01, min_ess = 400, min_ebfmi = 0.3)
draws_gate_failures <- function(dg) {
  # One finite number per diagnostic, or NA (fails): +Inf ESS or -Inf R-hat used to satisfy the
  # inequalities below, and a negative count is invalid, not better (review F02 of 2026-10-26).
  g <- function(...) {
    for (k in c(...)) if (!is.null(dg[[k]])) {
      x <- suppressWarnings(as.numeric(dg[[k]]))
      return(if (length(x) == 1 && !is.logical(dg[[k]]) && is.finite(x)) x else NA_real_)
    }
    NA_real_
  }
  v <- c(divergences = g("divergences"), treedepth_hits = g("treedepth_hits"), max_rhat = g("max_rhat"),
         ess_bulk = g("min_ess_bulk", "min_ess_bulk_per_fit"), ess_tail = g("min_ess_tail", "min_ess_tail_per_fit"),
         ebfmi = g("min_ebfmi"))
  ok <- c(divergences = v[["divergences"]] == 0, treedepth_hits = v[["treedepth_hits"]] == 0,
          max_rhat = v[["max_rhat"]] > 0 && v[["max_rhat"]] <= DRAWS_GATE$max_rhat, ess_bulk = v[["ess_bulk"]] >= DRAWS_GATE$min_ess,
          ess_tail = v[["ess_tail"]] >= DRAWS_GATE$min_ess, ebfmi = v[["ebfmi"]] >= DRAWS_GATE$min_ebfmi)
  ok[is.na(ok)] <- FALSE
  sprintf("%s %s", names(ok)[!ok], format(v[!ok], digits = 4))
}

# Unit of a tumor BED difference over draws: Gy at alpha/beta 10 only when every
# draw has alpha/beta 10; otherwise the threshold applies at each draw's alpha/beta.
draws_bed_unit <- function(dr) {
  if (all(abs(dr$ab - 10) < 1e-9)) "Gy\u2081\u2080" else "Gy, \u03b1/\u03b2 varies by draw"
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
                note = "Mucosal onset (elapsed days), analyst sensitivity prior: Fowler 2003 day 7; Pardo-Montero 2021 refitted Tk 19.44 inclusive days (about 18.4 elapsed) jointly with other parameters. Neither is a direct biological range."),
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
toxicity_draws <- function(n, seed = TOXICITY_SEED, prior = toxicity_prior()) {
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
# (used for the "fixed toxicity" option and for tests). organ_ab: the organ
# alpha/beta values (cord and brainstem, other organs) of the Organs at risk
# tab, so a fixed draw gives the tab's own organ check.
fixed_toxicity_draws <- function(p, n, organ_ab = c(cord = 2, oar = 3)) {
  data.frame(ab_e = rep(p$ab_e, n), alpha_e = p$alpha_e, Tp_e = log(2) / p$lambda_e,
             Tk_e = p$Tk_e, ab_l = p$ab_l, ab_cord = organ_ab[["cord"]], ab_oar = organ_ab[["oar"]],
             lambda_e = p$lambda_e)
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

# Lower quantile of the FULL ensemble, with not-evaluable draws (NA) counted
# as failures below every evaluable value. Type 7 interpolates between the
# sorted values at floor(h) and floor(h) + 1, h = (n - 1) p + 1. The result is
# defined only when both lie among the evaluable draws; otherwise it is NA
# ("unavailable"). Raw values, negative ones included, are never floored.
q_full_ensemble <- function(x, p) {
  n <- length(x); k <- sum(is.na(x))
  h <- (n - 1) * p + 1
  if (n == 0 || floor(h) <= k) return(NA_real_)
  v <- sort(x[!is.na(x)])
  lo <- floor(h) - k; hi <- min(lo + 1, length(v))
  v[lo] + (h - floor(h)) * (v[hi] - v[lo])
}

# Summary per regimen. threshold: minimum tumor BED change (Gy) under each
# model. Shares are over ALL draws: a draw that is not evaluable does not
# meet any criterion and is counted in n_invalid_*. min_q10, the conservative
# score, is the 10th percentile of min(dDI, dDD) over ALL draws, with a draw
# not evaluable under either model counted as a failure (q_full_ensemble);
# n_valid_both counts the draws evaluable under both models.
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
    min_q10 = apply(dmin, 1, q_full_ensemble, p = .10),
    share_di = rowMeans(meet_di), share_dd = rowMeans(meet_dd), share_both = rowMeans(both),
    n_invalid_di = rowSums(!td$valid_di), n_invalid_dd = rowSums(!td$valid_dd),
    n_valid_both = rowSums(td$valid_di & td$valid_dd))
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
# Share within the toxicity limits for the search table. At the ends, and within 5 points of the
# required share, it also gives the count; near the requirement it gives one decimal and says which
# side of it the row is on (review M03 of 2026-10-26: 88.7% and 90.3% both displayed as "90%").
fmt_tox_share <- function(x, n, req) {
  cnt <- round(x * n)
  near <- !is.na(x) & abs(x - req) < 0.05
  ends <- !is.na(x) & (x < 0.05 | x > 0.95) & !near
  side <- ifelse(x >= req - 1e-12, "meets", "below")
  ifelse(near, sprintf("%.1f%% (%d of %d; %s the required %g%%)", 100 * x, cnt, n, side, 100 * req),
         ifelse(ends, sprintf("%s (%d of %d)", fmt_share(x), cnt, n), fmt_share(x)))
}
# BEDs in summaries to 0.5 Gy.
round_half <- function(x) round(x * 2) / 2
fmt_half <- function(x, plus = TRUE) ifelse(is.na(x), "n/a", sprintf(if (plus) "%+.1f" else "%.1f", round_half(x)))

# Monte Carlo standard error of a share with effective sample size ess.
share_mcse <- function(p, ess) sqrt(p * (1 - p) / ess)

# --- search under uncertainty ------------------------------------------------

# Candidate grid like suggest_regimens(), on the regular candidate calendar
# (candidate_calendar(): skipped days are never applied to candidates).
regimen_grid <- function(ref, m_range, d_step, d_max, min_days_saved, calendar = default_calendar()) {
  calendar <- candidate_calendar(calendar)
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
# The interactive search screens with a fixed subset of the draws (every 4th,
# so 1000 of 4000); the screen's choice is then re-checked on every draw
# (run_uncertainty_search()).
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

# The search as the app runs it. 1) Screen the grid with SEARCH_DRAWS draws
# (uncertain_search). 2) For each number of fractions, re-check the
# toxicity-share requirement of the screen's choice on ALL draws; when it
# fails, the next-ranked dose of the screen is checked, and so on. 3) Summarize
# the listed doses, their grid neighbors and the failed choices with all draws.
# Returns list(res = listed rows (screen values), full = summary rows with
# column role "listed", "neighbor" or "failed_full_check").
run_uncertainty_search <- function(ref, di, dd, tx, limits, threshold = 0, rule = "min_q10",
                                   m_range = 5:34, d_step = 0.05, d_max = 5, min_days_saved = 1,
                                   calendar = default_calendar(), tox_share_min = 0.9, oar_fun = NULL) {
  screen <- uncertain_search(ref, search_subset(di), search_subset(dd), search_subset(tx), limits, threshold, rule,
                             m_range = m_range, d_step = d_step, d_max = d_max, min_days_saved = min_days_saved,
                             calendar = calendar, tox_share_min = tox_share_min, oar_fun = oar_fun)
  grid <- attr(screen, "grid")
  if (is.null(grid) || nrow(screen) == 0) return(list(res = screen, full = NULL))
  tox_all <- function(r) {
    oar <- if (is.null(oar_fun)) NULL else oar_fun(r, tx)
    rowMeans(toxicity_ok_draws(r, ref, tx, limits, oar)$ok)
  }
  cand <- grid[grid$eligible, , drop = FALSE]
  cand <- cand[order(cand$m, -cand$score, cand$d), , drop = FALSE]
  chosen <- list(); failed <- list()
  for (m in unique(cand$m)) {
    cm <- cand[cand$m == m, , drop = FALSE]
    for (j in seq_len(nrow(cm))) {
      if (tox_all(cm[j, c("m", "d", "T")]) >= tox_share_min - 1e-12) { chosen[[length(chosen) + 1]] <- cm[j, ]; break }
      failed[[length(failed) + 1]] <- cm[j, c("m", "d", "T")]
    }
  }
  res <- if (length(chosen)) do.call(rbind, chosen) else grid[0, , drop = FALSE]
  res <- res[order(-res$score), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "grid") <- grid; attr(res, "rule") <- rule
  if (nrow(res) == 0) { attr(res, "reason") <- "no_regimen_passes_all_draws" }
  rows <- list()
  for (i in seq_len(nrow(res))) {
    g <- grid_neighbors(grid, res$m[i], res$d[i], d_step)[, c("m", "d", "T")]
    g$role <- ifelse(abs(g$d - res$d[i]) < 1e-9, "listed", "neighbor")
    rows[[length(rows) + 1]] <- g
  }
  if (length(failed)) rows[[length(rows) + 1]] <- cbind(do.call(rbind, failed), role = "failed_full_check")
  nb <- do.call(rbind, rows)
  rank <- c(listed = 1, failed_full_check = 2, neighbor = 3)[nb$role]
  nb <- nb[order(nb$m, nb$d, rank), , drop = FALSE]
  nb <- nb[!duplicated(paste(nb$m, round(nb$d, 6))), , drop = FALSE]
  regs <- nb[, c("m", "d", "T")]
  rownames(regs) <- NULL
  td <- tumor_deltas(regs, ref, di, dd)
  tox <- toxicity_ok_draws(regs, ref, tx, limits, if (is.null(oar_fun)) NULL else oar_fun(regs, tx))
  full <- cbind(summarize_draws(regs, td, tox, threshold), role = nb$role, listed = nb$role == "listed",
                outside_data = nb$d > tumor_data_support()$d_max + 1e-9, stringsAsFactors = FALSE)
  list(res = res, full = full, n_failed = length(failed))
}
