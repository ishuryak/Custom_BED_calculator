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
V4_MARGINS <- list(tumor = 5, early = 5, late = 5)   # defaults; the main page lets the user change each
V4_MARGIN_RANGE <- c(0, 20)
# The three main-page limits (Gy of BED), each between V4_MARGIN_RANGE[1] and [2]. Stops with
# a message naming the first one that is missing or out of range (version 4.4.0).
v4_check_margins <- function(tumor, early, late) {
  x <- list(tumor = tumor, early = early, late = late)
  words <- c(tumor = "tumor BED loss", early = "early BED increase", late = "late BED increase")
  for (k in names(x)) {
    v <- x[[k]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < V4_MARGIN_RANGE[1] || v > V4_MARGIN_RANGE[2])
      stop(sprintf("The largest allowed %s must be a number from %g to %g Gy.", words[[k]],
                   V4_MARGIN_RANGE[1], V4_MARGIN_RANGE[2]), call. = FALSE)
    x[[k]] <- as.numeric(v)
  }
  x
}
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

# Contiguous runs of a sorted grid: list of c(lo, hi). Two values belong to one
# run when they are one grid step apart (review F01 of 2026-10-31: a future
# setting with a gap must not be shown as one unbroken range).
v4_runs <- function(x, step) {
  x <- sort(unique(round(x, 9)))
  if (length(x) == 0) return(list())
  brk <- which(diff(x) > step * 1.5)
  starts <- c(1, brk + 1); ends <- c(brk, length(x))
  lapply(seq_along(starts), function(k) c(x[starts[k]], x[ends[k]]))
}

# "2.65–2.80", or one value when the range has a single dose. With step, each
# contiguous run is shown separately ("2.10–2.20, 2.40–2.50").
v4_range <- function(x, digits = 2, step = NULL) {
  f <- paste0("%.", digits, "f")
  one <- function(lo, hi) if (abs(hi - lo) < 1e-9) sprintf(f, lo) else sprintf(paste0(f, "–", f), lo, hi)
  if (is.null(step)) return(one(min(x), max(x)))
  paste(vapply(v4_runs(x, step), function(r) one(r[1], r[2]), ""), collapse = ", ")
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

# A tumor limit as a signed threshold for display: "−3", and "0" rather than "−0".
v4_neg <- function(x) if (isTRUE(x == 0)) "0" else sprintf("\u2212%g", x)

# The plot's dashed limit lines, one per endpoint, each in its bars' color (version 4.4.0:
# the early and late limits can differ, so one shared "toxicity" line is no longer right).
V4_PLOT_COLOURS <- c("Tumor (DI model)" = "#0072B2", "Tumor (DD model)" = "#56B4E9",
                     "Early toxicity" = "#D62728", "Late toxicity" = "#009E73")
v4_plot_limits <- function(margins) {
  data.frame(limit = c("tumor", "early", "late"), y = c(0 - margins$tumor, margins$early, margins$late),
             colour = unname(V4_PLOT_COLOURS[c("Tumor (DI model)", "Early toxicity", "Late toxicity")]))
}
# units: the BED units of tumor, early and late, for example c("Gy₁₀", "Gy₁₀", "Gy₃").
v4_plot_caption <- function(units, margins) sprintf(paste(
  "Bars span the listed doses per fraction. Tumor in %s, early toxicity in %s, late toxicity in %s.",
  "Dashed lines, in the color of their bars: tumor limit %s, early limit +%g, late limit +%g."),
  units[1], units[2], units[3], v4_neg(margins$tumor), margins$early, margins$late)

# Settings copied to the Report card (and so to Results) for one regimen: comparison
# 35 x 2 Gy, tumor threshold -margins$tumor, early and late limits equal to the comparison
# plus margins$early and margins$late (5, 5 and 5 by default),
# 5 fractions a week from a Monday. Search, the Both-model screen, Uncertainty and
# Organs at risk keep their own criteria (review F04 of 2026-10-31; V4_TRANSFER_SCOPE).
v4_report_inputs <- function(m, d, margins = V4_MARGINS) {
  list(m = m, d = d, fx_week = "5", start_day = "1", skip_days = "", auto_T = TRUE,
       comparator = "custom", ref_m = V4_STANDARD$m, ref_d = V4_STANDARD$d, ref_fx_week = "5", ref_auto_T = TRUE,
       min_gain = 0 - margins$tumor, early_mode = "reference", tol_early = margins$early,
       late_mode = "reference", tol_late = margins$late)
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
    step <- V4_GRID$d_step
    models <- if (all(both)) "both models" else if (any(both))
      sprintf("both at %s Gy, else %s only", v4_range(x$d[both], step = step), one) else sprintf("%s only", one)
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
      if (any(low)) sprintf("support < %d%% at %s Gy", round(100 * support_min), v4_range(x$d[low], step = step)),
      if (nzchar(trial)) trial)
    data.frame(m = x$m[1], dose = v4_range(x$d, step = step), total = v4_range(x$m * x$d, 1, step = x$m[1] * step),
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

# Refit draws for the support column. A successful load is kept for the life of
# the R process. A failed load returns NULL (the column then says "not
# available") and is retried at most once every V4_DRAWS_RETRY_SEC seconds, so a
# transient failure does not last until a restart and a persistent one is not
# retried on every recalculation (review F12 of 2026-10-31). If the deployed
# draw files are ever replaced, restart the app so the new files are read.
V4_DRAWS_CACHE <- new.env()
V4_DRAWS_RETRY_SEC <- 300
v4_draws <- function(loader = function() list(di = load_draws("di"), dd = load_draws("dd")),
                     now = Sys.time(), retry_sec = V4_DRAWS_RETRY_SEC, cache = V4_DRAWS_CACHE) {
  if (isTRUE(cache$ok)) return(cache$x)
  if (!is.null(cache$failed_at) &&
      as.numeric(difftime(now, cache$failed_at, units = "secs")) < retry_sec) return(NULL)
  x <- tryCatch(loader(), error = function(e) { cache$error <- conditionMessage(e); NULL })
  if (is.null(x)) { cache$failed_at <- now; return(NULL) }
  cache$x <- x; cache$ok <- TRUE; cache$failed_at <- NULL; cache$error <- NULL
  x
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

# Tumor endpoints of one or more schedules against the standard, under the
# same validity rule as the search (tumor_endpoint(): the model must be
# available and both BEDs positive). x: compute_all() rows.
v4_tumor_endpoints <- function(x, p, anchor, r0) {
  dd_available <- dd_setup(p, anchor)$status != "infeasible"
  list(di = tumor_endpoint(x$bed_di, r0$bed_di), dd = tumor_endpoint(x$bed_dd, r0$bed_dd, dd_available))
}

# Why one tumor model cannot credit a schedule, in words.
v4_endpoint_words <- function(e, name, k = 1) {
  switch(e$reason[k],
    dd_unavailable = sprintf("%s model not available", name),
    negative_reference = sprintf("%s BED of the standard not positive, so %s cannot be used", name, name),
    negative_candidate = sprintf("%s BED not positive", name),
    sprintf("%s change %s Gy", name, v4_signed_range(e$delta[k])))
}

# Which main-page conditions a schedule fails, compared with the standard
# (used to explain where a row's dose range starts and stops). Empty when it
# meets all three. r0: compute_all() of the standard. Uses the search's own
# validity rules (review F01 of 2026-10-31): a model whose standard BED is not
# positive, or that is unavailable, never credits a schedule.
v4_failed_conditions <- function(m, d, p, anchor, r0, margins = V4_MARGINS) {
  x <- compute_all(m, d, overall_time(m), p, anchor)
  e <- v4_tumor_endpoints(x, p, anchor, r0)
  meets <- function(en) isTRUE(en$valid[1]) && en$delta[1] >= -margins$tumor - 1e-9
  c(if (!(meets(e$di) || meets(e$dd)))
      sprintf("tumor BED not within %g Gy of the standard under either model: %s; %s", margins$tumor,
              v4_endpoint_words(e$di, "DI"), v4_endpoint_words(e$dd, "DD")),
    if (!(is.finite(x$bed_early) && x$bed_early > 0)) "early BED not positive"
    else if (x$bed_early > r0$bed_early + margins$early + 1e-9)
      sprintf("early BED more than %g Gy above the standard", margins$early),
    if (x$bed_late > r0$bed_late + margins$late + 1e-9)
      sprintf("late BED more than %g Gy above the standard", margins$late))
}

# Lowest and highest value of a BED change over a row's doses, with the dose at
# which each occurs. Under the DD model the change need not rise with dose, so
# the extremes are not assumed to sit at the row's lowest and highest dose.
v4_extremes_text <- function(v, d, name) {
  ok <- is.finite(v)
  if (!any(ok)) return(sprintf("The %s change is not available for these doses.", name))
  v <- v[ok]; d <- d[ok]
  lo <- which.min(v); hi <- which.max(v)
  if (abs(v[hi] - v[lo]) < 1e-9)
    return(sprintf("The %s change is %s Gy at every dose.", name, v4_signed_range(v[lo])))
  sprintf("The %s change runs from %s Gy (at %.2f Gy) to %s Gy (at %.2f Gy).", name,
          v4_signed_range(v[lo]), d[lo], v4_signed_range(v[hi]), d[hi])
}

# Plain-language worked example for one row of the main list: the doses that
# qualify, why the next dose below and above does not, and where each tumor
# model's lowest and highest change occur.
v4_range_example <- function(x, p, anchor, r0, margins = V4_MARGINS, grid = V4_GRID) {
  m <- x$m[1]; lo <- min(x$d); hi <- max(x$d)
  below <- round(lo - grid$d_step, 6); above <- round(hi + grid$d_step, 6)
  why <- function(dd) {
    f <- v4_failed_conditions(m, dd, p, anchor, r0, margins)
    if (length(f) == 0) sprintf("%.2f Gy also meets them (it belongs to another listed range)", dd)
    else sprintf("%.2f Gy fails (%s)", dd, paste(f, collapse = "; "))
  }
  why_lo <- if (below <= 2 + 1e-9) "the search starts above 2 Gy per fraction" else why(below)
  why_hi <- if (above > grid$d_max + 1e-9) sprintf("the search stops at %g Gy per fraction", grid$d_max) else why(above)
  paste(sprintf("For example, with %d fractions the doses %s Gy meet all three conditions. Below them, %s; above them, %s.",
                m, v4_range(x$d, step = grid$d_step), why_lo, why_hi),
        if (any(is.finite(x$d_di))) sub("\\.$", sprintf(", so the row shows %s.", v4_signed_range(x$d_di)),
                                        v4_extremes_text(x$d_di, x$d, "tumor DI"))
        else v4_extremes_text(x$d_di, x$d, "tumor DI"),
        v4_extremes_text(x$d_dd, x$d, "tumor DD"))
}

# Units and the 2-Gy-fraction equivalent of the tumor margin for the current
# tumor alpha/beta (review F02 of 2026-10-31): five BED units are 4.2 Gy EQD2 at
# alpha/beta 10 Gy but 2.8 Gy at 2.5 Gy.
v4_margin_eqd2 <- function(margin, ab) margin / (1 + 2 / ab)

# What the two tumor models are in the current settings: calibrated DD at the
# actual anchor, directly entered DD, zero DI rate, or DD unavailable.
v4_model_relation <- function(p, anchor) {
  cal <- dd_setup(p, anchor)
  anc <- sprintf("%g × %g Gy in %g days", anchor$m, anchor$d, anchor$T)
  is_std <- isTRUE(all.equal(unlist(anchor[c("m", "d", "T")]), unlist(V4_STANDARD[c("m", "d", "T")])))
  switch(cal$status,
    ok = paste0("The DD model shares α and the background growth rate with DI, and its onset constant is",
                " solved so that both models give the same tumor BED at the anchor, ", anc,
                if (is_std) " (the standard)" else " (not the standard, so the two models differ at 35 × 2 Gy)",
                ". Their agreement is therefore not an independent confirmation. This calibrated DD model is a",
                " variant of the published dose-dependent model, not its independent fit."),
    no_repopulation = paste("The DI repopulation rate is zero, so the calibrated DD model has no accelerated",
                            "repopulation either and the two models give the same tumor BED."),
    fitted = sprintf(paste("The DD model uses its own entered parameters (α %g per Gy, C/α %g Gy, growth %g per day)",
                           "and is not calibrated to DI, so the two models need not agree at the standard or any other schedule."),
                     cal$alpha, cal$c_over_alpha, cal$g),
    "The DD model is not available with these settings, so only the DI model is used.")
}

# Which tumor models can be used against the standard (same rule as the search).
v4_usable_models <- function(p, anchor, r0) {
  e <- v4_tumor_endpoints(data.frame(bed_di = 1, bed_dd = 1), p, anchor, r0)
  list(di = e$di$ref_valid, dd = e$dd$ref_valid, dd_available = e$dd$reason[1] != "dd_unavailable")
}

# Row index and dose sent by the browser, checked against the current list
# (review F05 of 2026-10-31). Returns list(m, d) of a listed point, or NULL.
v4_valid_row <- function(i, n) {
  i <- suppressWarnings(as.numeric(i))
  length(i) == 1 && is.finite(i) && i >= 1 && i == round(i) && i <= n
}
v4_validate_transfer <- function(i, dose, s, points, tol = 1e-6) {
  if (!v4_valid_row(i, nrow(s))) return(NULL)
  dose <- suppressWarnings(as.numeric(dose))
  if (length(dose) != 1 || !is.finite(dose)) return(NULL)
  m <- s$m[i]
  d <- points$d[points$m == m]
  k <- which(abs(d - dose) < tol)
  if (length(k) != 1) return(NULL)
  list(m = m, d = d[k])
}

# What "Open on the Report card" copies, and what other pages keep (review F04).
v4_transfer_scope <- function(margins = V4_MARGINS) paste(
  "This copies the regimen, the comparison 35 × 2 Gy in 47 days (5 fractions a week from a Monday), the tumor",
  sprintf("threshold %s and the early and late limits (standard + %g and standard + %g) to Check a regimen,",
          v4_neg(margins$tumor), margins$early, margins$late), "which its All BED values section also uses.",
  "The custom search with its both-model screen, Uncertainty and Organs at risk keep their own additional criteria (for example",
  "both tumor models, a stricter ranking or organ limits), so they can be stricter and need not list this regimen.")

# The trials the Bayesian refit of the DI and DD tumor models was fitted to (version 4.4.0; checks and
# published-outcome columns added 2026-11-09). Doses, fraction numbers, times and the published outcome
# come from the fitting file named in the draws manifest, so the list moves with the data; references and
# drugs come from data/fit_trial_list.csv (tools/build_fit_trial_list.py). Before joining, both files must
# have unique nonempty (ref_no, arm) keys, the same key set and the same PMID for every key, and the fitting
# file must match the SHA-256 the manifest records. Any failure is an error of class "v4_fit_trials_error",
# so the dialog can say the list is unavailable (v4_fit_trials_dialog) instead of showing a wrong citation.
V4_FIT_TRIAL_LIST <- "data/fit_trial_list.csv"
v4_fit_trials_fail <- function(...)
  stop(structure(class = c("v4_fit_trials_error", "error", "condition"),
                 list(message = paste0(...), call = NULL)))
v4_fit_trials <- function(dir = DRAWS_DIR, list_path = V4_FIT_TRIAL_LIST, data_dir = "data") {
  fail <- v4_fit_trials_fail
  man <- draws_manifest(dir)
  ds <- man$dataset
  if (is.null(ds)) fail("The draws manifest names no fitting data file.")
  path <- file.path(data_dir, ds)
  if (!file.exists(path)) fail("The fitting data file ", ds, " named in the draws manifest is missing.")
  if (!is.null(man$data_sha256)) {
    sha <- digest::digest(readBin(path, "raw", n = file.size(path)), algo = "sha256", serialize = FALSE)
    if (!identical(sha, man$data_sha256))
      fail("The fitting data file ", ds, " differs from the one the draws manifest records (SHA-256 ", sha,
           ", manifest ", man$data_sha256, ").")
  }
  read <- function(f, cc) tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, colClasses = cc),
                                   error = function(e) fail("Cannot read ", f, ": ", conditionMessage(e)),
                                   warning = function(w) fail("Cannot read ", f, ": ", conditionMessage(w)))
  x <- read(path, c(pmid = "character", source_value_pct = "character"))
  ref <- read(list_path, c(pmid = "character"))
  for (s in list(list(d = x, name = path, need = c("ref_no", "arm", "pmid", "source_value_pct", "tcp_actual_pct")),
                 list(d = ref, name = list_path, need = c("ref_no", "arm", "pmid", "trial", "citation")))) {
    gone <- setdiff(s$need, names(s$d))
    if (length(gone)) fail("Missing column(s) ", paste(gone, collapse = ", "), " in ", s$name, ".")
  }
  key <- function(d) paste(trimws(d$ref_no), trimws(d$arm), sep = " / ")
  for (s in list(list(d = x, name = paste("the fitting data", ds)), list(d = ref, name = list_path))) {
    blank <- which(is.na(s$d$ref_no) | !nzchar(trimws(s$d$ref_no)) | is.na(s$d$arm) | !nzchar(trimws(s$d$arm)))
    if (length(blank)) fail("Empty ref_no or arm in row(s) ", paste(blank, collapse = ", "), " of ", s$name, ".")
    dup <- unique(key(s$d)[duplicated(key(s$d))])
    if (length(dup)) fail("Repeated arm in ", s$name, ": ", paste(dup, collapse = "; "), ".")
  }
  miss <- setdiff(key(x), key(ref)); extra <- setdiff(key(ref), key(x))
  if (length(miss) || length(extra))
    fail("The fitted trial list does not match the fitting data ", ds, ".",
         if (length(miss)) paste0(" Arms missing from ", list_path, ": ", paste(miss, collapse = "; "), "."),
         if (length(extra)) paste0(" Arms not in the fitting data: ", paste(extra, collapse = "; "), "."),
         " Rebuild ", list_path, ".")
  r <- ref[match(key(x), key(ref)), , drop = FALSE]
  norm_pmid <- function(p) sub("^0+", "", sub("\\.0+$", "", trimws(as.character(p))))
  bad <- which(is.na(r$pmid) | norm_pmid(x$pmid) != norm_pmid(r$pmid))
  if (length(bad))
    fail("PMID differs between the fitting data ", ds, " and ", list_path, " for: ",
         paste(sprintf("%s (data %s, list %s)", key(x)[bad], x$pmid[bad], r$pmid[bad]), collapse = "; "), ".")
  # The fitted rate must be the published value, or its complement for a reported failure rate.
  pub <- suppressWarnings(as.numeric(x$source_value_pct))
  is_fail <- toupper(as.character(x$source_value_is_failure)) == "TRUE"
  bad <- which(is.na(pub) | abs(ifelse(is_fail, 100 - pub, pub) - x$tcp_actual_pct) > 1e-9)
  if (length(bad))
    fail("Published and fitted outcome disagree in the fitting data ", ds, " for: ",
         paste(sprintf("%s (published %s%s, fitted %g)", key(x)[bad], x$source_value_pct[bad],
                       ifelse(is_fail[bad], " failure", ""), x$tcp_actual_pct[bad]), collapse = "; "), ".")
  fx <- x$n_fractions
  fc <- x$fraction_components; fc[is.na(fc)] <- ""
  y <- round(x$n_patients * x$tcp_actual_pct / 100)
  endpoint <- trimws(sub("\\s*--.*$", "", x$endpoint))
  endpoint <- ifelse(is_fail, paste("freedom from", trimws(sub("\\(cumulative incidence\\)", "",
                                                                 sub("FAILURE", "failure", endpoint, fixed = TRUE)))),
                     endpoint)
  out <- data.frame(
    ref_no = x$ref_no, arm_key = x$arm,
    trial = r$trial, pmid = r$pmid, doi = r$doi, citation = r$citation,
    arm = ifelse(is.na(r$arm_label) | !nzchar(r$arm_label), x$arm, r$arm_label),
    total = sprintf("%.1f", x$total_dose_gy),
    dose = ifelse(nzchar(fc), gsub(";", " + ", gsub("x", " × ", fc)),
                  sprintf("%g", round(x$dose_per_fraction_gy, 4))),
    fractions = ifelse(abs(fx - round(fx)) < 1e-9, sprintf("%d", as.integer(round(fx))), sprintf("≈ %.1f", fx)),
    time = sprintf("%g", round(x$total_time_days, 2)),
    time_basis = x$time_rule,
    n = x$n_patients,
    published = ifelse(is_fail, sprintf("%s (from %s%% failure)", sprintf("%g", round(100 - pub, 4)), x$source_value_pct),
                       x$source_value_pct),
    published_is_failure = is_fail,
    endpoint = endpoint,
    horizon = sprintf("%s y", x$timepoint_years),
    estimator = ifelse(x$estimate_type == "actuarial/KM", "actuarial (Kaplan–Meier)", x$estimate_type),
    figure_read = x$value_source == "figure_read",
    outcome_locator = x$location,
    fit_count = sprintf("%d / %d", as.integer(y), as.integer(x$n_patients)),
    fit_rate = sprintf("%.1f", 100 * y / x$n_patients),
    outcome_source = if (is.null(r$outcome_source)) "" else ifelse(is.na(r$outcome_source), "", r$outcome_source),
    schedule_source = if (is.null(r$schedule_source)) "" else ifelse(is.na(r$schedule_source), "", r$schedule_source),
    chemotherapy = r$chemotherapy, other_drugs = r$other_drugs,
    shuryak2018_ref = ifelse(is.na(r$shuryak2018_ref), "", as.character(r$shuryak2018_ref)),
    trial_note = ifelse(is.na(r$trial_note), "", r$trial_note),
    stringsAsFactors = FALSE)
  out[order(out$trial, seq_len(nrow(out))), , drop = FALSE]
}

# The dialog body for the fitted trials, or an explicit "unavailable" message when the metadata checks of
# v4_fit_trials fail, so the page stays usable and never shows a citation it could not verify.
v4_fit_trials_dialog <- function(explainer_href, ...) {
  t <- tryCatch(v4_fit_trials(...), v4_fit_trials_error = function(e) e)
  if (inherits(t, "v4_fit_trials_error"))
    return(list(ok = FALSE, title = "Trials the tumor models were refitted to (list unavailable)",
                body = tags$div(class = "alert alert-warning",
                                tags$p("The list of fitted trials is unavailable because its references could not be",
                                       "checked against the fitting data:"),
                                tags$p(class = "small", conditionMessage(t)))))
  list(ok = TRUE, t = t,
       title = sprintf("Trials the tumor models were refitted to (%d trials, %d arms or strata)",
                       length(unique(t$pmid)), nrow(t)),
       body = tagList(v4_fit_trials_intro(t, explainer_href), v4_fit_trials_table(t)))
}

# Footnotes of the fitted-trial list: a lettered note for each trial with a note on drugs or on
# the schedule values used in the fit (data/fit_trial_list.csv, trial_note), and * for trials
# that were not in Shuryak, Hall and Brenner 2018 and 2019. The 2019 paper used the 2018 data
# set (the same 16 publications, 7283 patients, 35 regimens; its Methods and references 5, 6 and
# 12-25, read 2026-11-06), so one list serves both. Returns list(marks = named by trial, notes).
V4_SHURYAK_2018 <- "Shuryak, Hall and Brenner, Radiother Oncol 2018;127:20\u201326"
V4_SHURYAK_2019 <- "Shuryak, Hall and Brenner, Int J Radiat Oncol Biol Phys 2019;104:272\u2013278"
v4_fit_trial_notes <- function(t) {
  trials <- unique(t$trial)
  note <- vapply(trials, function(tr) {
    x <- unique(t$trial_note[t$trial == tr & nzchar(t$trial_note)])
    if (length(x)) x[1] else ""
  }, character(1))
  new <- vapply(trials, function(tr) !any(nzchar(t$shuryak2018_ref[t$trial == tr])), logical(1))
  k <- 0L; marks <- setNames(character(length(trials)), trials); notes <- list()
  for (i in seq_along(trials)) if (nzchar(note[[i]])) {
    k <- k + 1L; marks[[i]] <- letters[k]
    notes[[length(notes) + 1]] <- list(mark = letters[k], text = paste0(trials[i], ". ", note[[i]]))
  }
  if (any(new)) {
    marks[new] <- paste0(marks[new], ifelse(nzchar(marks[new]), ",", ""), "*")
    notes[[length(notes) + 1]] <- list(mark = "*", text = paste0(
      "Added in the refit: not among the trials of the earlier fits (", V4_SHURYAK_2018, ", and ", V4_SHURYAK_2019,
      "). Both papers used the same 16 publications. Their other trials are all in the refit except five that",
      " were dropped (not randomized, full text not available, or control not reported by arm), as described",
      " in the refit explainer."))
  }
  list(marks = marks, notes = notes)
}

# The dialog's opening text: what the trials are, what they were used for, and how to read
# the schedule columns.
v4_fit_trials_intro <- function(t, explainer_href) {
  n_trials <- length(unique(t$pmid))
  no_chemo <- all(startsWith(t$chemotherapy, "none"))
  tagList(
    tags$p(sprintf(paste(
      "The Bayesian refit of the DI and DD tumor models used these %d arms or strata of %d randomized trials of",
      "radiotherapy for head and neck cancer.%s Its parameter samples give the Refit support column and the refit",
      "tumor preset (Model settings). The default tumor values (scenario D18) were chosen by the developer and were",
      "not fitted to outcome data."), nrow(t), n_trials,
      if (no_chemo) " No arm included chemotherapy with radiotherapy." else "")),
    tags$p(class = "small text-muted", HTML(paste0(
      "Doses, fraction numbers and overall times are the values entered in the fit. A fraction number marked \u2248 ",
      "is an average over the arm (mean delivered dose, or a mix of prescriptions). Time basis says whether the ",
      "overall time is planned or delivered, as reported by each trial, and trials count planned days in different ",
      "ways. Published outcome is the rate each trial reported for the arm, with its own endpoint, time point, ",
      "estimator and page, table or figure. Count in fit is the rounded number of patients counted as controlled ",
      "that the fit used in place of that rate, an approximation for Kaplan\u2013Meier and cumulative-incidence estimates. ",
      "Where the schedule values were checked in a companion report, it is named under the reference. ",
      "Lettered notes give the schedule details that the columns cannot show. How the trials were chosen and ",
      "read is in the <a href=\"", htmltools::htmlEscape(explainer_href, attribute = TRUE),
      "\" target=\"_blank\" rel=\"noopener\">refit explainer</a>."))))
}

# HTML table for the fitted-trial list: one row per arm, the trial and its reference
# spanning that trial's arms, with footnote marks (v4_fit_trial_notes) and the notes below.
v4_fit_trials_table <- function(t) {
  fn <- v4_fit_trial_notes(t)
  head <- c("Trial and reference", "Arm", "Total dose (Gy)", "Dose per fraction (Gy)", "Fractions",
            "Overall time (days)", "Time basis", "Patients", "Published outcome (%)",
            "Endpoint, time point, estimator", "Outcome source (page, table or figure)", "Count in fit")
  rows <- list()
  for (tr in unique(t$trial)) {
    a <- t[t$trial == tr, , drop = FALSE]
    for (i in seq_len(nrow(a))) {
      first <- if (i == 1) tags$td(rowspan = nrow(a), class = "v4-trial-ref",
        tags$strong(tr), if (nzchar(fn$marks[[tr]])) tags$sup(fn$marks[[tr]]), tags$br(),
        tags$span(class = "small", a$citation[1]), tags$br(),
        tags$a(href = sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", a$pmid[1]), target = "_blank",
               rel = "noopener", class = "small", sprintf("PubMed %s", a$pmid[1])),
        if (nzchar(a$doi[1])) tagList(" · ", tags$a(href = paste0("https://doi.org/", a$doi[1]), target = "_blank",
                                                  rel = "noopener", class = "small", "DOI")),
        if (nzchar(a$schedule_source[1])) tagList(tags$br(), tags$span(class = "small v4-schedule-source",
                                                  paste0("Schedule: ", a$schedule_source[1]))))
      pub <- tagList(a$published[i], if (a$published_is_failure[i]) tags$sup("\u2020"),
                     if (a$figure_read[i]) tags$span(class = "small", " (read from figure)"))
      rows[[length(rows) + 1]] <- tags$tr(first, tags$td(a$arm[i]), tags$td(a$total[i]), tags$td(a$dose[i]),
        tags$td(a$fractions[i]), tags$td(a$time[i]), tags$td(a$time_basis[i]), tags$td(a$n[i]), tags$td(pub),
        tags$td(sprintf("%s, %s, %s", a$endpoint[i], a$horizon[i], a$estimator[i])),
        tags$td(class = "small", a$outcome_locator[i]), tags$td(a$fit_count[i]))
    }
  }
  tagList(
    div(class = "table-responsive",
        tags$table(class = "table table-sm table-bordered v4-trials",
                   tags$thead(tags$tr(lapply(head, tags$th))), tags$tbody(rows))),
    div(class = "small v4-trial-notes",
        lapply(fn$notes, function(n) tags$p(tags$sup(n$mark), " ", n$text)),
        if (any(t$published_is_failure)) tags$p(tags$sup("\u2020"), " ", paste(
          "Reported as a cumulative incidence of failure. The fit uses its complement, freedom from failure,",
          "which also counts patients who had a competing event such as death without failure.")),
        if (any(t$figure_read)) tags$p(paste(
          "Read from figure: the rate was read from the trial's published curve at the time point shown."))))
}
V4_TRANSFER_SCOPE <- v4_transfer_scope()
