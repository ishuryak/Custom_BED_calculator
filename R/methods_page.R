# The page "How the models were fitted (Bayesian)" (version 4.6.1, review of 2026-11-12).
# Narrative: how_fitted.md. Numbers: draws/model_methods.json, built offline by
# tools/build_model_methods.py from the retained fit identities. The page restates no fitted
# number by hand: priors, settings and diagnostics come from the record, and the result
# summaries are the two existing briefs (refit_in_brief.md, early_fit_in_brief.md).

FITTED_PAGE <- "How the models were fitted (Bayesian)"
MODEL_METHODS_FILE <- "model_methods.json"
# Every link that opens the page, and the section it opens at.
FITTED_LINKS <- c(how_fitted_home = "fitted-tumor", how_fitted_unc = "fitted-tumor",
                  how_fitted_settings = "fitted-tumor", early_how = "fitted-early",
                  home_how_fitted = "fitted-bayes")

# The record, checked against the draws the app actually loads. ok = FALSE (with the reason)
# when it is missing, of another schema, or describes other draws than those deployed.
load_model_methods <- function(dir = DRAWS_DIR) {
  fail <- function(why) list(ok = FALSE, why = why)
  path <- file.path(dir, MODEL_METHODS_FILE)
  if (!file.exists(path) || file.size(path) == 0) return(fail("the methods record is missing"))
  m <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(m) || !identical(m$schema, "model-methods-1")) return(fail("the methods record is unreadable"))
  tm <- tryCatch(draws_manifest(dir)$files, error = function(e) NULL)
  em <- tryCatch(jsonlite::fromJSON(file.path(dir, "early_fit_manifest.json"), simplifyVector = FALSE)$files,
                 error = function(e) NULL)
  same <- is.list(tm) && is.list(em) &&
    identical(m$tumor$di$draws_sha256, tm$tumor_di_primary.csv$sha256) &&
    identical(m$tumor$dd$draws_sha256, tm$tumor_dd_primary.csv$sha256) &&
    identical(m$early$draws$sha256, em$early_fit_draws.csv$sha256)
  if (!same) return(fail("the methods record describes other fitted draws than those the calculator uses"))
  list(ok = TRUE, m = m)
}

md_escape <- function(x) gsub("|", "\\|", x, fixed = TRUE)

prior_table_md <- function(priors) {
  rows <- vapply(priors, function(p) sprintf("| %s | %s | %s | %s |", md_escape(p$parameter),
                                             md_escape(if (nzchar(p$units)) p$units else "–"),
                                             md_escape(p$prior), md_escape(p$status)), "")
  paste(c("| Parameter | Units | Prior or setting | Basis |", "|---|---|---|---|", rows), collapse = "\n")
}

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

sampler_text <- function(s, what) {
  sprintf("%s: %d chains, %s warm-up and %s sampling iterations per chain, target acceptance %g%s, seed %s.",
          what, s$chains, fmt_n(s$warmup), fmt_n(s$samples), s$adapt_delta,
          if (!is.null(s$max_treedepth)) sprintf(", maximum tree depth %d", s$max_treedepth) else "", s$seed)
}

brief_body <- function(file) paste(readLines(file, encoding = "UTF-8", warn = FALSE)[-1], collapse = "\n")

# The page's markdown with every placeholder filled. Stops if one is left, so a renamed
# placeholder cannot reach a reader as literal braces.
methods_page_text <- function(rec, page = "how_fitted.md", tumor_brief = "refit_in_brief.md",
                              early_brief = "early_fit_in_brief.md", support = tumor_data_support(),
                              sources = "") {
  txt <- paste(readLines(page, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  if (!isTRUE(rec$ok)) {
    check <- sprintf(paste0("<div class=\"alert alert-warning\">The technical details of the fits cannot be shown: %s. ",
                            "The summaries below are unaffected.</div>"), rec$why)
    m <- NULL
  } else { check <- ""; m <- rec$m }
  na <- "not available"
  t <- m$tumor; e <- m$early; l <- m$late
  dd_diag <- t$dd$diagnostics; di_diag <- t$di$diagnostics; ed <- e$diagnostics
  gate <- if (!is.null(m)) jsonlite::fromJSON(file.path(DRAWS_DIR, "draws_manifest.json"))$gate_policy else NULL
  vals <- if (is.null(m)) list() else list(
    TUMOR_ARMS = t$n_arms, TUMOR_UNITS = t$n_fitting_units,
    EARLY_ARMS = e$n_arms, EARLY_FAMILIES = e$n_families, EARLY_N = fmt_n(e$n_evaluable),
    LATE_AB = sprintf("%g", l$point_ab), LATE_MED = sprintf("%g", l$uncertainty$median),
    LATE_LO = sprintf("%g", l$uncertainty$central_95[[1]]), LATE_HI = sprintf("%g", l$uncertainty$central_95[[2]]),
    AB_LO = sprintf("%g", t$di$ab_prior_95[[1]]), AB_HI = sprintf("%g", t$di$ab_prior_95[[2]]),
    DD_ONSET = sprintf("%g", t$dd$onset_min), N_DI = fmt_n(t$draws$distinct_di_draws), PER_DI = t$draws$dd_per_di_draw,
    TUMOR_ROWS = fmt_n(t$draws$rows), EARLY_ROWS = fmt_n(e$draws$rows),
    TUMOR_DI_PRIORS = prior_table_md(t$di$priors), TUMOR_DD_PRIORS = prior_table_md(t$dd$priors),
    EARLY_PRIORS = prior_table_md(e$priors),
    TUMOR_SAMPLER = paste(sampler_text(t$di$sampler, "DI fit"),
      sprintf(paste("Calibrated DD: %s DI samples, each with %d chains of %d warm-up and %d sampling iterations,",
                    "lengthened (up to %d attempts, then targeted repairs) until the checks passed; %d fits were",
                    "repaired for R-hat and %d for effective sample size."),
              fmt_n(t$dd$n_di_draws), t$dd$sampler$chains, t$dd$sampler$warmup, t$dd$sampler$samples,
              t$dd$sampler$max_attempts, dd_diag$outer_fits_repaired, dd_diag$outer_fits_repaired_ess)),
    TUMOR_DIAG = sprintf(paste("DI: %d divergences, %d tree-depth hits, largest R-hat %.3f, smallest bulk and tail ESS",
                               "%s and %s, smallest E-BFMI %.2f (E-BFMI is an energy diagnostic; 0.3 or more is",
                               "acceptable). Calibrated DD: no fit failing the checks, largest R-hat %.3f, smallest",
                               "bulk and tail ESS per fit %s and %s. Both use the checks of policy %s, which include",
                               "tree-depth hits. Every parameter and every regimen contrast is checked."),
                         di_diag$divergences, di_diag$treedepth_hits, di_diag$max_rhat, fmt_n(round(di_diag$min_ess_bulk)),
                         fmt_n(round(di_diag$min_ess_tail)), di_diag$min_ebfmi, dd_diag$max_rhat,
                         fmt_n(round(dd_diag$min_ess_bulk_per_fit)), fmt_n(round(dd_diag$min_ess_tail_per_fit)),
                         t$draws$gate_policy),
    EARLY_TK = sprintf("%g", e$fixed$onset_elapsed_days), EARLY_H = sprintf("%g", e$fixed$repair_halftime_h),
    EARLY_SAMPLER = sampler_text(e$sampler, "Early fit"),
    EARLY_DIAG = sprintf(paste("%d divergences, largest R-hat %.3f, smallest bulk and tail ESS %s and %s,",
                               "smallest E-BFMI %.2f, over %d quantities. Tree-depth hits (%d) are recorded but,",
                               "unlike the tumor checks, are not part of this fit's pass rule."),
                         ed$divergences, ed$max_rhat, fmt_n(round(ed$min_ess_bulk)), fmt_n(round(ed$min_ess_tail)),
                         ed$min_ebfmi, ed$n_quantities, ed$treedepth_hits),
    RHAT_MAX = sprintf("%g", gate$max_rhat), ESS_MIN = sprintf("%g", gate$min_ess_bulk),
    PPC_IN = di_diag$ppc_inside_90)
  vals$CHECK <- check
  vals$TUMOR_BRIEF <- brief_body(tumor_brief)
  vals$EARLY_BRIEF <- brief_body(early_brief)
  vals$TUMOR_DMAX <- sprintf("%g", support$d_max)
  vals$SOURCES <- sources
  for (k in names(vals)) txt <- gsub(paste0("{{", k, "}}"), as.character(vals[[k]]), txt, fixed = TRUE)
  left <- regmatches(txt, gregexpr("\\{\\{[A-Z_0-9]+\\}\\}", txt))[[1]]
  if (length(left)) {
    if (!is.null(m)) stop("unfilled placeholders on the methods page: ", paste(unique(left), collapse = ", "))
    txt <- gsub("\\{\\{[A-Z_0-9]+\\}\\}", na, txt)
  }
  txt
}
