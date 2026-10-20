# Freezes version 2.6 scenario files and their full results BEFORE the 2.7
# interface changes, so tests/testthat/test_v27_colleague_review.R can check that
# compatible 2.6 files still replay to the same numbers. Run once, from the app
# folder, on the 2.6 code (commit 0c51ae6):
#   Rscript docs/development/2026_10_20_v27/make_v26_replay_fixtures.R
source("tools/env.R")
invisible(use_bed_library(required = TRUE))
source("R/bed_models.R")
source("R/export.R")
stopifnot(APP_VERSION == "2.6.0")

p <- default_params()
std <- list(m = 35, d = 2, T = 47)
hypno <- list(m = 20, d = 2.75, T = 26)
cap <- compute_all(20, 2.75, 26, p, std)$bed_early
base <- list(m_range = 5:34, d_step = 0.05, d_max = 5, tol_early = 0, tol_late = 0,
             early_cap = cap, late_cap = NULL, min_days_saved = 7, fenwick = FALSE)
cand <- list(
  either_best_p1 = c(base, list(criterion = "either", rank_by = "best", min_gain = 1)),
  both_worst_m3  = c(base, list(criterion = "both", rank_by = "worst", min_gain = -3)),
  both_worst_p1  = c(base, list(criterion = "both", rank_by = "worst", min_gain = 1)))
out <- file.path("tests", "testthat", "fixtures")
save_case <- function(name, sc, res) {
  writeLines(scenario_json(sc), file.path(out, paste0("v26_", name, ".json")), useBytes = TRUE)
  feas <- attr(res, "feasible")
  writeLines(jsonlite::toJSON(list(rows = as.data.frame(res), feasible = as.data.frame(feas)),
                              digits = NA, na = "string", dataframe = "rows"),
             file.path(out, paste0("v26_", name, "_results.json")), useBytes = TRUE)
  cat(name, ": ", nrow(res), " rows\n", sep = "")
}
for (nm in names(cand)) {
  sc <- build_scenario(p, std, std, search = cand[[nm]], context = list(), analysis = "candidate")
  sc$scenario_id <- scenario_id(sc)
  res <- rerun_scenario(jsonlite::fromJSON(scenario_json(sc)))
  save_case(paste0("candidate_", nm), sc, res)
}
agree_args <- list(m_range = 1:19, min_days_saved = 1, d_step = 0.05, d_max = 5, tol_early = 0,
                   tol_late = 0, early_cap = cap, late_cap = NULL, fenwick = FALSE,
                   criterion = "both", rank_by = "worst", min_gain = -3)
sc <- build_scenario(p, hypno, std, search = agree_args, context = list(agree_cmp = "hypno_test"),
                     analysis = "agreement", agreement = list(margin = -3, ab_values = SCENARIO_AB_VALUES))
sc$scenario_id <- scenario_id(sc)
res <- rerun_scenario(jsonlite::fromJSON(scenario_json(sc)))
save_case("agreement_default", sc, res)
