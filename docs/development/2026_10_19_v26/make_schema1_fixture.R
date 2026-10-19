# Builds a genuine schema-1 scenario file with the version 2.4 engine
# (commit ad16504) for the R1 ranking case, and its original results.
# Run from a worktree of ad16504:  Rscript make_schema1_fixture.R <out_dir>
args <- commandArgs(TRUE); out <- args[1]
source("R/bed_models.R"); source("R/export.R")
p <- modifyList(default_params(), list(alpha = 0.2, ab = 10, lambda_di = 0.4, lambda_dd = 5, Tk = 28, g = 0))
ref <- list(m = 35, d = 2, T = 47)
search <- list(m_range = 32L, d_step = 0.05, d_max = 2.5, criterion = "either", rank_by = "worst",
               tol_early = 30, tol_late = 30, early_cap = NULL, late_cap = NULL,
               min_gain = 0, min_days_saved = 1, fenwick = FALSE)
sc <- build_scenario(p, ref, ref, search, context = list(note = "R1 ranking case, built by v2.4"))
res <- do.call(suggest_regimens, c(list(p = p, ref = ref, anchor = ref), search))
sc$scenario_id <- scenario_id(sc)
writeLines(as.character(scenario_json(sc)), file.path(out, "schema1_v24_r1case.json"))
write.csv(res, file.path(out, "schema1_v24_r1case_results.csv"), row.names = FALSE)
print(res[, c("m", "d", "d_di", "d_dd", "score")]); cat("app", APP_VERSION, "equation", EQUATION_VERSION, "\n")
