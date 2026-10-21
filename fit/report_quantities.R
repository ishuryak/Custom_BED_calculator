# Clinical comparisons under the refit draws, with the app's own engine.
# Writes fit/output/report_quantities.json. Run from the repository root.
source("tools/env.R"); invisible(use_bed_library(TRUE))
for (f in c("bed_models.R", "schedule.R", "uncertainty.R")) source(file.path("R", f))
ref <- list(m = 35, d = 2, T = 47)
regs <- data.frame(name = c("HYPNO 55 Gy/20 (26 d)", "17 x 3.05 Gy (23 d)", "18 x 2.90 Gy (24 d)", "19 x 2.80 Gy (25 d)",
                            "HYHOPE 46.5 Gy/15 (19 d)", "30 x 2 Gy (40 d)", "HYPNO control 66 Gy/33, 6/week (38 d)"),
                   m = c(20, 17, 18, 19, 15, 30, 33), d = c(2.75, 3.05, 2.9, 2.8, 3.1, 2, 2), T = c(26, 23, 24, 25, 19, 40, 38))
out <- list()
for (v in c("primary", "paper_restrictions")) {
  di <- load_draws("di", v); dd <- load_draws("dd", v)
  td <- tumor_deltas(regs, ref, di, dd)
  s <- summarize_draws(regs, td, threshold = 0)
  out[[v]] <- data.frame(name = regs$name, s[, c("di_median", "di_q10", "di_q90", "di_q025", "di_q975", "dd_median",
                                                 "dd_q10", "dd_q90", "share_di", "share_dd", "n_invalid_di", "n_invalid_dd")])
}
# deterministic D18 and refit-preset values for the same regimens
p18 <- default_params()
pts <- jsonlite::fromJSON("draws/refit_point_values.json")$values
pref <- modifyList(p18, list(alpha = pts$di$alpha, lambda_di = pts$di$lambda, Tk = pts$di$Tk, g = pts$di$g))
det <- function(p) { a <- compute_all(regs$m, regs$d, regs$T, p, default_anchor()); r <- compute_all(35, 2, 47, p, default_anchor())
  data.frame(name = regs$name, d_di = a$bed_di - r$bed_di, d_dd = a$bed_dd - r$bed_dd) }
out$d18 <- det(p18); out$refit_point_di <- det(pref)
jsonlite::write_json(out, "fit/output/report_quantities.json", digits = NA, pretty = TRUE)
print(out$primary[, 1:4]); print(out$d18); print(out$refit_point_di)
