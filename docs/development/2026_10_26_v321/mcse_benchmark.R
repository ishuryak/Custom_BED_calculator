# Monte Carlo standard errors of the Uncertainty tab's 4000-draw summaries, for named scenarios
# (review F07 of the 2026-10-26 v3.2 review). Resamples the 1000 distinct DI draws with their 4
# paired DD rows (cluster bootstrap, B = 2000, seed 1). Run from the app root:
#   Rscript docs/development/2026_10_26_v321/mcse_benchmark.R
source("R/bed_models.R"); source("R/uncertainty.R")
di <- load_draws("di", "primary"); dd <- load_draws("dd", "primary")
ref <- list(m = 35, d = 2, T = 47)
scen <- list(default_18x2.9_24d = list(m = 18, d = 2.9, T = 24), HYPNO_20x2.75_26d = list(m = 20, d = 2.75, T = 26))
cl <- utils::read.csv("draws/tumor_di_primary.csv")$di_draw; stopifnot(length(cl) == nrow(di), length(unique(cl)) == 1000); ids <- unique(cl); B <- 2000; set.seed(1)
boot_idx <- lapply(seq_len(B), function(b) unlist(split(seq_along(cl), cl)[as.character(sample(ids, replace = TRUE))], use.names = FALSE))
out <- list()
for (nm in names(scen)) {
  td <- tumor_deltas(scen[[nm]], ref, di, dd)
  for (mod in c("di", "dd")) {
    x <- as.numeric(td[[paste0("delta_", mod)]][1, ])
    stat <- function(v) c(median = median(v, na.rm = TRUE), q025 = unname(quantile(v, .025, na.rm = TRUE)),
                          q975 = unname(quantile(v, .975, na.rm = TRUE)), share_ge_p1 = mean(v >= 1, na.rm = TRUE),
                          share_ge_m3 = mean(v >= -3, na.rm = TRUE))
    est <- stat(x); se <- apply(sapply(boot_idx, function(i) stat(x[i])), 1, sd)
    out[[length(out) + 1]] <- data.frame(scenario = nm, model = mod, estimand = names(est), estimate = est, mcse = se)
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
res$estimate <- signif(res$estimate, 5); res$mcse <- signif(res$mcse, 3)
print(res, row.names = FALSE)
write.csv(res, "docs/development/2026_10_26_v321/mcse_benchmark.csv", row.names = FALSE)
