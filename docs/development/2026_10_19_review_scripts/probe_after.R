# Probe outputs after the version 2.5 fixes (compare with probe_v24.R on ad16504).
source("R/bed_models.R"); source("R/export.R")
p <- default_params(); an <- default_anchor(); q <- modifyList(p, list(lambda_di = 0))
cap <- compute_all(20, 2.75, 26, q, an)$bed_early
a <- agreement_regimens(q, list(m = 35, d = 2, T = 150), an, margin = -3, early_cap = cap,
                        min_days_saved = 1, d_max = 5, m_range = 1:34)
print(a[, c("m", "d", "T", "n_scenarios", "n_valid_both", "n_pass_both", "n_invalid", "d_di_min", "d_di_max")])
hy <- list(m = 20, d = 2.75, T = 26)
a2 <- agreement_regimens(p, hy, an, margin = -3, early_cap = compute_all(20, 2.75, 26, p, an)$bed_early,
                         min_days_saved = 1, d_max = 5, m_range = 1:19)
print(a2[, c("m", "d", "T", "d_di", "d_dd", "n_valid_both", "n_pass_both")], digits = 9)
s <- suggest_regimens(p, default_ref(), an, d_max = 5, min_gain = 1, min_days_saved = 7,
                      early_cap = compute_all(20, 2.75, 26, p, an)$bed_early)
print(s[, c("m", "d", "T", "d_di", "d_dd", "credited_by")], digits = 9)
