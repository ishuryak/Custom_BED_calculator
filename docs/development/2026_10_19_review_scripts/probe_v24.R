# Independent probes of version 2.4 (commit ad16504), run from the app folder:
#   Rscript docs/development/2026_10_19_review_scripts/probe_v24.R
# Reproduces the external review's counterexamples and adds checks of my own.
source("R/bed_models.R"); source("R/export.R")
p <- default_params(); an <- default_anchor()
hy <- list(m = 20, d = 2.75, T = 26); ref70 <- default_ref()
cat("== Fixtures ==\n")
f <- compute_all(c(35, 20, 33, 15), c(2, 2.75, 2, 3.1), c(47, 26, 38, 19), p, an)
print(f[, c("m", "d", "T", "bed_di", "bed_dd", "bed_early", "bed_late")], digits = 10)
L <- legacy_params(); fl <- compute_all(c(35, 20), c(2, 2.75), c(47, 26), L, an)
print(fl[, c("bed_di", "bed_dd", "bed_early", "bed_late")], digits = 10)
cat("Fenwick T47/26/19:", fenwick_limit(c(47, 26, 19)), "\n")

cat("\n== F01: negative reference credited in agreement (150-day comparator) ==\n")
q <- modifyList(p, list(lambda_di = 0, g = 0))
ref150 <- list(m = 35, d = 2, T = 150)
cap <- compute_all(20, 2.75, 26, q, an)$bed_early
a <- agreement_regimens(q, ref150, an, margin = -3, early_cap = cap, min_days_saved = 1,
                        d_max = 5, m_range = 1:34)
print(a[, c("m", "d", "T", "d_di", "d_dd", "n_both_within", "n_scenarios")])
st <- sensitivity_table(list(m = 32, d = 2.10, T = 44), ref150, an, q)
r0 <- do.call(rbind, lapply(seq_len(nrow(st)), function(i) {
  qq <- modifyList(q, resolve_preset(st$preset[i], st$ab[i])); qq$ab <- st$ab[i]
  compute_all(c(32, 35), c(2.1, 2), c(44, 150), qq, an)[2, c("bed_di", "bed_dd")]
}))
print(cbind(st[, c("preset", "ab", "d_di", "d_dd")], ref_di = r0$bed_di, ref_dd = r0$bed_dd))
cat("valid both (ref>0):", sum(r0$bed_di > 0 & !is.na(r0$bed_dd) & r0$bed_dd > 0), "\n")

cat("\n== F02: pre-onset course still gains from g; r=2 counterexample ==\n")
r17 <- list(m = 17, d = 3.05, T = 23)
for (gg in c(0, 0.0154, 0.1)) {
  x <- compute_all(c(17, 20), c(3.05, 2.75), c(23, 26), modifyList(p, list(g = gg)), an)
  cat(sprintf("g=%.4f dDI=%.8f dDD=%.8f\n", gg, x$bed_di[1] - x$bed_di[2], x$bed_dd[1] - x$bed_dd[2]))
}
x <- compute_all(c(17, 20), c(3.05, 2.75), c(23, 26), modifyList(p, list(ab = 2)), an)
cat(sprintf("ab=2: dDI=%.8f, acute %.8f vs HYPNO %.8f\n", x$bed_di[1] - x$bed_di[2], x$bed_early[1], x$bed_early[2]))

cat("\n== F05: replay accepts unsupported metadata; ID depends on key order ==\n")
sc <- build_scenario(p, ref70, an, search = list(m_range = 5:34, d_step = 0.05, d_max = 5,
  criterion = "either", rank_by = "best", tol_early = 0, tol_late = 0, early_cap = cap,
  late_cap = NULL, min_gain = 1, min_days_saved = 7, fenwick = FALSE), context = list())
sc2 <- jsonlite::fromJSON(scenario_json(sc))
bad <- sc2; bad$equation_version <- "9.9"; bad$hinge$type <- "hard"; bad$clocks$tumor <- "elapsed"
cat("fabricated replay rows:", nrow(rerun_scenario(bad)), "(same as real:", nrow(rerun_scenario(sc2)), ")\n")
rev <- sc; rev$params <- rev(rev$params)
cat("ID original", scenario_id(sc), "reordered", scenario_id(rev), "\n")

cat("\n== F06: acute unit at ab_e = 5 ==\n")
cat("HYPNO acute cap at ab_e=5:", compute_all(20, 2.75, 26, modifyList(p, list(ab_e = 5)), an)$bed_early, "\n")

cat("\n== Own checks ==\n")
# DD onset reported as Inf when lambda_di = 0 (agreement table prints 'Inf')
q0 <- modifyList(p, list(lambda_di = 0))
a0 <- agreement_regimens(q0, hy, an, margin = -3, early_cap = compute_all(20, 2.75, 26, q0, an)$bed_early,
                         min_days_saved = 1, d_max = 5, m_range = 1:19)
cat("lambda_di=0 agreement rows:", nrow(a0), " dd_onset:", head(a0$dd_onset, 3), "\n")
# HYHOPE 12-fraction level, 44.4 Gy at 4 fx/week: protocol's 51.3 Gy10 on inclusive clock
for (T in 17:20) cat(sprintf("12 x 3.7, T=%d inclusive-clock acute = %.2f\n", T,
  bed_early(12, 3.7, T, 10, 0.35, log(2) / 2.5, 7)))
cat(sprintf("15 x 3.1, T=19 inclusive-clock acute = %.2f, elapsed = %.2f\n",
  bed_early(15, 3.1, 19, 10, 0.35, log(2) / 2.5, 7), bed_early(15, 3.1, 18, 10, 0.35, log(2) / 2.5, 7)))
# 2018 fitted DD (own alpha, max rate, C, g) for GORTEC/TROG/HYPNO comparisons
fit <- list(alpha = 0.224, lam = 1.17, C = 14.5, g = 0.0154, ab = 10)
dd_fit <- function(m, d, T) {
  K <- bed_lq(m, d, 10); on <- fit$C * T / (fit$alpha * K)
  K - fit$lam * p_kill(K, T, fit$alpha) * smooth_max0(T - on) / fit$alpha - fit$g * T / fit$alpha
}
cat(sprintf("fitted DD: GORTEC 63/31.5 in 22d vs 70/35 48d: %+.2f ; HYPNO 55/20 vs 66/33@38d: %+.2f\n",
  dd_fit(31.5, 2, 22) - dd_fit(35, 2, 48), dd_fit(20, 2.75, 26) - dd_fit(33, 2, 38)))
cat(sprintf("fitted DD: TROG 59.4/33 in 24d vs 70/35 in 50d: %+.2f\n", dd_fit(33, 1.8, 24) - dd_fit(35, 2, 50)))
# Illustrative clinical translation of the HYPNO test-minus-control BED difference.
# gamma50 (percentage points of control per 1% dose change) is an ASSUMED range,
# not a verified source value; LRC near 50% is at the steepest part of the curve.
eq <- function(b) eqd2(b, 10)
for (m in c("DI", "DD")) {
  dB <- if (m == "DI") f$bed_di[2] - f$bed_di[3] else f$bed_dd[2] - f$bed_dd[3]
  dpct <- 100 * eq(dB) / eq(if (m == "DI") f$bed_di[3] else f$bed_dd[3])
  cat(sprintf("%s: dBED %.2f Gy10 = %.1f%% of control EQD2 -> LRC change %.1f to %.1f pp at gamma50 1.5-2.5\n",
              m, dB, dpct, 1.5 * dpct, 2.5 * dpct))
}
cat(sprintf("HYPNO upper 90%% HR 1.30 at control LRC 51.2%%: test LRC %.1f%% (%.1f pp)\n",
            100 * 0.512^1.30, 100 * (0.512^1.30 - 0.512)))
