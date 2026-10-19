# Conditional illustration on the "Model and sources" page: what HYPNO's
# test-minus-control tumor BED deficit under scenario D18 would mean in
# locoregional control IF a published dose-response slope applied to it.
#
# Extra assumptions, none supplied by the calculator or the sources:
#  A1. The normalized slope gamma50 of Stuschke & Thames 1999 (logit fit,
#      gamma50 3.2, 95% CI 1.6-4.8, Table 1 p115 and Appendix A eqs 1-3 p119;
#      RTOG 83-13 0.8 +/- 0.4, p117) applies unchanged on the calculator's
#      time-corrected EQD2 scale. Stuschke defines gamma50 = D dP/dD at TCD50
#      with D the physical or isoeffective conventional dose; a normalized
#      slope is not invariant under an additive time correction
#      (gamma_E = [E50 / (D50 f'(D50))] gamma_D), so A1 is a transfer, not a
#      derivation.
#  A2. A logistic curve P(E) = 1 / (1 + exp(4 gamma (1 - E/E50))) passes
#      through HYPNO's control-arm 3-year locoregional control, 0.512
#      (Bentzen 2026 p391-392, Table 3), at the control arm's EQD2.
#  A3. The fractional EQD2 change equals the fractional BED change (EQD2 is
#      BED / (1 + 2/ab), so the factor cancels in the ratio).
#  A4. The HYPNO upper 90% confidence limit of the hazard ratio, 1.30, is
#      converted to a control-probability difference under proportional
#      hazards: 0.512 - 0.512^1.30. This is not an estimated bound on the
#      absolute difference.
# Population, chemotherapy (75% weekly cisplatin in HYPNO, none in most
# hyperfractionation trials) and endpoint timing are assumed transportable.
#
# Run from the app folder: Rscript docs/development/2026_10_19_v26/dose_response_illustration.R

illustration_table <- function(gammas = c(0.8, 1.2, 1.6, 3.2, 4.8), p0 = 0.512, hr_upper = 1.30) {
  p <- default_params(); an <- default_anchor()
  x <- compute_all(c(20, 33), c(2.75, 2), c(26, 38), p, an)
  rows <- list()
  for (model in c("DI", "DD")) {
    b <- if (model == "DI") x$bed_di else x$bed_dd
    ratio <- (b[1] - b[2]) / b[2]
    for (g in gammas) {
      l0 <- qlogis(p0)
      p1 <- plogis(l0 + (4 * g + l0) * ratio)   # exact logistic through (E0, p0)
      rows[[length(rows) + 1]] <- data.frame(model = model, control_bed = b[2], delta_bed = b[1] - b[2],
        fractional_change = ratio, gamma50 = g, loss_pp_linear = -100 * g * ratio,
        loss_pp_logistic = 100 * (p0 - p1))
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "ph_upper_loss_pp") <- 100 * (p0 - p0^hr_upper)
  out
}

if (sys.nframe() == 0) {
  source("R/bed_models.R")
  t <- illustration_table()
  print(t, digits = 6)
  cat("HR 1.30 under proportional hazards at control 0.512:", attr(t, "ph_upper_loss_pp"), "points\n")
}
