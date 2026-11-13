# Independent R cross-checks for the likelihood-only comparison (plan revision 3, test T6).
# Run by fit/frequentist/r_check.py from the folder holding r_fixture.json (bare file names:
# Windows R cannot resolve /mnt/c paths). Writes r_result.json. Needs lme4 (user library);
# never part of the app's release library.
fx <- jsonlite::fromJSON("r_fixture.json", simplifyVector = TRUE)
out <- list(R = R.version.string, lme4 = as.character(utils::packageVersion("lme4")))

# Tumor: cloglog on failures, unit fixed effects, one arm random effect, 25-point adaptive quadrature,
# at the fixed (alpha/beta, Tk) of the Python fit. Coefficients: B -> -alpha, hinge -> lambda, T -> g.
t <- as.data.frame(fx$tumor[c("unit", "B", "hinge", "T", "fail", "ctrl")])
t$unit <- factor(t$unit); t$arm <- factor(seq_len(nrow(t)))
ft <- lme4::glmer(cbind(fail, ctrl) ~ 0 + unit + B + hinge + T + (1 | arm), data = t,
                  family = binomial("cloglog"), nAGQ = 25,
                  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5),
                                               check.conv.singular = "ignore"))
fe <- lme4::fixef(ft)
out$tumor <- list(alpha = -unname(fe["B"]), lam = unname(fe["hinge"]), g = unname(fe["T"]),
                  sigma = unname(sqrt(lme4::VarCorr(ft)$arm[1])),
                  intercepts = unname(fe[grep("^unit", names(fe))]), unit_levels = levels(t$unit),
                  logLik = as.numeric(stats::logLik(ft)))

# Early, zero random effects: binomial GLM on synthetic counts from known interior parameters.
# Coefficients: dD -> kappa, dQ -> kappa / ab, dTe -> -kappa r. Profile interval of the dD coefficient.
e <- as.data.frame(fx$early_glm[c("y", "N", "grp", "dD", "dQ", "dTe", "chemo")])
fg <- stats::glm(cbind(y, N - y) ~ grp + dD + dQ + dTe + chemo, data = e, family = binomial("logit"),
                 control = stats::glm.control(epsilon = 1e-14, maxit = 100))
ci <- suppressMessages(stats::confint(fg, parm = "dD", level = 0.95))
out$early_glm <- list(coef = as.list(stats::coef(fg)), dD_ci = unname(ci), logLik = as.numeric(stats::logLik(fg)))

# Early, both random effects, Laplace (descriptive only: Laplace can move fixed effects and variances).
r <- as.data.frame(fx$early_real[c("y", "N", "grp", "dD", "dQ", "dTe", "chemo", "fam")])
r$fam <- factor(r$fam); r$arm <- factor(seq_len(nrow(r)))
fl <- lme4::glmer(cbind(y, N - y) ~ grp + dD + dQ + dTe + chemo + (1 | fam) + (1 | arm), data = r,
                  family = binomial("logit"),
                  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
vc <- lme4::VarCorr(fl)
out$early_laplace <- list(coef = as.list(lme4::fixef(fl)), tau = unname(sqrt(vc$fam[1])),
                          sigma_arm = unname(sqrt(vc$arm[1])), logLik = as.numeric(stats::logLik(fl)))
jsonlite::write_json(out, "r_result.json", auto_unbox = TRUE, digits = NA)
cat("r_result.json written\n")
