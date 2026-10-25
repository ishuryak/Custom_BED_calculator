// Calibrated dose-dependent (DD) repopulation model, conditional on one draw of
// the DI fit (as in the app's dd_calibration): alpha, alpha/beta and g are taken
// from the DI draw, and the onset constant is solved so that DD repopulation
// equals DI repopulation R0 = lambda_di * smax(T0 - Tk) at the anchor regimen.
// Only lambda_dd, the trial intercepts and the arm-level overdispersion are fitted.
// A smooth exponential penalty exp((onset_min - onset_i) / onset_scale) keeps every
// fitted arm's onset at or after onset_min days.
functions {
#include common_functions.stanfunctions
  real inv_smooth_max0(real y) {
    real z = 3 * y;
    if (z > 30) return y;
    return log(expm1(z)) / 3;
  }
}
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> trial;
  vector[N] d;
  vector[N] nfx;
  vector[N] T;
  array[N] int<lower=1> npat;
  array[N] int<lower=0> y;
  real<lower=0> alpha;        // from the DI draw
  real<lower=0> ab;           // from the DI draw
  real<lower=0> g;            // from the DI draw
  real<lower=0> R0;           // DI repopulation at the anchor (log cell units)
  real<lower=0> K0;           // anchor BED (Gy, at ab)
  real<lower=0> T0;           // anchor overall time (days)
  real onset_min;             // days
  real<lower=0> onset_scale;  // days
  real prior_log_lambda_mean; real<lower=0> prior_log_lambda_sd;
}
transformed data {
  vector[N] bed;
  vector[J] bedbar = rep_vector(0, J);
  vector[J] Tbar = rep_vector(0, J);
  vector[J] cnt = rep_vector(0, J);
  real pkill0 = -expm1(-alpha * K0 / T0);
  for (i in 1:N) {
    bed[i] = nfx[i] * d[i] * (1 + d[i] / ab);
    bedbar[trial[i]] += bed[i]; Tbar[trial[i]] += T[i]; cnt[trial[i]] += 1;
  }
  bedbar = bedbar ./ cnt; Tbar = Tbar ./ cnt;
}
parameters {
  real<lower=0> lambda_dd;
  vector[J] c;
  real<lower=0> sigma_raw;
  vector[N] z;
}
transformed parameters {
  real onset_ref = T0 - inv_smooth_max0(R0 / (lambda_dd * pkill0));
  real c_over_alpha = K0 * onset_ref / T0;
  vector[N] onset;
  vector[N] eta;
  for (i in 1:N) {
    real pkill = -expm1(-alpha * bed[i] / T[i]);
    onset[i] = c_over_alpha * T[i] / bed[i];
    eta[i] = c[trial[i]] - alpha * (bed[i] - bedbar[trial[i]])
             + lambda_dd * pkill * smooth_max0(T[i] - onset[i])
             + g * (T[i] - Tbar[trial[i]]) + sigma_raw * z[i];
  }
}
model {
  lambda_dd ~ lognormal(prior_log_lambda_mean, prior_log_lambda_sd);
  c ~ normal(0, 3);
  sigma_raw ~ normal(0, 0.5);
  z ~ std_normal();
  target += -sum(exp((onset_min - onset) / onset_scale));
  for (i in 1:N) target += binomial_lpmf(y[i] | npat[i], exp(-exp(eta[i])));
}
generated quantities {
  // Calibration check: DD repopulation at the anchor, which must equal R0.
  real rep_anchor = lambda_dd * pkill0 * smooth_max0(T0 - onset_ref);
  real onset_min_fitted = min(onset);
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = binomial_lpmf(y[i] | npat[i], exp(-exp(eta[i])));
}
