// Dose-independent (DI) repopulation model, refit of Shuryak, Hall and Brenner
// 2018 Table 1 with per-trial intercepts and arm-level overdispersion.
// log(-log TCP_i) = c[trial] - alpha*(BED_LQ_i - mean_j) + lambda*smax(T_i - Tk) + g*(T_i - mean_j) + sigma*z_i
functions {
#include common_functions.stanfunctions
}
data {
  int<lower=1> N;                 // arms
  int<lower=1> J;                 // trials
  array[N] int<lower=1, upper=J> trial;
  vector[N] d;                    // dose per fraction (Gy)
  vector[N] nfx;                  // number of fractions (trial averages may be non-integer)
  vector[N] T;                    // overall time (days)
  array[N] int<lower=1> npat;
  array[N] int<lower=0> y;        // controlled = round(npat * TCP)
  real prior_log_ab_mean; real<lower=0> prior_log_ab_sd;   // tumor alpha/beta
  int<lower=0, upper=1> fix_ab; real ab_fixed;
  int<lower=0, upper=1> fix_g;  real g_fixed;
  int<lower=0, upper=1> overdispersion;
}
parameters {
  real log_ab_raw;
  real<lower=0> alpha;
  real<lower=0> lambda;
  real<lower=0, upper=60> Tk;
  real<lower=0.0046, upper=0.021> g_free;
  vector[J] c;
  real<lower=0> sigma_raw;
  vector[N] z;
}
transformed parameters {
  real ab = fix_ab ? ab_fixed : exp(prior_log_ab_mean + prior_log_ab_sd * log_ab_raw);
  real g = fix_g ? g_fixed : g_free;
  real sigma = overdispersion ? sigma_raw : 0;
  vector[N] eta;
  vector[N] bed;
  vector[J] bedbar = rep_vector(0, J);
  vector[J] Tbar = rep_vector(0, J);
  vector[J] cnt = rep_vector(0, J);
  for (i in 1:N) {
    bed[i] = nfx[i] * d[i] * (1 + d[i] / ab);
    bedbar[trial[i]] += bed[i]; Tbar[trial[i]] += T[i]; cnt[trial[i]] += 1;
  }
  bedbar = bedbar ./ cnt; Tbar = Tbar ./ cnt;
  // c[j] is log(-log TCP) of a hypothetical arm at trial j's mean BED and mean
  // time without repopulation: centering removes the c-alpha correlation.
  for (i in 1:N) {
    eta[i] = c[trial[i]] - alpha * (bed[i] - bedbar[trial[i]]) + lambda * smooth_max0(T[i] - Tk) + g * (T[i] - Tbar[trial[i]]) + sigma * z[i];
  }
}
model {
  log_ab_raw ~ std_normal();
  alpha ~ lognormal(log(0.2), 1);
  lambda ~ lognormal(log(0.1), 1.5);
  Tk ~ normal(28, 10);
  c ~ normal(0, 3);
  sigma_raw ~ normal(0, 0.5);
  z ~ std_normal();
  for (i in 1:N) target += binomial_lpmf(y[i] | npat[i], exp(-exp(eta[i])));
}
generated quantities {
  vector[N] tcp_arm;            // expected TCP of each arm including its own random effect
  array[N] int y_rep;           // replicate with a NEW arm-level effect (posterior predictive)
  vector[N] log_lik;
  for (i in 1:N) {
    real e_new = eta[i] - sigma * z[i] + sigma * normal_rng(0, 1);
    tcp_arm[i] = exp(-exp(eta[i]));
    y_rep[i] = binomial_rng(npat[i], exp(-exp(e_new)));
    log_lik[i] = binomial_lpmf(y[i] | npat[i], exp(-exp(eta[i])));
  }
}
