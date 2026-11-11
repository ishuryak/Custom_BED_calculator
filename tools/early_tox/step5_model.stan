// Plan step 5: arm-level grade >=3 acute mucositis, within-trial contrasts of early BED.
//   y_i ~ Binomial(N_i, p_i)
//   logit p_i = eta_fam[i] + kappa * dE_i + gamma * chemo_i + sigma_arm * zarm_i
//   dE_i = dD_i + dQ_i / ab - r * dTe_i          (arm minus its family reference arm)
//   eta_j = m0 + delta' * scale_group_j + tau0 * z_j
// Switches: fixed_E (benchmark: ab and r fixed at the Fowler values), use_gamma (C1),
// use_od (C3), n_groups (M1 scale-family offsets), prior_only (prior predictive / SBC).
data {
  int<lower=1> I;
  int<lower=1> J;
  array[I] int<lower=1, upper=J> fam;
  array[I] int<lower=0> N;
  array[I] int<lower=0> y;
  vector[I] dD;
  vector[I] dQ;
  vector[I] dTe;
  vector[I] chemo;
  int<lower=0> G;                         // number of non-reference scale groups
  matrix[J, G] grp;                       // family j belongs to scale group g (G1 is the reference)
  int<lower=0, upper=1> fixed_E;
  real<lower=0> ab_fix;
  real<lower=0> r_fix;
  int<lower=0, upper=1> use_gamma;
  int<lower=0, upper=1> use_od;
  real<lower=0> Lstar_sd;
  real ab_mu;
  real<lower=0> ab_sd;
  real r_mu;
  real<lower=0> r_sd;
  real<lower=0> gamma_sd;
  real<lower=0> od_sd;
  int<lower=0> K;                         // decision schedules, each minus 70 Gy / 35
  vector[K] kD;
  vector[K] kQ;
  vector[K] kTe;
  int<lower=0, upper=1> prior_only;
}
parameters {
  real m0;
  real<lower=0> tau0;
  vector[J] z;
  vector[G] delta;
  real<lower=0> Lstar;
  array[1 - fixed_E] real log_ab;
  array[1 - fixed_E] real log_r;
  array[use_gamma] real gamma;
  array[use_od] real<lower=0> sigma_arm;
  vector[use_od ? I : 0] zarm;
}
transformed parameters {
  real kappa = Lstar / 10;
  real ab = ab_fix;
  real r = r_fix;
  if (!fixed_E) {
    ab = exp(log_ab[1]);
    r = exp(log_r[1]);
  }
  vector[J] eta = m0 + grp * delta + tau0 * z;
  vector[I] dE = dD + dQ / ab - r * dTe;
  vector[I] lin = eta[fam] + kappa * dE;
  if (use_gamma)
    lin += gamma[1] * chemo;
  if (use_od)
    lin += sigma_arm[1] * zarm;
}
model {
  m0 ~ normal(0, 2.5);
  tau0 ~ normal(0, 1);
  z ~ std_normal();
  delta ~ normal(0, 1);
  Lstar ~ normal(0, Lstar_sd);
  if (!fixed_E) {
    log_ab ~ normal(ab_mu, ab_sd);
    log_r ~ normal(r_mu, r_sd);
  }
  if (use_gamma)
    gamma ~ normal(0, gamma_sd);
  if (use_od) {
    sigma_arm ~ normal(0, od_sd);
    zarm ~ std_normal();
  }
  if (!prior_only)
    y ~ binomial_logit(N, lin);
}
generated quantities {
  vector[I] log_lik;
  array[I] int y_rep;
  vector[K] dE_dec = kD + kQ / ab - r * kTe;      // early BED of decision schedule minus 70/35
  vector[K] lor_dec = kappa * dE_dec;            // log odds ratio against 70/35
  real T_obs = 0;
  real T_rep = 0;
  for (i in 1:I) {
    real p = inv_logit(lin[i]);
    log_lik[i] = binomial_logit_lpmf(y[i] | N[i], lin[i]);
    y_rep[i] = binomial_rng(N[i], p);
    T_obs += square(y[i] - N[i] * p) / (N[i] * p * (1 - p));
    T_rep += square(y_rep[i] - N[i] * p) / (N[i] * p * (1 - p));
  }
}
