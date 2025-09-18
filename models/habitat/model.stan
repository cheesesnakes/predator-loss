// generated with brms 2.23.0
functions {
  /* dirichlet-logit log-PDF
   * Args:
   *   y: vector of real response values
   *   mu: vector of category logit probabilities
   *   phi: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
     return dirichlet_lpdf(y | softmax(mu) * phi);
   }
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=2> ncat;  // number of categories
  array[N] vector[ncat] Y;  // response array
  int<lower=1> K_mubiomass;  // number of population-level effects
  matrix[N, K_mubiomass] X_mubiomass;  // population-level design matrix
  int<lower=1> Kc_mubiomass;  // number of population-level effects after centering
  int<lower=1> K_musponge;  // number of population-level effects
  matrix[N, K_musponge] X_musponge;  // population-level design matrix
  int<lower=1> Kc_musponge;  // number of population-level effects after centering
  int<lower=1> K_musubstrate;  // number of population-level effects
  matrix[N, K_musubstrate] X_musubstrate;  // population-level design matrix
  int<lower=1> Kc_musubstrate;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_mubiomass_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_musponge_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_musubstrate_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc_mubiomass] Xc_mubiomass;  // centered version of X_mubiomass without an intercept
  vector[Kc_mubiomass] means_X_mubiomass;  // column means of X_mubiomass before centering
  matrix[N, Kc_musponge] Xc_musponge;  // centered version of X_musponge without an intercept
  vector[Kc_musponge] means_X_musponge;  // column means of X_musponge before centering
  matrix[N, Kc_musubstrate] Xc_musubstrate;  // centered version of X_musubstrate without an intercept
  vector[Kc_musubstrate] means_X_musubstrate;  // column means of X_musubstrate before centering
  for (i in 2:K_mubiomass) {
    means_X_mubiomass[i - 1] = mean(X_mubiomass[, i]);
    Xc_mubiomass[, i - 1] = X_mubiomass[, i] - means_X_mubiomass[i - 1];
  }
  for (i in 2:K_musponge) {
    means_X_musponge[i - 1] = mean(X_musponge[, i]);
    Xc_musponge[, i - 1] = X_musponge[, i] - means_X_musponge[i - 1];
  }
  for (i in 2:K_musubstrate) {
    means_X_musubstrate[i - 1] = mean(X_musubstrate[, i]);
    Xc_musubstrate[, i - 1] = X_musubstrate[, i] - means_X_musubstrate[i - 1];
  }
}
parameters {
  vector[Kc_mubiomass] b_mubiomass;  // regression coefficients
  real Intercept_mubiomass;  // temporary intercept for centered predictors
  vector[Kc_musponge] b_musponge;  // regression coefficients
  real Intercept_musponge;  // temporary intercept for centered predictors
  vector[Kc_musubstrate] b_musubstrate;  // regression coefficients
  real Intercept_musubstrate;  // temporary intercept for centered predictors
  real<lower=0> phi;  // precision parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  array[M_3] vector[N_3] z_3;  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_mubiomass_1;  // actual group-level effects
  vector[N_2] r_2_musponge_1;  // actual group-level effects
  vector[N_3] r_3_musubstrate_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_mubiomass_1 = (sd_1[1] * (z_1[1]));
  r_2_musponge_1 = (sd_2[1] * (z_2[1]));
  r_3_musubstrate_1 = (sd_3[1] * (z_3[1]));
  lprior += normal_lpdf(b_mubiomass | 0,0.5);
  lprior += normal_lpdf(Intercept_mubiomass | 0,0.5);
  lprior += normal_lpdf(b_musponge | 0,0.5);
  lprior += normal_lpdf(Intercept_musponge | 0,0.5);
  lprior += normal_lpdf(b_musubstrate | 0,0.5);
  lprior += normal_lpdf(Intercept_musubstrate | 0,0.5);
  lprior += exponential_lpdf(phi | 2);
  lprior += exponential_lpdf(sd_1 | 2);
  lprior += exponential_lpdf(sd_2 | 2);
  lprior += exponential_lpdf(sd_3 | 2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mubiomass = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] musponge = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] musubstrate = rep_vector(0.0, N);
    // linear predictor matrix
    array[N] vector[ncat] mu;
    mubiomass += Intercept_mubiomass + Xc_mubiomass * b_mubiomass;
    musponge += Intercept_musponge + Xc_musponge * b_musponge;
    musubstrate += Intercept_musubstrate + Xc_musubstrate * b_musubstrate;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mubiomass[n] += r_1_mubiomass_1[J_1[n]] * Z_1_mubiomass_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      musponge[n] += r_2_musponge_1[J_2[n]] * Z_2_musponge_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      musubstrate[n] += r_3_musubstrate_1[J_3[n]] * Z_3_musubstrate_1[n];
    }
    for (n in 1:N) {
      mu[n] = transpose([0, mubiomass[n], musponge[n], musubstrate[n]]);
    }
    for (n in 1:N) {
      target += dirichlet_logit_lpdf(Y[n] | mu[n], phi);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_3[1]);
}
generated quantities {
  // actual population-level intercept
  real b_mubiomass_Intercept = Intercept_mubiomass - dot_product(means_X_mubiomass, b_mubiomass);
  // actual population-level intercept
  real b_musponge_Intercept = Intercept_musponge - dot_product(means_X_musponge, b_musponge);
  // actual population-level intercept
  real b_musubstrate_Intercept = Intercept_musubstrate - dot_product(means_X_musubstrate, b_musubstrate);
}

