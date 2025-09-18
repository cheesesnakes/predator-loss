// generated with brms 2.23.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int<lower=1> Ksp;  // number of special effects terms
  // covariates of special effects terms
  vector[N] Csp_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for noise-free variables
  int<lower=1> Mme_1;  // number of groups
  vector[N] Xn_1;  // noisy values
  vector<lower=0>[N] noise_1;  // measurement noise
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector[Ksp] bsp;  // special effects coefficients
  real<lower=0> shape;  // shape parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  // parameters for noise free variables
  vector[Mme_1] meanme_1;  // latent means
  vector<lower=0>[Mme_1] sdme_1;  // latent SDs
  vector[N] zme_1;  // standardized latent values
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N] Xme_1;  // actual latent values
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_1 = (sd_1[1] * (z_1[1]));
  // compute actual latent values
  Xme_1 = meanme_1[1] + sdme_1[1] * zme_1;
  lprior += normal_lpdf(b | 0,0.5);
  lprior += normal_lpdf(Intercept | 0,0.5);
  lprior += normal_lpdf(bsp | 0,0.5);
  lprior += exponential_lpdf(shape | 2);
  lprior += exponential_lpdf(sd_1 | 2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += (bsp[1]) * Xme_1[n] + (bsp[2]) * Xme_1[n] * Csp_1[n] + r_1_1[J_1[n]] * Z_1_1[n];
    }
    target += neg_binomial_2_log_glm_lpmf(Y | Xc, mu, b, shape);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += normal_lpdf(Xn_1 | Xme_1, noise_1);
  target += std_normal_lpdf(zme_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

