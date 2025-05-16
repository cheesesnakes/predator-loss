data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of protection levels
  int<lower=0> B; // number of behavioral types
  int<lower=0> S; // number of plots
  array[N, B] real<lower=0> D; // observed durations
  array[N] int<lower=0> bites; // number of bites observed
  
  // Predictor variables
  array[S] int<lower=1, upper=2> predator; // predator presence (0 or 1)
  array[N] int<lower=1, upper=S> plot; // plot number
  array[S, 3] real<lower=0, upper=1> rugosity_raw; // rugosity (continuous variable)
  array[S] int<lower=1, upper=2> protection; // protection level (0 or 1)
  array[S] real<lower=0, upper=1> biomass; // biomass availability (continuous variable)
}
parameters {
  // Rugosity
  array[S] real<lower=1e-6, upper=(1 - 1e-6)> rugosity;
  array[S] real<lower=1e-6> phi_rug;
  
  // Top-level coefficients for duration
  real beta_risk_foraging;
  real beta_risk_movement;
  real beta_risk_vigilance;
  
  real alpha_pi_foraging;
  real alpha_pi_movement;
  real alpha_pi_vigilance;
  
  real beta_pi_risk_foraging;
  real beta_pi_risk_movement;
  real beta_pi_risk_vigilance;
  
  // Top-level coefficients for bites
  real beta_risk_bites;
  
  // Risk model coefficients
  ordered[2] alpha_risk_protection;
  ordered[2] beta_predator;
  real beta_rug;
  real<lower=0> sigma_risk;
  real beta_res;
  
  // Observation model
  real<lower=0> sigma_D_foraging;
  real<lower=0> sigma_D_movement;
  real<lower=0> sigma_D_vigilance;
  
  // Latent variables non-centered
  array[N] real z_risk;
}
transformed parameters {
  array[S] real mu_risk;
  array[N] real risk;
  array[N] real lambda_risk;
  array[N, B] real mu_D;
  array[N, B] real<lower=0, upper=1> pi;
  
  // Compute mean risk per plot
  for (s in 1 : S) {
    mu_risk[s] = alpha_risk_protection[protection[s]]
                 + beta_predator[predator[s]] + beta_rug * rugosity[s]
                 + beta_res * biomass[s];
  }
  
  // Compute risk and other parameters per observation
  for (n in 1 : N) {
    risk[n] = mu_risk[plot[n]] + sigma_risk * z_risk[n];
    lambda_risk[n] = exp(beta_risk_bites * risk[n]);
    
    // Foraging
    pi[n, 1] = inv_logit(alpha_pi_foraging + beta_pi_risk_foraging * risk[n]);
    mu_D[n, 1] = beta_risk_foraging * risk[n];
    
    // Movement
    pi[n, 2] = inv_logit(alpha_pi_movement + beta_pi_risk_movement * risk[n]);
    mu_D[n, 2] = beta_risk_movement * risk[n];
    
    // Vigilance
    pi[n, 3] = inv_logit(alpha_pi_vigilance
                         + beta_pi_risk_vigilance * risk[n]);
    mu_D[n, 3] = beta_risk_vigilance * risk[n];
  }
}
model {
  // Priors
  beta_risk_bites ~ normal(0, 1);
  
  beta_risk_foraging ~ normal(-1, 0.5);
  beta_pi_risk_foraging ~ normal(1, 0.5);
  alpha_pi_foraging ~ normal(-1, 2);
  
  beta_risk_movement ~ normal(-1, 0.5);
  beta_pi_risk_movement ~ normal(0, 0.5);
  alpha_pi_movement ~ normal(-1, 2);
  
  beta_risk_vigilance ~ normal(1, 0.5);
  beta_pi_risk_vigilance ~ normal(-1, 0.5);
  alpha_pi_vigilance ~ normal(-1, 2);
  
  beta_rug ~ normal(-1, 1);
  beta_res ~ normal(0, 1);
  
  beta_predator ~ normal(0, 1);
  alpha_risk_protection ~ normal(0, 1);
  
  z_risk ~ normal(0, 0.5);
  sigma_risk ~ normal(0, 1) T[1e-6, ];
  
  rugosity ~ beta(1, 1);
  phi_rug ~ gamma(2, 2);
  
  sigma_D_foraging ~ normal(0, 1) T[1e-6, ];
  sigma_D_movement ~ normal(0, 1) T[1e-6, ];
  sigma_D_vigilance ~ normal(0, 1) T[1e-6, ];
  
  // Rugosity observations
  for (s in 1 : S) {
    for (i in 1 : 3) {
      rugosity_raw[s, i] ~ beta(rugosity[s] * phi_rug[s],
                                (1 - rugosity[s]) * phi_rug[s]);
    }
  }
  
  // Observation model
  for (n in 1 : N) {
    for (b in 1 : B) {
      if (D[n, b] == 0) {
        target += bernoulli_lpmf(1 | pi[n, b]);
      } else {
        target += bernoulli_lpmf(0 | pi[n, b]);
        if (b == 1) {
          target += lognormal_lpdf(D[n, b] | mu_D[n, b], sigma_D_foraging);
          if (D[n, b] > 0) {
            target += poisson_lpmf(bites[n] | lambda_risk[n]);
          }
        } else if (b == 2) {
          target += lognormal_lpdf(D[n, b] | mu_D[n, b], sigma_D_movement);
        } else if (b == 3) {
          target += lognormal_lpdf(D[n, b] | mu_D[n, b], sigma_D_vigilance);
        }
      }
    }
  }
}
generated quantities {
  array[N, B] real D_pred;
  array[N] int bites_pred;
  
  for (n in 1 : N) {
    real risk_draw = risk[n];
    real lambda = exp(beta_risk_bites * risk_draw);
    
    // Foraging
    real mu_foraging = beta_risk_foraging * risk_draw;
    real pi_foraging = inv_logit(alpha_pi_foraging
                                 + beta_pi_risk_foraging * risk_draw);
    if (bernoulli_rng(pi_foraging)) {
      D_pred[n, 1] = 0;
      bites_pred[n] = 0;
    } else {
      D_pred[n, 1] = lognormal_rng(mu_foraging, sigma_D_foraging);
      bites_pred[n] = poisson_rng(lambda);
    }
    
    // Movement
    real mu_movement = beta_risk_movement * risk_draw;
    real pi_movement = inv_logit(alpha_pi_movement
                                 + beta_pi_risk_movement * risk_draw);
    if (bernoulli_rng(pi_movement)) {
      D_pred[n, 2] = 0;
    } else {
      D_pred[n, 2] = lognormal_rng(mu_movement, sigma_D_movement);
    }
    
    // Vigilance
    real mu_vigilance = beta_risk_vigilance * risk_draw;
    real pi_vigilance = inv_logit(alpha_pi_vigilance
                                  + beta_pi_risk_vigilance * risk_draw);
    if (bernoulli_rng(pi_vigilance)) {
      D_pred[n, 3] = 0;
    } else {
      D_pred[n, 3] = lognormal_rng(mu_vigilance, sigma_D_vigilance);
    }
  }
}
