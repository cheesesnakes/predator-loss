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
  
  array[S] real<lower=1e-6, upper=(1 - 1e-6)> rugosity; // accounting for multiple rugosity measurements
  array[S] real<lower=1e-6> phi_rug; // standard deviation of rugosity in each plot
  
  // Top-level coefficients for duration
  vector[B] beta_risk; // effect of risk on duration
  
  vector[B] alpha_pi; // baseline zero-inflation
  vector[B] beta_pi_risk; // effect of risk on zero-inflation
  
  // Top-level coefficients for bites
  real beta_risk_bites; // effect of risk on bites
  
  // Risk model coefficients
  ordered[2] alpha_risk_protection; // random risk per individual
  ordered[2] beta_predator; // effect of predator presence on risk
  real beta_rug; // effect of rugosity on risk
  real beta_res; // effect of biomass on risk
  real<lower=0> sigma_risk; // standard deviation of riskrgy
  
  // Observation model
  vector<lower=0>[B] sigma_D; // standard deviation of duration, non-zero
  
  // Latent variables non-centred
  array[N] real z_risk; // latent risk variable
}
transformed parameters {
  array[S] real mu_risk; // mean risk
  
  array[N, B] real mu_D; // mean duration
  
  array[N] real lambda_risk; // mean bites
  
  array[N, B] real<lower=0, upper=1> pi; // zero-inflation probability
  
  // Latent variables per observation
  array[N] real risk; // latent risk variable
  
  // mean risk per plot
  for (s in 1 : S) {
    mu_risk[s] = alpha_risk_protection[protection[s]]
                 + beta_predator[predator[s]] + beta_rug * rugosity[s]
                 + beta_res * biomass[s];
  }
  
  for (n in 1 : N) {
    // Latent risk variable
    risk[n] = exp(mu_risk[plot[n]] + sigma_risk * z_risk[n]); // log-normal non-centred
    
    for (b in 1 : B) {
      // zero-inflation probability
      real logit_input = alpha_pi[b] + beta_pi_risk[b] * risk[n];
      // zero-inflation probability
      pi[n, b] = inv_logit(logit_input);
      
      // Mean duration
      mu_D[n, b] = beta_risk[b] * risk[n];
    }
    
    // Mean bites
    
    lambda_risk[n] = exp(beta_risk_bites * risk[n]);
  }
}
model {
  // Priors
  
  // Top-level priors
  beta_risk_bites ~ normal(0, 1);
  
  // priors for foraging
  beta_risk[1] ~ normal(-1, 1);
  beta_pi_risk[1] ~ normal(1, 1);
  alpha_pi[1] ~ normal(-1, 2);
  
  // priors for vigilance
  beta_risk[2] ~ normal(1, 1);
  beta_pi_risk[2] ~ normal(-1, 1);
  alpha_pi[2] ~ normal(-1, 2);
  
  // priors for movement
  beta_risk[3] ~ normal(-1, 1);
  beta_pi_risk[3] ~ normal(0, 1);
  alpha_pi[3] ~ normal(-1, 2);
  
  // Risk model priors
  beta_rug ~ normal(-1, 1);
  beta_res ~ normal(0, 1);
  
  for (p in 1 : 2) {
    beta_predator[p] ~ normal(0, 1);
  }
  for (r in 1 : 2) {
    alpha_risk_protection[r] ~ normal(0, 1);
  }
  
  sigma_risk ~ normal(0, 1) T[1e-6, ]; // standard deviation of risk
  
  // rugosity priors
  rugosity ~ beta(1, 1); // beta distribution for rugosity
  phi_rug ~ gamma(2, 2); // precision of rugosity
  
  // Latent variables priors
  
  for (b in 1 : B) {
    sigma_D[b] ~ normal(0, 1) T[1e-6, ]; // standard deviation of duration
  }
  
  z_risk ~ normal(0, 0.5); // latent risk variable
  
  // rugosity
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
        target += bernoulli_lpmf(1 | pi[n, b]); // zero-inflation
      } else {
        target += bernoulli_lpmf(0 | pi[n, b]); // non-zero inflation
        target += lognormal_lpdf(D[n, b] | mu_D[n, b], sigma_D[b]); // duration model
        if (b == 1 && D[n, b] > 0) {
          target += poisson_lpmf(bites[n] | lambda_risk[n]); // bites model
        }
      }
    }
  }
}
generated quantities {
  // posterior predictive checks
  array[N, B] real D_pred; // predicted durations
  array[N] real bites_pred; // predicted bites
  
  for (n in 1 : N) {
    real risk_draw = risk[n]; // risk from current posterior draw
    real lambda = exp(beta_risk_bites * risk_draw);
    
    for (b in 1 : B) {
      real mu = beta_risk[b] * risk_draw;
      real pi_draw = inv_logit(alpha_pi[b] + beta_pi_risk[b] * risk_draw);
      
      if (bernoulli_rng(pi_draw)) {
        D_pred[n, b] = 0;
      } else {
        D_pred[n, b] = lognormal_rng(mu, sigma_D[b]);
      }
      
      if (b == 1 && D_pred[n, b] > 0) {
        bites_pred[n] = poisson_rng(lambda);
      } else if (b == 1) {
        bites_pred[n] = 0;
      }
    }
  }
}
