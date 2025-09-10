data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of protection levels
  int<lower=0> B; // number of behavioral types
  int<lower=0> S; // number of plots
  int<lower=0> C; // number of size classes
  int<lower=0> G; // number of guilds
  int<lower=0> K; // number of species
  int<lower=0> L; // number of locations sampled
  
  // response variables
  
  array[N, B] real<lower=0> D; // observed durations
  array[N] int<lower=0> bites; // number of bites observed
  
  // index variables
  array[N] int<lower=1, upper=S> plot; // plot index
  array[N] int<lower=1, upper=K> species; // species index
  array[S] int<lower=1, upper=L> location; // location index
  
  // plot-level predictors
  array[S] int<lower=1, upper=2> protection; // protection level (0 or 1)
  array[S] int<lower=1, upper=2> predator; // predator presence (0 or 1)
  array[S, 3] real<lower=0, upper=1> rugosity_raw; // rugosity (continuous variable)
  array[S] real<lower=0, upper=1> biomass; // biomass availability (continuous variable)
  
  // individual-level predictors
  array[N] int<lower=1, upper=C> size; // size class of each observation
  array[N] int<lower=1, upper=2> group; // group membership (0 or 1)
  array[K, G] int<lower=0, upper=1> guild; // guild of each observation
  
  // simulation parameters
  int<lower=0> R; // number of simulations for counterfactuals
  array[R] real sim_var; // standarised values for counterfactuals
}
parameters {
  // Parameters to estimate rugosity
  
  array[S] real<lower=1e-6, upper=(1 - 1e-6)> rugosity; // accounting for multiple rugosity measurements
  array[S] real<lower=1e-6> phi_rug; // standard deviation of rugosity in each plot
  
  // parameters for zero-inflation
  array[S, B] real alpha_pi_plot; // plot-level random intercepts for zero-inflation
  array[B] real<lower=1e-6> sd_alpha_pi_plot; // standard deviation of random intercepts for zero-inflation
  
  array[L, B] real alpha_pi_loc; // location-level random intercepts for zero-inflation
  array[B] real<lower=1e-6> sd_alpha_pi_loc; // standard deviation of random intercepts for zero-inflation
  
  array[P, B] real beta_pi; // effect of protection on zero-inflation
  array[B] real gamma_pi; //effect of rugosity on zero-inflation
  array[B] real zeta_pi; // effect of biomass on zero-inflation
  array[2, B] real eta_pi; // effect of predator presence on zero-inflation
  array[C, B] real omega_pi; // effect of size class on zero-inflation
  array[G, B] real epsilon_pi; // effect of guild on zero-inflation
  array[2, B] real delta_pi; // effect of group on zero-inflation
  
  // parameters for duration model
  array[S, B] real alpha_D_plot; // plot-level random intercepts for duration
  array[B] real<lower=1e-6> sd_alpha_D_plot; // standard deviation of random intercepts for duration
  
  array[L, B] real alpha_D_loc; // location-level random intercepts for duration
  array[B] real<lower=1e-6> sd_alpha_D_loc; // standard deviation of random intercepts for duration
  
  array[P, B] real beta_D; // effect of protection on duration
  array[B] real gamma_D; // effect of rugosity on duration
  array[B] real zeta_D; // effect of biomass on duration
  array[2, B] real eta_D; // effect of predator presence on duration
  array[C, B] real omega_D; // effect of size class on duration
  array[G, B] real epsilon_D; // effect of guild on duration
  array[2, B] real delta_D; // effect of group on duration
  array[B] real<lower=1e-6> sigma_D; // standard deviation of duration
  
  // parameters for bites model
  array[L] real alpha_bites_loc; // location-level random intercepts for bites
  real<lower=1e-6> sd_alpha_bites_loc; // standard deviation of random intercepts for bites
  
  array[S] real alpha_bites_plot; // plot-level random intercepts for bites
  real<lower=1e-6> sd_alpha_bites_plot; // standard deviation of random intercepts for bites
  
  array[P] real beta_bites; // effect of protection on bites
  real gamma_bites; // effect of rugosity on bites
  real zeta_bites; // effect of biomass on bites
  array[2] real eta_bites; // effect of predator presence on bites
  array[C] real omega_bites; // effect of size class on bites
  array[G] real epsilon_bites; // effect of guild on bites
  array[2] real delta_bites; // effect of group on bites
}
transformed parameters {
  // standardize rugosity and biomass
  array[S] real rugosity_std; // standardized rugosity
  array[S] real biomass_std; // standardized biomass
  
  array[N] real log_lambda; // expected number of bites
  array[N, B] real log_mu; // expected duration
  array[N, B] real pi; // probability of non-zero duration
  
  for (s in 1 : S) {
    rugosity_std[s] = (rugosity[s] - mean(rugosity)) / sd(rugosity);
    biomass_std[s] = (biomass[s] - mean(biomass)) / sd(biomass);
  }
  
  for (n in 1 : N) {
    for (b in 1 : B) {
      log_mu[n, b] = alpha_D_plot[plot[n], b]
                     + beta_D[protection[plot[n]], b]
                     + gamma_D[b] * rugosity_std[plot[n]]
                     + zeta_D[b] * biomass_std[plot[n]]
                     + eta_D[predator[plot[n]], b] + omega_D[size[n], b]
                     + dot_product(epsilon_D[ : , b], guild[species[n]])
                     + delta_D[group[n], b];
      
      // calculate probability of non-zero duration
      
      pi[n, b] = inv_logit(
                           alpha_pi_plot[plot[n], b]
                           + beta_pi[protection[plot[n]], b]
                           + gamma_pi[b] * rugosity_std[plot[n]]
                           + zeta_pi[b] * biomass_std[plot[n]]
                           + eta_pi[predator[plot[n]], b]
                           + omega_pi[size[n], b]
                           + dot_product(epsilon_pi[ : , b],
                                         guild[species[n]])
                           + delta_pi[group[n], b]);
    }
    
    // calculate expected number of bites
    
    log_lambda[n] = alpha_bites_plot[plot[n]]
                    + beta_bites[protection[plot[n]]]
                    + gamma_bites * rugosity_std[plot[n]]
                    + zeta_bites * biomass_std[plot[n]]
                    + eta_bites[predator[plot[n]]] + omega_bites[size[n]]
                    + dot_product(epsilon_bites, guild[species[n]])
                    + delta_bites[group[n]];
  }
}
model {
  // Priors
  for (s in 1 : S) {
    rugosity[s] ~ beta(1, 1);
    phi_rug[s] ~ exponential(1);
  }
  
  // Zero-inflation parameters
  for (b in 1 : B) {
    sd_alpha_pi_loc[b] ~ exponential(1);
    sd_alpha_pi_plot[b] ~ exponential(1);
    
    beta_pi[ : , b] ~ normal(0, 1);
    gamma_pi[b] ~ normal(0, 1);
    zeta_pi[b] ~ normal(0, 1);
    eta_pi[ : , b] ~ normal(0, 1);
    omega_pi[ : , b] ~ normal(0, 1);
    epsilon_pi[ : , b] ~ normal(0, 1);
    delta_pi[ : , b] ~ normal(0, 1);
  }
  
  for (l in 1 : L) {
    for (b in 1 : B) {
      alpha_pi_loc[l, b] ~ normal(0, sd_alpha_pi_loc[b]);
    }
  }
  
  for (s in 1 : S) {
    for (b in 1 : B) {
      alpha_pi_plot[s, b] ~ normal(alpha_pi_loc[location[s], b],
                                   sd_alpha_pi_plot[b]);
    }
  }
  
  // Duration parameters
  
  for (b in 1 : B) {
    sd_alpha_D_loc[b] ~ exponential(1);
    sd_alpha_D_plot[b] ~ exponential(1);
    
    beta_D[ : , b] ~ normal(0, 1);
    gamma_D[b] ~ normal(0, 1);
    zeta_D[b] ~ normal(0, 1);
    eta_D[ : , b] ~ normal(0, 1);
    omega_D[ : , b] ~ normal(0, 1);
    epsilon_D[ : , b] ~ normal(0, 1);
    delta_D[ : , b] ~ normal(0, 1);
    sigma_D[b] ~ cauchy(0, 1); // using Cauchy prior for scale
  }
  
  for (l in 1 : L) {
    for (b in 1 : B) {
      alpha_D_loc[l, b] ~ normal(0, sd_alpha_D_loc[b]);
    }
  }
  
  for (s in 1 : S) {
    for (b in 1 : B) {
      alpha_D_plot[s, b] ~ normal(alpha_D_loc[location[s], b],
                                  sd_alpha_D_plot[b]);
    }
  }
  
  // Bites parameters
  sd_alpha_bites_loc ~ exponential(1);
  sd_alpha_bites_plot ~ exponential(1);
  beta_bites ~ normal(0, 1);
  gamma_bites ~ normal(0, 1);
  zeta_bites ~ normal(0, 1);
  eta_bites ~ normal(0, 1);
  omega_bites ~ normal(0, 1);
  epsilon_bites ~ normal(0, 1);
  delta_bites ~ normal(0, 1);
  
  for (l in 1 : L) {
    alpha_bites_loc[l] ~ normal(0, sd_alpha_bites_loc);
  }
  
  for (s in 1 : S) {
    alpha_bites_plot[s] ~ normal(alpha_bites_loc[location[s]],
                                 sd_alpha_bites_plot);
  }
  
  // rugosity estimation
  
  for (s in 1 : S) {
    rugosity_raw[s, 1] ~ beta(rugosity[s] * phi_rug[s],
                              (1 - rugosity[s]) * phi_rug[s]);
    rugosity_raw[s, 2] ~ beta(rugosity[s] * phi_rug[s],
                              (1 - rugosity[s]) * phi_rug[s]);
    rugosity_raw[s, 3] ~ beta(rugosity[s] * phi_rug[s],
                              (1 - rugosity[s]) * phi_rug[s]);
  }
  
  // Likelihood
  
  for (n in 1 : N) {
    for (b in 1 : B) {
      if (D[n, b] == 0) {
        target += bernoulli_lpmf(1 | pi[n, b]); // zero-inflation
      } else {
        target += bernoulli_lpmf(0 | pi[n, b]); // non-zero duration
        target += lognormal_lpdf(D[n, b] | log_mu[n, b], sigma_D[b]); // duration likelihood      
      }
    }
    if (D[n, 1] > 0) {
      target += poisson_log_lpmf(bites[n] | log_lambda[n] + log(D[n, 1])); // bites likelihood for behavior 1
    } else {
      target += poisson_log_lpmf(bites[n] | log(1e-10)); // no bites for other behaviors
    }
  }
}
generated quantities {
  // posterior predictive checks
  array[N, B] real D_pred; // predicted durations
  array[N] real bites_pred; // predicted bites
  
  for (n in 1 : N) {
    real lambda_pred = exp(log_lambda[n]);
    for (b in 1 : B) {
      real mu_pred = log_mu[n, b];
      real pi_pred = pi[n, b];
      if (bernoulli_rng(pi_pred)) {
        D_pred[n, b] = 0;
      } else {
        D_pred[n, b] = lognormal_rng(mu_pred, sigma_D[b]);
      }
    }
    if (D_pred[n, 1] > 0) {
      bites_pred[n] = poisson_rng(lambda_pred * D_pred[n, 1]);
    } else {
      // if no duration, set bites to a small value
      bites_pred[n] = 0; // or some other small value
    }
  }
  // counterfactual protection 
  array[P, N, B] real D_pred_protection; // predicted durations with protection
  array[P, N] real bites_pred_protection; // predicted bites with protection
  
  for (p in 1 : P) {
    for (n in 1 : N) {
      real rug_std = mean(rugosity_std);
      real bio_std = mean(biomass_std);
      for (b in 1 : B) {
        real mu_cf = alpha_D_plot[plot[n], b] + beta_D[p, b]
                     + gamma_D[b] * rug_std + zeta_D[b] * bio_std
                     + eta_D[1, b]
                     + // assuming predator present
                     omega_D[size[n], b]
                     + dot_product(epsilon_D[ : , b], guild[species[n]])
                     + delta_D[group[n], b];
        
        real pi_cf = inv_logit(
                               alpha_pi_plot[plot[n], b] + beta_pi[p, b]
                               + gamma_pi[b] * rug_std + zeta_pi[b] * bio_std
                               + eta_pi[1, b] + omega_pi[size[n], b]
                               + dot_product(epsilon_pi[ : , b],
                                             guild[species[n]])
                               + delta_pi[group[n], b]);
        
        if (bernoulli_rng(pi_cf)) {
          D_pred_protection[p, n, b] = 0;
        } else {
          D_pred_protection[p, n, b] = lognormal_rng(mu_cf, sigma_D[b]);
        }
      }
      
      real lambda_cf = exp(
                           alpha_bites_plot[plot[n]] + beta_bites[p]
                           + gamma_bites * rug_std + zeta_bites * bio_std
                           + eta_bites[1] + omega_bites[size[n]]
                           + dot_product(epsilon_bites, guild[species[n]])
                           + delta_bites[group[n]]);
      
      bites_pred_protection[p, n] = lambda_cf;
    }
  }
  
  // counterfactual protection varying rugosity
  array[P, R, N, B] real D_pred_protection_rugosity;
  array[P, R, N] real bites_pred_protection_rugosity;
  
  for (p in 1 : P) {
    for (n in 1 : N) {
      real bio_std = mean(biomass_std);
      for (r in 1 : R) {
        real rug_std = sim_var[r];
        for (b in 1 : B) {
          real mu_cf = alpha_D_plot[plot[n], b] + beta_D[p, b]
                       + gamma_D[b] * rug_std + zeta_D[b] * bio_std
                       + eta_D[1, b] + omega_D[size[n], b]
                       + dot_product(epsilon_D[ : , b], guild[species[n]])
                       + delta_D[group[n], b];
          
          real pi_cf = inv_logit(
                                 alpha_pi_plot[plot[n], b] + beta_pi[p, b]
                                 + gamma_pi[b] * rug_std
                                 + zeta_pi[b] * bio_std + eta_pi[1, b]
                                 + omega_pi[size[n], b]
                                 + dot_product(epsilon_pi[ : , b],
                                               guild[species[n]])
                                 + delta_pi[group[n], b]);
          
          if (bernoulli_rng(pi_cf)) {
            D_pred_protection_rugosity[p, r, n, b] = 0;
          } else {
            D_pred_protection_rugosity[p, r, n, b] = lognormal_rng(mu_cf,
                                                                   sigma_D[b]);
          }
        }
        real lambda_cf = exp(
                             alpha_bites_plot[plot[n]] + beta_bites[p]
                             + gamma_bites * rug_std + zeta_bites * bio_std
                             + eta_bites[1] + omega_bites[size[n]]
                             + dot_product(epsilon_bites, guild[species[n]])
                             + delta_bites[group[n]]);
        
        bites_pred_protection_rugosity[p, r, n] = lambda_cf;
      }
    }
  }
  
  // counterfactual protection varying biomass
  array[P, R, N, B] real D_pred_protection_biomass;
  array[P, R, N] real bites_pred_protection_biomass;
  
  for (p in 1 : P) {
    for (n in 1 : N) {
      real rug_std = mean(rugosity_std);
      for (r in 1 : R) {
        real bio_std = sim_var[r];
        for (b in 1 : B) {
          real mu_cf = alpha_D_plot[plot[n], b] + beta_D[p, b]
                       + gamma_D[b] * rug_std + zeta_D[b] * bio_std
                       + eta_D[1, b] + omega_D[size[n], b]
                       + dot_product(epsilon_D[ : , b], guild[species[n]])
                       + delta_D[group[n], b];
          
          real pi_cf = inv_logit(
                                 alpha_pi_plot[plot[n], b] + beta_pi[p, b]
                                 + gamma_pi[b] * rug_std
                                 + zeta_pi[b] * bio_std + eta_pi[1, b]
                                 + omega_pi[size[n], b]
                                 + dot_product(epsilon_pi[ : , b],
                                               guild[species[n]])
                                 + delta_pi[group[n], b]);
          
          if (bernoulli_rng(pi_cf)) {
            D_pred_protection_biomass[p, r, n, b] = 0;
          } else {
            D_pred_protection_biomass[p, r, n, b] = lognormal_rng(mu_cf,
                                                                  sigma_D[b]);
          }
        }
        real lambda_cf = exp(
                             alpha_bites_plot[plot[n]] + beta_bites[p]
                             + gamma_bites * rug_std + zeta_bites * bio_std
                             + eta_bites[1] + omega_bites[size[n]]
                             + dot_product(epsilon_bites, guild[species[n]])
                             + delta_bites[group[n]]);
        bites_pred_protection_biomass[p, r, n] = lambda_cf;
      }
    }
  }
  
  // counterfactual protection varying predator presence
  array[P, 2, N, B] real D_pred_protection_predator;
  array[P, 2, N] real bites_pred_protection_predator;
  
  for (p in 1 : P) {
    for (n in 1 : N) {
      real rug_std = mean(rugosity_std);
      real bio_std = mean(biomass_std);
      for (g in 1 : 2) {
        for (b in 1 : B) {
          real mu_cf = alpha_D_plot[plot[n], b] + beta_D[p, b]
                       + gamma_D[b] * rug_std + zeta_D[b] * bio_std
                       + eta_D[g, b]
                       + // assuming predator absent
                       omega_D[size[n], b]
                       + dot_product(epsilon_D[ : , b], guild[species[n]])
                       + delta_D[group[n], b];
          
          real pi_cf = inv_logit(
                                 alpha_pi_plot[plot[n], b] + beta_pi[p, b]
                                 + gamma_pi[b] * rug_std
                                 + zeta_pi[b] * bio_std + eta_pi[g, b]
                                 + omega_pi[size[n], b]
                                 + dot_product(epsilon_pi[ : , b],
                                               guild[species[n]])
                                 + delta_pi[group[n], b]);
          
          if (bernoulli_rng(pi_cf)) {
            D_pred_protection_predator[p, g, n, b] = 0;
          } else {
            D_pred_protection_predator[p, g, n, b] = lognormal_rng(mu_cf,
                                                                   sigma_D[b]);
          }
        }
        real lambda_cf = exp(
                             alpha_bites_plot[plot[n]] + beta_bites[p]
                             + gamma_bites * rug_std + zeta_bites * bio_std
                             + eta_bites[g]
                             // assuming predator absent
                             + omega_bites[size[n]]
                             + dot_product(epsilon_bites, guild[species[n]])
                             + delta_bites[group[n]]);
        
        bites_pred_protection_predator[p, g, n] = lambda_cf;
      }
    }
  }
}
