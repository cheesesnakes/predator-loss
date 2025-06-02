data {
  int N; // number of observations
  int S; // number of species
  int L; // number of locations
  int T; // number of transects
  int P; // number of protection levels
  int G; // number of feeding guilds
  int C; // number of size classes
  int Q; // total number of quadrats sampled
  
  array[N] int<lower=1, upper=S> species; // species index
  array[T] int<lower=1, upper=L> location; // location index
  array[N] int<lower=1, upper=T> transect; // transect index
  array[T] int<lower=1, upper=P> protection; // protection level index
  array[S, G] real<lower=0, upper=1> guild; // feeding guild membership matrix
  array[N] int<lower=1, upper=C> size; // size class index
  array[N] int<lower=0> count; // count of individuals
  
  array[Q] int<lower=1, upper=T> transect_quadrat; // transect index for quadrats
  array[Q] real<lower=0, upper=1> coral; // coral cover
  array[Q] real<lower=0, upper=1> biomass; // biomass 
}
parameters {
  array[L] real alpha; // intercepts for locations
  array[G] real gamma; // effect of feeding guilds
  array[P] real delta; // effect of protection levels
  array[P, G] real beta; // interaction between protection and guild
  array[C] real theta; // size class effects
  array[P, C] real phi; // interaction between protection and size class
  array[C, G] real epsilon; // interaction between size class and guild
  real omega; // effect of habitat complexity
  array[G] real omega_G; // effect of habitat complexity by guild
  array[C] real omega_C; // effect of habitat complexity by size class
  real tau; // effect of resource availability
  array[G] real tau_G; // effect of resource availability by guild
  array[C] real tau_C; // effect of resource availability by size class
  array[L] real<lower=0> phi_nb; // location-specific dispersion parameters for negative binomial
  
  array[L] real alpha_psi; // zero-inflation intercepts for locations
  array[G] real gamma_psi; // zero-inflation effect of feeding guilds
  array[P] real delta_psi; // zero-inflation effect of protection levels
  array[P, G] real beta_psi; // zero-inflation interaction between protection and guild
  array[C] real theta_psi; // zero-inflation size class effects
  array[C, G] real epsilon_psi; // zero-inflation interaction between size class and guild
  array[P, C] real phi_psi; // zero-inflation interaction between protection and size class
  real omega_psi; // zero-inflation effect of habitat complexity
  array[G] real omega_psi_G; // zero-inflation effect of habitat complexity by guild
  array[C] real omega_psi_C; // zero-inflation effect of habitat complexity by size class
  real tau_psi; // zero-inflation effect of resource availability
  array[G] real tau_psi_G; // zero-inflation effect of resource availability by guild
  array[C] real tau_psi_C; // zero-inflation effect of resource availability by size class
  
  array[T] real<lower=1e-6, upper=1 - 1e-6> H; // habitat complexity estimate
  array[T] real<lower=1e-6> pi_H; //precision of habitat complexity
  array[T] real<lower=1e-6, upper=1 - 1e-6> R; // resource availability estimate
  array[T] real<lower=1e-6> pi_R; // precision of resource availability
}
model {
  // Priors for count
  alpha ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5);
  omega ~ normal(0, 0.5); // habitat complexity effect
  omega_G ~ normal(0, 0.5); // habitat complexity effect by guild
  omega_C ~ normal(0, 0.5); // habitat complexity effect by size class
  tau ~ normal(0, 0.5); // resource availability effect
  tau_G ~ normal(0, 0.5); // resource availability effect by guild
  tau_C ~ normal(0, 0.5); // resource availability effect by size class
  phi_nb ~ exponential(1); // dispersion prior
  
  for (p in 1 : P) {
    phi[p] ~ normal(0, 0.5);
    for (g in 1 : G) {
      beta[p, g] ~ normal(0, 0.5);
    }
  }
  
  for (c in 1 : C) {
    theta[c] ~ normal(0, 0.5);
    for (g in 1 : G) {
      epsilon[c, g] ~ normal(0, 0.5);
    }
  }
  
  // Priors for zero-inflation
  
  alpha_psi ~ normal(0, 0.5);
  gamma_psi ~ normal(0, 0.5);
  delta_psi ~ normal(0, 0.5);
  theta_psi ~ normal(0, 0.5);
  omega_psi ~ normal(0, 0.5); // habitat complexity effect
  omega_psi_G ~ normal(0, 0.5); // habitat complexity effect by guild
  omega_psi_C ~ normal(0, 0.5); // habitat complexity effect by size class
  tau_psi ~ normal(0, 0.5); // resource availability effect
  tau_psi_G ~ normal(0, 0.5); // resource availability effect by guild
  tau_psi_C ~ normal(0, 0.5); // resource availability effect by size class
  
  for (c in 1 : C) {
    for (g in 1 : G) {
      epsilon_psi[c, g] ~ normal(0, 0.5);
    }
  }
  
  for (p in 1 : P) {
    beta_psi[p] ~ normal(0, 0.5);
    phi_psi[p] ~ normal(0, 0.5);
  }
  
  // Priors for habitat complexity and resource availability
  H ~ beta(1, 1); // uniform prior for habitat complexity
  R ~ beta(1, 1); // uniform prior for resource availability
  pi_H ~ exponential(1);
  pi_R ~ exponential(1);
  
  // estimate habitat complexity and resource availability
  for (q in 1 : Q) {
    coral[q] ~ beta(H[transect_quadrat[q]] * pi_H[transect_quadrat[q]],
                    (1 - H[transect_quadrat[q]]) * pi_H[transect_quadrat[q]]);
    biomass[q] ~ beta(R[transect_quadrat[q]] * pi_R[transect_quadrat[q]],
                      (1 - R[transect_quadrat[q]])
                      * pi_R[transect_quadrat[q]]);
  }
  
  // Likelihood
  for (n in 1 : N) {
    int t = transect[n];
    int s = species[n];
    int c = size[n];
    int p = protection[t];
    int l = location[t];
    
    real log_lambda = alpha[l] + dot_product(gamma, guild[s]) + delta[p]
                      + dot_product(beta[p], guild[s]) + phi[p, c]
                      + dot_product(epsilon[c], guild[s]) + theta[c]
                      + omega * H[t] + tau * R[t]
                      + dot_product(omega_G, guild[s]) * H[t]
                      + omega_C[c] * H[t]
                      + dot_product(tau_G, guild[s]) * R[t] + tau_C[c] * R[t];
    
    real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                     + delta_psi[p] + dot_product(beta_psi[p], guild[s])
                     + dot_product(epsilon_psi[c], guild[s]) + phi_psi[p, c]
                     + theta_psi[c] + omega_psi * H[t] + tau_psi * R[t]
                     + dot_product(omega_psi_G, guild[s]) * H[t]
                     + omega_psi_C[c] * H[t]
                     + dot_product(tau_psi_G, guild[s]) * R[t]
                     + tau_psi_C[c] * R[t];
    
    real lambda_n = exp(log_lambda);
    real psi_n = inv_logit(logit_psi);
    
    if (count[n] == 0) {
      target += bernoulli_lpmf(1 | psi_n);
    } else {
      target += bernoulli_lpmf(0 | psi_n)
                + neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb[l]);
    }
  }
}
generated quantities {
  array[N] real log_lik;
  array[N] real pred_count;
  array[P, T, S, C] real abundance;
  
  // log likelihood and predicted counts
  for (n in 1 : N) {
    int t = transect[n];
    int s = species[n];
    int c = size[n];
    int p = protection[t];
    int l = location[t];
    
    real log_lambda = alpha[l] + dot_product(gamma, guild[s]) + delta[p]
                      + dot_product(beta[p], guild[s]) + phi[p, c]
                      + dot_product(epsilon[c], guild[s]) + theta[c]
                      + omega * H[t] + tau * R[t]
                      + dot_product(omega_G, guild[s]) * H[t]
                      + omega_C[c] * H[t]
                      + dot_product(tau_G, guild[s]) * R[t] + tau_C[c] * R[t];
    
    real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                     + delta_psi[p] + dot_product(beta_psi[p], guild[s])
                     + phi_psi[p, c] + dot_product(epsilon_psi[c], guild[s])
                     + theta_psi[c] + omega_psi * H[t] + tau_psi * R[t]
                     + dot_product(omega_psi_G, guild[s]) * H[t]
                     + omega_psi_C[c] * H[t]
                     + dot_product(tau_psi_G, guild[s]) * R[t]
                     + tau_psi_C[c] * R[t];
    
    real lambda_n = exp(log_lambda);
    real psi_n = inv_logit(logit_psi);
    
    if (count[n] == 0) {
      log_lik[n] = bernoulli_lpmf(1 | psi_n);
    } else {
      log_lik[n] = bernoulli_lpmf(0 | psi_n)
                   + neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb[l]);
    }
    
    if (count[n] == 0) {
      pred_count[n] = 0; // zero-inflated case
    } else {
      pred_count[n] = neg_binomial_2_rng(lambda_n, phi_nb[l]); // non-zero case
    }
  }
  
  // counterfactual abundance estimation
  
  real mean_H = mean(H);
  real mean_R = mean(R);
  
  for (pr in 1 : P) {
    for (t in 1 : T) {
      for (s in 1 : S) {
        for (c in 1 : C) {
          int p = pr;
          int l = location[t];
          
          real log_lambda = alpha[l] + dot_product(gamma, guild[s])
                            + delta[p] + dot_product(beta[p], guild[s])
                            + phi[p, c] + dot_product(epsilon[c], guild[s])
                            + theta[c] + omega * mean_H + tau * mean_R
                            + dot_product(omega_G, guild[s]) * mean_H
                            + omega_C[c] * mean_H
                            + dot_product(tau_G, guild[s]) * mean_R
                            + tau_C[c] * mean_R;
          
          real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                           + delta_psi[p]
                           + dot_product(beta_psi[p], guild[s])
                           + phi_psi[p, c]
                           + dot_product(epsilon_psi[c], guild[s])
                           + theta_psi[c] + omega_psi * mean_H
                           + tau_psi * mean_R
                           + dot_product(omega_psi_G, guild[s]) * mean_H
                           + omega_psi_C[c] * mean_H
                           + dot_product(tau_psi_G, guild[s]) * mean_R
                           + tau_psi_C[c] * mean_R;
          
          real psi_tsc = inv_logit(logit_psi);
          
          if (bernoulli_rng(psi_tsc)) {
            abundance[p, t, s, c] = 0; // zero-inflated case
          } else {
            // Non-zero expected abundance
            abundance[p, t, s, c] = neg_binomial_2_rng(exp(log_lambda),
                                                       phi_nb[l]); // non-zero case
          }
        }
      }
    }
  }
}
