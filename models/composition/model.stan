data {
  int N; // number of observations
  int S; // number of species
  int L; // number of locations
  int T; // number of transects
  int P; // number of protection levels
  int G; // number of feeding guilds
  int C; // number of size classes
  
  array[N] int<lower=1, upper=S> species; // species index
  array[T] int<lower=1, upper=L> location; // location index
  array[N] int<lower=1, upper=T> transect; // transect index
  array[T] int<lower=1, upper=P> protection; // protection level index
  array[S, G] int<lower=0, upper=1> guild; // feeding guild membership matrix
  array[N] int<lower=1, upper=C> size; // size class index
  array[N] int<lower=0> count; // count of individuals
}
parameters {
  array[L] real alpha; // intercepts for locations
  array[G] real gamma; // effect of feeding guilds
  array[P] real delta; // effect of protection levels
  array[P, G] real beta; // interaction between protection and guild
  array[C] real theta; // size class effects
  array[P, C] real phi; // interaction between protection and size class
  array[C, G] real epsilon; // interaction between size class and guild
  real<lower=0> phi_nb; // dispersion parameter for NB
  
  array[L] real alpha_psi; // zero-inflation intercepts for locations
  array[G] real gamma_psi; // zero-inflation effect of feeding guilds
  array[P] real delta_psi; // zero-inflation effect of protection levels
  array[P, G] real beta_psi; // zero-inflation interaction between protection and guild
  array[C] real theta_psi; // zero-inflation size class effects
  array[C, G] real epsilon_psi; // zero-inflation interaction between size class and guild
  array[P, C] real phi_psi; // zero-inflation interaction between protection and size class
}
model {
  // Priors
  alpha ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5);
  theta ~ normal(0, 0.5);
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
  
  alpha_psi ~ normal(0, 0.5);
  gamma_psi ~ normal(0, 0.5);
  delta_psi ~ normal(0, 0.5);
  theta_psi ~ normal(0, 0.5);
  
  for (c in 1 : C) {
    epsilon_psi[c] ~ normal(0, 0.5);
    for (g in 1 : G) {
      epsilon_psi[c, g] ~ normal(0, 0.5);
    }
  }
  
  for (p in 1 : P) {
    beta_psi[p] ~ normal(0, 0.5);
    phi_psi[p] ~ normal(0, 0.5);
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
                      + dot_product(epsilon[c], guild[s]) + theta[c];
    
    real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                     + delta_psi[p] + dot_product(beta_psi[p], guild[s])
                     + dot_product(epsilon_psi[c], guild[s]) + phi_psi[p, c]
                     + theta_psi[c];
    
    real lambda_n = exp(log_lambda);
    real psi_n = inv_logit(logit_psi);
    
    if (count[n] == 0) {
      target += bernoulli_lpmf(1 | psi_n);
    } else {
      target += bernoulli_lpmf(0 | psi_n)
                + neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb);
    }
  }
}
generated quantities {
  array[N] real log_lik;
  array[N] real pred_count;
  array[T, S, C] real abundance;
  
  for (n in 1 : N) {
    int t = transect[n];
    int s = species[n];
    int c = size[n];
    int p = protection[t];
    int l = location[t];
    
    real log_lambda = alpha[l] + dot_product(gamma, guild[s]) + delta[p]
                      + dot_product(beta[p], guild[s]) + phi[p, c]
                      + dot_product(epsilon[c], guild[s]) + theta[c];
    
    real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                     + delta_psi[p] + dot_product(beta_psi[p], guild[s])
                     + phi_psi[p, c] + dot_product(epsilon_psi[c], guild[s])
                     + theta_psi[c];
    
    real lambda_n = exp(log_lambda);
    real psi_n = inv_logit(logit_psi);
    
    log_lik[n] = bernoulli_lpmf(1 | psi_n) * (count[n] == 0)
                 + bernoulli_lpmf(0 | psi_n) * (count[n] > 0)
                 + neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb)
                   * (count[n] > 0);
    
    if (count[n] == 0) {
      pred_count[n] = 0; // zero-inflated case
    } else {
      pred_count[n] = neg_binomial_2_rng(lambda_n, phi_nb); // non-zero case
    }
  }
  
  for (t in 1 : T) {
    for (s in 1 : S) {
      for (c in 1 : C) {
        int p = protection[t];
        int l = location[t];
        real log_lambda = alpha[l] + dot_product(gamma, guild[s]) + delta[p]
                          + dot_product(beta[p], guild[s]) + phi[p, c]
                          + dot_product(epsilon[c], guild[s]) + theta[c];
        real logit_psi = alpha_psi[l] + dot_product(gamma_psi, guild[s])
                         + delta_psi[p] + dot_product(beta_psi[p], guild[s])
                         + phi_psi[p, c]
                         + dot_product(epsilon_psi[c], guild[s])
                         + theta_psi[c];
        real psi_tsc = inv_logit(logit_psi);
        
        if (bernoulli_rng(psi_tsc)) {
          abundance[t, s, c] = 0; // zero-inflated case
        } else {
          // Non-zero expected abundance
          abundance[t, s, c] = neg_binomial_2_rng(exp(log_lambda), phi_nb); // non-zero case
        }
      }
    }
  }
}
