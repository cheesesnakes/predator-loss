data {
  int N; // number of observations
  int S; // number of species
  int L; // number of locations
  int T; // number of deployments
  int P; // number of protection levels
  int G; // number of feeding guilds
  int C; // number of size classes
  
  array[N] int<lower=1, upper=S> species; // species index
  array[T] int<lower=1, upper=L> location; // location index
  array[N] int<lower=1, upper=T> transect; // transect index
  array[T] int<lower=1, upper=P> protection; // protection level index
  array[S, G] real<lower=0, upper=1> guild; // feeding guild membership matrix
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
  real<lower=0> phi_nb; // dispersion parameters for negative binomial
}
model {
  // Priors for count
  alpha ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5);
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
    
    real lambda_n = exp(log_lambda);
    
    target += neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb);
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
                      + dot_product(epsilon[c], guild[s]) + theta[c];
    
    real lambda_n = exp(log_lambda);
    
    log_lik[n] = neg_binomial_2_lpmf(count[n] | lambda_n, phi_nb);
    
    pred_count[n] = neg_binomial_2_rng(lambda_n, phi_nb); // non-zero case
  }
  
  // counterfactual abundance estimation
  
  for (pr in 1 : P) {
    for (t in 1 : T) {
      for (s in 1 : S) {
        for (c in 1 : C) {
          int p = pr;
          int l = location[t];
          
          real log_lambda = alpha[l] + dot_product(gamma, guild[s])
                            + delta[p] + dot_product(beta[p], guild[s])
                            + phi[p, c] + dot_product(epsilon[c], guild[s]);
          
          // Non-zero expected abundance
          abundance[p, t, s, c] = neg_binomial_2_rng(exp(log_lambda), phi_nb); // non-zero case
        }
      }
    }
  }
}
