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
  real<lower=0> phi_nb; // dispersion parameter for NB

  array[L] real alpha_psi; // zero-inflation intercepts for locations
  array[G] real gamma_psi; // zero-inflation effect of feeding guilds
  array[P] real delta_psi; // zero-inflation effect of protection levels
  array[P, G] real beta_psi; // zero-inflation interaction between protection and guild
  array[C] real theta_psi; // zero-inflation size class effects
  array[P, C] real phi_psi; // zero-inflation interaction between protection and size class
}

transformed parameters {
  array[T, S, C] real lambda; // expected counts
  array[T, S, C] real psi; // zero-inflation proabability


  for (t in 1 : T) {
    for (s in 1 : S) {
      for (c in 1 : C) {

        // Zero-inflation probability
        psi[t, s, c] = inv_logit(
          alpha_psi[location[t]] +
          dot_product(gamma_psi, guild[s]) +
          delta_psi[protection[t]] +
          dot_product(beta_psi[protection[t]], guild[s]) +
          phi_psi[protection[t], c] +
          theta_psi[c]
        );

        // Expected counts
        lambda[t, s, c] = exp(
          alpha[location[t]] +
          dot_product(gamma, guild[s]) +
          delta[protection[t]] +
          dot_product(beta[protection[t]], guild[s]) +
          phi[protection[t], c] +
          theta[c]
        );
      }
    }
  }
}

model {
  // Priors
  alpha ~ normal(0, 1);
  gamma ~ normal(0, 1);
  delta ~ normal(0, 1);
  theta ~ normal(0, 1);
  phi_nb ~ exponential(1); // dispersion prior

  for (p in 1 : P)
    phi[p] ~ normal(0, 1);

  for (p in 1 : P)
    for (g in 1 : G)
      beta[p, g] ~ normal(0, 1);

  // Likelihood
  for (n in 1 : N) {
    if (count[n] == 0) {
      target += bernoulli_lpmf(1 | psi[transect[n], species[n], size[n]]);
    } else {
      target += bernoulli_lpmf(0 | psi[transect[n], species[n], size[n]]) +
                neg_binomial_2_lpmf(count[n] | lambda[transect[n], species[n], size[n]], phi_nb);
    }
  }
}

generated quantities {
  array[N] real log_lik;
  array[N] real pred_count;
  array[T, S, C] real abundance;

  for (n in 1 : N) {
    if (count[n] == 0) {
      log_lik[n] = log(psi[transect[n], species[n], size[n]]);
    } else {
      log_lik[n] = log1m(psi[transect[n], species[n], size[n]]) +
                    neg_binomial_2_lpmf(count[n] | lambda[transect[n], species[n], size[n]], phi_nb);
    }
    
    real mu = lambda[transect[n], species[n], size[n]];
    real psi_value = psi[transect[n], species[n], size[n]];

    if (bernoulli_rng(psi_value)) {
      pred_count[n] = 0; // zero-inflated case
    } else {
      pred_count[n] = neg_binomial_2_rng(mu, phi_nb); // non-zero case
    }
}

  for (t in 1 : T) {
    for (s in 1 : S) {
      for (c in 1 : C) {
        if (bernoulli_rng(psi[t, s, c])) {
          abundance[t, s, c] = 0; // zero-inflated case
        } else {
          // Non-zero case
        abundance[t, s, c] = neg_binomial_2_rng(lambda[t, s, c], phi_nb);
      }
      }
    }
  }
}
