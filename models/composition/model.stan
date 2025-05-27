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
}

transformed parameters {
  array[T, S, C] real lambda; // expected counts

  for (t in 1 : T) {
    for (s in 1 : S) {
      for (c in 1 : C) {
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
    count[n] ~ neg_binomial_2(lambda[transect[n], species[n], size[n]], phi_nb);
  }
}

generated quantities {
  array[N] real log_lik;
  array[N] real pred_count;
  array[T, S, C] real abundance;

  for (n in 1 : N) {
    log_lik[n] = neg_binomial_2_lpmf(count[n] | lambda[transect[n], species[n], size[n]], phi_nb);
    pred_count[n] = neg_binomial_2_rng(lambda[transect[n], species[n], size[n]], phi_nb);
  }

  for (t in 1 : T) {
    for (s in 1 : S) {
      for (c in 1 : C) {
        abundance[t, s, c] = neg_binomial_2_rng(lambda[t, s, c], phi_nb);
      }
    }
  }
}
