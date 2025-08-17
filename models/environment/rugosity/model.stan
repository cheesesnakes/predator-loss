data {
  int<lower=0> P; // number of protection levels
  int<lower=0> S; // number of plots
  int<lower=0> L; // number of locations
  
  // response variables
  array[S, 3] real<lower=0, upper=1> rugosity_raw; // rugosity (continuous variable)
  
  // index variables
  array[S] int<lower=1, upper=L> location; // location index

  // plot-level predictors
  array[S] int<lower=1, upper=2> protection; // protection level (0 or 1)
  
}
parameters {
  // Parameters to estimate rugosity
  
  array[S] real<lower=1e-6, upper=(1 - 1e-6)> rugosity; // accounting for multiple rugosity measurements
  array[S] real<lower=1e-6> phi_rug; // standard deviation of rugosity in each plot

  // parameters for effect of protection on rugosity

  array[L] real alpha; // random intercepts for each location
  array[P] real beta; // effect of protection on rugosity
  real<lower=1e-6> psi; // standard deviation of rugosity in each protection level

}
transformed parameters {
  // transformed parameters for protection model

  array[S] real mu; // location specific mean rugosity

  for (s in 1:S){
    int l = location[s];

    mu[s] = inv_logit(alpha[l] + beta[protection[s]]); // mean rugosity for each location
  }

}
model {
  // Priors
  for (s in 1 : S) {
    rugosity[s] ~ beta(1, 4);
    phi_rug[s] ~ exponential(1);
  }

  alpha ~ normal(0, 1);
  beta[1] ~ normal(0, 1);
  beta[2] ~ normal(0.5, 1);
  psi ~ exponential(1);

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
  
  for (s in 1 : S) {
    rugosity[s] ~ beta(mu[s]*psi, (1-mu[s])*psi);
}
}

generated quantities {
  // posterior predictive checks
  array[S] real rugosity_pred;
  for (s in 1:S) {

    int l = location[s];
    real a = alpha[l];
    int p = protection[s];
    real b = beta[p];

    real m = inv_logit(a + b);

    rugosity_pred[s] = beta_rng(m * psi, (1 - m) * psi);
  }

  // counterfactual rugosity by protection

  array[P, S] real rugosity_cf;

  for (p in 1:P) {
    for (s in 1:S) {
      
    int l = location[s];
    real a = alpha[l];
    real b = beta[p];

    real m = inv_logit(a + b);

    rugosity_cf[p, s] = beta_rng(m * psi, (1 - m) * psi);
  }
}
}