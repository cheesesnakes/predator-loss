data {
  int<lower=0> P; // number of protection levels
  int<lower=0> S; // number of plots
  int<lower=0> L; // number of locations
  
  // response
  array[S] real<lower=0, upper=1> biomass; // biomass cover of plot
  
  // index variables
  array[S] int<lower=1, upper=L> location; // location index

  // plot-level predictors
  array[S] int<lower=1, upper=2> protection; // protection level (0 or 1)
  
}
parameters {
  // parameters for effect of protection on biomass cover

  array[L] real alpha; // random intercepts for each location
  array[P] real beta; // effect of protection on biomass cover
  real<lower=1e-6> psi; // standard deviation of biomass cover in each protection level

}
transformed parameters {
  // transformed parameters for protection model

  array[S] real mu; // location specific mean biomass cover

  for (s in 1:S){
    int l = location[s];

    mu[s] = inv_logit(alpha[l] + beta[protection[s]]); // mean biomass cover for each location
  }

}
model {
  // Priors
  
  alpha ~ normal(0, 1);
  beta[1] ~ normal(0, 1);
  beta[2] ~ normal(0.5, 1);
  psi ~ exponential(1);

  
  // Likelihood
  
  for (s in 1 : S) {
    biomass[s] ~ beta(mu[s]*psi, (1-mu[s])*psi);
}
}

generated quantities {
  // posterior predictive checks
  array[S] real biomass_pred;
  for (s in 1:S) {

    int l = location[s];
    real a = alpha[l];
    int p = protection[s];
    real b = beta[p];

    real m = inv_logit(a + b);

    biomass_pred[s] = beta_rng(m * psi, (1 - m) * psi);
  }

  // counterfactual biomass by protection

  array[P, S] real biomass_cf;

  for (p in 1:P) {
    for (s in 1:S) {
      
    int l = location[s];
    real a = alpha[l];
    real b = beta[p];

    real m = inv_logit(a + b);

    biomass_cf[p, s] = beta_rng(m * psi, (1 - m) * psi);
  }
}
}