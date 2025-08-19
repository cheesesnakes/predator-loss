data {
  int<lower=0> P; // number of protection levels
  int<lower=0> S; // number of transects
  int<lower=0> L; // number of locations
  int<lower=0> Q; // total number of quadrats sampled

  // response
  array[Q] real<lower=0, upper=1> biomass; // biomass cover of quadrat
  
  // index variables
  array[S] int<lower=1, upper=L> location; // location index
  array[Q] int<lower=1, upper=S> transect_quadrat; // transect index for quadrats
  
  // plot-level predictors
  array[S] int<lower=1, upper=2> protection; // protection level (0 or 1)
  
}
parameters {
  // parameters for effect of protection on biomass cover

  array[L] real alpha; // random intercepts for each location
  array[P] real beta; // effect of protection on biomass cover
  real<lower=1e-6> psi; // standard deviation of biomass cover in each protection level

  // sample estimation
  array[S] real<lower=1e-6, upper=1 - 1e-6> R; // resource availability estimate
  array[S] real<lower=1e-6> pi_R; // precision of resource availability
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

  // Priors for habitat complexity and resource availability
  R ~ beta(1, 1); // uniform prior for resource availability
  pi_R ~ exponential(1);
  
  // estimate habitat complexity and resource availability
  for (q in 1 : Q) {
    biomass[q] ~ beta(R[transect_quadrat[q]] * pi_R[transect_quadrat[q]],
                      (1 - R[transect_quadrat[q]])
                      * pi_R[transect_quadrat[q]]);
  }
  // Likelihood
  
  for (s in 1 : S) {
    R[s] ~ beta(mu[s]*psi, (1-mu[s])*psi);
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