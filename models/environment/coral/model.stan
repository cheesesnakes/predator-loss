data {
  int<lower=0> P; // number of protection levels
  int<lower=0> S; // number of transects
  int<lower=0> L; // number of locations
  int<lower=0> Q; // total number of quadrats sampled
  
  // response
  array[Q] real<lower=0, upper=1> coral; // coral cover of plot
  
  // index variables
  array[S] int<lower=1, upper=L> location; // site index
  array[Q] int<lower=1, upper=S> transect_quadrat; // transect index for quadrats
  
  // plot-level predictors
  array[S] int<lower=0, upper=1> protection; // protection level (0 or 1)
}
parameters {
  // parameters for effect of protection on coral cover
  
  array[S] real alpha_s; // random intercepts for each transect
  real<lower=0> sigma_s; // standard deviation of transect intercepts
  real beta; // effect of protection on coral cover
  real<lower=1e-6> psi; // standard deviation of coral cover in each protection level
}
transformed parameters {
  // transformed parameters for protection model
  
  array[Q] real mu; // location specific mean coral cover
  array[S] real a_s; // non-centered alpha_s
  
  for (q in 1 : Q) {
    int t = transect_quadrat[q];
    int l = location[t];
    
    a_s[t] = alpha_s[t] * sigma_s;
    
    mu[q] = inv_logit(alpha_s[t] + beta * protection[t]); // mean coral cover for each location
  }
}
model {
  // Priors
  
  alpha_s ~ normal(0, 1);
  sigma_s ~ exponential(1);
  beta ~ normal(0, 1);
  psi ~ exponential(1);
  
  // Likelihood
  
  for (q in 1 : Q) {
    coral[q] ~ beta(mu[q] * psi, (1 - mu[q]) * psi);
  }
}
generated quantities {
  // posterior predictive checks
  array[Q] real coral_pred;
  for (q in 1 : Q) {
    int s = transect_quadrat[q];
    int l = location[s];
    real as = a_s[s];
    int p = protection[s];
    
    real m = inv_logit(as + beta * p);
    
    coral_pred[q] = beta_rng(m * psi, (1 - m) * psi);
  }
  
  // counterfactual coral by protection
  
  array[P, Q] real coral_cf;
  
  for (p in 1 : P) {
    for (q in 1 : Q) {
      int s = transect_quadrat[q];
      int l = location[s];
      real as = a_s[s];
      
      real m = inv_logit(as + beta * (p - 1));
      
      coral_cf[p, q] = beta_rng(m * psi, (1 - m) * psi);
    }
  }
}
