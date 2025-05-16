# simulate predator presence based on protection
predator_presence <- function(protection) {
    return(rbinom(1, 1, beta_protection_predator[protection])+1)

}


# simulate rugosity

rugosity <- function(protection) {
    D_max = 190 # chain length

    D_chain = D_protection[protection] # chain length

    sd_chain = sd_chain[protection]  # standard deviation of the chain length

    a = rnorm(3, 0, 1) # random normal values
    
    a = D_chain + a * sd_chain # calculate chain length

    b = D_max - a # calculate chain length
    
    r = c(0, 0, 0)
    
    for (i in 1:3) {

        r[i] = rbeta(1, a[i], b[i]) # generate random values from beta distribution
    }

    return(1 - r) # return the random values
}

# simulate resource availability

biomass <- function(protection) {

    points = 100

    a = points_mean[protection] + rnorm(1, 0, points_sd) # calculate resource availability

    b = points - a # calculate resource availability

    return(rbeta(1, a, b)) # generate random values from beta distribution
}

# simulate risk


risk <- function(ind) {

    # calculate mean risk

    ## plot - level predictors
    mu = exp(alpha_risk_protected[individual_protected[ind]] 
    + beta_rugosity * mean_rugosity[individual_plot[ind]] 
    + beta_predator[predator[individual_plot[ind]]])

    sigma = 1

    # calculate gamma parameters

    a = mu^2 / sigma^2 # shape parameter
    b = sigma^2 / mu # scale parameter

    return(rgamma(n = 1, shape = a, scale = b)) # generate random values from gamma distribution
}

# simulate behaviour duration

expit <- function(x) {
    return(1 / (1 + exp(-x))) # calculate inverse logit
}

total_time <- function(behaviour, ind) {

    # calculate zero inflation probability

    pi_0 <- expit(-5 + beta_pi_risk[behaviour] * individual_risk[ind] + beta_pi_biomass[behaviour] * plot_biomass[individual_plot[ind]]) # calculate zero inflation probability

    if(rbinom(1, 1, pi_0) == 1) { # if zero inflation occurs
        return(0) # return zero
    } else {
        # calculate total time for foraging, vigilance and movement
        mu = beta_risk[behaviour] * individual_risk[ind] + beta_biomass[behaviour] * plot_biomass[individual_plot[ind]] # calculate mean time
        sigma = 1 # standard deviation
        return(rlnorm(1, mu, sigma)) # generate random values from normal distribution
    }


}

# simulate number of bites

bites <- function(foraging, ind){

    if(foraging > 0){

        mu = beta_risk_bites * individual_risk[ind] + beta_biomass_bites * plot_biomass[individual_plot[ind]]  # calculate mean bites

        lambda = exp(mu) # calculate lambda for Poisson distribution
        return(rpois(1, lambda)) # generate random values from Poisson distribution
    }
    else {
        return(0) # return zero if foraging is zero
    }
}
