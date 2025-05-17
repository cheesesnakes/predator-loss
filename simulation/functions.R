# simulate predator presence based on protection
predator_presence <- function(protection) {
    return(rbinom(1, 1, beta_protection_predator[protection])+1)

}


# simulate rugosity

rugosity <- function(protection) {
    D_max = 190 # chain length

    D_chain = D_protection[protection] # chain length

    sd_chain = sd_chain[protection]  # standard deviation of the chain length

    a = rnorm(n = 3, mean = 0, sd = 1) # random normal values
    
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

    return(rbeta(n = 1, shape1 = a, shape2 = b)) # generate random values from beta distribution
}

# simulate risk


risk <- function(ind) {

    # calculate mean risk

    ## plot - level predictors
    mu = alpha_risk_protected[individual_protected[ind]] 
    + beta_rugosity * mean_rugosity[individual_plot[ind]] 
    + beta_predator[predator[individual_plot[ind]]]
    + beta_biomass * plot_biomass[individual_plot[ind]] # calculate mean risk

    sigma = 1

    return(rlnorm(n = 1, meanlog =  mu, sdlog = sigma)) # generate random values from gamma distribution
}

# simulate behaviour duration

expit <- function(x) {
    return(1 / (1 + exp(-x))) # calculate inverse logit
}

total_time <- function(behaviour, ind) {

    # calculate zero inflation probability

    size = size_individual[ind] # get size class
    guild = foraging_individual[ind] # get foraging guild

    pi_0 <- expit(-5 + beta_pi_risk[behaviour] * individual_risk[ind] 
    + beta_size_pi[behaviour, size] 
    + beta_guild_pi[behaviour, guild]
    + beta_size_risk_pi[behaviour, size] * individual_risk[ind]
    + beta_guild_risk_pi[behaviour, guild] * individual_risk[ind]) # calculate zero inflation probability

    if(rbinom(n = 1, size= 1, prob = pi_0) == 1) { # if zero inflation occurs
        return(0) # return zero
    } else {
        # calculate total time for foraging, vigilance and movement
        mu = beta_risk[behaviour] * individual_risk[ind] 
        + beta_size[behaviour, size] 
        + beta_guild[behaviour, guild] 
        + beta_size_risk[behaviour, size] * individual_risk[ind]
        + beta_guild_risk[behaviour, guild] * individual_risk[ind] # calculate mean time

        sigma = 1 # standard deviation
        return(rlnorm(n = 1, meanlog =  mu, sdlog = sigma)) # generate random values from normal distribution
    }


}

# simulate number of bites

bites <- function(foraging, ind){

    if(foraging > 0){

        size = size_individual[ind] # get size class
        guild = foraging_individual[ind] # get foraging guild

        mu = beta_risk_bites * individual_risk[ind] 
        + beta_size_bites[size] 
        + beta_guild_bites[guild]
        + beta_size_risk_bites[size] * individual_risk[ind]
        + beta_guild_risk_bites[guild] * individual_risk[ind] # calculate mean bites

        lambda = exp(mu) # calculate lambda for Poisson distribution
        
        return(rpois(n = 1, lambda = lambda)) # generate random values from Poisson distribution
    }
    else {
        return(0) # return zero if foraging is zero
    }
}
