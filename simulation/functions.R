# simulate predator presence based on protection
predator_presence <- function(protection) {
    return(rbinom(1, 1, beta_protection_predator[protection]) + 1)
}


# simulate rugosity

rugosity <- function(protection) {
    D_max <- 190 # chain length

    D_chain <- D_protection[protection] # chain length

    sd_chain <- sd_chain[protection] # standard deviation of the chain length

    a <- rnorm(n = 3, mean = 0, sd = 1) # random normal values

    a <- D_chain + a * sd_chain # calculate chain length

    b <- D_max - a # calculate chain length

    r <- c(0, 0, 0)

    for (i in 1:3) {
        r[i] <- rbeta(1, a[i], b[i]) # generate random values from beta distribution
    }

    return(1 - r) # return the random values
}

# simulate resource availability

biomass <- function(protection) {
    points <- 100

    a <- points_mean[protection] + rnorm(1, 0, points_sd) # calculate resource availability

    b <- points - a # calculate resource availability

    return(rbeta(n = 1, shape1 = a, shape2 = b)) # generate random values from beta distribution
}

# simulate behaviour duration

expit <- function(x) {
    return(1 / (1 + exp(-x))) # calculate inverse logit
}

total_time <- function(behaviour, ind) {
    # calculate zero inflation probability

    size <- size_individual[ind] # get size class
    guild <- guild_species[individual_spp[ind]] # get foraging guild

    mu <- 0 # initial

    intercept = rnorm(1, alpha_pi[behaviour], 0.5) # plot-level random intercept, pooled at site

    pi_0 <- expit(
        intercept+
            beta_pi_protection[behaviour] * protected[individual_plot[ind]]
            + beta_pi_rugosity[behaviour] * mean_rugosity[individual_plot[ind]]
            + beta_pi_biomass[behaviour] * plot_biomass[individual_plot[ind]]
            + beta_pi_predator[behaviour] * predator[individual_plot[ind]]
            + beta_size_pi[behaviour, size]
            + beta_guild_pi[behaviour, guild]
            + beta_group_pi[behaviour] * group[ind]
    ) # calculate zero inflation probability

    if (rbinom(n = 1, size = 1, prob = pi_0) == 1) { # if zero inflation occurs
        return(0) # return zero
    } else {
        # calculate total time for foraging, vigilance and movement
        intercept = rnorm(1, alpha_duration[behaviour], 0.5) # plot-level random intercept, pooled at site
        
        mu <- (intercept + beta_protection[behaviour] * protected[individual_plot[ind]]
            + beta_rugosity[behaviour] * mean_rugosity[individual_plot[ind]]
            + beta_biomass[behaviour] * plot_biomass[individual_plot[ind]]
            + beta_predator[behaviour] * predator[individual_plot[ind]]
            + beta_size[behaviour, size]
            + beta_guild[behaviour, guild]
            + beta_group[behaviour] * group[ind])

        sigma <- 0.5 # standard deviation
        return(rlnorm(n = 1, meanlog = mu, sdlog = sigma)) # generate random values from normal distribution
    }
}

# simulate number of bites

bites <- function(foraging, ind) {
    if (foraging > 0) {
        size <- size_individual[ind] # get size class
        guild <- guild_species[individual_spp[ind]] # get foraging guild

        intercept = rnorm(1, alpha_bites, 0.5) # plot-level random intercept, pooled at site
        mu <- (intercept +
            beta_protection_bites * protected[individual_plot[ind]]
            + beta_rugosity_bites * mean_rugosity[individual_plot[ind]]
            + beta_biomass_bites * plot_biomass[individual_plot[ind]]
            + beta_predator_bites * predator[individual_plot[ind]]
            + beta_size_bites[size]
            + beta_guild_bites[guild]
            + beta_group_bites * group[ind])

        lambda <- exp(mu) # calculate lambda for Poisson distribution

        return(rpois(n = 1, lambda = lambda)) # generate random values from Poisson distribution
    } else {
        return(0) # return zero if foraging is zero
    }
}
