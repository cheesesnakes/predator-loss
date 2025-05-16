P <- 2 # number of protection levels
M <- 60 # number of plots
N <- 500 # number of individuals

beta_protection_predator <- c(0.1, 0.7) # effect of protection on predator presence

D_protection <- c(100, 150) # effect of protection on rugosity
sd_chain <- c(20, 10) # standard deviation of rugosity

points_mean <- c(30, 70) # mean resource availability for protected and unprotected areas
points_sd <- 5 # standard deviation of resource availability

alpha_size <- c(10, 7, 5, 3, 3, 3) # effect of size class on abundance
alpha_size <- alpha_size / sum(alpha_size) # normalize proportions

size_classes <- c(1, 2, 3, 4, 5, 6) # size classes

beta_protection_size <- matrix(0, nrow = length(size_classes), ncol = P) # initialize matrix for size class effects
beta_protection_size[, 1] <- c(0.5, 0.5, 0.5, 0.2, 0.1, 0.1) # effect of size class on abundance inside PA
beta_protection_size[, 2] <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5) # effect of size class on abundance outside PA

alpha_foraging <- c(10, 7, 5, 3) # effect of foraging guild on abundance
alpha_foraging <- alpha_foraging / sum(alpha_foraging) # normalize proportions

foraging_guilds <- c(1, 2, 3, 4) # foraging guilds

beta_protection_foraging <- matrix(0, nrow = length(foraging_guilds), ncol = P) # initialize matrix for foraging guild effects
beta_protection_foraging[, 1] <- c(0.5, 0.5, 0.5, 0.1) # effect of foraging guild on abundance inside PA
beta_protection_foraging[, 2] <- c(0.5, 0.5, 0.5, 0.5) # effect of foraging guild on abundance outside PA

alpha_risk_protected <- c(0, 1) # effect of protection on risk

# effect of rugosity on risk
beta_rugosity <- -0.5

# effect of predator presence on risk
beta_predator <- c(0, 0.5)

beta_risk <- c(-0.5, 1, -0.1) # effect of risk on foraging, vigilance and movement
beta_biomass <- c(0.5, 0, 0.01) # effect of resource availability on foraging, vigilance and movement

beta_pi_risk <- c(0.5, -0.5, 0.5) # zero inflation effect
beta_pi_biomass <- c(-0.5, 0, 0) # zero inflation effect


beta_risk_bites <- 0.1 # effect of risk on bites
beta_biomass_bites <- 0.5 # effect of resource availability on bites

