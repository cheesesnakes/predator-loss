P <- 2 # number of protection levels
L <- 12 # number of sites
M <- 60 # number of plots
N <- 500 # number of individuals
R <- 150 # number of species

beta_protection_predator <- c(0.1, 0.7) # effect of protection on predator presence

D_protection <- c(100, 150) # effect of protection on rugosity
sd_chain <- c(20, 10) # standard deviation of rugosity

points_mean <- c(30, 70) # mean resource availability for protected and unprotected areas
points_sd <- 5 # standard deviation of resource availability

alpha_size <- c(10, 7, 5, 3, 3, 3) # effect of size class on abundance
alpha_size <- alpha_size / sum(alpha_size) # normalize proportions

size_classes <- c(1, 2, 3, 4, 5, 6) # size classes
size_levels <- size_classes # for counterfactuals

beta_protection_size <- matrix(0, nrow = length(size_classes), ncol = P) # initialize matrix for size class effects
beta_protection_size[, 1] <- c(0.5, 0.5, 0.5, 0.2, 0.1, 0.1) # effect of size class on abundance inside PA
beta_protection_size[, 2] <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5) # effect of size class on abundance outside PA

alpha_foraging <- c(10, 7, 5, 3) # effect of foraging guild on abundance
alpha_foraging <- alpha_foraging / sum(alpha_foraging) # normalize proportions

foraging_guilds <- c(1, 2, 3, 4) # foraging guilds
guild_levels <- foraging_guilds # for counterfactuals

beta_protection_foraging <- matrix(0, nrow = length(foraging_guilds), ncol = P) # initialize matrix for foraging guild effects
beta_protection_foraging[, 1] <- c(0.5, 0.5, 0.5, 0.1) # effect of foraging guild on abundance inside PA
beta_protection_foraging[, 2] <- c(0.5, 0.5, 0.5, 0.5) # effect of foraging guild on abundance outside PA


# outcome model parameters

## Zero inflation parameters
alpha_pi <- c(0.25, 0.25, 0.25) # zero-inflation intercepts
beta_pi_protection <- c(0.25, -0.25, 0.25) # protection effect on zero inflation
beta_pi_rugosity <- c(0.25, -0.25, -0.25) # rugosity effect on zero inflation
beta_pi_biomass <- c(-0.25, 0, 0.25) # biomass effect on zero inflation
beta_pi_predator <- c(0.25, -0.25, 0.25) # predator effect on zero inflation
beta_size_pi <- rep(0.1, length(size_classes)) # size class effect on zero inflation
beta_size_pi <- matrix(rep(beta_size_pi, each = 3), nrow = length(size_classes), ncol = 6) # make matrix for size class effect on zero inflation
beta_guild_pi <- rep(0.1, length(foraging_guilds)) # foraging guild effect on zero inflation
beta_guild_pi <- matrix(rep(beta_guild_pi, each = 3), nrow = length(foraging_guilds), ncol = 6) # make matrix for foraging guild effect on zero inflation
beta_group_pi <- c(0, 0, -0.5) # group size effect on zero inflation

# Duration model parameters

alpha_duration <- c(3, 2, 1.5) # intercepts for foraging, vigilance and movement
beta_protection <- c(-0.5, 0.5, -0.5) # protection effect on foraging, vigilance and movement
beta_rugosity <- c(0.15, -0.15, -0.1) # effect of rugosity behaviour (foraging, movement, vigilance)
beta_biomass <- c(0.1, 0, 0.1) # effect of biomass behaviour
beta_predator <- c(-1, 0.5, 0.5) # effect of predator presence behaviour
beta_size <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1) # effect of size class on foraging, vigilance and movement
beta_size <- matrix(rep(beta_size, each = 3), nrow = length(size_classes), ncol = 6) # make matrix for size class effect on foraging, vigilance and movement
beta_guild <- c(0.1, 0.1, 0.1, 0.1) # effect of foraging guild on foraging, vigilance and movement
beta_guild <- matrix(rep(beta_guild, each = 3), nrow = length(foraging_guilds), ncol = 6) # make matrix for foraging guild effect on foraging, vigilance and movement
beta_group <- c(0, 0, 0.1) # effect of group size on foraging, vigilance and movement

# bite rate parameters

alpha_bites <- 1 # intercept for bites
beta_protection_bites <- 0.5 # effect of protection on bites
beta_rugosity_bites <- 0.1 # effect of rugosity on bites
beta_biomass_bites <- 0 # effect of biomass on bites
beta_predator_bites <- 0.1 # effect of predator presence on bites
beta_size_bites <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1) # effect of size class on bites
beta_guild_bites <- c(0.1, 0.1, 0.1, 0.1) # effect of foraging guild on bites
beta_group_bites <- 1 # effect of group size on bites