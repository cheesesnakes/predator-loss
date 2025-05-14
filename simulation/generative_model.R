pacman::p_load(here)

# Generative model for simulation

## define species and traits
regional_pool <- 1:150 # regional pool of species
feeding_guild <- 1:4 # feeding guilds
regional_pool_fg <- sample(feeding_guild, length(regional_pool), replace = TRUE) # assign feeding guilds to species in the regional pool
size <- 1:4 # size classes


## determine abundance in the regional pool

beta_abn_fg <- c(5, 3, 3, 2) # effect of feeding guild on abundance

lambda_abn <- exp(beta_abn_fg[regional_pool_fg] + rnorm(length(regional_pool), 0, 2)) # calculate lambda abundance for each species based on its feeding guild
regional_abundance <- rpois(length(regional_pool), lambda_abn) # generate random abundances
total_regional_abundance <- sum(regional_abundance) # total abundance of all species in the regional pool
regional_relative_abundance <- regional_abundance / total_regional_abundance # relative abundance of each species in the regional pool

## define protection levels
P <- 2 # number of protection levels
beta_abn_inside_fg <- c(2, 2, 3, 2) # proportions of abundances for inside PA
beta_abn_inside_fg <- beta_abn_inside_fg / sum(beta_abn_inside_fg) # normalize proportions
beta_abn_outside <- c(3, 3, 2, 1) # proportions of abundances for outside PA
beta_abn_outside <- beta_abn_outside / sum(beta_abn_outside) # normalize proportions

# Calculate the probability of each species being in the local community
prob_local <- matrix(0, nrow = length(regional_pool), ncol = P) # initialize probability matrix for local communities
for (i in 1:length(regional_pool)) {
  prob_local[i, 1] <- regional_relative_abundance[i] * beta_abn_inside_fg[regional_pool_fg[i]] # calculate probability for inside PA
  prob_local[i, 2] <- regional_relative_abundance[i] * beta_abn_outside[regional_pool_fg[i]] # calculate probability for outside PA
}

# Generate local communities
N <- 2000 # number of individuals in the local community
community_inside_PA <- sample(regional_pool, size = N, replace = TRUE, prob = prob_local[,1]) # sample individuals from the regional pool for inside PA
community_outside_PA <- sample(regional_pool, size = N, replace = TRUE, prob = prob_local[,2]) # sample individuals from the regional pool for outside PA

length(unique(community_inside_PA)) # number of species in the local community inside PA
length(unique(community_outside_PA)) # number of species in the local community outside PA

## Individual traits
p_size <- c(4, 4, 3, 3, 2, 1) # proportions of size classes
p_size <- p_size / sum(p_size) # normalize proportions

beta_size_inside <- c(2, 2, 1, 1, 1, 1) # effect of size class on abundance inside PA
beta_size_outside <- c(2, 2, 1, 1, 0.5, 0.5) # effect of size class on abundance outside PA
beta_size_inside <- beta_size_inside / sum(beta_size_inside) # normalize proportions
beta_size_outside <- beta_size_outside / sum(beta_size_outside) # normalize proportions

prob_size <- matrix(0, nrow = length(size), ncol = P) # initialize probability matrix for size classes
for (s in 1:length(size)) {
  prob_size[s, 1] <- p_size[s] * beta_size_inside[s] # calculate probability for inside PA
  prob_size[s, 2] <- p_size[s] * beta_size_outside[s] # calculate probability for outside PA
}

# Generate size classes for local communities
size_inside_PA <- sample(size, size = N, replace = TRUE, prob = prob_size[,1]) # sample size classes for inside PA
size_outside_PA <- sample(size, size = N, replace = TRUE, prob = prob_size[,2]) # sample size classes for outside PA

# Create data frames for local communities

fish_inside_PA <- data.frame(species = paste("Fish", 1:N, sep = "_"), size = size_inside_PA)
fish_outside_PA <- data.frame(species = paste("Fish", 1:N, sep = "_"), size = size_outside_PA)

# Add feeding guilds to the data frames
fish_inside_PA$feeding_guild <- regional_pool_fg[community_inside_PA]
fish_outside_PA$feeding_guild <- regional_pool_fg[community_outside_PA]

fish_inside_PA$protection <- rep("2", N)
fish_outside_PA$protection <- rep("1", N)
fish <- rbind(fish_inside_PA, fish_outside_PA)

head(fish)
