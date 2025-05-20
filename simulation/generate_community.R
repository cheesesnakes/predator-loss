pacman::p_load(dplyr, tidyr, purrr, here, stringr)

# Generative model for simulation

## define species and traits
regional_pool <- 1:150 # regional pool of species
feeding_guild <- 1:4 # feeding guilds
size <- 1:4 # size classes

regional_pool_fg <- sample(feeding_guild, length(regional_pool), replace = TRUE) # assign feeding guilds to species in the regional pool

regional_pool_df <- data.frame(species = regional_pool, feeding_guild = regional_pool_fg)

## determine abundance in the regional pool

beta_abn_fg <- c(5, 3, 3, 2) # effect of feeding guild on abundance

lambda_abn <- exp(beta_abn_fg[regional_pool_fg] + rnorm(length(regional_pool), 0, 2)) # calculate lambda abundance for each species based on its feeding guild
regional_pool_df$abundance <- rpois(length(regional_pool), lambda_abn) # generate random abundances
regional_pool_df$p <- regional_pool_df$abundance / sum(regional_pool_df$abundance) # calculate proportions of abundances

## define protection levels
P <- 2 # number of protection levels
beta_abn_inside_fg <- c(2, 2, 3, 2) # proportions of abundances for inside PA
beta_abn_inside_fg <- beta_abn_inside_fg / sum(beta_abn_inside_fg) # normalize proportions
beta_abn_outside <- c(3, 3, 2, 1) # proportions of abundances for outside PA
beta_abn_outside <- beta_abn_outside / sum(beta_abn_outside) # normalize proportions

# Calculate the probability of each species being in the local community
prob_local <- matrix(0, nrow = length(regional_pool), ncol = P) # initialize probability matrix for local communities
for (i in 1:length(regional_pool)) {
  prob_local[i, 1] <- regional_pool_df$p[i] * beta_abn_inside_fg[regional_pool_fg[i]] # calculate probability for inside PA
  prob_local[i, 2] <- regional_pool_df$p[i] * beta_abn_outside[regional_pool_fg[i]] # calculate probability for outside PA
}

# Generate local communities
Tr <- 20 # number of transsects
N_t <- sample(1:500, Tr, replace = TRUE) # number of individuals in each transsect

fish <- data.frame(Transect = 1:Tr, Protection = rep(1:P, each = Tr / P)) # initialize empty data frame for fish data

sampling <- function(Transect, protection) {
  N <- N_t[Transect] # get the number of individuals for the current transect
  p <- prob_local[, protection] # get the probabilities for the specified protection level
  return(sample(regional_pool, size = N, replace = TRUE, prob = p))
}

fish <- fish %>%
  mutate(Species = map2(Transect, Protection, ~ sampling(.x, .y))) %>%
  unnest(Species)

fish %>%
  group_by(Protection) %>%
  summarise(
    richness = n_distinct(Species),
    abundance = n()
  )

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

sample_size <- function(Protection) {
  p <- prob_size[, Protection] # get the probabilities for the specified protection level
  return(sample(size, size = 1, replace = TRUE, prob = p))
}

fish <- fish %>%
  mutate(Size = map(Protection, ~ sample_size(.x))) %>%
  unnest(Size) %>%
  select(-Protection) %>%
  group_by(Transect, Species, Size) %>%
  summarise(abundance = n()) %>%
  ungroup() %>%
  complete(Transect, Species, Size, fill = list(abundance = 0))

transects <- data.frame(Transect = 1:Tr, Location = rep(1:10, each = Tr / 10), Protection = rep(1:P, length.out = Tr)) # create a data frame for transects

species <- data.frame(Species = regional_pool, Feeding_guild = regional_pool_fg) # create a data frame for species

colnames(transect) <- str_to_lower(colnames(transect))
colnames(fish) <- str_to_lower(colnames(fish))
colnames(species) <- str_to_lower(colnames(species))

write.csv(fish, here("outputs/simulated_data/fish.csv"), row.names = F) # write fish data to CSV file
write.csv(transects, here("outputs/simulated_data/transects.csv"), row.names = F) # write transects data to CSV file
write.csv(species, here("outputs/simulated_data/species.csv"), row.names = F) # write species data to CSV file
