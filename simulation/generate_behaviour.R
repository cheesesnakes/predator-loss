pacman::p_load(dplyr, tidyr, ggplot2, patchwork, here)

# set save directories

data_dir <- "outputs/simulated_data"
plot_dir <- "figures/simulated_data"

if (!dir.exists(here(data_dir))) {
    dir.create(here(data_dir), recursive = TRUE)
}
if (!dir.exists(here(plot_dir))) {
    dir.create(here(plot_dir), recursive = TRUE)
}

# set plotting theme
theme_set(theme_bw()+ theme(text = element_text(size = 20), legend.position = "top")) # set theme for ggplot2

# load parameters and functions
source(here("simulation", "params.R")) # load parameters
source(here("simulation", "functions.R")) # load functions
source(here("simulation", "plots.R")) # load functions for behaviour simulation

# set seed for reproducibility
set.seed(123) # set seed for random number generation

# Generate plot level predictors

plot_id <- 1:M # plot IDs
protected <- rep(c(1, 2), each = M/2) # protected or unprotected

## simulate predator presence

# fewer predators outside of protected areas

predator <- sapply(protected, predator_presence) # generate predator presence for each plot

## simulate rugosity

# higher rugosity inside protected areas

plot_rugosity <- sapply(protected, rugosity) # generate rugosity for each plot
plot_rugosity <- t(plot_rugosity) # transpose the matrix to have plots as rows and rugosity as columns

mean_rugosity <- apply(plot_rugosity, 1, mean) # calculate mean rugosity for each plot

plot_rugosity <- data.frame(plot_rugosity)
colnames(plot_rugosity) <- c("sample_1", "sample_2", "sample_3") # rename columns
plot_rugosity$plot_id <- plot_id # add plot IDs to the data frame
write.csv(plot_rugosity, here(data_dir, "rugosity.csv"), row.names = F) # write rugosity to CSV file

## simulate resource availability

plot_biomass <- sapply(protected, biomass) # generate resource availability for each plot

# make data frame of predictors

predictors <- data.frame(plot_id, protected, predator, mean_rugosity, plot_biomass) # create data frame of predictors
colnames(predictors) <- c("plot_id", "protected", "predator", "mean_rugosity", "biomass") # rename columns
write.csv(predictors, here(data_dir, "predictors.csv"), row.names = F) # write predictors to CSV file

# invidiual level predictors

individual_id <- 1:N # individual IDs
individual_plot <- sample(plot_id, size = N, replace = TRUE) # assign individuals to plots
individual_protected <- rep(0, N) # initialize vector for protected status of individuals
for (i in 1:N) {
    individual_protected[i] <- predictors$protected[predictors$plot_id == individual_plot[i]] # assign protected status to individuals based on their plot
}

# size class

size_classes <- c(1, 2, 3, 4, 5, 6) # size classes

# calculate size class probabilities

prob_size <- matrix(0, nrow = length(size_classes), ncol = P) # initialize probability matrix for size classes

for (s in 1:length(size_classes)) {
    prob_size[s, 1] <- alpha_size[s] * beta_protection_size[s, 1] # calculate probability for inside PA
    prob_size[s, 2] <- alpha_size[s] * beta_protection_size[s, 2] # calculate probability for outside PA
}

# generate size class for individuals

size_individual <- rep(0, N) # sample size classes for inside PA and outside PA

for (i in 1:N) {
    size_individual[i] <- sample(size_classes, size = 1, replace = TRUE, prob = prob_size[, individual_protected[i]]) # sample size class for each individual
}

# foraging guilds

foraging_guilds <- c(1, 2, 3, 4) # foraging guilds

# calculate foraging guild probabilities

prob_foraging <- matrix(0, nrow = length(foraging_guilds), ncol = P) # initialize probability matrix for foraging guilds

for (f in 1:length(foraging_guilds)) {
    prob_foraging[f, 1] <- alpha_foraging[f] * beta_protection_foraging[f, 1] # calculate probability for inside PA
    prob_foraging[f, 2] <- alpha_foraging[f] * beta_protection_foraging[f, 2] # calculate probability for outside PA
}

# generate foraging guild for individuals

foraging_individual <- rep(0, N) # sample foraging guilds for inside PA and outside PA

for (i in 1:N) {
    foraging_individual[i] <- sample(foraging_guilds, size = 1, replace = TRUE, prob = prob_foraging[, individual_protected[i]]) # sample foraging guild for each individual
}

# create data frame of individual level predictors

individuals <- data.frame(individual_id, individual_plot, size_individual, foraging_individual) # create data frame of individual level predictors

colnames(individuals) <- c("ind_id", "plot_id", "size_class", "guild") # rename columns

write.csv(individuals, here(data_dir, "individuals.csv"), row.names = F) # write individual level predictors to CSV file

# Calculate risk per individual

# baseline risk by protection level

individual_risk <- sapply(1:N, function(x) risk(x)) # generate risk for each individual

plot_risk()

ggsave(here(plot_dir, "plots_risk.png"), width = 12, height = 9) # save plots to file

# calculate total duration for foraging, vigilance and movement

foraging <- sapply(1:N, function(x) total_time(1)) # generate foraging time for each individual
vigilance <- sapply(1:N, function(x) total_time(2)) # generate vigilance time for each individual
movement <- sapply(1:N, function(x) total_time(3)) # generate movement time for each individual

# calculate number of bites

bites <- sapply(foraging, function(x) bites(x)) # generate number of bites for each individual

# create data frame of individual level behaviours
response <- data.frame(individual_id, foraging, vigilance, movement, bites) # create data frame of individual level behaviours
colnames(response) <- c("ind_id", "foraging", "vigilance", "movement", "bites") # rename columns
write.csv(response, here(data_dir, "behaviour.csv"), row.names = F) # write individual level behaviours to CSV file

plot_beh()

ggsave(here(plot_dir, "plots_behaviour.png"), width = 24, height = 18) # save plots to file

beh_dis()

ggsave(here(plot_dir, "plots_behaviour_dis.png"), width = 24, height = 18) # save plots to file