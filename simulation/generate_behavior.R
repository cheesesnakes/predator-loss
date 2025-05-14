pacman::p_load(dplyr, tidyr, ggplot2, patchwork)

theme_set(theme_bw()+ theme(text = element_text(size = 20), legend.position = "top")) # set theme for ggplot2

P <- 2 # number of protection levels
M <- 60 # number of plots
N <- 500 # number of individuals

# Generate plot level predictors

plot_id <- 1:M # plot IDs
protected <- rep(c(1, 2), each = M/2) # protected or unprotected

## simulate predator presence

beta_protection_predator <- c(0.1, 0.5) # effect of protection on predator presence

# fewer predators outside of protected areas

predator_presence <- function(protection) {
    return(rbinom(1, 1, beta_protection_predator[protection])+1)

}

predator <- sapply(protected, predator_presence) # generate predator presence for each plot

## simulate rugosity

D_protection <- c(100, 150) # effect of protection on rugosity

# higher rugosity inside protected areas

rugosity <- function(protection) {
    D_max = 190 # chain length

    D_chain = D_protection[protection] # chain length

    sd_chain = 10  # standard deviation of the chain length

    a = rnorm(3, 0, 1) # random normal values
    
    a = D_chain + a * sd_chain # calculate chain length

    b = D_max - a # calculate chain length
    
    r = c(0, 0, 0)
    for (i in 1:3) {

        r[i] = rbeta(1, a[i], b[i]) # generate random values from beta distribution
    }

    return(r) # return the random values
}
plot_rugosity <- sapply(protected, rugosity) # generate rugosity for each plot
plot_rugosity <- t(plot_rugosity) # transpose the matrix to have plots as rows and rugosity as columns

mean_rugosity <- apply(plot_rugosity, 1, mean) # calculate mean rugosity for each plot

plot_rugosity <- data.frame(plot_rugosity)
colnames(plot_rugosity) <- c("sample_1", "sample_2", "sample_3") # rename columns
plot_rugosity$plot_id <- plot_id # add plot IDs to the data frame
write.csv(plot_rugosity, "simulation/generated_data_rugosity.csv", row.names = F) # write rugosity to CSV file

## simulate resource availability

points_mean <- c(50, 70) # mean resource availability for protected and unprotected areas
points_sd <- 5

biomass <- function(protection) {

    points = 100

    a = points_mean[protection] + rnorm(1, 0, points_sd) # calculate resource availability

    b = points - a # calculate resource availability

    return(rbeta(1, a, b)) # generate random values from beta distribution
}

plot_biomass <- sapply(protected, biomass) # generate resource availability for each plot

# make data frame of predictors

predictors <- data.frame(plot_id, protected, predator, mean_rugosity, plot_biomass) # create data frame of predictors
colnames(predictors) <- c("plot_id", "protected", "predator", "mean_rugosity", "biomass") # rename columns
write.csv(predictors, "simulation/generated_data_predictors.csv", row.names = F) # write predictors to CSV file

# invidiual level predictors

individual_id <- 1:N # individual IDs
individual_plot <- sample(plot_id, size = N, replace = TRUE) # assign individuals to plots
individual_protected <- rep(0, N) # initialize vector for protected status of individuals
for (i in 1:N) {
    individual_protected[i] <- predictors$protected[predictors$plot_id == individual_plot[i]] # assign protected status to individuals based on their plot
}
size_classes <- c(1, 2, 3, 4, 5, 6) # size classes
alpha_size <- c(10, 7, 5, 3, 3, 3) # effect of size class on abundance
alpha_size <- alpha_size / sum(alpha_size) # normalize proportions
beta_protection_size <- matrix(0, nrow = length(size_classes), ncol = P) # initialize matrix for size class effects
beta_protection_size[, 1] <- c(0.5, 0.5, 0.5, 0.2, 0.1, 0.1) # effect of size class on abundance inside PA
beta_protection_size[, 2] <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5) # effect of size class on abundance outside PA

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
foraging_guilds <- c(1, 2, 3, 4) # foraging guilds
alpha_foraging <- c(10, 7, 5, 3) # effect of foraging guild on abundance
alpha_foraging <- alpha_foraging / sum(alpha_foraging) # normalize proportions
beta_protection_foraging <- matrix(0, nrow = length(foraging_guilds), ncol = P) # initialize matrix for foraging guild effects
beta_protection_foraging[, 1] <- c(0.5, 0.5, 0.5, 0.1) # effect of foraging guild on abundance inside PA
beta_protection_foraging[, 2] <- c(0.5, 0.5, 0.5, 0.5) # effect of foraging guild on abundance outside PA

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
write.csv(individuals, "simulation/generated_data_individuals.csv", row.names = F) # write individual level predictors to CSV file

# Calculate risk per individual
# baseline risk by protection level

alpha_risk_protected <- c(0.5, 0.5) # effect of protection on risk

# effect of rugosity on risk
beta_rugosity <- -0.1

# effect of biomass on risk
beta_biomass <- -0.1

# effect of predator presence on risk
beta_predator <- c(0, 0.5)

# effect of size class on risk

beta_size <- c(0.5, 0.5, 0.5, 0.2, 0.1, 0.1)

# effect of foraging guild on risk

beta_foraging <- c(0.5, 0.5, 0.5, 0.1)

risk <- function(ind) {

    # individual random effect

    i = rnorm(1, 0, 1) # individual random effect

    # calculate mean risk

    ## plot - level predictors
    mu = alpha_risk_protected[individual_protected[ind]] + beta_rugosity * mean_rugosity[individual_plot[ind]] + beta_biomass * plot_biomass[individual_plot[ind]] + beta_predator[predator[individual_plot[ind]]]
    
    ## individual - level predictors

    mu = mu + beta_size[size_individual[ind]] + beta_foraging[foraging_individual[ind]] + i # calculate mean risk

    sigma = 1

    return(rnorm(1, mu, sigma)) # generate random values from normal distribution
}

individual_risk <- sapply(1:N, function(x) risk(x)) # generate risk for each individual

# box plot risk by protection level
p1 <- ggplot(individuals, aes(x = as.factor(individual_protected), y = individual_risk)) + # create box plot of risk by protection level
    geom_boxplot() +
    labs(x = "Protection Level", y = "Risk") 

# box plot risk by size class
p2 <- ggplot(individuals, aes(x = as.factor(size_individual), y = individual_risk)) + # create box plot of risk by size class
    geom_boxplot() +
    labs(x = "Size Class", y = "Risk") 

# box plot risk by foraging guild
p3 <- ggplot(individuals, aes(x = as.factor(foraging_individual), y = individual_risk)) + # create box plot of risk by foraging guild
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Risk") 

# historgram of risk

p4 <- ggplot(individuals, aes(x = individual_risk)) + # create histogram of risk
    geom_density() +
    labs(x = "Risk", y = "Count") 

(p1 + p2 + p3 )/ p4# combine plots

ggsave("simulation/plots_risk.png", width = 12, height = 9) # save plots to file

# calculate total duration for foraging, vigilance and movement

beta_risk <- c(-0.5, 0.5, -0.5) # effect of risk on foraging, vigilance and movement
beta_pi_risk <- c(0.5, -0.5, 0.5) # zero inflation effect

expit <- function(x) {
    return(1 / (1 + exp(-x))) # calculate inverse logit
}

total_time <- function(behaviour) {

    # calculate zero inflation probability

    pi_0 <- expit(-5 + beta_pi_risk[behaviour] * individual_risk) # calculate zero inflation probability

    if(rbinom(1, 1, pi_0) == 1) { # if zero inflation occurs
        return(0) # return zero
    } else {
        # calculate total time for foraging, vigilance and movement
        mu = beta_risk[behaviour] * individual_risk # calculate mean time
        sigma = 1 # standard deviation
        return(rlnorm(1, mu, sigma)) # generate random values from normal distribution
    }


}

foraging <- sapply(1:N, function(x) total_time(1)) # generate foraging time for each individual
vigilance <- sapply(1:N, function(x) total_time(2)) # generate vigilance time for each individual
movement <- sapply(1:N, function(x) total_time(3)) # generate movement time for each individual
# calculate number of bites

beta_risk_bites <- 0.5 # effect of risk on bites
bites <- function(foraging){

    if(foraging > 0){

        mu = beta_risk_bites * individual_risk # calculate mean bites

        lambda = exp(mu) # calculate lambda for Poisson distribution
        return(rpois(1, lambda)) # generate random values from Poisson distributionq
    }
    else {
        return(0) # return zero if foraging is zero
    }
}

bites <- sapply(foraging, function(x) bites(x)) # generate number of bites for each individual

# box plot behaviours by protection level
p1 <- ggplot(response, aes(x = as.factor(individual_protected), y = foraging)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Foraging")

p2 <- ggplot(response, aes(x = as.factor(individual_protected), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Vigilance")

p3 <- ggplot(response, aes(x = as.factor(individual_protected), y = movement)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Movement")

p4 <- ggplot(response, aes(x = as.factor(individual_protected), y = bites)) +
    geom_boxplot() +
    labs(x = "Protection Level", y = "Bites")

# box plot behaviours by size class
p5 <- ggplot(response, aes(x = as.factor(size_individual), y = foraging)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Foraging")

p6 <- ggplot(response, aes(x = as.factor(size_individual), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Vigilance")

p7 <- ggplot(response, aes(x = as.factor(size_individual), y = movement)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Movement")

p8 <- ggplot(response, aes(x = as.factor(size_individual), y = bites)) +
    geom_boxplot() +
    labs(x = "Size Class", y = "Bites")

# box plot behaviours by foraging guild
p9 <- ggplot(response, aes(x = as.factor(foraging_individual), y = foraging)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Foraging")

p10 <- ggplot(response, aes(x = as.factor(foraging_individual), y = vigilance)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Vigilance")

p11 <- ggplot(response, aes(x = as.factor(foraging_individual), y = movement)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Movement")

p12 <- ggplot(response, aes(x = as.factor(foraging_individual), y = bites)) +
    geom_boxplot() +
    labs(x = "Foraging Guild", y = "Bites")

# Combine plots grouped by predictor
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8  + p9 + p10 + p11 + p12 + plot_layout(ncol = 4, nrow = 3, byrow = TRUE) 
ggsave("simulation/plots_behaviour.png", width = 24, height = 18) # save plots to file

# create data frame of individual level behaviours
response <- data.frame(individual_id, foraging, vigilance, movement, bites) # create data frame of individual level behaviours
colnames(response) <- c("ind_id", "foraging", "vigilance", "movement", "bites") # rename columns
write.csv(response, "simulation/generated_data_behaviour.csv", row.names = F) # write individual level behaviours to CSV file