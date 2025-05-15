library(here)

source(here("models", "functions.R"))

# import data

run <- "validation"

# test
data <- load_data(run)

response <- data$response
individuals <- data$individuals
predictors <- data$predictors
rugosity <- data$rugosity

colnames(response)
colnames(individuals)
colnames(predictors)
colnames(rugosity)

stan_data <- prep_stan(response, individuals, predictors, rugosity)

fit_result <- fit_model(stan_data, here("models", "model.stan"), run = run)

