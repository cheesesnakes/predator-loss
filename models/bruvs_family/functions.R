pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  if (run == "validation") {
    data_dir <- "outputs/simulated_data"
  } else {
    data_dir <- "outputs/standardised_data"
  }
  fish <- read.csv(here(data_dir, "family_counts_bruvs.csv"))
  predictors <- read.csv(here(data_dir, "predictors_bruvs.csv"))

  return(list(fish = fish, predictors = predictors))
}

# prepare the data for Stan

prep_stan <- function(fish, predictors) {
  # Prepare the data for Stan
  # fish: data frame with fish data
  # predictors: data frame with transect data
  # family: data frame with feeding family data
  stan_data <- list(
    N = nrow(fish), # number of fish
    T = nrow(predictors), # number of deployments
    P = length(unique(predictors$protection)), # number of protection levels
    L = length(unique(predictors$site)), # number of sites
    G = length(unique(fish$family)), # number of fish families
    C = length(unique(fish$size)), # number of size classes
    transect = as.integer(fish$deployment_id), # transect IDs
    size = as.integer(fish$size), # size class IDs
    family = as.integer(fish$family), # feeding family IDs
    count = fish$abundance, # abundance data
    protection = as.integer(predictors$protection), # protection level IDs
    location = as.integer(predictors$site) # location IDs
  )

  return(stan_data)
}

# function to fit the model
fit_model <- function(stan_data, model_file, run = "validation") {
  # Set output directory
  output_dir <- here("outputs")

  # Compile model
  fit <- tryCatch(
    {
      cmdstanr::cmdstan_model(model_file)
    },
    error = function(e) {
      stop("Error compiling model: ", e$message)
    }
  )

  # Fit model
  fit_result <- fit$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 50,
    # adapt_delta = 0.95,
    # max_treedepth = 15,
    init = 0.01
  )

  ## optimize the model
  # fit_mle <- fit$optimize(
  #  data = stan_data,
  #  seed = 123,
  #  iter = 1000
  # )
  # Save the fit result
  fit_result$save_object(file = here(output_dir, paste0(run, "_bruvs-family_fit.rds")))

  # fit_mle$save_object(file = here(output_dir, paste0(run, "_fit_mle.rds")))
  return(fit_result)
}
