pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  if (run == "validation") {
    data_dir <- "outputs/simulated_data"
  } else {
    data_dir <- "outputs/standardised_data"
  }
  response <- read.csv(here(data_dir, "behaviour.csv"))
  individuals <- read.csv(here(data_dir, "individuals.csv"))
  predictors <- read.csv(here(data_dir, "predictors.csv"))
  rugosity <- read.csv(here(data_dir, "rugosity.csv"))
  return(list(response = response, individuals = individuals, predictors = predictors, rugosity = rugosity))
}

# prepare the data for Stan

prep_stan <- function(response, individuals, predictors, rugosity) {
  # Prepare the data for Stan
  # response: data frame with response variables
  # individuals: data frame with individual IDs
  # predictors: data frame with predictor variables
  # rugosity: data frame with rugosity values
  stan_data <- list(
    N = nrow(response),
    P = length(unique(predictors$protected)),
    B = 3,
    S = nrow(predictors),
    D = matrix(c(response$foraging, response$movement, response$vigilance), ncol = 3),
    individual = individuals$ind_id,
    C = length(unique(individuals$size_class)),
    size = individuals$size_class,
    G = length(unique(individuals$guild)),
    guild = individuals$guild,
    bites = response$bites,
    predator = predictors$predator,
    plot = individuals$plot_id,
    rugosity_raw = matrix(c(rugosity$sample_1, rugosity$sample_2, rugosity$sample_3), ncol = 3),
    protection = predictors$protected,
    biomass = predictors$biomass
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
    refresh = 100,
    adapt_delta = 0.95,
    max_treedepth = 15
  )

  ## optimize the model
  #fit_mle <- fit$optimize(
  #  data = stan_data,
  #  seed = 123,
  #  iter = 1000
  #)
  # Save the fit result
  fit_result$save_object(file = here(output_dir, paste0(run, "_fit.rds")))

  #fit_mle$save_object(file = here(output_dir, paste0(run, "_fit_mle.rds")))
  return(fit_result)
}
