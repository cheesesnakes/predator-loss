pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  data_dir <- "outputs/standardised_data"
  transects <- read.csv(here(data_dir, "transects.csv"))
  benthic <- read.csv(here(data_dir, "transect_benthic.csv"))
  return(list(benthic = benthic, transects = transects))
}

# prepare the data for Stan

prep_stan <- function(benthic, transects) {
  # Prepare the data for Stan
  # response: data frame with response variables
  # individuals: data frame with individual IDs
  # benthic: data frame with predictor variables
  stan_data <- list(
    P = length(unique(transects$protection)),
    S = nrow(transects),
    L = length(unique(transects$site)),
    Q = nrow(benthic),
    location = transects$site,
    protection = transects$protection,
    coral = benthic$coral,
    transect_quadrat = benthic$quadrat
  )

  return(stan_data)
}

# function to fit the model
fit_model <- function(stan_data, model_file) {
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
    # adapt_delta = 0.95,
    # max_treedepth = 15
  )

  ## optimize the model
  # fit_mle <- fit$optimize(
  #  data = stan_data,
  #  seed = 123,
  #  iter = 1000
  # )
  # Save the fit result
  fit_result$save_object(file = here(output_dir, "standardised_coral_fit.rds"))

  return(fit_result)
}
