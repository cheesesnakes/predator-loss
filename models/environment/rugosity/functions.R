pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  data_dir <- "outputs/standardised_data"
  predictors <- read.csv(here(data_dir, "predictors.csv"))
  rugosity <- read.csv(here(data_dir, "rugosity.csv"))
  return(list(predictors = predictors, rugosity = rugosity))
}

# prepare the data for Stan

prep_stan <- function(predictors, rugosity) {
  # Prepare the data for Stan
  # response: data frame with response variables
  # individuals: data frame with individual IDs
  # predictors: data frame with predictor variables
  # rugosity: data frame with rugosity values
  stan_data <- list(
    P = length(unique(predictors$protected)),
    S = nrow(predictors),
    L = length(unique(predictors$site)),
    location = predictors$site,
    plot = predictors$plot_id,
    rugosity_raw = matrix(c(rugosity$sample_1, rugosity$sample_2, rugosity$sample_3), ncol = 3),
    protection = predictors$protected,
    R = 25
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
  fit_result$save_object(file = here(output_dir, "standardised_rugosity_fit.rds"))

  return(fit_result)
}

# function to extract posterior samples

extract <- function(p) {
# Load the fitted model
fit <- readRDS(here("outputs/standardised_rugosity_fit.rds"))

# Extract draws for each response
data = data.table()

for (r in responses) {
    param <- paste0(p,"_", r)
    message("Extracting ", param)

    draw <- fit$draws(param)
    df <- as.data.table(draw)

    # Create draw index
    df[, draw := as.numeric(iteration) * as.numeric(chain)]
    df[, c("iteration", "chain") := NULL]

    # Clean variable names
    df[, variable := str_remove(variable, paste0(param, "\\["))]
    df[, variable := str_remove(variable, "\\]")]

    if (r == "pi"){
        df[, model := "Zero-inflated"]
    } else {
        df[, model := "Behaviour"]
    }

    data <- rbind(data, df, use.names = TRUE, fill = TRUE)
    rm(df)
    gc()
}

rm(fit)
gc()

return(data)

}
