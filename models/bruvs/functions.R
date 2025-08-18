pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  if (run == "validation") {
    data_dir <- "outputs/simulated_data"
  } else {
    data_dir <- "outputs/standardised_data"
  }
  fish <- read.csv(here(data_dir, "individuals_bruvs.csv"))
  predictors <- read.csv(here(data_dir, "predictors_bruvs.csv"))
  guild <- read.csv(here(data_dir, "traits_bruvs.csv"))

  return(list(fish = fish, predictors = predictors, guild = guild))
}

# prepare the data for Stan

prep_stan <- function(fish, predictors, guild) {
  # Prepare the data for Stan
  # fish: data frame with fish data
  # predictors: data frame with transect data
  # guild: data frame with feeding guild data
  stan_data <- list(
    N = nrow(fish), # number of fish
    T = nrow(predictors), # number of deployments
    S = length(unique(guild$fish_spp)), # number of species
    P = length(unique(predictors$protection)), # number of protection levels
    L = length(unique(predictors$site)), # number of sites
    G = ncol(guild) - 1, # number of feeding guilds
    C = length(unique(fish$size)), # number of size classes
    transect = as.integer(fish$deployment_id), # transect IDs
    species = as.integer(fish$fish_spp), # species IDs
    size = as.integer(fish$size), # size class IDs
    guild = as.matrix(guild[-1]), # feeding guild IDs
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
  fit_result$save_object(file = here(output_dir, paste0(run, "_bruvs_fit.rds")))

  # fit_mle$save_object(file = here(output_dir, paste0(run, "_fit_mle.rds")))
  return(fit_result)
}


# function to extract posterior samples

extract <- function(p) {
  # Load the fitted model
  fit <- readRDS(here("outputs/standardised_bruvs_fit.rds"))

  interactions <- c("", "_G", "_C")

  data <- data.table()

  for (i in interactions)
  {
    param <- paste0(p, i)
    # Extract draws
    message("Extracting ", param)

    draw <- fit$draws(param)
    df <- as.data.table(draw)

    # Create draw index
    df[, draw := as.numeric(iteration) * as.numeric(chain)]
    df[, c("iteration", "chain") := NULL]

    # strip variable name

    if (i != "") {
      df[, variable := str_remove_all(variable, paste0(param, "\\["))]
      df[, variable := str_remove_all(variable, "\\]")]
    } else {
      df[, variable := NA]
    }

    df[, model := "abundance"]

    # add interaction type

    if (i == "_G") {
      df[, interaction := "Foraging guild"]
    } else if (i == "_C") {
      df[, interaction := "Size class"]
    } else {
      df[, interaction := "Overall"]
    }

    # save
    data <- rbind(data, df, use.names = TRUE, fill = TRUE)
    rm(df)
    gc()
  }

  # Clean up
  rm(fit)
  gc()

  return(data)
}
