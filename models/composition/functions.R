pacman::p_load(cmdstanr, here)

# Load the data

load_data <- function(run) {
  if (run == "validation") {
    data_dir <- "outputs/simulated_data"
  } else {
    data_dir <- "outputs/standardised_data"
  }
  fish <- read.csv(here(data_dir, "fish.csv"))
  transects <- read.csv(here(data_dir, "transects.csv"))
  species <- read.csv(here(data_dir, "species.csv"))
  guild <- read.csv(here(data_dir, "transect_guild.csv"))

  return(list(fish = fish, transects = transects, species = species, guild = guild))
}

# prepare the data for Stan

prep_stan <- function(fish, transects, species, guild) {
  # Prepare the data for Stan
  # fish: data frame with fish data
  # transects: data frame with transect data
  # species: data frame with species data
  stan_data <- list(
    N = nrow(fish), # number of fish
    T = nrow(transects), # number of transects
    S = length(unique(species$species)), # number of species
    P = length(unique(transects$protection)), # number of protection levels
    L = length(unique(transects$location)), # number of locations
    G = ncol(guild) - 1, # number of feeding guilds
    C = length(unique(fish$size)), # number of size classes
    transect = as.integer(fish$transect), # transect IDs
    species = as.integer(fish$species), # species IDs
    size = as.integer(fish$size), # size class IDs
    guild = as.matrix(guild[-1]), # feeding guild IDs
    count = fish$abundance, # abundance data
    protection = as.integer(transects$protection), # protection level IDs
    location = as.integer(transects$location) # location IDs
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
    # max_treedepth = 15
  )

  ## optimize the model
  # fit_mle <- fit$optimize(
  #  data = stan_data,
  #  seed = 123,
  #  iter = 1000
  # )
  # Save the fit result
  fit_result$save_object(file = here(output_dir, paste0(run, "_community_fit.rds")))

  # fit_mle$save_object(file = here(output_dir, paste0(run, "_fit_mle.rds")))
  return(fit_result)
}
