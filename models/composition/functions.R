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
  benthic <- read.csv(here(data_dir, "transect_benthic.csv"))

  return(list(fish = fish, transects = transects, species = species, guild = guild, benthic = benthic))
}

# prepare the data for Stan

prep_stan <- function(fish, transects, species, guild, benthic) {
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
    Q = length(benthic$quadrat), # number of quadrats
    transect = as.integer(fish$transect), # transect IDs
    species = as.integer(fish$species), # species IDs
    size = as.integer(fish$size), # size class IDs
    guild = as.matrix(guild[-1]), # feeding guild IDs
    count = fish$abundance, # abundance data
    protection = as.integer(transects$protection), # protection level IDs
    location = as.integer(transects$location), # location IDs

    transect_quadrat = benthic$transect_id, # transect index for benthic data
    coral = benthic$coral, # coral cover data
    biomass = benthic$resource # benthic biomass data
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
  fit_result$save_object(file = here(output_dir, paste0(run, "_community_fit.rds")))

  # fit_mle$save_object(file = here(output_dir, paste0(run, "_fit_mle.rds")))
  return(fit_result)
}


# function to extract posterior samples

extract <- function(p) {
# Load the fitted model
fit <- readRDS(here("outputs/standardised_community_fit.rds"))

models <- c("", "_psi")

data <- data.table()

for (m in models){

param <- paste0(p, m)
# Extract draws 
message("Extracting ", param)

draw <- fit$draws(param)
df <- as.data.table(draw)

# Create draw index
df[, draw := as.numeric(iteration) * as.numeric(chain)]
df[, c("iteration", "chain") := NULL]

# add model name
if (m == "_psi") {
  df[, model := "zero-inflation"]
} else {
  df[, model := "abundance"]
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
