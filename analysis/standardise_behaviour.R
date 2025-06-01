pacman::p_load(here)

source(here("analysis", "cleaning_behaviour.R"))

# standardise data
species_levels <- sort(unique(fish$spp))
size_levels <- sort(unique(fish$size))
individual_levels <- sort(unique(fish$ind_id))
plot_levels <- sort(unique(predictors$plot_id))
location_levels <- sort(unique(predictors$location))
site_levels <- sort(unique(predictors$site))
protection_levels <- c("Unprotected", "Protected")
guild_levels <- sort(colnames(guild)[-1]) # Exclude 'species' column

responses <- responses %>%
  mutate(
    ind_id = as.numeric(factor(ind_id, levels = individual_levels)),
    foraging_time = as.numeric(foraging_time),
    vigilance_time = as.numeric(vigilance_time),
    move_time = as.numeric(move_time),
    bites = as.numeric(bites),
    escape = as.numeric(escape),
    aggression = as.numeric(aggression)
  )

predictors <- predictors %>%
  mutate(
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    location = as.numeric(factor(location, levels = location_levels)),
    site = as.numeric(factor(site, levels = site_levels)),
    protection = as.numeric(factor(protection, levels = protection_levels)),
    damsel = as.numeric(damsel),
    predator = as.numeric(predator + 1),
    rugosity = as.numeric(rugosity),
    resource = as.numeric(resource),
    substrate = as.numeric(substrate),
    coral = as.numeric(coral),
    sponge = as.numeric(sponge)
  )

fish <- fish %>%
  mutate(
    ind_id = as.numeric(factor(ind_id, levels = individual_levels)),
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    spp = as.numeric(factor(spp, levels = species_levels)),
    group = as.numeric(group + 1),
    size = as.numeric(factor(size, levels = size_levels))
  )

guild <- guild %>%
  ungroup() %>%
  filter(species %in% species_levels) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  complete(species) %>%
  mutate(across(-species, ~ replace_na(., 0))) %>% # Fill NAs with 0
  mutate(species = as.numeric(species)) # Convert spp to numeric index


rugosity <- rugosity %>%
  mutate(plot_id = as.numeric(factor(plot_id, levels = plot_levels))) %>%
  filter(!is.na(plot_id))

# rename columns


responses <- responses %>%
  rename(
    foraging = foraging_time,
    vigilance = vigilance_time,
    movement = move_time,
  )

fish <- fish %>%
  rename(
    species = spp,
    size_class = size
  )

predictors <- predictors %>%
  rename(
    biomass = resource,
    protected = protection
  )

head(responses)
head(predictors)
head(fish)
head(guild)
head(rugosity)

rugosity %>%
  filter(is.na(plot_id))
# save standardised data

write.csv(responses, here("outputs/standardised_data/behaviour.csv"), row.names = FALSE)
write.csv(predictors, here("outputs/standardised_data/predictors.csv"), row.names = FALSE)
write.csv(fish, here("outputs/standardised_data/individuals.csv"), row.names = FALSE)
write.csv(guild, here("outputs/standardised_data/guild.csv"), row.names = FALSE)
write.csv(rugosity, here("outputs/standardised_data/rugosity.csv"), row.names = FALSE)
