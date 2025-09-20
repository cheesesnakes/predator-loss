pacman::p_load(here)

source(here("analysis", "cleaning_behaviour.R"))

# Create directory

dir.create(here("outputs/standardised_data"), showWarnings = FALSE)

# standardise data
species_levels <- sort(unique(plot_abundance$spp))
size_levels <- sort(unique(plot_size$size))
individual_levels <- sort(unique(fish$ind_id))
plot_levels <- sort(unique(predictors$plot_id))
location_levels <- sort(unique(predictors$location))
site_levels <- sort(unique(predictors$site))
protection_levels <- c("Unprotected", "Protected")
guild_levels <- sort(colnames(guild)[-1]) # Exclude 'species' column
family_levels <- sort(unique(fish$family))

responses <- responses %>%
  mutate(
    ind_id = as.numeric(factor(ind_id, levels = individual_levels)),
    foraging = as.numeric(foraging),
    vigilance = as.numeric(vigilance),
    move = as.numeric(move),
    bites = as.numeric(bites),
    escape = as.numeric(escape),
    aggression = as.numeric(aggression)
  )

moments_biomass <- c(mean(predictors$biomass, na.rm = TRUE), sd(predictors$biomass, na.rm = TRUE))

predictors <- predictors %>%
  mutate(
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    location = as.numeric(factor(location, levels = location_levels)),
    site = as.numeric(factor(site, levels = site_levels)),
    protection = as.numeric(factor(protection, levels = protection_levels)) - 1,
    damsel = as.numeric(damsel),
    predator = as.numeric(predator + 1),
    biomass = as.numeric(biomass),
    substrate = as.numeric(substrate),
    coral = as.numeric(coral),
    sponge = as.numeric(sponge),
    biomass = (biomass - mean(biomass, na.rm = TRUE)) / sd(biomass, na.rm = TRUE)
  )

fish <- fish %>%
  mutate(
    ind_id = as.numeric(factor(ind_id, levels = individual_levels)),
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    spp = as.numeric(factor(spp, levels = species_levels)),
    group = as.numeric(group + 1),
    size = as.numeric(factor(size, levels = size_levels)),
    family = as.numeric(factor(family, levels = family_levels))
  )

plot_abundance <- plot_abundance %>%
  mutate(
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    spp = as.numeric(factor(spp, levels = species_levels))
  )

plot_size <- plot_size %>%
  mutate(
    plot_id = as.numeric(factor(plot_id, levels = plot_levels)),
    size = as.numeric(factor(size, levels = size_levels))
  )

guild <- guild %>%
  ungroup() %>%
  filter(species %in% species_levels) %>%
  mutate(species = as.numeric(factor(species, levels = species_levels)))

scans <- scans %>%
  ungroup() %>%
  mutate(ind_id = as.numeric(factor(ind_id, levels = individual_levels)))

# rename columns

fish <- fish %>%
  rename(
    species = spp,
    size_class = size
  )

plot_abundance <- plot_abundance %>%
  rename(
    species = spp
  )

predictors <- predictors %>%
  rename(
    biomass = biomass,
    protected = protection
  )

head(responses)
head(predictors)
head(fish)
head(guild)

# save standardised data

write.csv(responses, here("outputs/standardised_data/behaviour.csv"), row.names = FALSE)
write.csv(predictors, here("outputs/standardised_data/predictors.csv"), row.names = FALSE)
write.csv(fish, here("outputs/standardised_data/individuals.csv"), row.names = FALSE)
write.csv(guild, here("outputs/standardised_data/guild.csv"), row.names = FALSE)
write.csv(rugosity, here("outputs/standardised_data/rugosity.csv"), row.names = FALSE)

# make data frame

data <- responses %>%
  select(ind_id, foraging, vigilance, move, bites, codingduration) %>%
  left_join(fish, by = c("ind_id" = "ind_id")) %>%
  left_join(predictors, by = c("plot_id" = "plot_id")) %>%
  left_join(guild, by = c("species" = "species")) %>%
  select(-c(sponge)) %>%
  mutate(
    predator = predator - 1,
    group = group - 1,
    size_class = as.factor(size_class)
  ) %>%
  filter(guild != "Piscivore")

# Calculate proportion of time in each behavioural state

data_ts <- data %>%
  filter(codingduration > 45) %>%
  mutate(
    bite_rate = bites/foraging,
    foraging = foraging / codingduration,
    vigilance = vigilance / codingduration,
    move = move / codingduration
  )%>%
  select(-c(bites, codingduration))
