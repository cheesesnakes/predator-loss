pacman::p_load(
       here
)

source(here("analysis", "cleaning_bruvs.R"))

# code factors

deployment_levels <- sort(unique(predictors$deployment_id))
location_levels <- sort(unique(predictors$location))
protection_levels <- c("Unprotected", "Protected")
site_levels <- sort(unique(predictors$site))
species_levels <- sort(unique(guild$fish_spp))


size_levels <- sort(unique(individuals$size))
starts <- as.numeric(sub("-.*", "", size_levels))
size_levels <- size_levels[order(starts)]

# make colnames lower case
colnames(predictors) <- tolower(colnames(predictors))
colnames(guild) <- tolower(colnames(guild))
colnames(individuals) <- tolower(colnames(individuals))

# standardise predictors
predictors <- predictors %>%
    mutate(
        deployment_id = as.numeric(factor(deployment_id, levels = deployment_levels)),
        location = as.numeric(factor(location, levels = location_levels)),
        protection = as.numeric(factor(protection, levels = protection_levels)),
        site = as.numeric(factor(site, levels = site_levels))
    )

# standardise guild
guild <- guild %>%
    mutate(
        fish_spp = as.numeric(factor(fish_spp, levels = species_levels))
    )

# standardise individuals
individuals <- individuals %>%
    filter(fish_spp %in% species_levels) %>%
    mutate(
        deployment_id = as.numeric(factor(deployment_id, levels = deployment_levels)),
        fish_spp = as.numeric(factor(fish_spp, levels = species_levels)),
        size = as.numeric(factor(size, levels = size_levels))
    )%>%
    rename(
        abundance = count,
    )

# summary

summary(predictors)
summary(guild)
summary(individuals)


# save data
write.csv(predictors, here("outputs/standardised_data/predictors_bruvs.csv"), row.names = FALSE)
write.csv(guild, here("outputs/standardised_data/traits_bruvs.csv"), row.names = FALSE)
write.csv(individuals, here("outputs/standardised_data/individuals_bruvs.csv"), row.names = FALSE)