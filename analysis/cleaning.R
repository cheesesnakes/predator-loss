pacman::p_load(dplyr, tidyr, here, stringr)

# Transect data

transects <- read.csv(here("data", "transect-metadata.csv"))

transects <- transects %>%
       select(c(transect_id, location, aspect, site, status)) %>%
       rename(
              transect = transect_id,
              protection = status
       )

# Fish data

fish <- read.csv(here("data", "transect-data_fish.csv"))
small_fish <- read.csv(here("data", "transect-data_damsel.csv"))


fish <- fish %>%
       select(-c(remarks, media.code)) %>%
       rename(
              transect = transect.id,
              species = fish.spp
       ) %>%
       pivot_longer(-c(transect, species), names_to = "size", values_to = "count")

small_fish <- small_fish %>%
       select(-c(x, territory.no)) %>%
       rename(
              transect = transect.id,
              species = fish.spp
       ) %>%
       pivot_longer(-c(transect, species), names_to = "size", values_to = "count") %>%
       group_by(transect, species) %>%
       summarise(count = sum(count)) %>%
       ungroup() %>%
       mutate(size = "0-10")

fish <- fish %>%
       bind_rows(small_fish)

fish <- fish %>%
       group_by(transect, species, size) %>%
       summarise(count = sum(count, na.rm = TRUE)) %>%
       ungroup() %>%
       replace_na(list(count = 0)) %>%
       complete(transect, species, size, fill = list(count = 0))


# species data

traits <- read.csv(here("data", "transect-fish-traits.csv"))

phylo <- traits %>%
       select(c(Family, Genus)) %>%
       distinct()

traits <- traits %>%
       select(c(Genus, Species, Length..cm., Height, Feeding.Guild, Habitat, Water.Column)) %>%
       mutate(species = paste(Genus, Species)) %>%
       rename(
              length = Length..cm.,
              height = Height, guild = Feeding.Guild, habitat = Habitat, column = Water.Column
       ) %>%
       select(-c(Species))

trait_mode <- function(x) {
       ux <- table(x)
       return(as.character(ux[which.max(ux)]))
}

species_df <- data.frame(species = unique(fish$species)) %>%
       left_join(traits, by = c("species" = "species")) %>%
       mutate(Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus)) %>%
       group_by(Genus) %>%
       mutate_at(c("length", "height", "guild", "habitat", "column"), ~ ifelse(is.na(.), trait_mode(.), as.character(.)))

species_noinfo <- species_df %>%
       filter(is.na(guild)) %>%
       select(species) %>%
       distinct() %>%
       write.csv(here("outputs/standardised_data/species_noinfo.csv"), row.names = FALSE)

species_df <- species_df %>%
       filter(!is.na(guild)) %>%
       distinct(species, .keep_all = TRUE)

species_list <- unique(species_df$species)

# code factors

transect_levels <- sort(unique(transects$transect))
location_levels <- sort(unique(transects$location))
aspect_levels <- sort(unique(transects$aspect))
site_levels <- sort(unique(transects$site))
species_levels <- sort(unique(species_df$species))
size_levels <- sort(unique(fish$size))
guild_levels <- sort(unique(species_df$guild))
column_levels <- sort(unique(species_df$column))
habitat_levels <- sort(unique(species_df$habitat))
protection_levels <- c("unprotected", "protected")

transects <- transects %>%
       mutate(
              transect = as.numeric(factor(transect, levels = transect_levels)),
              protection = as.numeric(factor(protection, levels = protection_levels)),
              location = as.numeric(factor(location, levels = location_levels)),
              aspect = as.numeric(factor(aspect, levels = aspect_levels)),
              site = as.numeric(factor(site, levels = site_levels))
       )

fish_noinfo <- fish %>%
       filter(!species %in% species_list) %>%
       select(transect, species, size) %>%
       distinct() %>%
       write.csv(here("outputs/standardised_data/fish_noinfo.csv"), row.names = FALSE)

fish <- fish %>%
       filter(species %in% species_list) %>%
       mutate(
              transect = as.numeric(factor(transect, levels = transect_levels)),
              species = as.numeric(factor(species, levels = species_levels)),
              size = as.numeric(factor(size, levels = size_levels))
       )

species_df <- species_df %>%
       mutate(
              species = as.numeric(factor(species, levels = species_levels)),
              guild = as.numeric(factor(guild, levels = guild_levels)),
              habitat = as.numeric(factor(habitat, levels = habitat_levels)),
              column = as.numeric(factor(column, levels = column_levels))
       )

colnames(transects) <- str_to_lower(colnames(transects))
colnames(fish) <- str_to_lower(colnames(fish))
colnames(species_df) <- str_to_lower(colnames(species_df))

summary(transects)
summary(fish)
summary(species_df)

# Save cleaned data
write.csv(transects, here("outputs/standardised_data/transects.csv"), row.names = FALSE)
write.csv(fish, here("outputs/standardised_data/fish.csv"), row.names = FALSE)
write.csv(species_df, here("outputs/standardised_data/species.csv"), row.names = FALSE)
