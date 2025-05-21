pacman::p_load(dplyr, tidyr, here, stringr, purrr)

# Transect data

transects <- read.csv(here("data", "transect-metadata.csv"))

transects <- transects %>%
       select(c(transect_id, location, aspect, site, status)) %>%
       rename(
              transect = transect_id,
              protection = status
       )%>%
       mutate_all(~str_trim(.))

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
       complete(transect, species, size, fill = list(count = 0))%>%
       mutate_if(~is.character(.), ~str_trim(.))


# species data

traits <- read.csv(here("data", "transect-fish-traits.csv"))

traits <- traits %>%
       select(c(Family, Genus, Species, Length..cm., Height, Feeding.Guild, Habitat, Water.Column)) %>%
       mutate(species = paste(Genus, Species)) %>%
       rename(
              length = Length..cm.,
              height = Height, guild = Feeding.Guild, habitat = Habitat, column = Water.Column
       ) %>%
       select(-c(Species)) %>%
       mutate_if(~is.character(.), ~str_trim(.))%>%
       mutate(Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus))

phylo <- traits %>%
       select(c(Family, Genus)) %>%
       distinct()

# Helper function to get mode
get_mode <- function(x) {
       x <- na.omit(x)
       if (length(x) == 0) {
              return(NA_character_)
       }
       tbl <- table(x)
       max_freq <- max(tbl)
       mode_candidates <- names(tbl)[tbl == max_freq]
       return(mode_candidates[1]) # return the first if there's a tie
}

# Trait imputation function
trait_mode <- function(family = NA, genus = NA, trait = NA) {
       # Try genus-level
       trait_vals <- traits %>%
              filter(Genus == genus) %>%
              pull(!!sym(trait))

       mode_val <- get_mode(trait_vals)

       # If genus-level fails, try family-level
       if (is.na(mode_val)) {
              trait_vals <- traits %>%
                     filter(Family == family) %>%
                     pull(!!sym(trait))
              mode_val <- get_mode(trait_vals)
       }

       return(mode_val)
}

# Main species data processing
species_df <- data.frame(species = unique(fish$species)) %>%
       left_join(traits, by = "species") %>%
       mutate(
              Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus)
       ) %>%
       left_join(phylo, by = "Genus", relationship = "many-to-many") %>%
       mutate(
              Family = ifelse(is.na(Family.x), Family.y, Family.x),
              guild = ifelse(
                     is.na(guild),
                     map2_chr(Family, Genus, ~ trait_mode(family = .x, genus = .y, trait = "guild")),
                     guild
              ),
              guild = na_if(guild, "NULL") # clean up any residual "NULL"
       )

# Save species with no guild info
species_noinfo <- species_df %>%
       filter(is.na(guild)) %>%
       select(species) %>%
       distinct()

write.csv(species_noinfo, here("outputs/standardised_data/species_noinfo.csv"), row.names = FALSE)

species_df <- species_df %>%
       filter(!is.na(guild)) %>%
       distinct(species, .keep_all = TRUE)

species_list <- unique(species_df$species)
