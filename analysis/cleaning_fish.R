pacman::p_load(dplyr, tidyr, here, stringr, purrr, flextable)

# Benthic cover data

transects <- read.csv(here("data", "transect-metadata.csv"))
benthic <- read.csv("data/chapter-3_benthic_cover-transects.csv")

transects <- transects %>%
       select(transect_id, n_cover)

category_names <- read.csv(here("data/benthic_categories.csv"))

benthic <- benthic %>%
       left_join(category_names, by = "category") %>%
       select(transect_id, quadrat, names, count) %>%
       group_by(transect_id, quadrat, names) %>%
       summarise(cover = sum(count)) %>%
       left_join(transects, by = "transect_id") %>%
       mutate(p = cover / n_cover) %>%
       select(-c(cover, n_cover)) %>%
       pivot_wider(
              names_from = names,
              values_from = p,
              values_fill = 0
       )

# Clean transect metadata
transects <- read.csv(here("data", "transect-metadata.csv"))

transects <- transects %>%
       select(c(transect_id, location, aspect, site, status)) %>%
       rename(
              transect = transect_id,
              protection = status
       ) %>%
       mutate_all(~ str_trim(.))

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
       mutate_if(~ is.character(.), ~ str_trim(.))


fish <- fish %>%
       mutate(
              size = str_replace(size, "X", ""),
              size = str_replace(size, "cm", ""),
              size = str_replace(size, "\\.", "-")
       )

# species data

traits <- read.csv(here("data", "fish-traits.csv"))
head(count)
phylo <- traits %>%
       select(c(Order, Family, Genus)) %>%
       filter(Family != "") %>%
       distinct()

traits <- traits %>%
       select(c(Family, Genus, Species, Length..cm., Height, Feeding.Guild, Habitat, Water.Column)) %>%
       mutate(species = paste(Genus, Species)) %>%
       rename(
              length = Length..cm.,
              height = Height, guild = Feeding.Guild, habitat = Habitat, column = Water.Column
       ) %>%
       select(-c(Species)) %>%
       mutate_if(~ is.character(.), ~ str_trim(.)) %>%
       mutate(Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus)) %>%
       mutate(guild = str_to_sentence(guild))


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
                     is.na(guild) | guild == "Unknown",
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

dir.create(here("outputs/errors"), showWarnings = FALSE)
write.csv(species_noinfo, here("outputs/errors/species_noinfo.csv"), row.names = FALSE)

species_df <- species_df %>%
       filter(!is.na(guild)) %>%
       distinct(species, .keep_all = TRUE) %>%
       mutate(guild = str_to_sentence(guild))

species_list <- unique(species_df$species)

fish <- fish %>%
       filter(species %in% species_list)

# Guild membership

guild <- species_df %>%
       select(species, guild) %>%
       mutate(
              FG1 = str_split(guild, ", ", simplify = T)[, 1],
              FG2 = str_split(guild, ", ", simplify = T)[, 2]
       ) %>%
       select(-c(guild)) %>%
       pivot_longer(cols = c(FG1, FG2), names_to = "type", values_to = "guild") %>%
       select(-c(type)) %>%
       mutate(guild = str_to_sentence(str_trim(guild))) %>%
       distinct() %>%
       filter(guild != "") %>%
       mutate(value = 1) %>%
       complete(species, guild, fill = list(value = 0)) %>%
       group_by(species) %>%
       reframe(
              value = value / sum(value),
              guild = guild
       ) %>%
       pivot_wider(names_from = guild, values_from = value)

# Multiguild

multiguild <- guild[which(rowSums(guild[, -1] > 0) > 1), ] %>%
       select(species) %>%
       filter(species %in% fish$species) %>%
       rename(Species = species)

print(paste0("Number of species with multiple guilds: ", nrow(multiguild)))

multiguild <- multiguild %>%
       flextable() %>%
       autofit()

save_as_image(
       multiguild,
       path = here("figures", "transects_multiguild.png"),
       type = "png",
       width = 12,
       height = 8
)

# Use only dominant guild

guild <- guild %>%
       mutate(
              Herbivore = ifelse(Herbivore > 0, 1, 0),
              Piscivore = ifelse(Piscivore > 0, 1, 0),
              Invertivore = ifelse(Invertivore == 1, 1, 0),
       ) %>%
       select(-c(Planktivore))

# make counts data

guild <- guild %>%
       pivot_longer(-species, names_to = "guild", values_to = "value") %>%
       filter(value > 0) %>%
       select(-c(value))

counts <- fish %>%
       left_join(guild, by = "species") %>%
       group_by(transect, guild, size) %>%
       summarise(count = sum(count, na.rm = TRUE)) %>%
       ungroup()

head(fish)
head(species_df)
head(transects)
head(benthic)
head(guild)
head(counts)
