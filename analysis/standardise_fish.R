pacman::p_load(
       here
)

source(here("analysis", "cleaning_fish.R"))

# Guild membership

guild <- species_df%>%
  select(species, guild)%>%
  mutate(FG1 = str_split(guild, ", ", simplify = T)[,1],
         FG2 = str_split(guild, ", ", simplify = T)[,2])%>%
  select(-c(guild))%>%
  pivot_longer(cols = c(FG1, FG2), names_to = "type", values_to = "guild")%>%
  select(-c(type))%>%
  mutate(guild = str_to_sentence(str_trim(guild)))%>%
  distinct()%>%
  filter(guild != "")%>%
  mutate(value = 1)%>%
  complete(species, guild, fill = list(value = 0))%>%
  group_by(species)%>%
  reframe(value = value/sum(value),
          guild = guild)%>%
  pivot_wider(names_from = guild, values_from = value)
  

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
       select(species) %>%
       distinct() %>%
       write.csv(here("outputs/standardised_data/fish_noinfo.csv"), row.names = FALSE)

fish <- fish %>%
       filter(species %in% species_list) %>%
       mutate(
              transect = as.numeric(factor(transect, levels = transect_levels)),
              species = as.numeric(factor(species, levels = species_levels)),
              size = as.numeric(factor(size, levels = size_levels)),
       )

species_df <- species_df %>%
       mutate(
              species = as.numeric(factor(species, levels = species_levels)),
              guild = as.numeric(factor(guild, levels = guild_levels)),
              habitat = as.numeric(factor(habitat, levels = habitat_levels)),
              column = as.numeric(factor(column, levels = column_levels))
       )

guild <- guild %>%
       mutate(
              species = as.numeric(factor(species, levels = species_levels))
       )

colnames(transects) <- str_to_lower(colnames(transects))
colnames(fish) <- str_to_lower(colnames(fish))
colnames(species_df) <- str_to_lower(colnames(species_df))
colnames(guild) <- str_to_lower(colnames(guild))

fish <- rename(fish, abundance = count)

fish <- filter(fish, abundance > 0)

summary(transects)
summary(fish)
summary(species_df)
summary(guild)

# Save cleaned data
write.csv(transects, here("outputs/standardised_data/transects.csv"), row.names = FALSE)
write.csv(fish, here("outputs/standardised_data/fish.csv"), row.names = FALSE)
write.csv(species_df, here("outputs/standardised_data/species.csv"), row.names = FALSE)
write.csv(guild, here("outputs/standardised_data/transect_guild.csv"), row.names = FALSE)
