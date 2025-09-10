pacman::p_load(here, flextable)

source(here("analysis", "cleaning_bruvs.R"))

colnames(fish) <- str_to_sentence(colnames(fish))

# sampling

print(paste0("Number of locations: ", n_distinct(predictors$location)))

print(paste0("Locations", ": ", paste(unique(predictors$location), collapse = ", ")))
print(paste0("Number of sites: ", n_distinct(predictors$site)))
print(paste0("Sites", ": ", paste(str_to_sentence(unique(predictors$site)), collapse = ", ")))

print(paste0("Number of deployments: ", n_distinct(predictors$deployment_id)))

# treatment

print(paste0("Protection levels: ", paste(unique(predictors$protection), collapse = ", ")))

# species

print(paste0("Number of species: ", n_distinct(individuals$fish_spp)))
print(paste0("Number of individual sampled: ", sum(individuals$count)))


print(paste0("Size classes: ", paste(unique(individuals$size), collapse = ", ")))

assemblage <- individuals %>%
  group_by(fish_spp) %>%
  summarise(Abundance = sum(count)) %>%
  mutate(Species = str_to_sentence(fish_spp)) %>%
  arrange(desc(Abundance)) %>%
  distinct() %>%
  select(Species, Abundance)

colnames(assemblage) <- c("Species", "Abundance")

assemblage <- flextable(assemblage, col_keys = colnames(assemblage)) %>%
  italic(j = 1) %>%
  set_table_properties(layout = "fixed", align = "center", width = 0.95) %>%
  theme_box()

print(assemblage)

save_as_docx(assemblage, path = here("outputs/bruvs_summary.docx"))

# Calculate guild wise abundance

individuals %>%
  left_join(guild) %>%
  group_by(guild) %>%
  summarise(Richness = length(unique(fish_spp)), Abundance = sum(count), .groups = "drop") %>%
  mutate(guild = ifelse(is.na(guild), "Unknown", guild)) %>%
  print(.)
