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

# traits

guild_n <- guild %>%
  ungroup() %>%
  select(-c(fish_spp)) %>%
  summarise_all(~ sum(. == 1, na.rm = TRUE)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Trait",
    values_to = "Count"
  )

guild_n <- guild_n %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Count"),
    digits = 0, na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2, align = "right", part = "body") %>%
  bold(part = "header") %>%
  autofit() %>%
  theme_box()

print(guild_n)

save_as_docx(guild_n, path = here("outputs", "bruvs_guild_n.docx"))
save_as_image(
  guild_n,
  path = here("figures", "bruvs_guild_n.png"),
  type = "png",
  width = 12,
  height = 8
)

# Calculate guild wise abundance


guild_long <- guild %>%
  pivot_longer(
    cols = -fish_spp,
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  filter(Value == 1) %>%
  select(-Value)


individuals %>%
  left_join(guild_long) %>%
  group_by(Trait) %>%
  summarise(Richness = length(unique(fish_spp)), Abundance = sum(count), .groups = "drop") %>%
  mutate(Trait = ifelse(is.na(Trait), "Unknown", Trait)) %>%
  print(.)
