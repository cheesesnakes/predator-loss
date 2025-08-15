pacman::p_load(here, flextable)

source(here("analysis", "cleaning_fish.R"))

colnames(fish) <- str_to_sentence(colnames(fish))

# sampling

print(paste0("Number of locations: ", n_distinct(transects$location)))

print(paste0("Locations", ": ", paste(unique(transects$location), collapse = ", ")))
print(paste0("Number of sites: ", n_distinct(transects$site)))
print(paste0("Sites", ": ", paste(str_to_sentence(unique(transects$site)), collapse = ", ")))

print(paste0("Number of transects: ", n_distinct(transects$transect)))

# treatment

print(paste0("Protection levels: ", paste(unique(transects$protection), collapse = ", ")))

# species

print(paste0("Number of species: ", n_distinct(fish$Species)))
print(paste0("Number of orders: ", n_distinct(species_df$Order)))
print(paste0("Number of families: ", n_distinct(species_df$Family)))
print(paste0("Number of individual sampled: ", sum(fish$Count)))

print("Feeding guilds: ")
table(str_to_sentence(unlist(str_split(species_df$guild, ", "))))

print(paste0("Species with more than 1 feeding guild", ": ", sum(grepl(", ", species_df$guild))))
print(paste("Species with no feeding guild", ": ", sum(species_df$guild == "")))

print(paste0("Species with no information: ", nrow(species_noinfo)))

print(paste0("Size classes: ", paste(unique(fish$Size), collapse = ", ")))

assemblage <- fish %>%
  group_by(Species) %>%
  summarise(Abundance = sum(Count)) %>%
  mutate(Species = str_to_sentence(Species)) %>%
  left_join(species_df, by = c("Species" = "species")) %>%
  arrange(desc(Abundance)) %>%
  distinct() %>%
  select(Order, Family, Species, guild, column, habitat, Abundance)

colnames(assemblage) <- c("Order", "Family", "Species", "Guild", "Column", "Habitat", "Abundance")

assemblage <- flextable(assemblage, col_keys = colnames(assemblage)) %>%
  italic(j = 3) %>%
  set_table_properties(layout = "fixed", align = "center", width = 0.95) %>%
  theme_box()

print(assemblage)

save_as_docx(assemblage, path = here("outputs/assemblage_summary.docx"))

# traits

guild_n <- guild %>%
  ungroup() %>%
  select(-c(species)) %>%
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

save_as_docx(guild_n, path = here("outputs", "transects_guild_n.docx"))
save_as_image(
  guild_n,
  path = here("figures", "transects_guild_n.png"),
  type = "png",
  width = 12,
  height = 8
)
