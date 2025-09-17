pacman::p_load(here, ggplot2, GGally, purrr, flextable)

source(here("analysis", "cleaning_behaviour.R"))

# sampling

print(paste0("Number of plots: ", length(unique(predictors$plot_id))))
print(paste0("Number of locations: ", length(unique(predictors$location))))
print(paste0("Locations sampled: ", paste(unique(predictors$location), collapse = ", ")))
print(paste0("Number of sites: ", length(unique(predictors$site))))
print(paste0("Sites samples: ", paste(unique(predictors$site), collapse = ", ")))

print(paste0("Number of individuals: ", length(unique(fish$ind_id))))
print(paste0("Number of species: ", length(unique(fish$spp))))
print(paste0("Number of size classes: ", length(unique(fish$size))))
print(paste0("Size classes sampled: ", paste(unique(fish$size), collapse = ", ")))

print(paste0("Guilds samples: ", paste(colnames(guild)[-1], collapse = ", ")))

print(paste0("Number of grouping individuals: ", sum(fish$group)))

# Treatment

print(paste0("Number of treatments: ", length(unique(predictors$protection))))
print(paste0("Treatments: ", paste(unique(predictors$protection), collapse = ", ")))

# covariates

colnames(predictors) <- str_to_sentence(colnames(predictors))
str_to_sentence(colnames(predictors)[5:ncol(predictors)])

ggpairs(predictors[, 4:ncol(predictors)],
  lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5)),
  upper = list(continuous = wrap("cor", size = 3))
) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(here("figures", "plot-predictors.png"), width = 12, height = 12)

predictors_summary <- predictors %>%
  select(-c(Plot_id, Location, Site)) %>%
  group_by(Protection) %>%
  pivot_longer(
    cols = -c(Protection),
    names_to = "Predictor",
    values_to = "Value"
  ) %>%
  group_by(Protection, Predictor) %>%
  summarise(across(everything(), list(
    n = ~ length(.x),
    median = ~ median(.x, na.rm = TRUE),
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.fn}"), .groups = "drop")

colnames(predictors_summary) <- str_to_title(colnames(predictors_summary))

predictors_summary <- predictors_summary %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Sd", "Min", "Q1", "Median", "Q3", "Max"),
    digits = 2, na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 3:9, align = "right", part = "body") %>%
  bold(part = "header") %>%
  autofit() %>%
  merge_v(j = "Protection") %>%
  valign(j = "Protection", valign = "center") %>%
  theme_box()

print(predictors_summary)

save_as_docx(predictors_summary, path = here("outputs", "predictors_summary.docx"))
save_as_image(
  predictors_summary,
  path = here("figures", "predictors_summary.png"),
  type = "png",
  width = 12,
  height = 8
)

# responses

p <- ggpairs(responses[2:ncol(responses)],
  lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5)),
  upper = list(continuous = wrap("cor", size = 3))
) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(p, filename = here("figures", "behaviour-responses.png"), width = 12, height = 12)

responses_summary <- responses %>%
  select(-c(ind_id)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Response",
    values_to = "Value"
  ) %>%
  group_by(Response) %>%
  summarise(across(everything(), list(
    n = ~ sum(.x > 0, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.fn}"), .groups = "drop")

colnames(responses_summary) <- str_to_title(colnames(responses_summary))
responses_summary$Response <- str_to_sentence(str_replace(responses_summary$Response, "_", " "))

responses_summary <- responses_summary %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Sd", "Min", "Q1", "Median", "Q3", "Max"),
    digits = 2, na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 3:9, align = "right", part = "body") %>%
  bold(part = "header") %>%
  theme_box()

print(responses_summary)

save_as_docx(responses_summary, path = here("outputs", "responses_summary.docx"))

save_as_image(
  responses_summary,
  path = here("figures", "responses_summary.png"),
  type = "png",
  width = 12,
  height = 8
)

# traits

guild_n <- guild %>%
  ungroup() %>%
  group_by(guild) %>%
  tally()

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

save_as_docx(guild_n, path = here("outputs", "guild_n.docx"))
save_as_image(
  guild_n,
  path = here("figures", "guild_n.png"),
  type = "png",
  width = 12,
  height = 8
)

# Calculate guild wise abundance

fish %>%
  left_join(guild, by = c("spp" = "species")) %>%
  rename(Species = spp) %>%
  filter(is.na(guild)) %>%
  select(Species) %>%
  distinct() %>%
  arrange(Species) %>%
  flextable() %>%
  italic(j = 1) %>%
  print(.)

fish %>%
  left_join(guild, by = c("spp" = "species")) %>%
  group_by(guild) %>%
  summarise(Richness = length(unique(spp)), Abundance = n(), .groups = "drop") %>%
  mutate(guild = ifelse(is.na(guild), "Unknown", guild)) %>%
  print(.)

source(here("analysis", "standardise_behaviour.R"))

# number of observations

data %>%
  filter(codingduration > 100) %>%
  group_by(protected) %>%
  tally()


data_ts %>%
  select(species) %>%
  mutate(species = species_levels[species]) %>%
  group_by(species) %>%
  tally() %>%
  arrange(desc(n)) %>%
  flextable() %>%
  italic(j = 1) %>%
  set_table_properties(width = .5, layout = "fixed") %>%
  colformat_double(
    j = c("n"),
    digits = 0, na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2, align = "right", part = "body") %>%
  bold(part = "header") %>%
  autofit() %>%
  theme_box() %>%
  print(.)

print(paste0("Total time sampled: ", sum(data_ts$codingduration) / 60))
