pacman::p_load(
  here,
  dplyr,
  tidyr,
  flextable,
  patchwork,
  stringr,
  ggplot2,
  data.table,
  ggridges
)

theme_set(theme_bw() +
  theme(legend.position = "top", legend.text = element_text(size = 20)))

# set and create directories

run <- "validation"

figure_dir <- switch(run,
  validation = here("figures", "simulated_data"),
  standardised = here("figures")
)

output_dir <- switch(run,
  validation = here("outputs", "validation"),
  standardised = here("outputs")
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load raw data

if (run == "standardised"){
source(here("analysis", "standardise_behaviour.R"))
}else{
  source(here("simulation", "generate_behaviour.R"))
  fish <- individuals
  predictors <- predictors %>%
  rename(rugosity = mean_rugosity)
  
}


traits <- guild %>%
  pivot_longer(
    cols = -c(species),
    names_to = "guild",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  group_by(species) %>%
  summarise(guild = first(guild))


fish <- fish %>%
  left_join(traits, by = "species") %>%
  select(-c(group, species))

predictors <- predictors %>%
  select(c(plot_id, biomass, rugosity, predator))

rm(rugosity, responses, traits)
gc()

## counterfactual protection

fit <- readRDS(here("outputs", paste0(run, "_behaviour_fit.rds")))

protection <- fit$draws("D_pred_protection", format = "matrix")

protection_df <- as.data.table(protection) %>%
  mutate(Draw = seq_len(nrow(protection))) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection, fit)
gc()

setDT(protection_df)
protection_df[, variable := str_remove(variable, "D_pred_protection\\[")]
protection_df[, variable := str_remove(variable, "\\]")]
protection_df[, c("protection", "individual", "response") := tstrsplit(variable, ",")]
protection_df[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside", "Inside")
)]
protection_df[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
protection_df[, variable := NULL]

dcast(protection_df,
  Draw + individual + response ~ protection,
  value.var = "value"
) -> protection_df
protection_df[, Delta := log2(Inside + 1) - log2(Outside + 1)]
protection_df[, P := Inside > Outside] # P(Inside > Outside)

posterior_protection <- protection_df[, .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = .(Draw, response)]

# counterfactual protection bites
fit <- readRDS(here("outputs", paste0(run, "_behaviour_fit.rds")))
protection_bites <- fit$draws("bites_pred_protection", format = "matrix")
protection_bites <- as.data.table(protection_bites) %>%
  mutate(Draw = seq_len(nrow(protection_bites))) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(fit)
gc()

setDT(protection_bites)
protection_bites[, variable := str_remove(variable, "bites_pred_protection\\[")]
protection_bites[, variable := str_remove(variable, "\\]")]
protection_bites[, c("protection", "individual") := tstrsplit(variable, ",")]
protection_bites[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside", "Inside")
)]
protection_bites[, variable := NULL]
protection_bites <- protection_bites[value > 0] # filter out zero bites

dcast(protection_bites, Draw + individual ~ protection, value.var = "value") -> protection_bites
protection_bites[, Delta := log2(Inside + 1) - log2(Outside + 1)]
protection_bites[, P := Inside > Outside] # P(Inside > Outside)

posterior_protection_bites <- protection_bites[, .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = Draw]

posterior_protection_bites[, response := "Bites"] # add response column for bites

posterior_protection <- full_join(posterior_protection, posterior_protection_bites)


# Summarise posteriors
summary_protection <- posterior_protection %>%
  group_by(response) %>%
  summarise(
    Mean = mean(mean, na.rm = TRUE),
    Median = median(mean, na.rm = TRUE),
    Q025 = quantile(mean, 0.10, na.rm = TRUE),
    Q250 = quantile(mean, 0.25, na.rm = TRUE),
    Q750 = quantile(mean, 0.75, na.rm = TRUE),
    Q975 = quantile(mean, 0.90, na.rm = TRUE),
    P = mean(P, na.rm = TRUE)
  )

summary_protection_table <- summary_protection %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(Behaviour = response) %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Median", "Q025", "Q250", "Q750", "Q975"),
    digits = 3,
    na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2:8, align = "right", part = "body") %>%
  bold(part = "header") %>%
  theme_box()

print(summary_protection_table)

save_as_docx(
  summary_protection_table,
  path = here(output_dir, "counterfactual_protection_summary.docx")
)

# plot counterfactual protection

ggplot(summary_protection, aes(y = Median, x = response, col = response)) +
  geom_errorbar(aes(ymin = Q025, ymax = Q975), width = 0.1, linewidth = 1) +
  geom_errorbar(aes(ymin = Q250, ymax = Q750), width = 0, linewidth = 2) +
  geom_point(size = 4) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  labs(
    y = "Log-fold change",
    x = "Density"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
  ) -> p1

ggsave(
  here(figure_dir, "counterfactual_protection.png"),
  plot = p1,
  width = 5,
  height = 5
)

rm(
  p1,
  posterior_protection,
  posterior_protection_bites,
  summary_protection_table
)
gc()

# add individual information

protection_bites[, response := "Bites"]

protection <- full_join(protection_df, protection_bites)

rm(protection_df, protection_bites)
gc()

protection <- protection %>%
  mutate(individual = as.integer(individual)) %>%
  left_join(fish, by = c("individual" = "ind_id"))

setDT(protection)

protection[, size_class := factor(size_class,
  levels = seq_along(size_levels),
  labels = size_levels
)]

posterior_protection_individual <- protection[guild %in% c("Herbivore", "Invertivore"), .(
  mean = mean(Delta, na.rm = TRUE),
  P = mean(Delta > 0, na.rm = TRUE)
), by = .(Draw, response, guild, size_class)]

posterior_protection_individual[, .(
  mean = median(mean, na.rm = TRUE),
  q025 = quantile(mean, 0.11, na.rm = TRUE),
  q250 = quantile(mean, 0.25, na.rm = TRUE),
  q750 = quantile(mean, 0.75, na.rm = TRUE),
  q975 = quantile(mean, 0.89, na.rm = TRUE),
  P = mean(P, na.rm = TRUE)
), by = .(response, guild, size_class)] -> summary_protection_individual

p <- summary_protection_individual %>%
  complete() %>%
  ggplot(aes(x = size_class, y = mean, col = guild)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = q025, ymax = q975),
    width = 0,
    linewidth = 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    aes(ymin = q250, ymax = q750),
    width = 0,
    linewidth = 1.5, position = position_dodge(width = 0.5)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  geom_line(aes(group = guild), linewidth = 0.5, position = position_dodge(width = 0.5), linetype = "dashed") +
  facet_wrap(~ interaction(guild, response, sep = " | "), scales = "free", ncol = 2) +
  labs(
    x = "Size Class",
    y = "Log-fold change",
    col = "Foraging Guild"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    legend.position = "none",
  )

ggsave(
  here(figure_dir, "counterfactual_protection_summary.png"),
  plot = p,
  width = 6,
  height = 8
)

# summary across guilds only

posterior_protection_guild <- protection[guild %in% c("Herbivore", "Invertivore"), .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = mean(Delta > 0, na.rm = TRUE)
), by = .(Draw, response, guild)]

summary_protection_guild <- posterior_protection_guild %>%
  group_by(response, guild) %>%
  summarise(
    Mean = mean(mean, na.rm = TRUE),
    Median = median(mean, na.rm = TRUE),
    q025 = quantile(mean, 0.025, na.rm = TRUE),
    q250 = quantile(mean, 0.25, na.rm = TRUE),
    q750 = quantile(mean, 0.75, na.rm = TRUE),
    q975 = quantile(mean, 0.975, na.rm = TRUE),
    P = mean(P, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(Behaviour = response, Guild = guild) %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Median", "q025", "q250", "q750", "q975", "P"),
    digits = 3,
    na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2:3, align = "right", part = "body") %>%
  merge_v(j = 1) %>%
  bold(part = "header") %>%
  theme_box()

print(summary_protection_guild)

save_as_docx(
  summary_protection_guild,
  path = here(output_dir, "counterfactual_protection_guild_summary.docx")
)

# summary across size classes only

posterior_protection_size <- protection[guild %in% c("Herbivore", "Invertivore"), .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = mean(Delta > 0, na.rm = TRUE)
), by = .(Draw, response, size_class)]

summary_protection_size <- posterior_protection_size %>%
  group_by(response, size_class) %>%
  summarise(
    Mean = mean(mean, na.rm = TRUE),
    Median = median(mean, na.rm = TRUE),
    q025 = quantile(mean, 0.025, na.rm = TRUE),
    q250 = quantile(mean, 0.25, na.rm = TRUE),
    q750 = quantile(mean, 0.75, na.rm = TRUE),
    q975 = quantile(mean, 0.975, na.rm = TRUE),
    P = mean(P, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(Behaviour = response, Size_Class = size_class) %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Median", "q025", "q250", "q750", "q975", "P"),
    digits = 3,
    na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2:3, align = "right", part = "body") %>%
  bold(part = "header") %>%
  merge_v(j = 1) %>%
  theme_box()

print(summary_protection_size)

save_as_docx(
  summary_protection_size,
  path = here(output_dir, "counterfactual_protection_size_summary.docx")
)
