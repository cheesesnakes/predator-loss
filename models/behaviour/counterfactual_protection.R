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
  
# load data
source(here("analysis", "standardise_behaviour.R"))

traits <- traits %>%
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

fit <- readRDS(here("outputs", "standardised_behaviour_fit.rds"))

protection <- fit$draws("D_pred_protection", format = "matrix")

protection_df <- as.data.table(protection) %>%
  mutate(Draw = 1:nrow(protection)) %>%
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
protection_df[, Delta := Inside - Outside]
protection_df[, P := Inside > Outside] # P(Inside > Outside)

summary_protection <- protection_df[, .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = .(Draw, response)]

p1 <- summary_protection %>%
  pivot_longer(
    cols = -c(Draw, response),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = str_to_title(statistic)) %>%
  filter(statistic != "Mean" &
    statistic != "Median" & statistic != "P") %>%
  ggplot(aes(x = value, y = statistic, fill = response)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(-120, 120)) +
  labs(y = "Probability Density", x = expression(Delta ~ "Duration (s)")) +
  scale_fill_brewer(palette = "Set1", name = "Behaviour") +
  theme_bw() +
  theme(legend.position = "none")

p2 <- summary_protection %>%
  ggplot(aes(x = mean, fill = response)) +
  geom_density(alpha = 0.5) +
  labs(y = "Probability Density", x = expression("Mean" ~ Delta ~ "Duration (s)")) +
  scale_fill_brewer(palette = "Set1", name = "Behaviour") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(-120, 120)) +
  theme_bw() +
  theme(legend.position = "top")

p3 <- summary_protection %>%
  ggplot(aes(x = P, fill = response)) +
  geom_density(alpha = 0.5) +
  labs(y = "Probability Density", x = "P(Inside > Outside)") +
  scale_fill_brewer(palette = "Set1", name = "Behaviour") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none")

# counterfactual protection bites
fit <- readRDS(here("outputs", "standardised_behaviour_fit.rds"))
protection_bites <- fit$draws("bites_pred_protection", format = "matrix")
protection_bites <- as.data.table(protection_bites) %>%
  mutate(Draw = 1:nrow(protection_bites)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(fit, protection_bites, protection_df)
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
protection_bites[, Delta := Inside - Outside]
protection_bites[, P := Inside > Outside] # P(Inside > Outside)

summary_protection_bites <- protection_bites[, .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = Draw]

summary_protection_bites[, response := "Bites"] # add response column for bites

summary_protection_table <- summary_protection %>%
  filter(mean < 1.2 * 1e2, mean > -1.2 * 1e2) %>%
  full_join(summary_protection_bites) %>%
  group_by(response) %>%
  summarise(
    Mean = mean(mean, na.rm = TRUE),
    Median = median(mean, na.rm = TRUE),
    Q025 = quantile(mean, 0.025, na.rm = TRUE),
    Q250 = quantile(mean, 0.25, na.rm = TRUE),
    Q750 = quantile(mean, 0.75, na.rm = TRUE),
    Q975 = quantile(mean, 0.975, na.rm = TRUE),
    P = mean(P, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(Behaviour = response) %>%
  flextable() %>%
  set_table_properties(width = .95, layout = "autofit") %>%
  colformat_double(
    j = c("Mean", "Median", "Q025", "Q250", "Q750", "Q975"),
    digits = 2,
    na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2:8, align = "right", part = "body") %>%
  bold(part = "header") %>%
  autofit() %>%
  theme_box()

print(summary_protection_table)

save_as_docx(
  summary_protection_table,
  path = here("outputs", "counterfactual_protection_summary.docx")
)

p4 <- summary_protection_bites %>%
  pivot_longer(
    cols = -c(Draw, response),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = str_to_title(statistic)) %>%
  filter(statistic != "Mean" &
    statistic != "Median" & statistic != "P") %>%
  ggplot(aes(x = value, y = statistic)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(y = "Probability Density", x = expression(Delta ~ "Bite rate")) +
  theme_bw()

p5 <- summary_protection_bites %>%
  ggplot(aes(x = mean)) +
  geom_density(alpha = 0.5, fill = "black") +
  labs(y = "Probability Density", x = expression("Mean" ~ Delta ~ "Bite rate")) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  theme_bw()

p6 <- summary_protection_bites %>%
  ggplot(aes(x = P)) +
  geom_density(alpha = 0.5, fill = "black") +
  labs(y = "Probability Density", x = "P(Inside > Outside)") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()

(p2 + p3) / (p5 + p6) +
  plot_layout(guides = "collect")

ggsave(
  here("figures", "counterfactual_protection.png"),
  width = 6,
  height = 8
)

rm(
  p1,
  p2,
  p3,
  p4,
  p5,
  p6,
  summary_protection,
  summary_protection_bites,
  summary_protection_table,
  summary_protection_bites_table
)
gc()

# add individual information

protection_bites[, Outside := NULL]
protection_bites[, Inside := NULL]
protection_bites[, response := "Bites"]

protection_df[, Outside := NULL]
protection_df[, Inside := NULL]

protection <- full_join(protection_df, protection_bites)

rm(protection_df, protection_bites)
gc()

protection <- protection %>%
  mutate(individual = as.integer(individual)) %>%
  left_join(fish, by = c("individual" = "ind_id"))

setDT(protection)
protection[, size_class := factor(size_class,
  levels = 1:length(size_levels),
  labels = size_levels
)]

summary_protection <- protection[guild %in% c("Herbivore", "Invertivore"), .(
  mean = mean(Delta, na.rm = TRUE),
  P = mean(Delta > 0, na.rm = TRUE)
), by = .(Draw, response, guild, size_class)]

p1 <- summary_protection %>%
  filter(response != "Bites") %>%
  ggplot(aes(x = mean, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = expression("Mean" ~ Delta ~ "Duration (s)")) +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(-120, 120)) +
  facet_wrap(~response, ncol = 1) +
  theme_bw() +
  theme(legend.position = "top")

p2 <- summary_protection %>%
  filter(response != "Bites") %>%
  ggplot(aes(x = P, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = "P(Inside > Outside)") +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(~response, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none")

p3 <- summary_protection %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = mean, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = expression("Mean" ~ Delta ~ "Bite rate")) +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_bw() +
  theme(legend.position = "top")
p4 <- summary_protection %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = P, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = "P(Inside > Outside)") +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none")

p <- (p1 + p2) / (p3 + p4) +
  plot_layout(
    ncol = 1,
    heights = c(3, 1),
    # Top row gets 75%, bottom 25%
    guides = "collect"
  ) &
  theme(legend.position = "top")

ggsave(
  here("figures", "counterfactual_protection_individuals.png"),
  p,
  width = 6,
  height = 8
)
