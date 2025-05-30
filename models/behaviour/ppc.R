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

# extract counterfactual draws

# predictive checks
fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))

D_pred <- fit$draws("D_pred", format = "matrix")

D_pred_df <- as.data.table(D_pred) %>%
  mutate(Draw = 1:nrow(D_pred)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(D_pred, fit)
gc()

setDT(D_pred_df)
D_pred_df[, variable := str_remove(variable, "D_pred\\[")]
D_pred_df[, variable := str_remove(variable, "\\]")]
D_pred_df[, c("individual", "response") := tstrsplit(variable, ",")]
D_pred_df[, individual := as.integer(individual)]
D_pred_df[, response := as.integer(response)]
D_pred_df[, variable := NULL]

# set levels
D_pred_df[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]

summary_Dpred <- D_pred_df[, .(
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  q025 = quantile(value, 0.025, na.rm = TRUE),
  q250 = quantile(value, 0.25, na.rm = TRUE),
  q750 = quantile(value, 0.75, na.rm = TRUE),
  q975 = quantile(value, 0.975, na.rm = TRUE)
), by = .(Draw, response)]

p1 <- summary_Dpred %>%
  pivot_longer(
    cols = -c(Draw, response),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = str_to_title(statistic)) %>%
  filter(statistic != "Mean" & statistic != "Median") %>%
  ggplot(aes(x = value, y = statistic, fill = response)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Probability Density", x = "Duration (s)") +
  scale_x_continuous(limits = c(0, 60)) +
  scale_fill_brewer(palette = "Set1", name = "Behaviour") +
  theme_bw() +
  theme(legend.position = "top")

p2 <- summary_Dpred %>%
  ggplot(aes(x = mean, fill = response)) +
  geom_density(alpha = 0.5) +
  labs(y = "Probability Density", x = "Mean Duration (s)") +
  scale_x_continuous(limits = c(0, 60)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none")
p1 / p2

ggsave(here("figures", "D_pred.png"),
  width = 6,
  height = 6
)

# extract bite rate draws
fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))
bites_pred <- fit$draws("bites_pred", format = "matrix")

bites_pred_df <- as.data.table(bites_pred) %>%
  mutate(Draw = 1:nrow(bites_pred)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(bites_pred, fit)
gc()

setDT(bites_pred_df)
bites_pred_df[, variable := str_remove(variable, "bites_pred\\[")]
bites_pred_df[, variable := str_remove(variable, "\\]")]
bites_pred_df[, individual := as.integer(variable)]
bites_pred_df[, variable := NULL]

summary_bites <- bites_pred_df[value > 0, .(
  Mean = mean(value, na.rm = TRUE),
  Median = median(value, na.rm = TRUE),
  Q025 = quantile(value, 0.025, na.rm = TRUE),
  Q250 = quantile(value, 0.25, na.rm = TRUE),
  Q750 = quantile(value, 0.75, na.rm = TRUE),
  Q975 = quantile(value, 0.975, na.rm = TRUE)
), by = Draw]

p1 <- summary_bites %>%
  pivot_longer(
    cols = -Draw,
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = str_to_title(statistic)) %>%
  filter(statistic != "Mean" & statistic != "Median") %>%
  ggplot(aes(x = value, y = statistic)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(y = "Probability Density", x = "Bite rate") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()

p2 <- summary_bites %>%
  ggplot(aes(x = Mean)) +
  geom_density(alpha = 0.5) +
  labs(y = "Probability Density", x = "Mean Bite rate") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0, 5)) +
  theme_bw()

p1 / p2

ggsave(here("figures", "bites_pred.png"),
  width = 6,
  height = 6
)

rm(
  D_pred,
  bites_pred,
  D_pred_df,
  bites_pred_df,
  summary_Dpred,
  summary_bites,
  p1,
  p2
)
gc()


