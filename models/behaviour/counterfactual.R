pacman::p_load(here,
               dplyr,
               tidyr,
               flextable,
               patchwork,
               stringr,
               ggplot2,
               data.table,
               ggridges)

theme_set(theme_bw() +
            theme(legend.position = "top", legend.text = element_text(size = 20)))

# load data
source(here("analysis", "standardise_behaviour.R"))

traits <- traits %>%
  pivot_longer(cols = -c(species),
               names_to = "guild",
               values_to = "value") %>%
  filter(value > 0) %>%
  group_by(species) %>%
  summarise(guild = first(guild))


fish <- fish %>%
  left_join(traits, by = "species") %>%
  select(-c(group, species))

predictors <- predictors %>%
  select(c(plot_id, biomass, rugosity, predator))

rm(rugosity, responses, traits)

# extract counterfactual draws

# predictive checks
fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))

D_pred <- fit$draws("D_pred", format = "matrix")

D_pred_df <- as.data.table(D_pred) %>%
  mutate(Draw = 1:nrow(D_pred)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")
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
  scale_x_continuous(limits = c(0, 120)) +
  scale_fill_brewer(palette = "Set1", name = "Behaviour") +
  theme_bw() +
  theme(legend.position = "top")

p2 <- summary_Dpred %>%
  ggplot(aes(x = mean, fill = response)) +
  geom_density(alpha = 0.5) +
  labs(y = "Probability Density", x = "Mean Duration (s)") +
  scale_x_continuous(limits = c(0, 120)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none")
p1 / p2

ggsave(here("figures", "D_pred.png"),
       width = 6,
       height = 6)

bites_pred <- fit$draws("bites_pred", format = "matrix")

bites_pred_df <- as.data.table(bites_pred) %>%
  mutate(Draw = 1:nrow(bites_pred)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

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
  pivot_longer(cols = -Draw,
               names_to = "statistic",
               values_to = "value") %>%
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
       height = 6)

rm(D_pred,
   bites_pred,
   D_pred_df,
   bites_pred_df,
   summary_Dpred,
   summary_bites,
   p1,
   p2)
gc()

## counterfactual protection

protection <- fit$draws("D_pred_protection", format = "matrix")

protection_df <- as.data.table(protection) %>%
  mutate(Draw = 1:nrow(protection)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

rm(protection)
gc()

setDT(protection_df)
protection_df[, variable := str_remove(variable, "D_pred_protection\\[")]
protection_df[, variable := str_remove(variable, "\\]")]
protection_df[, c("protection", "individual", "response") := tstrsplit(variable, ",")]
protection_df[, protection := factor(protection,
                                     levels = 1:2,
                                     labels = c("Outside", "Inside"))]
protection_df[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
protection_df[, variable := NULL]

dcast(protection_df,
      Draw + individual + response ~ protection,
      value.var = "value") -> protection_df
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
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
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
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
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
protection_bites <- fit$draws("bites_pred_protection", format = "matrix")
protection_bites <- as.data.table(protection_bites) %>%
  mutate(Draw = 1:nrow(protection_bites)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")
setDT(protection_bites)
protection_bites[, variable := str_remove(variable, "bites_pred_protection\\[")]
protection_bites[, variable := str_remove(variable, "\\]")]
protection_bites[, c("protection", "individual") := tstrsplit(variable, ",")]
protection_bites[, protection := factor(protection,
                                        levels = 1:2,
                                        labels = c("Outside", "Inside"))]
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
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(y = "Probability Density", x = expression(Delta ~ "Bite rate")) +
  theme_bw()

p5 <- summary_protection_bites %>%
  ggplot(aes(x = mean)) +
  geom_density(alpha = 0.5, fill = "black") +
  labs(y = "Probability Density", x = expression("Mean" ~ Delta ~ "Bite rate")) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  theme_bw()

p6 <- summary_protection_bites %>%
  ggplot(aes(x = P)) +
  geom_density(alpha = 0.5, fill = "black") +
  labs(y = "Probability Density", x = "P(Inside > Outside)") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()

(p1 + p2 + p3) / (p4 + p5 + p6) +
  plot_layout(guides = "collect")

ggsave(
  here("figures", "counterfactual_protection.png"),
  width = 18,
  height = 12
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
                                  labels = size_levels)]

summary_protection <- protection[guild %in% c("Herbivore", "Invertivore"), .(mean = mean(Delta, na.rm = TRUE),
                                                                             P = mean(Delta > 0, na.rm = TRUE)), by = .(Draw, response, guild, size_class)]

p1 <- summary_protection %>%
  filter(response != "Bites") %>%
  ggplot(aes(x = mean, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = expression("Mean" ~ Delta ~ "Duration (s)")) +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_x_continuous(limits = c(-120, 120)) +
  facet_wrap( ~ response, ncol = 1) +
  theme_bw() +
  theme(legend.position = "top")

p2 <- summary_protection %>%
  filter(response != "Bites") %>%
  ggplot(aes(x = P, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = "P(Inside > Outside)") +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap( ~ response, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none")

p3 <- summary_protection %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = mean, y = size_class, fill = guild)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  labs(y = "Size Class", x = expression("Mean" ~ Delta ~ "Bite rate")) +
  scale_fill_brewer(palette = "Set1", name = "Feeding Guild") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
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
  width = 12,
  height = 12
)

rm(protection, summary_protection, p1, p2, p3, p4)
gc()

# counterfactual protection varying biomass

protection_biomass <- fit$draws("D_pred_protection_biomass", format = "matrix")
protection_biomass_bites <- fit$draws("bites_pred_protection_biomass", format = "matrix")

sample <- sample(1:4000, 100)
rm(fit)
gc()

biomass <- as.data.table(protection_biomass) %>%
  mutate(Draw = 1:nrow(protection_biomass)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")
rm(protection_biomass)
gc()

setDT(biomass)
biomass <- biomass[Draw %in% sample]
biomass[, variable := str_remove(variable, "D_pred_protection_biomass\\[")]
biomass[, variable := str_remove(variable, "\\]")]
biomass[, c("protection", "sample", "individual", "response") := tstrsplit(variable, ",")]

biomass[, protection := factor(protection,
                               levels = 1:2,
                               labels = c("Outside PA", "Inside PA"))]
biomass[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
biomass[, variable := NULL]
#dcast(biomass, Draw + individual + response ~ protection, value.var = "value") -> biomass
#biomass[, Delta := Inside - Outside]

biomass_bites <- as.data.table(protection_biomass_bites) %>%
  mutate(Draw = 1:nrow(protection_biomass_bites)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

rm(protection_biomass_bites)
gc()

setDT(biomass_bites)
biomass_bites <- biomass_bites[Draw %in% sample]
biomass_bites[, variable := str_remove(variable, "bites_pred_protection_biomass\\[")]
biomass_bites[, variable := str_remove(variable, "\\]")]
biomass_bites[, c("protection", "sample", "individual") := tstrsplit(variable, ",")]
biomass_bites[, protection := factor(protection,
                                     levels = 1:2,
                                     labels = c("Outside", "Inside"))]
biomass_bites[, variable := NULL]
biomass_bites[, response := "Bites"]
biomass_bites <- biomass_bites[value > 0] # filter out zero bites

biomass <- full_join(biomass, biomass_bites)

mean <- mean(predictors$biomass)
sd <- sd(predictors$biomass)
std_var <- seq(-1, 1, length.out = max(as.integer(biomass$sample)))

biomass <- biomass %>%
  mutate(individual = as.integer(individual),
  sample = as.integer(sample),
  biomass = (mean + std_var[sample] * sd))

setDT(biomass)
biomass <- biomass[!is.na(biomass)]

p1 <- biomass %>%
  filter(response == "Foraging") %>%
  ggplot(aes(x = biomass * 100, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), )

p2 <- biomass %>%
  filter(response == "Vigilance") %>%
  ggplot(aes(x = biomass * 100, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p3 <- biomass %>%
  filter(response == "Movement") %>%
  ggplot(aes(x = biomass * 100, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p4 <- biomass %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = biomass * 100, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Bite rate")) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  plot = p,
  file = here("figures", "counterfactual_protection_biomass.png"),
  width = 6,
  height = 8
)

rm(biomass, p)
gc()

# counterfactual protection varying rugosity
fit <- readRDS(here("outputs/standardised_behaviour_fit.rds")) 

protection_rugosity <- fit$draws("D_pred_protection_rugosity", format = "matrix")
protection_rugosity_bites <- fit$draws("bites_pred_protection_rugosity", format = "matrix")

rm(fit)
gc()

rugosity <- as.data.table(protection_rugosity) %>%
  mutate(Draw = 1:nrow(protection_rugosity)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

rm(protection_rugosity)
gc()
setDT(rugosity)
rugosity <- rugosity[Draw %in% sample]
rugosity[, variable := str_remove(variable, "D_pred_protection_rugosity\\[")]
rugosity[, variable := str_remove(variable, "\\]")]
rugosity[, c("protection", "sample", "individual", "response") := tstrsplit(variable, ",")]
rugosity[, protection := factor(protection,
                                levels = 1:2,
                                labels = c("Outside", "Inside"))]
rugosity[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
rugosity[, variable := NULL]
#dcast(rugosity, Draw + individual + response ~ protection, value.var = "value") -> rugosity
#rugosity[, Delta := Inside - Outside]
rugosity_bites <- as.data.table(protection_rugosity_bites) %>%
  mutate(Draw = 1:nrow(protection_rugosity_bites)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")
rm(protection_rugosity_bites)
gc()
setDT(rugosity_bites)
rugosity_bites <- rugosity_bites[Draw %in% sample]
rugosity_bites[, variable := str_remove(variable, "bites_pred_protection_rugosity\\[")]
rugosity_bites[, variable := str_remove(variable, "\\]")]
rugosity_bites[, c("protection", "sample", "individual") := tstrsplit(variable, ",")]
rugosity_bites[, protection := factor(protection,
                                      levels = 1:2,
                                      labels = c("Outside", "Inside"))]
rugosity_bites[, variable := NULL]
rugosity_bites[, response := "Bites"]
rugosity_bites <- rugosity_bites[value > 0] # filter out zero bites

rugosity <- full_join(rugosity, rugosity_bites)

rm(rugosity_bites)
gc()

mean <- mean(predictors$rugosity)
sd <- sd(predictors$rugosity)

rugosity <- rugosity %>%
  mutate(individual = as.integer(individual),
  sample = as.integer(sample),
  rugosity = (mean + std_var[sample] * sd))
  

setDT(rugosity)
rugosity <- rugosity[!is.na(rugosity)]
sample <- sample(1:4000, 100)
rugosity <- rugosity[Draw %in% sample]

p1 <- rugosity %>%
  filter(response == "Foraging") %>%
  ggplot(aes(x = rugosity, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), )

p2 <- rugosity %>%
  filter(response == "Vigilance") %>%
  ggplot(aes(x = rugosity, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p3 <- rugosity %>%
  filter(response == "Movement") %>%
  ggplot(aes(x = rugosity, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p4 <- rugosity %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = rugosity, y = value, col = protection)) +
  geom_smooth(
    aes(group = Draw),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    alpha = 0.05,
    color = "gray"
  ) + # raw posterior draws
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) + # mean trend
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Bite rate")) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  plot = p,
  file = here("figures", "counterfactual_protection_rugosity.png"),
  width = 6,
  height = 8
)

rm(
  protection_rugosity,
  protection_rugosity_bites,
  rugosity,
  p1,
  p2,
  p3,
  p4
)

gc()

# counterfactual protection predator presence
fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))
protection_predator <- fit$draws("D_pred_protection_predator", format = "matrix")
protection_predator_bites <- fit$draws("bites_pred_protection_predator", format = "matrix")

rm(fit)
gc()

predator <- as.data.table(protection_predator) %>%
  mutate(Draw = 1:nrow(protection_predator)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

rm(protection_predator)
gc()

setDT(predator)
predator[, variable := str_remove(variable, "D_pred_protection_predator\\[")]
predator[, variable := str_remove(variable, "\\]")]
predator[, c("protection", "predator", "individual", "response") := tstrsplit(variable, ",")]
predator[, protection := factor(protection,
                                levels = 1:2,
                                labels = c("Outside", "Inside"))]
predator[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
predator[, variable := NULL]
#dcast(predator, Draw + individual + response ~ protection, value.var = "value") -> predator
#predator[, Delta := Inside - Outside]

predator_bites <- as.data.table(protection_predator_bites) %>%
  mutate(Draw = 1:nrow(protection_predator_bites)) %>%
  pivot_longer(cols = -c(Draw),
               names_to = "variable",
               values_to = "value")

rm(protection_predator_bites)
gc()

setDT(predator_bites)
predator_bites[, variable := str_remove(variable, "bites_pred_protection_predator\\[")]
predator_bites[, variable := str_remove(variable, "\\]")]
predator_bites[, c("protection", "predator", "individual") := tstrsplit(variable, ",")]
predator_bites[, protection := factor(protection,
                                      levels = 1:2,
                                      labels = c("Outside", "Inside"))]
predator_bites[, variable := NULL]
predator_bites[, response := "Bites"]
predator_bites <- predator_bites[value > 0] # filter out zero bites

predator <- full_join(predator, predator_bites)

rm(predator_bites)
gc()

setDT(predator)

predator <- predator[!is.na(predator)]
predator[, predator := factor(predator,
                              levels = c(1, 2),
                              labels = c("Absent", "Present"))]

summary_predator <- predator[,.(
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  q025 = quantile(value, 0.025, na.rm = TRUE),
  q250 = quantile(value, 0.25, na.rm = TRUE),
  q750 = quantile(value, 0.75, na.rm = TRUE),
  q975 = quantile(value, 0.975, na.rm = TRUE)
), by = .(Draw, response, predator, protection)]

p1 <- summary_predator %>%
  filter(response == "Foraging") %>%
  ggplot(aes(y = predator, x = mean, fill = protection)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Predator Presence", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p2 <- summary_predator %>%
  filter(response == "Vigilance") %>%
  ggplot(aes(y = predator, x = mean, fill = protection)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Predator Presence", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p3 <- summary_predator %>%
  filter(response == "Movement") %>%
  ggplot(aes(y = predator, x = mean, fill = protection)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Predator Presence", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0, 120)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p4 <- summary_predator %>%
  filter(response == "Bites") %>%
  ggplot(aes(y = predator, x = mean, fill = protection)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  facet_wrap( ~ response, ncol = 1, scale = "free_y") +
  labs(x = "Predator Presence", y = expression("Bite rate")) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank(), )

p <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  plot = p,
  file = here("figures", "counterfactual_protection_predator.png"),
  width = 6,
  height = 8
)
