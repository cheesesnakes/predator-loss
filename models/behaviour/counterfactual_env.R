pacman::p_load(
  here,
  dplyr,
  tidyr,
  flextable,
  patchwork,
  stringr,
  ggplot2,
  data.table,
  ggridges,
  tictoc
)

theme_set(theme_bw() +
  theme(legend.position = "top", legend.text = element_text(size = 20)))

# load data
print("Loading data\n\n")
tic()

source(here("analysis", "standardise_behaviour.R"))

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
  select(-c(group, species, size_class, plot_id))

predictors <- predictors %>%
  select(c(plot_id, biomass, rugosity, predator))

rugosity_raw <- rugosity

rm(rugosity, responses, traits)
gc()

toc()

# counterfactual protection varying biomass

print("Counterfactual protection varying biomass")
print("============================================\n")
print("Loading fit...\n")
tic()

fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))
protection_biomass <- fit$draws("D_pred_protection_biomass", format = "matrix")
protection_biomass_bites <- fit$draws("bites_pred_protection_biomass", format = "matrix")

sample <- sample(1:4000, 100)

rm(fit)
gc()

biomass <- as.data.table(protection_biomass) %>%
  mutate(Draw = 1:nrow(protection_biomass)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )
rm(protection_biomass)
gc()

setDT(biomass)
biomass <- biomass[Draw %in% sample]
biomass <- biomass[value < 1000]
biomass[, variable := str_remove(variable, "D_pred_protection_biomass\\[")]
biomass[, variable := str_remove(variable, "\\]")]
biomass[, c("protection", "sample", "individual", "response") := tstrsplit(variable, ",")]

biomass[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
biomass[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
biomass[, variable := NULL]

biomass_bites <- as.data.table(protection_biomass_bites) %>%
  mutate(Draw = 1:nrow(protection_biomass_bites)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection_biomass_bites)
gc()

setDT(biomass_bites)
biomass_bites <- biomass_bites[Draw %in% sample]
biomass_bites[, variable := str_remove(variable, "bites_pred_protection_biomass\\[")]
biomass_bites[, variable := str_remove(variable, "\\]")]
biomass_bites[, c("protection", "sample", "individual") := tstrsplit(variable, ",")]
biomass_bites[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
biomass_bites[, variable := NULL]
biomass_bites[, response := "Bites"]
biomass_bites <- biomass_bites[value > 0] # filter out zero bites

biomass <- full_join(biomass, biomass_bites)

mean <- mean(predictors$biomass)
sd <- sd(predictors$biomass)
std_var <- seq(-2, 2, length.out = max(as.integer(biomass$sample))) # standarised variables 

biomass <- biomass %>%
  mutate(
    individual = as.integer(individual),
    sample = as.integer(sample),
    biomass = (mean + std_var[sample] * sd)
  )%>%
  left_join(fish, by = c("individual" = "ind_id"))

setDT(biomass)
biomass <- biomass[!is.na(biomass)]
biomass <- biomass[!is.na(guild)]
biomass <- biomass[guild != "Piscivore"]

toc()

print("Plotting counterfactual protection varying biomass...\n")
tic()

p1 <- biomass %>%
  filter(response == "Foraging") %>%
  ggplot(aes(x = biomass * 100, y = value, col = guild, group = interaction(Draw, guild))) +
  # Raw posterior draws, colored by guild
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15,
  ) +
  # Mean trend for each guild (one line per guild)
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p2 <- biomass %>%
  filter(response == "Vigilance") %>%
  ggplot(aes(x = biomass * 100, y = value, col = guild, group = interaction(Draw, guild))) +
  # Raw posterior draws, colored by guild
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15,
  ) +
  # Mean trend for each guild (one line per guild)
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p3 <- biomass %>%
  filter(response == "Movement") %>%
  ggplot(aes(x = biomass * 100, y = value, col = guild, group = interaction(Draw, guild))) +
  # Raw posterior draws, colored by guild
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15,
  ) +
  # Mean trend for each guild (one line per guild)
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p4 <- biomass %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = biomass * 100, y = value, col = guild, group = interaction(Draw, guild))) +
  # Raw posterior draws, colored by guild
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15,
  ) +
  # Mean trend for each guild (one line per guild)
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Biomass Cover (%)", y = expression("Bite rate")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

g1 <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  plot = g1,
  file = here("figures", "counterfactual_protection_biomass.png"),
  width = 8,
  height = 8
)

rm(biomass, p1, p2, p3, p4, protection_biomass, protection_biomass_bites)
gc()

toc()

# counterfactual protection varying rugosity
print("Counterfactual protection varying rugosity")
print("============================================\n")
print("Loading fit\n")
tic()

fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))

protection_rugosity <- fit$draws("D_pred_protection_rugosity", format = "matrix")
protection_rugosity_bites <- fit$draws("bites_pred_protection_rugosity", format = "matrix")

rm(fit)
gc()

rugosity <- as.data.table(protection_rugosity) %>%
  mutate(Draw = 1:nrow(protection_rugosity)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection_rugosity)
gc()

setDT(rugosity)
rugosity <- rugosity[Draw %in% sample]
rugosity <- rugosity[value < 1000]
rugosity[, variable := str_remove(variable, "D_pred_protection_rugosity\\[")]
rugosity[, variable := str_remove(variable, "\\]")]
rugosity[, c("protection", "sample", "individual", "response") := tstrsplit(variable, ",")]
rugosity[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
rugosity[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
rugosity[, variable := NULL]

rugosity_bites <- as.data.table(protection_rugosity_bites) %>%
  mutate(Draw = 1:nrow(protection_rugosity_bites)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )
rm(protection_rugosity_bites)
gc()
setDT(rugosity_bites)

rugosity_bites <- rugosity_bites[Draw %in% sample]
rugosity_bites[, variable := str_remove(variable, "bites_pred_protection_rugosity\\[")]
rugosity_bites[, variable := str_remove(variable, "\\]")]
rugosity_bites[, c("protection", "sample", "individual") := tstrsplit(variable, ",")]
rugosity_bites[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
rugosity_bites[, variable := NULL]
rugosity_bites[, response := "Bites"]
rugosity_bites <- rugosity_bites[value > 0] # filter out zero bites

rugosity <- full_join(rugosity, rugosity_bites)

rm(rugosity_bites)
gc()

head(rugosity)

mean <- mean(c(rugosity_raw$sample_1, rugosity_raw$sample_2, rugosity_raw$sample_3))
sd <- sd(c(rugosity_raw$sample_1, rugosity_raw$sample_2, rugosity_raw$sample_3))

rugosity <- rugosity %>%
  mutate(
    individual = as.integer(individual),
    sample = as.integer(sample),
    rugosity = (mean + std_var[sample] * sd)
  ) %>%
  left_join(fish, by = c("individual" = "ind_id"))


setDT(rugosity)
rugosity <- rugosity[!is.na(rugosity)]
rugosity <- rugosity[!is.na(guild)]
rugosity <- rugosity[guild != "Piscivore"]
sample <- sample(1:4000, 100)

toc()

print("Plotting counterfactual protection varying rugosity...\n")
tic()

p1 <- rugosity %>%
  filter(response == "Foraging") %>%
  ggplot(aes(x = rugosity, y = value, col = guild, group = interaction(Draw, guild))) +
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15
  ) +
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p2 <- rugosity %>%
  filter(response == "Vigilance") %>%
  ggplot(aes(x = rugosity, y = value, col = guild, group = interaction(Draw, guild))) +
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15
  ) +
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p3 <- rugosity %>%
  filter(response == "Movement") %>%
  ggplot(aes(x = rugosity, y = value, col = guild, group = interaction(Draw, guild))) +
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15
  ) +
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Duration (s)")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

p4 <- rugosity %>%
  filter(response == "Bites") %>%
  ggplot(aes(x = rugosity, y = value, col = guild, group = interaction(Draw, guild))) +
  geom_line(
    stat = "smooth",
    method = "gam",
    formula = y ~ s(x),
    alpha = 0.15
  ) +
  geom_smooth(
    aes(group = guild),
    method = "gam",
    formula = y ~ s(x),
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~response, ncol = 1, scale = "free_y") +
  labs(x = "Rugosity", y = expression("Bite rate")) +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())

g2 <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  plot = g2,
  file = here("figures", "counterfactual_protection_rugosity.png"),
  width = 8,
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

toc()

# counterfactual protection predator presence
print("Counterfactual protection predator presence")
print("============================================\n")
print("Loading fit...\n")
tic()

fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))
protection_predator <- fit$draws("D_pred_protection_predator", format = "matrix")
protection_predator_bites <- fit$draws("bites_pred_protection_predator", format = "matrix")

rm(fit)
gc()

predator <- as.data.table(protection_predator) %>%
  mutate(Draw = 1:nrow(protection_predator)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection_predator)
gc()

setDT(predator)
predator[, variable := str_remove(variable, "D_pred_protection_predator\\[")]
predator[, variable := str_remove(variable, "\\]")]
predator[, c("protection", "predator", "individual", "response") := tstrsplit(variable, ",")]
predator[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
predator[, response := factor(
  response,
  levels = 1:3,
  labels = c("Foraging", "Vigilance", "Movement")
)]
predator[, variable := NULL]

predator_bites <- as.data.table(protection_predator_bites) %>%
  mutate(Draw = 1:nrow(protection_predator_bites)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection_predator_bites)
gc()

setDT(predator_bites)
predator_bites[, variable := str_remove(variable, "bites_pred_protection_predator\\[")]
predator_bites[, variable := str_remove(variable, "\\]")]
predator_bites[, c("protection", "predator", "individual") := tstrsplit(variable, ",")]
predator_bites[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside PA", "Inside PA")
)]
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
  labels = c("Absent", "Present")
)]

dcast(predator, Draw + individual + response + protection ~ predator, value.var = "value") -> predator
predator[, value := log2(Present + 1) - log2(Absent + 1)]
predator[, P := Present > Absent] # P(Present > Absent)

predator <- predator %>%
  mutate(
    individual = as.integer(individual),
    Draw = as.integer(Draw)
  ) %>%
  left_join(fish, by = c("individual" = "ind_id"))

setDT(predator)
predator <- predator[!is.na(guild)]
predator <- predator[guild != "Piscivore"]

toc()

print("Generating summary table...\n")
tic()

summary_predator <- predator[, .(
  mean = mean(value, na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  q025 = quantile(value, 0.025, na.rm = TRUE),
  q250 = quantile(value, 0.25, na.rm = TRUE),
  q750 = quantile(value, 0.75, na.rm = TRUE),
  q975 = quantile(value, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = .(Draw, response, guild)]

summary_predator_table <- summary_predator %>%
  group_by(response, guild) %>%
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
  set_table_properties(width = .95, layout = "fixed") %>%
  colformat_double(
    j = c("Mean", "Median", "Q025", "Q250", "Q750", "Q975"),
    digits = 2,
    na_str = "-"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 2:8, align = "right", part = "body") %>%
  bold(part = "header") %>%
  theme_box()

save_as_docx(
  summary_predator_table,
  path = here("outputs", "counterfactual_predator_summary.docx")
)

toc()

print("Plotting counterfactual protection predator presence...\n")
tic()

g3 <- summary_predator %>%
  ggplot(aes(y = mean, x = response)) +
  geom_boxplot(
    aes(fill = guild),
    alpha = 0.5,
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  labs(x = "Behaviours", y = "Log-fold change") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), )

ggsave(
  plot = g3,
  file = here("figures", "counterfactual_protection_predator.png"),
  width = 5,
  height = 5
)

gc()

toc()

# Combine all plots into one figure

print("Combining all plots into one figure...\n")
tic()

g <- (g1 / g2) + plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.title = element_blank())

ggsave(
  here("figures", "figure5.png"),
  g,
  width = 12,
  height = 16
)

toc()
