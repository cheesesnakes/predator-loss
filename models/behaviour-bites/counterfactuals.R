pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))
source(here("functions", "summaries.R"))

# Set default ggplot theme
theme_set(
    theme_bw() +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 14),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )
)

dir.create(here("outputs", "behaviour-bites"), showWarnings = FALSE)
dir.create(here("figures", "behaviour-bites"), showWarnings = FALSE)

# analysis

# load models

model <- readRDS(here("outputs", "behaviour-bites", "model.rds"))

# plot posterior proportions by protection

counterfactual_protection <- data_ts %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
    select(ind_id, plot_id, site, family, species, guild, size_class, predator, group, rugosity_mean, rugosity_se, biomass) %>%
    distinct() %>%
    expand_grid(protected = c(0, 1)) %>%
    add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(behaviour = "Bites")

# counterfactual rugosity

std_var <- seq(-2, 2, by = 0.2)

counterfactual_rugosity <- data_ts %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
    select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, group, biomass) %>%
    distinct() %>%
    expand_grid(rugosity_mean = std_var, rugosity_se = 0.1) %>%
    add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(behaviour = "Bites")


# counterfactual biomass
std_var <- seq(-2, 2, by = 0.2)

counterfactual_biomass <- data_ts %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
    select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, group, rugosity_mean, rugosity_se) %>%
    distinct() %>%
    expand_grid(biomass = std_var) %>%
    add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(behaviour = "Bites")

# counterfactual predator

counterfactual_predator <- data_ts %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
    select(ind_id, plot_id, family, species, site, guild, size_class, protected, group, rugosity_mean, rugosity_se, biomass) %>%
    distinct() %>%
    expand_grid(predator = c(0, 1)) %>%
    add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(behaviour = "Bites")

# counterfactual group

counterfactual_group <- data_ts %>%
    ungroup() %>%
    mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
    select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, rugosity_mean, rugosity_se, biomass) %>%
    distinct() %>%
    expand_grid(group = c(0, 1)) %>%
    add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
    ungroup() %>%
    mutate(behaviour = "Bites")

# Counterfactual effect of protection on time budget

effect_protection <- summarise_cf(posterior = counterfactual_protection, dist = "Abundance", do = "protected", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_protection.docx"))

effect_protection %>%
    ggplot(aes(x = effect)) +
    geom_density(fill = "blue", alpha = 0.5, col = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = median(effect_protection$effect), color = "red") +
    labs(x = "Log response ratio", y = "Density")

ggsave(here("figures", "behaviour-bites", "effect_plot.png"), plot = last_plot())

## Across guilds
effect_protection <- summarise_cf(posterior = counterfactual_protection, dist = "Abundance", do = "protected", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_protection_guild.docx"))

effect_protection %>%
    group_by(behaviour, guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = guild, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Foraging", y = "Log response ratio", color = "Guild")

ggsave(here("figures", "behaviour-bites", "effect_plot_guild.png"), plot = last_plot(), height = 6, width = 8)

## Across guilds and size classes
effect_protection <- summarise_cf(posterior = counterfactual_protection, dist = "Abundance", do = "protected", confound = c("behaviour", "guild", "size_class"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_protection_size.docx"))

effect_protection %>%
    group_by(behaviour, guild, size_class) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = size_class, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Size class", y = "Log response ratio", color = "Guild")

ggsave(here("figures", "behaviour-bites", "effect_plot_guild_size.png"), plot = last_plot(), height = 6, width = 8)

# Effect of rugosity

counterfactual_rugosity %>%
    group_by(behaviour, guild, rugosity_mean) %>%
    median_qi(prop, .width = c(0.5)) %>%
    mutate(rugosity_mean = moments_rugosity[2] * rugosity_mean + moments_rugosity[1]) %>%
    ggplot(aes(x = rugosity_mean, y = prop, col = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Rugosity", y = "Bite Rate") +
    facet_wrap(~guild, ncol = 2) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-bites", "counterfactual_rugosity.png"), plot = last_plot(), height = 6, width = 12)

effect_rugosity <- summarise_cf(posterior = counterfactual_rugosity, dist = "Abundance", do = "rugosity_mean", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_rugosity.docx"))

# Effect of biomass

counterfactual_biomass %>%
    group_by(behaviour, guild, biomass) %>%
    median_qi(prop, .width = c(0.5)) %>%
    mutate(biomass = (moments_biomass[2] * biomass + moments_biomass[1]) * 100) %>%
    ggplot(aes(x = biomass, y = prop, col = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Biomass Cover (%)", y = "Bite Rate") +
    facet_wrap(~guild, ncol = 2) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-bites", "counterfactual_biomass.png"), plot = last_plot(), height = 6, width = 12)

effect_biomass <- summarise_cf(posterior = counterfactual_biomass, dist = "Abundance", do = "biomass", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_biomass.docx"))

# Effect of predator presence
effect_predator <- summarise_cf(posterior = counterfactual_predator, dist = "Abundance", do = "predator", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_predator.docx"))

effect_predator %>%
    ggplot(aes(x = effect)) +
    geom_density(fill = "blue", alpha = 0.5, col = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = median(effect_predator$effect), color = "red") +
    labs(x = "Log response ratio", y = "Density")

ggsave(here("figures", "behaviour-bites", "counterfactual_predator.png"), plot = last_plot(), height = 12, width = 8)

# Effect of grouping
effect_group <- summ a rise_cf(posterior = counterfactual_group, dist = "Abundance", do = "group", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-bites", file = paste0("counterfactual_group.docx"))

effect_group %>%
    ggplot(aes(x = effect)) +
    geom_density(fill = "blue", alpha = 0.5, col = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = median(effect_group$effect), color = "red") +
    labs(x = "Log response ratio", y = "Density")

ggsave(here("figures", "behaviour-bites", "counterfactual_group.png"), plot = last_plot(), height = 12, width = 8)
