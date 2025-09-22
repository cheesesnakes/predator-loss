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

dir.create(here("outputs", "behaviour-prob"), showWarnings = FALSE)
dir.create(here("figures", "behaviour-prob"), showWarnings = FALSE)

# analysis

behaviours <- c("foraging", "vigilance", "move")
counterfactual_protection <- data.frame()
counterfactual_rugosity <- data.frame()
counterfactual_biomass <- data.frame()
counterfactual_predator <- data.frame()
counterfactual_group <- data.frame()

for (b in behaviours) {
    # load models

    model <- readRDS(here("outputs", "behaviour-prob", paste0("model_", b, ".rds")))

    # plot posterior proportions by protection

    posterior_protection <- data %>%
        ungroup() %>%
        mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
        select(ind_id, plot_id, site, family, species, guild, size_class, predator, group, rugosity_mean, rugosity_se, biomass) %>%
        distinct() %>%
        expand_grid(protected = c(0, 1)) %>%
        add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_protection <- rbind(counterfactual_protection, posterior_protection)


    # counterfactual rugosity

    std_var <- seq(-2, 2, by = 0.2)

    posterior_rugosity <- data %>%
        ungroup() %>%
        mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
        select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, group, biomass) %>%
        distinct() %>%
        expand_grid(rugosity_mean = std_var, rugosity_se = 0.1) %>%
        add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))


    counterfactual_rugosity <- rbind(counterfactual_rugosity, posterior_rugosity)

    # counterfactual biomass
    std_var <- seq(-2, 2, by = 0.2)

    posterior_biomass <- data %>%
        ungroup() %>%
        mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
        select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, group, rugosity_mean, rugosity_se) %>%
        distinct() %>%
        expand_grid(biomass = std_var) %>%
        add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_biomass <- rbind(counterfactual_biomass, posterior_biomass)

    # counterfactual predator

    posterior_predator <- data %>%
        ungroup() %>%
        mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
        select(ind_id, plot_id, family, species, site, guild, size_class, protected, group, rugosity_mean, rugosity_se, biomass) %>%
        distinct() %>%
        expand_grid(predator = c(0, 1)) %>%
        add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_predator <- rbind(counterfactual_predator, posterior_predator)

    # counterfactual group

    posterior_group <- data %>%
        ungroup() %>%
        mutate(size_class = factor(size_class, levels = c(1, 2, 3, 4))) %>%
        select(ind_id, plot_id, site, family, species, guild, size_class, protected, predator, rugosity_mean, rugosity_se, biomass) %>%
        distinct() %>%
        expand_grid(group = c(0, 1)) %>%
        add_epred_draws(model, value = "prop", allow_new_levels = TRUE) %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_group <- rbind(counterfactual_group, posterior_group)
}

# Counterfactual effect of protection on time budget

effect_protection <- summarise_cf(posterior = counterfactual_protection, do = "protected", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_protection.docx"))

effect_protection %>%
    group_by(behaviour) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Behaviour", y = "Log odds ratio")

ggsave(here("figures", "behaviour-prob", "effect_plot.png"), plot = last_plot())

## Across guilds

effect_protection <- summarise_cf(posterior = counterfactual_protection, do = "protected", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_protection_guild.docx"))

p1 <- effect_protection %>%
    group_by(behaviour, guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Behaviour", y = "Log odds ratio", color = "Guild")

ggsave(here("figures", "behaviour-prob", "effect_plot_guild.png"), plot = p, height = 6, width = 8)

## Across guilds and size classes
effect_protection <- summarise_cf(posterior = counterfactual_protection, do = "protected", confound = c("behaviour", "guild", "size_class"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_protection_size.docx"))

effect_protection %>%
    group_by(behaviour, guild, size_class) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = size_class, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Size class", y = "Log odds ratio", color = "Guild") +
    facet_wrap(~behaviour)

ggsave(here("figures", "behaviour-prob", "effect_plot_guild_size.png"), plot = last_plot(), height = 6, width = 8)

# Effect of rugosity

p2 <- counterfactual_rugosity %>%
    group_by(behaviour, guild, rugosity_mean) %>%
    median_qi(prop, .width = c(0.5)) %>%
    mutate(rugosity_mean = moments_rugosity[2] * rugosity_mean + moments_rugosity[1]) %>%
    ggplot(aes(x = rugosity_mean, y = prop, col = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Rugosity", y = "Probability") +
    facet_grid(guild ~ behaviour) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-prob", "counterfactual_rugosity.png"), plot = p2, height = 8, width = 12)

effect_rugosity <- summarise_cf(posterior = counterfactual_rugosity, do = "rugosity_mean", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_rugosity.docx"))

# Effect of biomass

p3 <- counterfactual_biomass %>%
    group_by(behaviour, guild, biomass) %>%
    median_qi(prop, .width = c(0.5)) %>%
    mutate(biomass = (moments_biomass[2] * biomass + moments_biomass[1]) * 100) %>%
    ggplot(aes(x = biomass, y = prop, col = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Biomass Cover (%)", y = "Probability") +
    facet_grid(guild ~ behaviour) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-prob", "counterfactual_biomass.png"), plot = p3, height = 8, width = 12)

effect_biomass <- summarise_cf(posterior = counterfactual_biomass, do = "biomass", confound = c("behaviour", "guild"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_biomass.docx"))

# Effect of predator presence

effect_predator <- summarise_cf(posterior = counterfactual_predator, do = "predator", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_predator.docx"))

effect_predator %>%
    group_by(behaviour) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect, col = behaviour)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Behaviour", y = expression(Delta ~ " P(time)")) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-prob", "counterfactual_predator.png"), plot = last_plot(), height = 6, width = 8)

# Effect of grouping

effect_group <- summarise_cf(posterior = counterfactual_group, do = "group", confound = c("behaviour"), groups = c("ind_id", "guild", "behaviour", "size_class"), path = "behaviour-prob", file = paste0("counterfactual_group.docx"))

effect_group %>%
    group_by(behaviour) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect, col = behaviour)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Behaviour", y = expression(Delta ~ " P(time)")) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-prob", "counterfactual_group.png"), plot = last_plot(), height = 6, width = 8)
