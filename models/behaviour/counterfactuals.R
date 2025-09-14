pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

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

# analysis

behaviours <- c("foraging", "vigilance", "move")
counterfactual_protection <- data.frame()
counterfactual_rugosity <- data.frame()
counterfactual_biomass <- data.frame()
counterfactual_predator <- data.frame()
counterfactual_group <- data.frame()

for (b in behaviours) {
    # load models

    model <- readRDS(here("models", "behaviour", paste0("model_", b, ".rds")))

    # plot posterior proportions by protection

    posterior_protection <- data %>%
        ungroup() %>%
        select(ind_id, plot_id, site, guild, size_class, predator, group) %>%
        distinct() %>%
        expand_grid(protected = c(0, 1), biomass = mean(data$biomass), rugosity = mean(data$rugosity)) %>%
        add_epred_draws(model, value = "prop") %>%
        ungroup()


    effect <- posterior_protection %>%
        pivot_wider(names_from = protected, values_from = prop, id_cols = c(ind_id, .draw, guild, size_class)) %>%
        rename("Unprotected" = `0`, "Protected" = `1`) %>%
        mutate(
            effect = Protected - Unprotected,
            behaviour = str_to_title(b)
        )

    counterfactual_protection <- rbind(counterfactual_protection, effect)


    # counterfactual rugosity

    std_var <- seq(-2, 2, by = 0.2)

    rugosity_std <- std_var * sd(data$rugosity) + mean(data$rugosity)

    posterior_rugosity <- data %>%
        ungroup() %>%
        select(ind_id, plot_id, site, guild, size_class, protected, predator, group) %>%
        distinct() %>%
        expand_grid(rugosity = rugosity_std, biomass = mean(data$biomass)) %>%
        add_epred_draws(model, value = "prop") %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))


    counterfactual_rugosity <- rbind(counterfactual_rugosity, posterior_rugosity)

    # counterfactual biomass
    std_var <- seq(-2, 2, by = 0.2)

    biomass_std <- std_var * sd(data$biomass) + mean(data$biomass)

    posterior_biomass <- data %>%
        ungroup() %>%
        select(ind_id, plot_id, site, guild, size_class, protected, predator, group) %>%
        distinct() %>%
        expand_grid(biomass = biomass_std, rugosity = mean(data$rugosity)) %>%
        add_epred_draws(model, value = "prop") %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_biomass <- rbind(counterfactual_biomass, posterior_biomass)

    # counterfactual predator

    posterior_predator <- data %>%
        ungroup() %>%
        select(ind_id, plot_id, site, guild, size_class, protected, group) %>%
        distinct() %>%
        expand_grid(predator = c(0, 1), biomass = mean(data$biomass), rugosity = mean(data$rugosity)) %>%
        add_epred_draws(model, value = "prop") %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_predator <- rbind(counterfactual_predator, posterior_predator)

    # counterfactual group

    posterior_group <- data %>%
        ungroup() %>%
        select(ind_id, plot_id, site, guild, size_class, protected, predator) %>%
        distinct() %>%
        expand_grid(group = c(0, 1), biomass = mean(data$biomass), rugosity = mean(data$rugosity)) %>%
        add_epred_draws(model, value = "prop") %>%
        ungroup() %>%
        mutate(behaviour = str_to_title(b))

    counterfactual_group <- rbind(counterfactual_group, posterior_group)
}

# Counterfactual effect of protection on time budget

counterfactual_protection %>%
    filter(guild != "Piscivore", !is.na(guild)) %>%
    group_by(behaviour) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Behaviour", y = "Difference in proportion")

ggsave(here("models", "behaviour", "effect_plot.png"), plot = last_plot())

## Across guilds

counterfactual_protection %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Behaviour", y = "Difference in proportion", color = "Guild")

ggsave(here("models", "behaviour", "effect_plot_guild.png"), plot = last_plot(), height = 6, width = 8)

## Across guilds and size classes

counterfactual_protection %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, guild, size_class) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = size_class, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Size class", y = "Difference in proportion", color = "Guild") +
    facet_wrap(~behaviour)

ggsave(here("models", "behaviour", "effect_plot_guild_size.png"), plot = last_plot(), height = 6, width = 8)
