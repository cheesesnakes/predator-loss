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

dir.create(here("outputs", "behaviour-time"), showWarnings = FALSE)
dir.create(here("figures", "behaviour-time"), showWarnings = FALSE)

# analysis

behaviours <- c("foraging", "vigilance", "move")
counterfactual_protection <- data.frame()
counterfactual_rugosity <- data.frame()
counterfactual_biomass <- data.frame()
counterfactual_predator <- data.frame()
counterfactual_group <- data.frame()

for (b in behaviours) {
    # load models

    model <- readRDS(here("outputs", "behaviour-time", paste0("model_", b, ".rds")))

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

ggsave(here("figures", "behaviour-time", "effect_plot.png"), plot = last_plot())

## Make table

table1 <- counterfactual_protection %>%
    filter(guild != "Piscivore", !is.na(guild)) %>%
    group_by(behaviour) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(behaviour, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        behaviour = "Behaviour",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed") 

print(table1)

save_as_docx(table1, path = here("outputs", "behaviour-time", paste0("effect_protection.docx")))

## Across guilds

counterfactual_protection %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = behaviour, y = effect, color = guild)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5)) +
    labs(x = "Behaviour", y = "Difference in proportion", color = "Guild")

ggsave(here("figures", "behaviour-time", "effect_plot_guild.png"), plot = last_plot(), height = 6, width = 8)

## Table 

table2 <- counterfactual_protection %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(behaviour, guild, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        guild = "Guild",
        behaviour = "Behaviour",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    merge_v(~behaviour) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table2)

save_as_docx(table2, path = here("outputs", "behaviour-time", paste0("effect_protection_guild.docx")))

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

ggsave(here("figures", "behaviour-time", "effect_plot_guild_size.png"), plot = last_plot(), height = 6, width = 8)

## Table

table3 <- counterfactual_protection %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, guild, size_class) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(behaviour, guild, size_class, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    mutate(size_class =size_levels[as.numeric(size_class)]) %>%
    flextable() %>%
    set_header_labels(
        guild = "Guild",
        size_class = "Size class",
        behaviour = "Behaviour",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    merge_v(~behaviour + guild) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "behaviour-time", paste0("effect_protection_guild_size.docx")))
# Effect of rugosity

counterfactual_rugosity %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, rugosity) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = rugosity, y = prop)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    labs(x = "Rugosity", y = "Proportion of time") +
    facet_wrap(~behaviour, ncol = 1) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-time", "counterfactual_rugosity.png"), plot = last_plot(), height = 12, width = 8)

# Effect of biomass

counterfactual_biomass %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, biomass) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = biomass, y = prop)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    labs(x = "Biomass Cover (%)", y = "Proportion of time") +
    facet_wrap(~behaviour, ncol = 1)

ggsave(here("figures", "behaviour-time", "counterfactual_biomass.png"), plot = last_plot(), height = 12, width = 8)

# Effect of predator presence

counterfactual_predator %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, predator) %>%
    median_qi(prop, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = as.factor(predator), y = prop)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Predator presence", y = "Proportion of time") +
    facet_wrap(~behaviour, ncol = 1) +
    scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-time", "counterfactual_predator.png"), plot = last_plot(), height = 12, width = 8)

## Table

table4 <- counterfactual_predator %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, predator) %>%
    median_qi(prop, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(behaviour, predator, prop, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        predator = "Predator",
        behaviour = "Behaviour",
        prop = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    merge_v(~behaviour) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table4)

save_as_docx(table4, path = here("outputs", "behaviour-time", paste0("prop_predator.docx")))

# Effect of grouping 

counterfactual_group %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, group) %>%
    median_qi(prop, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = as.factor(group), y = prop)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Grouping", y = "Proportion of time") +
    facet_wrap(~behaviour, ncol = 1) +
    scale_x_discrete(labels = c("0" = "Solitary", "1" = "Grouped")) +
    theme(legend.position = "none")

ggsave(here("figures", "behaviour-time", "counterfactual_group.png"), plot = last_plot(), height = 12, width = 8)

## Table

table5 <- counterfactual_group %>%
    filter(!is.na(guild), guild != "Piscivore") %>%
    group_by(behaviour, group) %>%
    median_qi(prop, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(behaviour, group, prop, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        group = "Grouping",
        behaviour = "Behaviour",
        prop = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    merge_v(~behaviour) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table5)

save_as_docx(table5, path = here("outputs", "behaviour-time", paste0("prop_group.docx")))