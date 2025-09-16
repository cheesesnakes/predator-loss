pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_behaviour.R"))

# make data frame

data <- plot_abundance %>%
  left_join(guild, by = "species")%>%
  filter(guild != "Piscivore") %>%
  group_by(plot_id, guild) %>%
  summarise(abundance = sum(abundance)) %>%
  left_join(predictors, by = "plot_id")%>%
  ungroup() 

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

dir.create(here("outputs", "abundance-plot"), showWarnings = FALSE)
dir.create(here("figures", "abundance-plot"), showWarnings = FALSE)

# analysis

# load models

model <- readRDS(here("outputs", "abundance-plot", paste0("model.rds")))

# plot posterior proportions by protection

posterior_protection <- data %>%
    ungroup() %>%
    select(plot_id, site, guild) %>%
    distinct() %>%
    expand_grid(protected = c(0, 1), biomass = 0, rugosity = 0) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()


effect <- posterior_protection %>%
    pivot_wider(names_from = protected, values_from = prop, id_cols = c(plot_id, .draw, guild)) %>%
    rename("Unprotected" = `0`, "Protected" = `1`) %>%
    mutate(
        effect = Protected - Unprotected
    )


# Counterfactual effect of protection on number individuals in plot

effect %>%
    group_by(guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = guild, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Foraging Guild", y = expression(Delta ~ "Expected Count"))

ggsave(here("figures", "abundance-plot", "effect_plot.png"), plot = last_plot())


# counterfactual rugosity

std_var <- seq(-2, 2, by = 0.2)

posterior_rugosity <- data %>%
    ungroup() %>%
    select(plot_id, site, guild, protected, biomass) %>%
    distinct() %>%
    expand_grid(rugosity = std_var) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

# Effect of rugosity

posterior_rugosity %>%
    mutate(rugosity = rugosity * moments_rugosity[2] + moments_rugosity[1]) %>%
    group_by(guild, rugosity) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = rugosity, y = prop, color = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Rugosity", y = "Count", col = "Guild", fill = "Guild") 

ggsave(here("figures", "abundance-plot", "counterfactual_rugosity.png"), plot = last_plot(), height = 6, width = 8)

# Standard effect of rugosity

effect_rugosity <- posterior_rugosity %>%
    filter(rugosity == 0 | rugosity == 1) %>%
    pivot_wider(names_from = rugosity, values_from = prop, id_cols = c(plot_id, .draw, guild)) %>%
    rename("Mean Rugosity" = `0`, "High Rugosity" = `1`) %>%
    mutate(
        effect = `High Rugosity` - `Mean Rugosity`
    ) %>%
    group_by(guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(guild, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9)

table6 <- effect_rugosity %>%
    flextable() %>%
    set_header_labels(
        guild = "Foraging Guild",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table6)

save_as_docx(table6, path = here("outputs", "abundance-plot", paste0("effect_rugosity.docx")))

# counterfactual biomass
std_var <- seq(-2, 2, by = 0.2)

posterior_biomass <- data %>%
    ungroup() %>%
    select(plot_id, site, guild, protected, rugosity) %>%
    distinct() %>%
    expand_grid(biomass = std_var) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

# plot

posterior_biomass %>%
    mutate(biomass = (biomass * moments_biomass[2] + moments_biomass[1])*100) %>%
    group_by(guild, biomass) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = biomass, y = prop, color = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Biomass (%)", y = "Count", col = "Guild", fill = "Guild")

ggsave(here("figures", "abundance-plot", "counterfactual_biomass.png"), plot = last_plot(), height = 6, width = 8)

# effect of biomass

effect_biomass <- posterior_biomass %>%
    filter(biomass == 0 | biomass == 1) %>%
    pivot_wider(names_from = biomass, values_from = prop, id_cols = c(plot_id, .draw, guild)) %>%
    rename("Mean Biomass" = `0`, "High Biomass" = `1`) %>%
    mutate(
        effect = `High Biomass` - `Mean Biomass`
    ) %>%
    group_by(guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(guild, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9)

table5 <- effect_biomass %>%
    flextable() %>%
    set_header_labels(
        guild = "Foraging Guild",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table5)

save_as_docx(table5, path = here("outputs", "abundance-plot", paste0("effect_biomass.docx")))