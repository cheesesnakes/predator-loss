pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable, patchwork)

source(here("analysis", "standardise_behaviour.R"))
source(here("functions", "summaries.R"))

# make data frame

data <- plot_abundance %>%
    left_join(guild, by = "species") %>%
    filter(!is.na(guild)) %>%
    group_by(plot_id, guild) %>%
    summarise(abundance = sum(abundance)) %>%
    left_join(predictors, by = "plot_id") %>%
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
    select(plot_id, site, guild, damsel, rugosity_mean, biomass) %>%
    distinct() %>%
    expand_grid(protected = c(0, 1)) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

# Calculate and summarise effect

effect <- summarise_cf(posterior = posterior_protection, do = "protected", dist = "Abundance", confound = "guild", groups = c("plot_id", "guild"), path = "abundance-plot", file = paste0("counterfactual_protection.docx"))

# Counterfactual effect of protection on number individuals in plot

effect %>%
    group_by(guild) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = guild, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
    labs(x = "Foraging Guild", y = "Log response ratio")

ggsave(here("figures", "abundance-plot", "effect_plot.png"), plot = last_plot())


# counterfactual rugosity

std_var <- seq(-2, 2, by = 0.2)

posterior_rugosity <- data %>%
    ungroup() %>%
    select(plot_id, site, guild, protected, biomass, damsel) %>%
    distinct() %>%
    expand_grid(rugosity_mean = std_var, rugosity_se = 0.1) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

# Effect of rugosity

p1 <- posterior_rugosity %>%
    mutate(rugosity = rugosity_mean * moments_rugosity[2] + moments_rugosity[1]) %>%
    group_by(guild, rugosity) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = rugosity, y = prop, color = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Rugosity", y = "Count", col = "Guild", fill = "Guild") +
    facet_wrap(~guild, ncol = 3) +
    theme(legend.position = "none")


ggsave(here("figures", "abundance-plot", "counterfactual_rugosity.png"), plot = p1, height = 6, width = 14)

# Standard effect of rugosity

effect <- summarise_cf(posterior = posterior_rugosity, do = "rugosity_mean", dist = "Abundance", confound = "guild", groups = c("plot_id", "guild"), path = "abundance-plot", file = paste0("counterfactual_rugosity.docx"))

# counterfactual biomass
std_var <- seq(-2, 2, by = 0.2)

posterior_biomass <- data %>%
    ungroup() %>%
    select(plot_id, site, guild, protected, rugosity_mean, rugosity_se) %>%
    distinct() %>%
    expand_grid(biomass = std_var, damsel = c(0, 1)) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

# plot

p2 <- posterior_biomass %>%
    mutate(biomass = (biomass * moments_biomass[2] + moments_biomass[1]) * 100) %>%
    group_by(guild, biomass) %>%
    median_qi(prop, .width = c(0.5)) %>%
    ggplot(aes(x = biomass, y = prop, color = guild, fill = guild)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Biomass (%)", y = "Count", col = "Guild", fill = "Guild") +
    facet_wrap(~guild, ncol = 3) +
    theme(legend.position = "none")

ggsave(here("figures", "abundance-plot", "counterfactual_biomass.png"), plot = p2, height = 6, width = 14)

# effect of biomass

effect <- summarise_cf(posterior = posterior_biomass, do = "biomass", confound = "guild", dist = "Abundance", groups = c("plot_id", "guild", "damsel"), path = "abundance-plot", file = paste0("counterfactual_biomass.docx"))

# Effect of damsels

effect <- summarise_cf(posterior = posterior_biomass, do = "damsel", confound = "guild", dist = "Abundance", groups = c("plot_id", "guild", "biomass"), path = "abundance-plot", file = paste0("counterfactual_damsel.docx"))

# Make plots

p <- p1 / p2

ggsave(here("figures", "figure2.png"), plot = p, height = 12, width = 14)
