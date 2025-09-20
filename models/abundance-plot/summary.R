pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_behaviour.R"))
source(here("functions", "summaries.R"))

dir.create(here("outputs", "abundance-plot"), showWarnings = FALSE)

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

# load models

model <- readRDS(here("outputs", "abundance-plot", paste0("model.rds")))

# Random effects

random_effects <- model %>%
    spread_draws(`r_site`[Site, ]) %>%
    summarise_draws() %>%
    ungroup() %>%
    select(-c(variable)) %>%
    select(Site, mean, median, sd, q5, q95)

# make table for random intercepts

table_re(df = random_effects, var_names = "Site", path = here("outputs", "abundance-plot", "intercepts.docx"))

# Parameters

summarise_params(model, type = "Fixed", vars = c("b_Intercept", "b_protected", "b_damsel", "bsp_merugosity_meanrugosity_se", "b_biomass", "b_guildInvertivore", "bsp_merugosity_meanrugosity_se:guildInvertivore", "b_guildInvertivore:biomass", "shape"), path = here("outputs", "abundance-plot", "fixed_effects.docx"))

# Make tables for conditional and marginal R2

summarise_r2(model, path = here("outputs", "abundance-plot", "R2.docx"))