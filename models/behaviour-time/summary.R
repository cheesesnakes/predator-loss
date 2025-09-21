pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))
source(here("functions", "summaries.R"))

dir.create(here("outputs", "behaviour-time"), showWarnings = FALSE)

# analysis

behaviours <- c("foraging", "vigilance", "move")

for (b in behaviours) {
    print(str_glue("\nSummarising model for {b}"))
    print("====================================\n")

    # load models

    model <- readRDS(here("outputs", "behaviour-time", paste0("model_", b, ".rds")))

    # summarise model

    summarise_params(model = model, type = "Hyperparameter", vars = c("cutone", "cutzero"), path = here("outputs", "behaviour-time", paste0(b, "_hyperparameters", ".docx")))

    # Random effects

    random_effects <- model %>%
        spread_draws(`r_site:plot_id`[l_s, ]) %>%
        summarise_draws() %>%
        ungroup() %>%
        mutate(
            Site = str_split(l_s, "_", simplify = TRUE)[, 1],
            Plot = str_split(l_s, "_", simplify = TRUE)[, 2]
        ) %>%
        select(-c(l_s, variable)) %>%
        select(Site, Plot, mean, median, sd, q5, q95) %>%
        mutate(Site = site_levels[as.numeric(Site)]) %>%
        arrange(as.numeric(Plot))

    # make table for random intercepts

    table_re(df = random_effects, var_names = "Site", path = here("outputs", "behaviour-time", paste0(b, "_intercepts", ".docx")))

    # random effects across species

    random_effects_species <- model %>%
        spread_draws(`r_family:species`[l_s, ]) %>%
        summarise_draws() %>%
        ungroup() %>%
        mutate(
            Family = str_split(l_s, "_", simplify = TRUE)[, 1],
            Species = str_split(l_s, "_", simplify = TRUE)[, 2]
        ) %>%
        mutate(
            Species = species_levels[as.numeric(Species)],
            Family = family_levels[as.numeric(Family)]
        ) %>%
        select(-c(l_s, variable)) %>%
        select(Family, Species, mean, median, sd, q5, q95) %>%
        arrange(Family, Species)

    # tables

    table_re(df = random_effects_species, var_names = "Family", path = here("outputs", "behaviour-time", paste0(b, "_family", ".docx")))

    # random effects across species (phi)

    random_effects_species_phi <- model %>%
        spread_draws(`r_family:species__phi`[l_s, ]) %>%
        summarise_draws() %>%
        ungroup() %>%
        mutate(
            Family = str_split(l_s, "_", simplify = TRUE)[, 1],
            Species = str_split(l_s, "_", simplify = TRUE)[, 2]
        ) %>%
        mutate(
            Species = species_levels[as.numeric(Species)],
            Family = family_levels[as.numeric(Family)]
        ) %>%
        select(-c(l_s, variable)) %>%
        select(Family, Species, mean, median, sd, q5, q95) %>%
        arrange(Family, Species)

    # tables

    table_re(df = random_effects_species_phi, var_names = "Family", path = here("outputs", "behaviour-time", paste0(b, "_family_phi", ".docx")))

    # Parameters

    summarise_params(model = model, type = "Parameter", vars = c("b_Intercept", "b_protected", "b_rugosity_mean", "b_biomass", "b_predator", "b_group", "b_rugosity_mean:guildInvertivore", "b_biomass:guildInvertivore", "b_protected:guildInvertivore"), path = here("outputs", "behaviour-prob", paste0(b, "_parameters", ".docx")))

    # Make tables for conditional and marginal R2

    summarise_r2(model = model, path = here("outputs", "behaviour-time", paste0(b, "_R2", ".docx")))
}
