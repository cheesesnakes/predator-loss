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

for (b in behaviours) {
    print(str_glue("\nSummarising model for {b}"))
    print("====================================\n")

    # load models

    model <- readRDS(here("models", "behaviour", paste0("model_", b, ".rds")))

    # summarise model

    hyperparameters <- model %>%
        spread_draws(phi, cutone, cutzero) %>%
        summarise_draws() %>%
        mutate(type = "Hyperparameter")

    print(hyperparameters)

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
        select(Site, Plot, mean, median, sd, q5, q95)

    print(random_effects)

    # Parameters

    parameters <- model %>%
        spread_draws(b_Intercept, b_protected, b_rugosity, b_biomass, b_predator, b_group) %>%
        summarise_draws() %>%
        mutate(
            type = "Parameter",
            variable = str_replace(variable, "b_", "")
        )

    print(parameters)
}
