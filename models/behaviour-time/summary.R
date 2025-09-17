pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

dir.create(here("outputs", "behaviour-time"), showWarnings = FALSE)

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

    model <- readRDS(here("outputs", "behaviour-time", paste0("model_", b, ".rds")))

    # summarise model

    hyperparameters <- model %>%
        spread_draws(phi, cutone, cutzero) %>%
        summarise_draws() %>%
        mutate(type = "Hyperparameter")

    table3 <- hyperparameters %>%
        select(c(variable, mean, median, q5, q95)) %>%
        mutate(variable = str_to_title(variable)) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        flextable() %>%
        set_header_labels(
            variable = "Variable",
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI"
        ) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        theme_box()

    print(table3)

    save_as_docx(table3, path = here("outputs", "behaviour-time", paste0(b, "_hyperparameters", ".docx")))

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

    # make table for random intercepts

    table1 <- random_effects %>%
        select(-c(sd)) %>%
        mutate(Site = site_levels[as.numeric(Site)]) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        arrange(as.numeric(Plot)) %>%
        flextable() %>%
        set_header_labels(
            Site = "Site",
            Plot = "Plot",
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI"
        ) %>%
        merge_v(~Site) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        theme_box()

    print(table1)

    save_as_docx(table1, path = here("outputs", "behaviour-time", paste0(b, "_intercepts", ".docx")))

    # random effects across species

    random_effects_species <- model %>%
        spread_draws(`r_species`[l_s, ]) %>%
        summarise_draws() %>%
        ungroup() %>%
        mutate(Species = species_levels[as.numeric(l_s)]) %>%
        select(-c(l_s, variable)) %>%
        select(Species, mean, median, sd, q5, q95)

    # tables

    table4 <- random_effects_species %>%
        select(-c(sd)) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        flextable() %>%
        set_header_labels(
            Species = "Species",
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI"
        ) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        italic(j = 1) %>%
        theme_box()

    print(table4)

    save_as_docx(table4, path = here("outputs", "behaviour-time", paste0(b, "_species", ".docx")))

    # Parameters

    parameters <- model %>%
        spread_draws(b_Intercept, b_protected, b_rugosity, b_biomass, b_predator, b_group) %>%
        summarise_draws() %>%
        mutate(
            type = "Parameter",
            variable = str_replace(variable, "b_", "")
        )

    # make table

    table2 <- parameters %>%
        select(c(variable, mean, median, q5, q95)) %>%
        mutate(variable = str_to_title(variable)) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        flextable() %>%
        set_header_labels(
            variable = "Variable",
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI"
        ) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        theme_box()

    print(table2)

    save_as_docx(table2, path = here("outputs", "behaviour-time", paste0(b, "_parameters", ".docx")))

    # Make tables for conditional and marginal R2

    conditional_r2 <- data.frame(bayes_R2(model, re.form = NULL), row.names = NULL) %>%
        mutate(type = "Conditional")

    marginal_r2 <- data.frame(bayes_R2(model, re.form = NA)) %>%
        mutate(type = "Marginal")

    R2 <- rbind(conditional_r2, marginal_r2)

    table3 <- R2 %>%
        flextable() %>%
        set_header_labels(
            type = "Type",
            R2 = "R^2"
        ) %>%
        theme_box()

    print(table3)

    save_as_docx(table3, path = here("outputs", "behaviour-time", "R2.docx"))
}
