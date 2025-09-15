pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_behaviour.R"))

dir.create(here("outputs", "behaviour-prob"), showWarnings = FALSE)

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

    model <- readRDS(here("outputs", "behaviour-prob", paste0("model_", b, ".rds")))

    # Random effects

    random_effects <- model %>%
        spread_draws(`r_site:plot_id`[l_s, ]) %>%
        summarise_draws() %>%
        ungroup() %>%
        mutate(
            Site = str_split(l_s, "_", simplify = TRUE)[, 1],
            Plot = str_split(l_s, "_", simplify = TRUE)[, 2],
        ) %>%
        select(-c(l_s, variable)) %>%
        select(Site, Plot, mean, median, sd, q5, q95)

    # make table for random intercepts

    table1 <- random_effects %>%
    select(-c(sd)) %>%
    mutate(Site = site_levels[as.numeric(Site)])%>%
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
        merge_v(~ Site) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        theme_box()
    
    print(table1)

    save_as_docx(table1, path = here("outputs", "behaviour-prob", paste0(b, "_intercepts", ".docx")))

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

    save_as_docx(table2, path = here("outputs", "behaviour-prob", paste0(b, "_parameters", ".docx")))
}
