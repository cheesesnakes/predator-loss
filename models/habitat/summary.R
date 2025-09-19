pacman::p_load(brms, here, tidybayes, ggplot2, marginaleffects, broom.mixed)

source(here("analysis", "standardise_fish.R"))

dir.create(here("outputs", "habitat"), showWarnings = FALSE)

# load models

model <- readRDS(here("outputs", "habitat", paste0("model.rds")))

# summarise model

hyperparameters <- model %>%
    spread_draws(phi) %>%
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

save_as_docx(table3, path = here("outputs", "habitat", paste0("_hyperparameters", ".docx")))

# Random effects

random_effects <- tidy(model, effects = "ran_vals") %>%
    select(level, estimate, conf.low, conf.high) %>%
    ungroup() %>%
    mutate(
        Site = str_split(level, "_", simplify = TRUE)[,1],
        Transect = str_split(level, "_", simplify = TRUE)[,2]
    ) %>%
    arrange(as.numeric(Transect)) %>%
    mutate(Response = rep(c("Biomass", "Sponge", "Substrate"), times = 24)) %>%
    select(c(Site, Transect, Response, estimate, conf.low, conf.high))


# make table for random intercepts

table1 <- random_effects %>%
    mutate(Site = site_levels[as.numeric(Site)],Site = str_to_title(Site)) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    arrange(as.numeric(Transect)) %>%
    flextable() %>%
    set_header_labels(
        Site = "Site",
        Transect = "Transect",
        Response = "Response",
        estimate = "Median",
        conf.low = "Lower CI",
        conf.high = "Upper CI"
    ) %>%
    merge_v(~Site) %>%
    set_table_properties(width = 0.5, layout = "fixed") %>%
    theme_box()

print(table1)

save_as_docx(table1, path = here("outputs", "habitat", paste0("_intercepts", ".docx")))

# Parameters

parameters <- tidy(model, effect = "fixed")%>%
select(-c(effect, component, std.error)) %>%
mutate(Response = str_split(term, "_", simplify = TRUE)[,1],
Variable = str_split(term, "_", simplify = TRUE)[,2],
Response = str_to_title(str_remove(Response, "mu")),
Variable = str_to_title(str_remove_all(Variable, "[()]"))
)%>%
arrange(Response) %>%
select(c(Response, Variable, estimate, conf.low, conf.high))

# make table

table2 <- parameters %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        estimate = "Median",
        conf.low = "Lower CI",
        conf.high = "Upper CI"
    ) %>%
    merge_v(~Response) %>%
    set_table_properties(width = 0.5, layout = "fixed") %>%
    theme_box()

print(table2)

save_as_docx(table2, path = here("outputs", "habitat", paste0("_parameters", ".docx")))
