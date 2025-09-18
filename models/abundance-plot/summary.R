pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_behaviour.R"))

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

table1 <- random_effects %>%
select(-c(sd)) %>%
mutate(Site = site_levels[as.numeric(Site)])%>%
mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        Site = "Site",
        mean = "Mean",
        median = "Median",
        sd = "SD",
        q5 = "Lower CI",
        q95 = "Upper CI"
    ) %>%
    set_table_properties(width = 0.5, layout = "fixed") %>%
    theme_box()

print(table1)

save_as_docx(table1, path = here("outputs", "abundance-plot", "intercepts.docx"))

# Parameters

get_variables(model)

parameters <- model %>%
    spread_draws(b_Intercept, b_protected, bsp_merugosity_meanrugosity_se, b_biomass, b_guildInvertivore, `bsp_merugosity_meanrugosity_se:guildInvertivore`, `b_guildInvertivore:biomass`, shape) %>%
    summarise_draws() %>%
    mutate(
        type = "Parameter",
        variable = str_replace(variable, "b_", "")
    )%>%
    mutate(variable = str_replace(variable, "bsp_merugosity_meanrugosity_se", "Rugosity"),
    variable = str_replace(variable, "guildInvertivore", "Guild (Invertivore)"))

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

save_as_docx(table2, path = here("outputs", "abundance-plot", paste0("parameters", ".docx")))

# Make tables for conditional and marginal R2

conditional_r2 <- data.frame(bayes_R2(model, re.form = NULL), row.names = NULL)%>%
mutate(type = "Conditional")

marginal_r2 <- data.frame(bayes_R2(model, re.form = NA))%>%
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

save_as_docx(table3, path = here("outputs", "abundance-plot", "R2.docx"))
