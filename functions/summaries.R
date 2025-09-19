pacman::p_load(here, dplyr, tidyr)

# Function for summarising parameters and hyperparameters

summarise_params <- function(model, type, vars, path) {
    params <- model %>%
        spread_draws(!!!syms(vars)) %>%
        summarise_draws() %>%
        mutate(type = type)

    posterior_probs <- model %>%
        spread_draws(!!!syms(vars)) %>%
        group_by(.draw) %>%
        pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
        group_by(variable) %>%
        summarise(P = mean(value > 0), .groups = "drop")

    params <- left_join(params, posterior_probs, by = "variable")

    params <- params %>%
            mutate(variable = case_when(
                str_detect(variable, "b_") ~ str_remove(variable, "b_"),
                str_detect(variable, "bsp_merugosity_meanrugosity_se") ~ "Rugosity",
                TRUE ~ variable
            ))
    
    table <- params %>%
        select(c(variable, mean, median, q5, q95, P)) %>%
        mutate(variable = str_to_title(variable)) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        flextable() %>%
        set_header_labels(
            variable = "Variable",
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI",
            P = "P(>0)"
        ) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        theme_box()

    print(table)

    save_as_docx(table, path = path)
}

# Function to summarise random effects

var_names <- "Site"

table_re <- function(df, var_names, path) {
    # tables

    table <- df %>%
        select(-c(sd)) %>%
        mutate_if(is.numeric, ~ round(., 3)) %>%
        flextable() %>%
        set_header_labels(
            mean = "Mean",
            median = "Median",
            sd = "SD",
            q5 = "Lower CI",
            q95 = "Upper CI"
        ) %>%
        set_table_properties(width = 0.5, layout = "fixed") %>%
        merge_v(j = var_names) %>%
        italic(j = 2) %>%
        theme_box()

    print(table)

    save_as_docx(table, path = path)
}

# Function for extracting and saving R2

summarise_r2 <- function(model, path) {
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

    save_as_docx(table3, path = path)
}
