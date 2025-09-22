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
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
        flextable() %>%
        set_header_labels(
            type = "Type",
            R2 = "R^2"
        ) %>%
        theme_box()

    print(table3)

    save_as_docx(table3, path = path)
}

# summaries of counterfactuals

summarise_cf <- function(posterior, do, confound = NULL, groups, path, file, dist = "proportion"){

    if (do  %in% c("rugosity_mean", "biomass")){
        posterior <- posterior %>%
            filter(!!sym(do) == 0 | !!sym(do) == 1)
    } 
    
    effect <- posterior %>%
        ungroup() %>%
        pivot_wider(names_from = !!sym(do), values_from = prop, id_cols = c(!!!syms(c(".draw", groups)))) %>%
        rename("Control" = `0`, "Treatment" = `1`) 

    if (dist == "Abundance") {
        effect <- effect %>%
            mutate(effect = log2(Treatment) - log2(Control))
    } else {
        effect <- effect %>%
            mutate(effect = log2(Treatment / (1 - Treatment)) - log2(Control / (1 - Control)))
    }

    P <- effect %>%
        group_by(!!!syms(confound)) %>%
        summarise(P = mean(effect > 0), .groups = "drop")
    
    effect_summary <- effect %>%
        group_by(!!!syms(confound)) %>%
        median_qi(effect, .width = c(0.5, 0.9)) %>%
        pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
        left_join(P) %>%
        select(!!!syms(confound), effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, P)

    table <- effect_summary %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
        flextable() %>%
        set_header_labels(
            behaviour = "Behaviour",
            guild = "Foraging Guild",
            effect = "Median",
            .lower_0.9 = "90% Lower CI",
            .lower_0.5 = "50% Lower CI",
            .upper_0.5 = "50% Upper CI",
            .upper_0.9 = "90% Upper CI",
            P = "P(>0)"
        ) %>%
        merge_v(j = confound) %>%
        theme_box() %>%
        set_table_properties(width = 0.5, layout = "fixed")

    print(table)

    save_as_docx(table, path = here("outputs", path, file))

    return(effect)
}
