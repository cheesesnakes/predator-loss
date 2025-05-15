pacman::p_load(here, ggplot2, dplyr, tidyr, flextable)

theme_set(theme_bw() + theme(text = element_text(size = 20), legend.position = "top"))

model <- "behaviour"

if (model == "behaviour") {
    # Load the behaviour model
    source(here("simulation/generate_behaviour.R"))
} else {
    # Load the health model
    source(here("simulation/generate_community.R"))
}

# Load the model data

if (model == "behaviour") {
    fit <- readRDS(here("outputs/validation_fit.rds"))
    fit_mle <- readRDS(here("outputs/validation_fit_mle.rds"))
} else {
    fit <- readRDS(here("outputs/validation_fit.rds"))
    fit_mle <- readRDS(here("outputs/validation_fit_mle.rds"))
}

# summary of the model fit
vars <- c("beta_risk", "beta_risk_bites", "beta_res", "beta_rug", "beta_predator")

fit$summary(variables = vars)

# Extract draws for specific parameters
draws <- fit$draws(vars, format = "df")

draws <- draws %>%
    select(-.chain, -.iteration, -.draw) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

# Maximum likelihood estimates

draws_mle <- fit_mle$draws(variables = vars, format = "df")

draw_mle <- draws_mle %>%
    select(-.chain, -.iteration, -.draw) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

# True values
true_values <- data.frame(
    parameter = unique(draws$parameter),
    value = c(beta_risk, beta_risk_bites, beta_biomass, beta_rugosity, beta_predator)
)

# Plot the draws
draws %>%
    ggplot(aes(x = value)) +
    geom_density() +
    geom_vline(data = draw_mle, aes(xintercept = value), color = "red", linetype = "dashed") +
    geom_vline(data = true_values, aes(xintercept = value), color = "blue") +
    facet_wrap(~parameter, scales = "free")

# Difference table

true_values %>%
    left_join(draw_mle, by = "parameter") %>%
    mutate(
        diff = value.x - value.y,
        abs_diff = abs(diff),
        rel_diff = (abs_diff * 100 /( value.x+1e-6)) %>% round(2), 
    ) %>%
    select(parameter, value.x, value.y, diff, abs_diff, rel_diff)%>%
    flextable() %>%
    set_header_labels(
        parameter = "Parameter",
        value.x = "True Value",
        value.y = "MLE Value",
        diff = "Difference",
        abs_diff = "Absolute Difference",
        rel_diff = "Relative Difference (%)"
    ) %>%
    set_table_properties(width = 0.5, layout = "autofit") %>%
    set_caption("Difference between true values and MLE estimates") %>%
    set_table_properties(
        width = 0.5,
        layout = "autofit"
    )
