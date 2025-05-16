pacman::p_load(here, ggplot2, dplyr, tidyr, flextable, bayesplot)

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
} else {
    fit <- readRDS(here("outputs/validation_fit.rds"))
}

# summary of the model fit
vars <- c("beta_risk",
          "beta_risk_bites", "beta_rug", "beta_res", "beta_predator")
print(fit$summary(variables = vars))

# Extract draws for specific parameters
draws <- fit$draws(vars, format = "df")

draws <- draws %>%
    select(-.chain, -.iteration, -.draw) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

# Maximum likelihood estimates

dist_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
    return(ux[which.max(tabulate(match(x, ux)))])
}

draws_mle <- draws %>%
    group_by(parameter) %>%
    summarise(mean = mean(value),
              median = median(value),
              mode = dist_mode(value)) %>%
    ungroup()

draws_mle$value <- draws_mle$mean

# True values
true_values <- data.frame(
    parameter = unique(draws$parameter),
    value = c(beta_risk, beta_risk_bites, beta_rugosity, beta_biomass, beta_predator)
)

# Plot the draws
draws %>%
    ggplot(aes(x = value)) +
    geom_density() +
    geom_vline(data = draws_mle, aes(xintercept = value), color = "red", linetype = "dashed") +
    geom_vline(data = true_values, aes(xintercept = value), color = "blue") +
    facet_wrap(~parameter, scales = "free")

ggsave(here("figures/model_validation.png"), width = 12, height = 8)
# Difference table

diff <- true_values %>%
    left_join(draws_mle, by = "parameter") %>%
    mutate(
        diff = value.x - value.y,
        abs_diff = abs(diff),
        rel_diff = (abs_diff * 100 /( value.x)) %>% round(2), 
    ) %>%
    select(parameter, value.x, value.y, diff, abs_diff, rel_diff)%>%
    flextable()%>%
    set_header_labels(
        parameter = "Parameter",
        value.x = "True Value",
        value.y = "MLE",
        diff = "Difference",
        abs_diff = "Absolute Difference",
        rel_diff = "Relative Difference (%)"
    ) %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("Model validation results") 

print(diff)
save_as_docx(diff, path = here("outputs/model_validation.docx"))


