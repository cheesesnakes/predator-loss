pacman::p_load(here, bayesplot, ggplot2, dplyr, tidyr)

source(here("models/environment/coral", "functions.R"))

theme_set(
  theme_bw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
    )
)

# import data

figure_dir <- here("figures", "environment")

# test
data <- load_data()

stan_data <- prep_stan(data$benthic, data$transects)

fit_result <- fit_model(stan_data, here("models/environment/coral", "model.stan"))

print("Completed model fitting.")

# summary of the model fit
# Retrieve all variable names from the fitted model
vars <- c("beta")

# Use grep to filter variable names matching the regex pattern
print("Extracting summary for parameters")
print(fit_result$summary(variables = vars))

# Extract draws for specific parameters
draws <- fit_result$draws(vars, format = "df")

# traceplot

mcmc_trace(draws, regex_pars = vars, facet_args = list(ncol = 1)) +
  scale_color_viridis_d()
ggsave(here(figure_dir, "traceplot.png"), width = 10, height = 8)

# Extract posterior draws for D_pred and convert to matrix
ppc <- fit_result$draws("coral_pred", format = "draws_matrix")

ppc <- ppc[complete.cases(ppc), ] # Remove rows with NA values

# Create matrix of observed durations
raw <- data$benthic
coral <- raw$coral

# Posterior predictive check for foraging durations
ppc_dens_overlay(y = coral, yrep = ppc) +
  labs(x = "Coral cover", y = "Density")

ggsave(here(figure_dir, "ppc.png"), width = 10, height = 8)
