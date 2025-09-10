pacman::p_load(here, bayesplot, ggplot2)

source(here("models/behaviour", "functions.R"))

theme_set(
  theme_bw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
    )
)

# import data

run <- "standardised"

figure_dir <- switch(run,
  validation = here("figures", "simulated_data"),
  standardised = here("figures", "behaviour")
)

# test
data <- load_data(run)

stan_data <- prep_stan(data$response, data$individuals, data$predictors, data$rugosity, data$guild)

fit_result <- fit_model(stan_data, here("models/behaviour", "model.stan"), run = run)

print("Completed model fitting.")

# summary of the model fit
# Retrieve all variable names from the fitted model
vars <- c("beta_pi", "beta_D", "beta_bites")

# Use grep to filter variable names matching the regex pattern
print("Extracting summary for parameters")
print(fit_result$summary(variables = vars))

# Extract draws for specific parameters
draws <- fit_result$draws(vars, format = "df")

# traceplot

mcmc_trace(draws, regex_pars = vars, facet_args = list(ncol = 2)) +
  scale_color_viridis_d()
ggsave(here(figure_dir, "traceplot.png"), width = 10, height = 8)

# draw_size <- fit_result$draws("beta_size", format = "df")
#
# mcmc_trace(draw_size, regex_pars = "beta_size", facet_args = list(ncol = 2))+
# scale_color_viridis_d()

# posterior predictive checks

# Extract posterior draws for D_pred and convert to matrix
ppc_D <- fit_result$draws("D_pred", format = "draws_matrix")

ppc_D <- ppc_D[complete.cases(ppc_D), ] # Remove rows with NA values

# Create matrix of observed durations
response <- data$response
data_D <- as.matrix(data.frame(
  foraging = response$foraging,
  vigilance = response$vigilance,
  movement = response$movement
))

random_predictions <- sample(seq_len(nrow(ppc_D)), 100)

# Posterior predictive check for foraging durations
ppc_dens_overlay(y = data_D[, 1], yrep = ppc_D[random_predictions, seq_len(nrow(response))]) +
  scale_x_continuous(limits = c(0, 120)) +
  labs(x = "Foraging duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_foraging.png"), width = 10, height = 8)

ppc_dens_overlay(y = data_D[, 2], yrep = ppc_D[random_predictions, nrow(response) + seq_len(nrow(response))]) +
  scale_x_continuous(limits = c(0, 120)) +
  labs(x = "Vigilance duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_vigilance.png"), width = 10, height = 8)

ppc_dens_overlay(y = data_D[, 3], yrep = ppc_D[random_predictions, (nrow(response) * 2) + seq_len(nrow(response))]) +
  scale_x_continuous(limits = c(0, 120)) +
  labs(x = "Movement duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_movement.png"), width = 10, height = 8)

# Extract posterior predictive draws for bites
ppc_bites <- fit_result$draws("bites_pred", format = "draws_matrix")

ppc_bites <- ppc_bites[complete.cases(ppc_bites), ] # Remove rows with NA values
# Posterior predictive check for bites
ppc_dens_overlay(yrep = ppc_bites[random_predictions, ], y = response$bites) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(x = "Bites", y = "Density")

ggsave(here(figure_dir, "ppc_bites.png"), width = 10, height = 8)
