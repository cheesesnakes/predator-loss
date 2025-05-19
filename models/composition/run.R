pacman::p_load(here, bayesplot, ggplot2)

source(here("models", "composition/functions.R"))

theme_set(
  theme_bw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
    )
)
# import data

run <- "validation"

figure_dir <- switch(run,
  validation = here("figures", "simulated_data/composition"),
  standardised = here("figures/composition")
)

dir.create(figure_dir, recursive = TRUE)

# test
data <- load_data(run)

stan_data <- prep_stan(
  fish = data$fish,
  transects = data$transects,
  species = data$species
)

fit_result <- fit_model(stan_data, here("models", "composition/model.stan"), run = run)

# summary of the model fit
# Retrieve all variable names from the fitted model
vars <- c("gamma", "delta")
# Use grep to filter variable names matching the regex pattern
print(fit_result$summary(variables = vars))

# Extract draws for specific parameters
draws <- fit_result$draws(vars, format = "df")

# traceplot

mcmc_trace(draws, regex_pars = vars, facet_args = list(ncol = 2)) +
  scale_color_viridis_d()
ggsave(here(figure_dir, "traceplot.png"), width = 10, height = 8)

# Extract posterior draws for pred_count and convert to matrix
ppc <- fit_result$draws("pred_count", format = "draws_matrix")

ppc
# Create matrix of observed durations
response <- data$fish
counts <- response$abundance

# Posterior predictive check for foraging durations
ppc_dens_overlay(y = counts, yrep = ppc) +
  scale_x_continuous(limits = c(0, 120)) +
  labs(x = "Abunddance", y = "Density")

ggsave(here(figure_dir, "ppc_abundance.png"), width = 10, height = 8)
