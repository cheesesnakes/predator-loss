pacman::p_load(here, bayesplot, ggplot2)

source(here("models", "functions.R"))

theme_set(
  theme_bw()+
  theme(
  legend.position = "top",
  legend.text = element_text(size = 12),)
)
# import data

run <- "validation"

figure_dir <- switch(run,
    validation = here("figures", "simulated_data"),
    standardised = here("figures")
)
# test
data <- load_data(run)

stan_data <- prep_stan(data$response, data$individuals, data$predictors, data$rugosity)

fit_result <- fit_model(stan_data, here("models", "model.stan"), run = run)

# summary of the model fit
# Retrieve all variable names from the fitted model
vars <- c("beta_risk", "beta_risk_bites", "beta_res", "beta_rug", "beta_predator") 
# Use grep to filter variable names matching the regex pattern
print(fit_result$summary(variables = vars))

# Extract draws for specific parameters
draws <- fit_result$draws(vars, format = "df")

# traceplot

mcmc_trace(draws, regex_pars = vars, facet_args = list(ncol = 2))+
scale_color_viridis_d()
ggsave(here(figure_dir, "traceplot.png"), width = 10, height = 8)

#draw_size <- fit_result$draws("beta_size", format = "df")
#
#mcmc_trace(draw_size, regex_pars = "beta_size", facet_args = list(ncol = 2))+
#scale_color_viridis_d()

# posterior predictive checks

# Extract posterior draws for D_pred and convert to matrix
ppc_D <- fit_result$draws("D_pred", format = "draws_matrix")


# Create matrix of observed durations
response <- data$response
data_D <- as.matrix(data.frame(
  foraging = response$foraging,
  vigilance = response$vigilance,
  movement = response$movement
))

random_predictions <- sample(1:nrow(ppc_D), 100)

# Posterior predictive check for foraging durations
ppc_dens_overlay(y = data_D[,1], yrep = ppc_D[random_predictions,1:500])+
scale_x_continuous(limits = c(0, 120))+
  labs(x = "Foraging duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_foraging.png"), width = 10, height = 8)

ppc_dens_overlay(y = data_D[,2], yrep = ppc_D[random_predictions,501:1000]) +
  scale_x_continuous(limits = c(0, 120))+
  labs(x = "Vigilance duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_vigilance.png"), width = 10, height = 8)

ppc_dens_overlay(y = data_D[,3], yrep = ppc_D[random_predictions,1001:1500])+
  scale_x_continuous(limits = c(0, 120))+
  labs(x = "Movement duration (s)", y = "Density")

ggsave(here(figure_dir, "ppc_movement.png"), width = 10, height = 8)

# Extract posterior predictive draws for bites
ppc_bites <- fit_result$draws("bites_pred", format = "draws_matrix")

# Posterior predictive check for bites
ppc_dens_overlay(yrep = ppc_bites[random_predictions,], y = response$bites)+
  scale_x_continuous(limits = c(0, 20))+
  labs(x = "Bites", y = "Density")

ggsave(here(figure_dir, "ppc_bites.png"), width = 10, height = 8)
