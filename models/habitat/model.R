pacman::p_load(brms, here, tidybayes, ggplot2, marginaleffects, DirichletReg)

source(here("analysis", "standardise_fish.R"))

dir.create(here("outputs", "habitat"), showWarnings = FALSE)

# make data

data <- benthic %>%
    left_join(transects, by = c("transect_id" = "transect"))

# normalisation and transformation for 0, 1

Y <- matrix(c(data$coral, data$biomass, data$sponge, data$substrate), ncol = 4)

Y_star <- DR_data(Y, trafo = TRUE)

# final data
data <- data %>%
    ungroup() %>%
    mutate(
        coral = Y_star[, 1],
        biomass = Y_star[, 2],
        sponge = Y_star[, 3],
        substrate = Y_star[, 4]
    )

# define formula

formula <- bf(cbind(coral, biomass, sponge, substrate) ~ protection + aspect + (1 | site:transect_id), family = dirichlet())

# define priors

priors <- c(
    set_prior("normal(0,0.5)", class = "Intercept", dpar = "mubiomass"),
    set_prior("normal(0,0.5)", class = "Intercept", dpar = "musponge"),
    set_prior("normal(0,0.5)", class = "Intercept", dpar = "musubstrate"),
    set_prior("normal(0,0.5)", class = "b", dpar = "mubiomass"),
    set_prior("normal(0,0.5)", class = "b", dpar = "musponge"),
    set_prior("normal(0,0.5)", class = "b", dpar = "musubstrate"),
    set_prior("exponential(2)", class = "sd", dpar = "mubiomass"),
    set_prior("exponential(2)", class = "sd", dpar = "musponge"),
    set_prior("exponential(2)", class = "sd", dpar = "musubstrate"),
    set_prior("exponential(2)", class = "phi")
)

# fit model

model <- brm(formula, data = data, family = Dirichlet(), iter = 2000, chains = 4, cores = 4, prior = priors)

# summarise model

print(summary(model))

# save model

saveRDS(model, file = here("outputs", "habitat", "model.rds"))

stancode <- make_stancode(formula, data = data, family = dirichlet(), prior = priors)
writeLines(stancode, con = here("models", "habitat", "model.stan"))
