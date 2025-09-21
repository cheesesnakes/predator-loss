pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

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

# Convert observations to binary outcomes

data <- scans %>%
  left_join(fish, by = "ind_id") %>%
  left_join(predictors, by = "plot_id") %>%
  left_join(guild, by = "species") %>%
  select(-c(sponge)) %>%
  mutate(
    predator = predator - 1,
    group = group - 1,
    size_class = as.factor(size_class)
  ) %>%
  filter(guild != "Piscivore") %>%
  rename(move = moving)

print(paste0("N = ", length(unique(data$ind_id)), " individuals"))

behaviours <- c("foraging", "vigilance", "move")

for (b in behaviours) {
  # define model formula
  if (b == "foraging") {
    formula <- bf(foraging ~ protected + rugosity_mean + biomass + predator + group + guild + rugosity_mean:guild + biomass:guild + protected:guild + size_class + (1 | site:plot_id) + (1 | family:species))
  } else if (b == "vigilance") {
    formula <- bf(vigilance ~ protected + rugosity_mean + biomass + predator + group + guild + rugosity_mean:guild + biomass:guild + protected:guild + size_class + (1 | site:plot_id) + (1 | family:species))
  } else if (b == "move") {
    formula <- bf(move ~ protected + rugosity_mean + biomass + predator + group + guild + rugosity_mean:guild + biomass:guild + protected:guild + size_class + (1 | site:plot_id) + (1 | family:species))
  }

  # priors

  priors <- c(
    set_prior("normal(0,0.5)", class = "Intercept"),
    set_prior("normal(0,0.5)", class = "b"),
    set_prior("exponential(2)", class = "sd")
  )

  # run model
  model <- brm(
    formula = formula,
    prior = priors,
    data = data,
    chains = 4,
    cores = 4,
    iter = 2000,
    family = bernoulli(link = "logit")
  )

  # generate stan code from formula, not fitted model
  stan_code <- make_stancode(formula, data = data, prior = priors, family = bernoulli(link = "logit"))
  writeLines(stan_code, con = here("models", "behaviour-prob", paste0("model_", b, ".stan")))

  # save model
  saveRDS(model, here("outputs", "behaviour-prob", paste0("model_", b, ".rds")))

  # summarise model
  # p1 <- plot(model)
  # ggsave(here("models", "behaviour-prob", paste0("model_plot_", b, ".png")), plot = p1)

  # posterior predictive check
  p2 <- pp_check(model)
  ggsave(here("figures", "behaviour-prob", paste0("pp_check_", b, ".png")), plot = p2)

  print(summary(model))
}
