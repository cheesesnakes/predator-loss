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

behaviours <- c("foraging", "vigilance", "move")

for (b in behaviours) {
  # define model formula
  if (b == "foraging") {
    formula <- bf(foraging ~ protected + rugosity + biomass + predator + group + guild + size_class + (1 | site))
  } else if (b == "vigilance") {
    formula <- bf(vigilance ~ protected + rugosity + biomass + predator + group + guild + size_class + (1 | site))
  } else if (b == "move") {
    formula <- bf(move ~ protected + rugosity + biomass + predator + group + guild + size_class + (1 | site))
  }

  # priors

  priors <- c(set_prior("normal(0,0.5)", class = "Intercept"),
            set_prior("normal(0,0.5)", class = "b"),
            set_prior("exponential(1)", class = "sd"),
            set_prior("exponential(1)", class = "phi"))

  # run model
  model <- ordbetareg(
    formula = formula,
    manual_prior = priors,
    data = data,
    chains = 4,
    cores = 4,
    iter = 2000
  )

  # generate stan code from formula, not fitted model
  stan_code <- make_stancode(formula, data = data_ts)
  writeLines(stan_code, con = here("models", "behaviour", paste0("model_", b, ".stan")))

  # save model
  saveRDS(model, here("models", "behaviour", paste0("model_", b, ".rds")))

  # summarise model
  # p1 <- plot(model)
  # ggsave(here("models", "behaviour", paste0("model_plot_", b, ".png")), plot = p1)

  # posterior predictive check
  p2 <- pp_check(model)
  ggsave(here("models", "behaviour", paste0("pp_check_", b, ".png")), plot = p2)

  print(summary(model))
}
