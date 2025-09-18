pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

dir.create(here("outputs", "behaviour-time"), showWarnings = FALSE)
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

# make data

# Calculate proportion of time in each behavioural state

data_ts <- data %>%
  filter(codingduration > 45) %>%
  mutate(
    foraging = foraging / codingduration,
    vigilance = vigilance / codingduration,
    move = move / codingduration
  )

behaviours <- c("foraging", "vigilance", "move")

for (b in behaviours) {
  # define model formula
  if (b == "foraging") {
    formula <- bf(foraging ~ protected + me(rugosity_mean, rugosity_se) + biomass + predator + group + guild + size_class + (1 | site:plot_id) + (1 | species),
    phi ~ 1 + (1|species))
  } else if (b == "vigilance") {
    formula <- bf(
      vigilance ~ protected + me(rugosity_mean, rugosity_se) + biomass + predator + group + guild + size_class + (1 | site:plot_id) + (1 | species),
      phi ~ 1 + (1 | species)
    )
  } else if (b == "move") {
    formula <- bf(
      move ~ protected + me(rugosity_mean, rugosity_se) + biomass + predator + group + guild + size_class + (1 | site:plot_id) + (1 | species),
      phi ~ 1 + (1 | species)
    )
  }

  # priors

  priors <- c(
    set_prior("normal(0,0.5)", class = "Intercept"),
    set_prior("normal(0,0.5)", class = "b"),
    set_prior("exponential(2)", class = "sd")
  )

  # run model
  model <- ordbetareg(
    formula = formula,
    manual_prior = priors,
    data = data_ts,
    chains = 4,
    cores = 4,
    iter = 2000
  )


  # save model
  saveRDS(model, here("outputs", "behaviour-time", paste0("model_", b, ".rds")))

  # summarise model
  # p1 <- plot(model)
  # ggsave(here("figures", "behaviour-time", paste0("model_plot_", b, ".png")), plot = p1)

  print(summary(model))

  # posterior predictive check
  p2 <- pp_check(model)
  ggsave(here("figures", "behaviour-time", paste0("pp_check_", b, ".png")), plot = p2)
}
