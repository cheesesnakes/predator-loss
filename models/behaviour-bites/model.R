pacman::p_load(brms, here, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

dir.create(here("outputs", "behaviour-bites"), showWarnings = FALSE)
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

data <- data_ts %>%
  filter(foraging > 0, bite_rate < 20) %>%
  ungroup()

print(paste0("N = ", length(unique(data$ind_id)), " individuals"))

formula <- bf(bite_rate ~ protected + rugosity_mean + biomass + predator + group + guild + rugosity_mean:guild + biomass:guild + protected:guild + size_class + (1 | site:plot_id) + (1 | family:species))

# priors

priors <- c(
  set_prior("normal(0,1)", class = "Intercept"),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("exponential(2)", class = "sd"),
  set_prior("exponential(1)", class = "shape")
)

# run model
model <- brm(
  formula = formula,
  prior = priors,
  data = data,
  chains = 4,
  cores = 4,
  iter = 2000,
  family = Gamma(link = "log")
)


# save model
saveRDS(model, here("outputs", "behaviour-bites", "model.rds"))

# summarise model
# p1 <- plot(model)
# ggsave(here("figures", "behaviour-bites", paste0("model_plot_", b, ".png")), plot = p1)

print(summary(model))

# posterior predictive check
p2 <- pp_check(model)
ggsave(here("figures", "behaviour-bites", "pp_check.png"), plot = p2)
