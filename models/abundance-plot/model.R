pacman::p_load(brms, here, tidybayes, ggplot2, marginaleffects)

source(here("analysis", "standardise_behaviour.R"))

dir.create(here("outputs", "abundance-plot"))
dir.create(here("figures", "abundance-plot"))

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

# make data frame

data <- plot_abundance %>%
  left_join(guild, by = "species") %>%
  filter(!is.na(guild)) %>%
  group_by(plot_id, guild) %>%
  summarise(abundance = sum(abundance)) %>%
  left_join(predictors, by = "plot_id") %>%
  ungroup()

print(paste0("N = ", length(unique(data$plot_id)), " plots"))

# define formula

formula <- bf(abundance ~ protected + guild + rugosity_mean + damsel + biomass + protected:guild + rugosity_mean:guild + biomass:guild + (1 | site))

# priors

priors <- c(
  set_prior("normal(0,0.5)", class = "Intercept"),
  set_prior("normal(0,0.5)", class = "b"),
  set_prior("exponential(2)", class = "sd"),
  set_prior("exponential(2)", class = "shape")
)

# run the model

model <- brm(
  formula = formula,
  prior = priors,
  data = data,
  chains = 4,
  cores = 4,
  iter = 2000,
  family = negbinomial(link = "log")
)

# generate stan code

stand_code <- make_stancode(
  formula = formula,
  data = data,
  family = negbinomial(link = "log"),
  prior = priors
)

writeLines(stand_code, con = here("models", "abundance-plot", "model.stan"))

# save model
saveRDS(model, here("outputs", "abundance-plot", "model.rds"))

# summarise model
print(summary(model))

# p1 <- plot(model)
# ggsave(here("models", "behaviour-prob", paste0("model_plot_", b, ".png")), plot = p1)

# posterior predictive check
p2 <- pp_check(model)
ggsave(here("figures", "abundance-plot", paste0("pp_check.png")), plot = p2)
