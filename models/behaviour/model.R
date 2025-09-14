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

# make data frame

data <- responses%>%
select(ind_id, foraging, vigilance, move)%>%
left_join(fish, by = c("ind_id" = "ind_id")) %>%
left_join(predictors, by = c("plot_id" = "plot_id"))%>%
left_join(guild, by =c("species" = "species"))%>%
select(-c(sponge, damsel))%>%
mutate(predator = predator - 1,
group = group - 1,
size_class = as.factor(size_class))

behaviours <- c("foraging", "vigilance", "move")
b <- behaviours[1]
for (b in behaviours) {

  # define model formula
  if (b == "foraging"){
  formula <- bf(foraging ~ protected + rugosity + biomass + predator + group + guild + size_class + (1|site:plot_id))
  }
  else if (b == "vigilance") {
    formula <- bf(vigilance ~ protected + rugosity + biomass + predator + group + guild + size_class + (1|site:plot_id))
  }
  else if (b == "move") {
    formula <- bf(move ~ protected + rugosity + biomass + predator + group + guild + size_class + (1|site:plot_id))
  }

  # run model
  model <- ordbetareg(
    formula = formula,
    data = data,
    chains = 1,
    cores = 1,
    iter = 1000
  )

  # generate stan code from formula, not fitted model
  stan_code <- make_stancode(formula, data = data)
  writeLines(stan_code, con = here("models", "behaviour", paste0("model_", b, ".stan")))

  write_rds(model, here("models", "behaviour", paste0("model_", b, ".rds")))

  # summarise model
  p1 <- plot(model)
  ggsave(here("models", "behaviour", paste0("model_plot_", b, ".png")), plot = p1)

  # posterior predictive check
  p2 <- pp_check(model)
  ggsave(here("models", "behaviour", paste0("pp_check_", b, ".png")), plot = p2)

  print(summary(model))
}

# plot posterior proportions by protection

posterior_props <- data %>%
  ungroup() %>%
  select(ind_id, plot_id, site) %>%
  distinct() %>%
  expand_grid(protected = c(0, 1), rugosity = mean(data$rugosity)) %>%
  add_epred_draws(model, value = "prop") %>%
  ungroup() %>%
  select(ind_id, .draw, protected, .category, prop)%>%
  pivot_wider(id_cols = c(ind_id, .draw, .category), names_from = protected, values_from = prop)%>%
  mutate(effect = `1` - `0`)



posterior_props%>%
group_by(.category)%>%
median_qi(effect, .width = c(0.9, 0.5)) %>%
ggplot(aes(x = .category, y = effect)) +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.5), ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  labs(x = "Behaviour", y = "Effect")

ggsave(here("models", "behaviour", "posterior_props.png"), plot = last_plot())
