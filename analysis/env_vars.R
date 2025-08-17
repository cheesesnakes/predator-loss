pacman::p_load(
  here,
  dplyr,
  tidyr,
  ggplot2
)

source(here("analysis", "cleaning_behaviour.R"))

rm(fish, fish_traits, guild, multiguild, phylo, responses, traits, get_mode, trait_mode)

head(predictors)

predictors_long <- predictors %>%
  group_by(plot_id, location, site, protection) %>%
  pivot_longer(cols = c(resource, rugosity, predator, damsel, coral, sponge, substrate), names_to = "variable", values_to = "value")


ggplot(predictors_long, aes(y = value, x = variable)) +
  geom_boxplot() +
  facet_wrap(~ protection) +
  labs(title = "Distribution of Environmental Variables",
       x = "Variable",
       y = "Value")

ggsave(here("figures", "behaviour_env.png"), width = 10, height = 6)
