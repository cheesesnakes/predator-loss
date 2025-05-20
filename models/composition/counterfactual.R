# Load required packages
pacman::p_load(
  here, dplyr, tidyr, flextable, patchwork, stringr,
  ggplot2, readr
)

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

# Load fitted model
fit <- readRDS(here("outputs/standardised_community_fit.rds"))

# Load covariates

source(here("analysis/cleaning.R"))

# Extract predicted counts (lambda)
predicted_counts_raw <- fit$draws("abundance", format = "matrix")

# Set dimensions
n_draws <- 4000
n_transects <- nrow(transects)
n_species <- nrow(species_df)
n_sizes <- n_distinct(fish$size)

# Sanity check
stopifnot(length(predicted_counts_raw) == n_transects * n_species * n_sizes * n_draws)

# Reshape array
predicted_counts <- array(
  predicted_counts_raw,
  dim = c(n_transects, n_species, n_sizes, n_draws),
  dimnames = list(
    Transect = as.character(1:n_transects),
    Species = as.character(1:n_species),
    Size = as.character(1:n_sizes),
    Draw = paste0("draw_", 1:n_draws)
  )
)

# Convert to long data frame
predicted_counts_df <- as.data.frame.table(predicted_counts, responseName = "count") %>%
  mutate(across(c(Transect, Species, Size), as.integer)) %>%
  mutate(Draw = as.integer(str_remove(Draw, "draw_")))

# Join metadata
predicted_counts_df <- predicted_counts_df %>%
  left_join(transects, by = c("Transect" = "transect")) %>%
  left_join(species_df, by = c("Species" = "species")) %>%
  mutate(
    guild = factor(guild,
      levels = length(guild_levels),
      labels = guild_levels
    ),
    protection = factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected")),
    size = factor(size, levels = length(size_levels), labels = size_levels),
  )

# Summarize posterior distribution
summary_df <- predicted_counts_df %>%
  group_by(guild, protection, size) %>%
  summarise(
    mean_count = mean(count),
    CI5 = quantile(count, 0.025),
    CI25 = quantile(count, 0.25),
    CI75 = quantile(count, 0.75),
    CI95 = quantile(count, 0.975),
    .groups = "drop"
  )

summary_df <- summary_df %>%
  flextable() %>%
  set_header_labels(
    guild = "Feeding guild",
    protection = "Protection",
    size = "Size class",
    mean_count = "Predicted Abundance",
    CI5 = "5% CI",
    CI25 = "25% CI",
    CI75 = "75% CI",
    CI95 = "95% CI"
  ) %>%
  set_caption("Predicted abundance of fish species by feeding guild, protection status, and size class") %>%
  set_table_properties(width = 0.5) %>%
  align(j = c(1:4), align = "center") %>%
  align(j = c(2:4), align = "right") %>%
  set_table_properties(layout = "autofit")

print(summary_df)

save_as_docx(summary_df, path = here("outputs/counterfactual_abundance_summary.docx"))

# Plot: Mean ± 95% credible interval
p <- ggplot(summary_df, aes(x = size, y = mean_count, shape = protection, color = guild, group = guild)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI5, ymax = CI95),
    width = 0, size= 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(aes(ymin = CI25, ymax = CI75),
    width = 0, position = position_dodge(width = 0.5), size = 1,
  ) +
  labs(y = "Predicted count (mean ± 95% CI)", x = "Size class")+
  coord_flip()

# Save figure
ggsave(here("figures/counterfactual_abundance.png"), p, width = 6, height = 10)
