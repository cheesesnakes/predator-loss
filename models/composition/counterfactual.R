# Load required packages
pacman::p_load(
  here, dplyr, tidyr, flextable, patchwork, stringr,
  ggplot2, readr, data.table
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

source(here("analysis/standardise_fish.R"))

# Extract predicted counts (lambda)
predicted_counts_raw <- fit$draws("abundance", format = "matrix")

rm(fit)
gc()

predicted_counts_df <- as.data.table(predicted_counts_raw)
predicted_counts_df[, Draw := 1:nrow(predicted_counts_df)]
predicted_counts_df <- predicted_counts_df%>%
  pivot_longer(
    cols = -Draw,
    names_to = "variable",
    values_to = "count"
  )%>%
  mutate(Draw = as.integer(Draw)) 

predicted_counts_df <- setDT(predicted_counts_df)
predicted_counts_df[, variable := str_remove(variable, "abundance\\[")]
predicted_counts_df[, variable := str_remove(variable, "\\]")]
predicted_counts_df[, Transect := as.integer(str_split(variable, ",", simplify=T)[,1])]
predicted_counts_df[, Species := as.integer(str_split(variable, ",", simplify=T)[,2])]
predicted_counts_df[, Size := as.integer(str_split(variable, ",", simplify=T)[,3])]
predicted_counts_df[, variable := NULL]

rm(predicted_counts_raw)
gc()

# Join metadata

transects <- select(transects, c(transect, protection))
species_df <- select(species_df, c(species, guild))

transects <- mutate_all(transects, ~as.integer(.))
species_df <- mutate_all(species_df, ~as.integer(.))

colnames(transects) <- c("Transect", "Protection")
colnames(species_df) <- c("Species", "Guild")

transects <- setDT(transects)
species_df <- setDT(species_df)

predicted_counts_df[transects, on="Transect", Protection := i.Protection]
predicted_counts_df[species_df, on = "Species", Guild := i.Guild]

colnames(predicted_counts_df) <- str_to_lower(colnames(predicted_counts_df))

gc()

# set factors

guild_levels <- str_split(guild_levels, ",", simplify = T)[,1]
predicted_counts_df[, guild := factor(guild,
                                      levels = 1:length(guild_levels),
                                      labels = guild_levels
)]

predicted_counts_df[, protection := factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected"))]

predicted_counts_df[, size := factor(size, levels = 1:length(size_levels), labels = size_levels)]

# sum count

predicted_counts_df[, species:=NULL]

predicted_counts_df[, count := sum(count), by = .(draw, protection, transect, guild, size)]

# Summarize posterior distribution

summary_df <- predicted_counts_df[, .(
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(guild, protection, size)]

summary_table <- summary_df %>%
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

print(summary_table)

save_as_docx(summary_table, path = here("outputs/counterfactual_abundance_summary.docx"))

gc()

# Plot: Mean ± 95% credible interval
p <- ggplot(summary_df, aes(x = size, y = mean_count, color = protection)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI5, ymax = CI95),
    width = 0, size= 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(aes(ymin = CI25, ymax = CI75),
    width = 0, position = position_dodge(width = 0.5), size = 1,
  ) +
  labs(y = "Predicted count (mean ± 50,95% CrI)", x = "Size class")+
  scale_color_brewer(palette = "Set1", name = "Protection") +
  scale_shape(name = "Protection") +
  facet_wrap(~guild, ncol = 1) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Save figure
ggsave(here("figures/counterfactual_abundance.png"), p, width = 6, height = 12)
