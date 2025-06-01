# Load required packages
pacman::p_load(
  here, dplyr, tidyr, flextable, patchwork, stringr,
  ggplot2, readr, data.table, ggridges
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
predicted_counts_df <- predicted_counts_df %>%
  pivot_longer(
    cols = -Draw,
    names_to = "variable",
    values_to = "count"
  ) %>%
  mutate(Draw = as.integer(Draw))

predicted_counts_df <- setDT(predicted_counts_df)
predicted_counts_df[, variable := str_remove(variable, "abundance\\[")]
predicted_counts_df[, variable := str_remove(variable, "\\]")]
predicted_counts_df[, protection := as.integer(str_split(variable, ",", simplify = T)[, 1])]
predicted_counts_df[, Transect := as.integer(str_split(variable, ",", simplify = T)[, 2])]
predicted_counts_df[, Species := as.integer(str_split(variable, ",", simplify = T)[, 3])]
predicted_counts_df[, Size := as.integer(str_split(variable, ",", simplify = T)[, 4])]
predicted_counts_df[, variable := NULL]

rm(predicted_counts_raw)
gc()

# Join metadata

guild <- guild %>%
  pivot_longer(
    cols = -species,
    names_to = "guild",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  select(-value) %>%
  mutate(guild = str_to_sentence(str_trim(guild))) %>%
  group_by(species) %>%
  summarise(guild = ifelse("Piscivore" %in% guild, "Piscivore", first(guild)))

species_df <- guild

colnames(species_df) <- c("Species", "Guild")

species_df <- setDT(species_df)

predicted_counts_df[species_df, on = "Species", Guild := i.Guild]

colnames(predicted_counts_df) <- str_to_lower(colnames(predicted_counts_df))

gc()

# set factors

guild_levels <- str_split(guild_levels, ",", simplify = T)[, 1]
predicted_counts_df[, guild := factor(guild,
  levels = 1:length(guild_levels),
  labels = guild_levels
)]

predicted_counts_df[, protection := factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected"))]

predicted_counts_df[, size := factor(size, levels = 1:length(size_levels), labels = size_levels)]

# sum count

predicted_counts_df[, `:=`(
  count = sum(count),
  N = uniqueN(species)
), by = .(Draw, protection, transect, guild, size)]

# summarise for across foraging guilds

summary_df1 <- predicted_counts_df[, .(
  N = max(N),
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(guild, protection)]

# Summarize posterior distribution guild + size

summary_df <- predicted_counts_df[, .(
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(guild, protection, size)]

dcast(predicted_counts_df, draw + transect + guild + size ~ protection, value.var = "count", fun.aggregate = sum) -> predicted_counts

predicted_counts[, P := Protected > Unprotected]

# P > Unprotected for guilds
predicted_counts[, .(P = mean(P, na.rm = T)), by = .(draw, guild)] -> protection_summary1

# P > Unprotected for guilds and size
predicted_counts[, .(P = mean(P, na.rm = T)), by = .(draw, guild, size)] -> protection_summary

head(protection_summary)

p <- protection_summary %>%
  filter(guild != "Unknown") %>%
  ggplot(aes(y = size, x = P, fill = guild)) +
  geom_density_ridges(scale = 1, alpha = 0.5) +
  labs(x = "P(Protected > Unprotected)", y = "Size class") +
  scale_fill_brewer(palette = "Set1", name = "Feeding guild") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  facet_wrap(~guild, ncol = 2)

ggsave(here("figures/counterfactual_abundance_probability.png"), p, width = 8, height = 8)

protection_summary[, .(P = mean(P, na.rm = T)), by = .(guild, size)] -> protection_summary

# Create summary table

# summary table guild + size
summary_df <- summary_df %>%
  left_join(protection_summary, by = c("guild", "size")) %>%
  arrange(guild, protection, size)

summary_table <- summary_df %>%
  select(guild, size, protection, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(guild, size) %>%
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
  merge_v(j = ~ guild + size + P) %>% # merge rows with same guild and size
  set_caption("Predicted abundance of fish species by feeding guild, protection status, and size class") %>%
  set_table_properties(width = 0.95, layout = "fixed") %>%
  align(j = c(1:4), align = "center") %>%
  align(j = c(2:4), align = "right") %>%
  theme_box()

print(summary_table)

save_as_docx(summary_table, path = here("outputs/counterfactual_abundance_summary.docx"))

gc()

# summary table for herbivores only

summary_df <- summary_df %>%
  filter(guild == "Herbivore") %>%
  # make size <10, 10-30, >30
  mutate(size = case_when(
    size == "<10" ~ "<10 cm",
    size == "10-20" ~ "10-30 cm",
    size == "20-30" ~ "10-30 cm",
    size == "30-40" ~ ">30 cm",
    size == "40-50" ~ ">30 cm"
  ))

# Create summary table for herbivores only

summary_table <- summary_df %>%
  select(guild, size, protection, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(guild, size) %>%
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
  merge_v(j = ~ guild + size + P) %>% # merge rows with same guild and size
  set_caption("Predicted abundance of herbivorous fish species by protection status and size class") %>%
  set_table_properties(width = 0.95, layout = "fixed") %>%
  align(j = c(1:4), align = "center") %>%
  align(j = c(2:4), align = "right") %>%
  theme_box()

print(summary_table)

save_as_docx(summary_table, path = here("outputs/counterfactual_abundance_summary_herbivore.docx"))

# summary table guild only
summary_df1 <- summary_df1 %>%
  left_join(protection_summary1, by = "guild") %>%
  arrange(guild, protection)

summary_table1 <- summary_df1 %>%
  select(guild, N, protection, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(guild) %>%
  flextable() %>%
  set_header_labels(
    guild = "Feeding guild",
    protection = "Protection",
    mean_count = "Predicted Abundance",
    CI5 = "5% CI",
    CI25 = "25% CI",
    CI75 = "75% CI",
    CI95 = "95% CI"
  ) %>%
  merge_v(j = ~ guild + N + P) %>% # merge rows with same guild
  set_caption("Predicted abundance of fish species by feeding guild and protection status") %>%
  set_table_properties(width = 0.95, layout = "fixed") %>%
  align(j = c(1:3), align = "center") %>%
  align(j = c(2:4), align = "right") %>%
  theme_box()

print(summary_table1)
save_as_docx(summary_table1, path = here("outputs/counterfactual_abundance_summary_guild.docx"))

# Plot: Mean ± 95% credible interval
p <- summary_df %>%
  filter(guild != "Unknown") %>%
  ggplot(aes(x = size, y = mean_count, color = protection)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI5, ymax = CI95),
    width = 0, size = 0.5, position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(aes(ymin = CI25, ymax = CI75),
    width = 0, position = position_dodge(width = 0.5), size = 1,
  ) +
  labs(y = "Predicted count (mean ± 50,95% CrI)", x = "Size class") +
  scale_color_brewer(palette = "Set1", name = "Protection") +
  facet_wrap(~guild, ncol = 2, scale = "free_y") +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Save figure
ggsave(here("figures/counterfactual_abundance.png"), p, width = 8, height = 8)
