# Load required packages
library(tictoc)

tic("Initialising...")

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

toc()

# Load fitted model
tic("Loading counterfactual abundance model...")

fit <- readRDS(here("outputs/standardised_community_fit.rds"))

toc()

# Load covariates

tic("Loading covariates...")

source(here("analysis/standardise_fish.R"))

toc()

# Extract predicted counts (lambda)
tic("Extracting predicted counts...")

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

toc()

# Join metadata

tic("Joining metadata...")

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
  summarise(guild = ifelse("Piscivore" %in% guild, "Piscivore", first(guild))) %>%
  mutate(specie = as.integer(species))

species_df <- guild

species_df <- species_df %>%
  select(c(species, guild))

species_df <- setDT(species_df)

colnames(predicted_counts_df) <- str_to_lower(colnames(predicted_counts_df))

predicted_counts_df[species_df, on = "species", guild := i.guild]

toc()
# set factors

tic("Setting factors...")

predicted_counts_df[, guild := as.factor(guild)]

predicted_counts_df[, protection := factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected"))]

predicted_counts_df[, size := factor(size, levels = 1:length(size_levels), labels = size_levels)]

toc()

# sum count

tic("Summarising predicted counts...")

predicted_counts_df <- predicted_counts_df[, .(
  count = sum(count),
  N = uniqueN(species)
), by = .(draw, protection, transect, guild, size)]

tic("...for guilds..")

summary_df1 <- predicted_counts_df[, .(
  N = max(N),
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(guild, protection)]

toc()

# Summarize posterior distribution guild + size

tic("...for guilds and size...")

summary_df <- predicted_counts_df[, .(
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(guild, protection, size)]

toc()

tic("Calculating probability P(Protected > Unprotected)...")

dcast(predicted_counts_df, draw + transect + guild + size ~ protection, value.var = "count", fun.aggregate = sum) -> predicted_counts

predicted_counts[, P := Protected > Unprotected]

# P > Unprotected for guilds
predicted_counts[, .(P = mean(P, na.rm = T)), by = .(draw, guild)] -> protection_summary_guild

# P > Unprotected for guilds and size
predicted_counts[, .(P = mean(P, na.rm = T)), by = .(draw, guild, size)] -> protection_summary

toc()

tic("Plotting probability P(Protected > Unprotected)...")

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

toc()

tic("Summarising P...")

protection_summary[, .(P = mean(P, na.rm = T)), by = .(guild, size)] -> protection_summary

toc()

# Create summary table

# summary table guild + size
tic("Creating summary table for guilds and size...")

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

toc()

# summary table for herbivores only

tic("Creating summary table for herbivores only...")

herbivore_count <- predicted_counts_df[guild == "Herbivore"]
size_levels <- c("<10", "10-30", ">30")
herbivore_count[, size := as.numeric(size)]
herbivore_count[, size := case_when(
  size < 2 ~ 1,
  size == 2 | size == 3 ~ 2,
  size > 3 ~ 3
)]
herbivore_count[, size := factor(size, levels = 1:3, labels = size_levels)]
herbivore_count[, guild := NULL]

herbivore_P <- protection_summary[guild == "Herbivore"]
herbivore_P[, size := as.numeric(size)]
herbivore_P[, size := case_when(
  size < 2 ~ 1,
  size == 2 | size == 3 ~ 2,
  size > 3 ~ 3
)]
herbivore_P[, size := factor(size, levels = 1:3, labels = size_levels)]
herbivore_P[, guild := NULL]
herbivore_P <- herbivore_P[, .(P = mean(P, na.rm = TRUE)), by = size]

summary_df <- herbivore_count[, .(
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(protection, size)]

summary_df <- summary_df %>%
  left_join(herbivore_P, by = "size") %>%
  arrange(protection, size)


# Create summary table for herbivores only

summary_table <- summary_df %>%
  select(size, protection, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(size) %>%
  flextable() %>%
  set_header_labels(
    protection = "Protection",
    size = "Size class",
    mean_count = "Predicted Abundance",
    CI5 = "5% CI",
    CI25 = "25% CI",
    CI75 = "75% CI",
    CI95 = "95% CI"
  ) %>%
  merge_v(j = ~ size + P) %>% # merge rows with same size
  set_caption("Predicted abundance of herbivorous fish species by protection status and size class") %>%
  set_table_properties(width = 0.95, layout = "fixed") %>%
  align(j = c(1:4), align = "center") %>%
  align(j = c(2:4), align = "right") %>%
  theme_box()

print(summary_table)

save_as_docx(summary_table, path = here("outputs/counterfactual_abundance_summary_herbivore.docx"))

toc()

# summary table guild only
protection_summary_guild[, .(P = mean(P, na.rm = T)), by = guild]

summary_df1 <- summary_df1 %>%
  left_join(protection_summary_guild, by = "guild") %>%
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

toc()
