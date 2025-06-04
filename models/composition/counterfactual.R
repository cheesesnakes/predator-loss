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
rm(predicted_counts_raw)
gc()

predicted_counts_df[, Draw := .I]

# Melt to long format using data.table's melt (faster than pivot_longer)
predicted_counts_df <- melt(predicted_counts_df,
  id.vars = "Draw",
  variable.name = "variable", value.name = "count"
)

# Remove prefix/suffix and split into multiple columns in one go
predicted_counts_df[, variable := str_remove_all(variable, "abundance\\[|\\]")]
predicted_counts_df[, c("protection", "Transect", "Species", "Size") :=
  tstrsplit(variable, ",", type.convert = TRUE)]

# Remove the original 'variable' column
predicted_counts_df[, variable := NULL]

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

toc()

# summary of abundance by protection status

head(predicted_counts_df)

summary_dist <- predicted_counts_df[guild != "Unknown", .(
  mean_count = mean(count)
), by = .(draw, protection, guild)]

summary_count <- predicted_counts_df[, .(
  mean_count = mean(count),
  CI5 = quantile(count, 0.025),
  CI25 = quantile(count, 0.25),
  CI75 = quantile(count, 0.75),
  CI95 = quantile(count, 0.975)
), by = .(protection, guild, size)]

p1 <- summary_count %>%
  filter(guild != "Unknown") %>%
  ggplot(aes(x = size, y = mean_count, color = protection)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = CI5, ymax = CI95),
    width = 0,
    position = position_dodge(width = 0.5),
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = CI25, ymax = CI75),
    width = 0,
    position = position_dodge(width = 0.5),
    linewidth = 1.5
  ) +
  geom_line(aes(group = protection), position = position_dodge(width = 0.5), linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Protection") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "Size class", y = "Abundance") +
  facet_wrap(~guild, ncol = 3, scale = "free_y") +
  theme_bw() +
  theme(
    legend.position = "top"
  )

p2 <- ggplot(summary_dist, aes(x = mean_count, fill = protection)) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean abundance", y = "Density") +
  scale_fill_brewer(palette = "Set1", name = "Protection") +
  facet_wrap(~guild, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# Reshape to wide format for protection status
tic("Reshaping to wide format...")
dcast(predicted_counts_df, draw + transect + guild + size + N ~ protection, value.var = "count", fun.aggregate = sum) -> predicted_counts_df

predicted_counts_df[, P := Protected > Unprotected]

predicted_counts_df[, Delta := log(Protected + 1) - log(Unprotected + 1)]

toc()

tic("...for guilds..")

summary_df_guild <- predicted_counts_df[, .(
  N = max(N),
  mean_count = mean(Delta),
  P = mean(P, na.rm = TRUE)
), by = .(draw, guild)]

toc()

# Summarize posterior distribution guild + size

tic("...for guilds and size...")

summary_df <- predicted_counts_df[, .(
  mean_count = mean(Delta),
  P = mean(P, na.rm = TRUE)
), by = .(draw, guild, size)]

toc()

tic("Plotting probability P(Protected > Unprotected)...")

p <- summary_df %>%
  filter(guild != "Unknown") %>%
  ggplot(aes(y = size, x = P, fill = guild)) +
  geom_density_ridges(scale = 1, alpha = 0.5) +
  labs(x = "P(Protected > Unprotected)", y = "Size class") +
  scale_fill_brewer(palette = "Set1", name = "Feeding guild") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
  ) +
  facet_wrap(~guild, ncol = 3)

ggsave(here("figures/counterfactual_abundance_probability.png"), p, width = 8, height = 4)

# Plot: Mean ± 95% credible interval

summary_df[, .(
  mean_count = mean(mean_count),
  CI5 = quantile(mean_count, 0.025, na.rm = TRUE),
  CI25 = quantile(mean_count, 0.25, na.rm = TRUE),
  CI75 = quantile(mean_count, 0.75, na.rm = TRUE),
  CI95 = quantile(mean_count, 0.975, na.rm = TRUE),
  P = mean(P, na.rm = TRUE)
), by = .(guild, size)] -> summary_df

p3 <- summary_df %>%
  filter(guild != "Unknown") %>%
  ggplot(aes(x = size, y = mean_count)) +
  geom_point(, size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = CI5, ymax = CI95),
    width = 0,
    position = position_dodge(width = 0.5),
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = CI25, ymax = CI75),
    width = 0,
    position = position_dodge(width = 0.5),
    linewidth = 1.5
  ) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.5
  ) +
  geom_line(aes(group = guild), position = position_dodge(width = 0.5), linetype = "dashed") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "Size class", y = "Log-fold change") +
  facet_wrap(~guild, ncol = 3, scale = "free_y") +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# Save figure
ggsave(here("figures/counterfactual_abundance.png"), p3, width = 9, height = 4)

(p2 / p1 / p3) +
  plot_layout(guides = "collect") +
  theme(legend.position = "top")
# Save figure
ggsave(here("figures/figure2.png"), width = 8, height = 8)

toc()

# Create summary table

# summary table guild + size
tic("Creating summary table for guilds and size...")

summary_table <- summary_df %>%
  select(guild, size, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(guild, size) %>%
  flextable() %>%
  set_header_labels(
    guild = "Feeding guild",
    P = "P(Protected > Unprotected)",
    protection = "Protection",
    size = "Size class",
    mean_count = "LFC",
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

summary_df <- herbivore_count[, .(
  mean_count = mean(Delta),
  CI5 = quantile(Delta, 0.025),
  CI25 = quantile(Delta, 0.25),
  CI75 = quantile(Delta, 0.75),
  CI95 = quantile(Delta, 0.975),
  P = mean(P, na.rm = TRUE)
), by = .(size)]

# Create summary table for herbivores only

summary_table <- summary_df %>%
  select(size, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(size) %>%
  flextable() %>%
  set_header_labels(
    size = "Size class",
    P = "P(Protected > Unprotected)",
    mean_count = "LFC",
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
summary_df_guild <- summary_df_guild[, .(
  N = max(N),
  mean_count = mean(mean_count),
  CI5 = quantile(mean_count, 0.025),
  CI25 = quantile(mean_count, 0.25),
  CI75 = quantile(mean_count, 0.75),
  CI95 = quantile(mean_count, 0.975),
  P = mean(P, na.rm = TRUE)
), by = guild]

summary_table1 <- summary_df_guild %>%
  select(guild, N, P, mean_count, CI5, CI25, CI75, CI95) %>%
  arrange(guild) %>%
  flextable() %>%
  set_header_labels(
    guild = "Feeding guild",
    P = "P(Protected > Unprotected)",
    mean_count = "LFC",
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
