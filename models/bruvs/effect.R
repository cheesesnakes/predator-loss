# Load libraries
pacman::p_load(
  here, dplyr, tidyr, flextable, patchwork, stringr,
  ggplot2, readr, data.table, ggridges, tictoc, fst
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

# Load posteriors

cat("Loading posteriors...")
tic()

predicted_counts_df <- read.fst(
  here("outputs/bruvs_posteriors.fst")
)

setDT(predicted_counts_df)

toc()

# Effect of protection on abundance across guilds

cat("\rEffect of protection on abundance across guilds")

tic()

## Calculate posterior distribution of counts in each deployment

cat("\rCalculating posterior distribution of counts in each deployment...")

posterior_guild <- predicted_counts_df[guild %in% c("Piscivore", "Invertivore"), .(
  mean_count = sum(count)
), by = .(draw, protection, transect, guild)]

posterior_guild <- dcast(
  posterior_guild,
  draw + guild + transect ~ protection,
  value.var = "mean_count",
  fill = 0
)

posterior_guild[, effect := (log2(Protected + 1) - log2(Unprotected + 1))]

head(posterior_guild)

## Summarize posterior distribution

cat("\rSummarizing posterior distribution...")

summary_guild <- posterior_guild[, .(
  mean_effect = mean(effect),
  CI5 = quantile(effect, 0.025),
  CI25 = quantile(effect, 0.25),
  CI75 = quantile(effect, 0.75),
  CI95 = quantile(effect, 0.975)
), by = .(guild)]

## Plot 1: posterior distribution of counts across protection

cat("\rPlotting posterior distribution of counts across protection...")

p1 <- posterior_guild %>%
  filter(guild == "Piscivore") %>%
  ggplot(aes(x = effect)) +
  geom_density(fill = "black", col = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = mean_effect), data = filter(summary_guild, guild == "Piscivore"), linetype = "dashed", linewidth = 1, col = "darkred") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", linewidth = 1) +
  labs(x = "Log fold change", y = "Density") +
  theme(legend.position = "none")


ggsave(here("figures/bruvs_effect_guild.png"), p1, width = 6, height = 4)

## Table 1: summary statistics

cat("\rCreating summary table...")

table_1 <- summary_guild %>%
  mutate(
    mean_effect = round(mean_effect, 3),
    CI5 = round(CI5, 3),
    CI25 = round(CI25, 3),
    CI75 = round(CI75, 3),
    CI95 = round(CI95, 3)
  ) %>%
  arrange(guild) %>%
  flextable() %>%
  set_header_labels(
    guild = "Feeding guild",
    mean_effect = "Effect",
    CI5 = "2.5%",
    CI25 = "25%",
    CI75 = "75%",
    CI95 = "97.5%"
  ) %>%
  bold(part = "header") %>%
  set_table_properties(width = 0.5, layout = "fixed") %>%
  theme_box()

print(table_1)

save_as_docx(table_1, path = here("figures/bruvs_effect_guild.docx"))

toc()

# Effect of protection on abundance with size classes across guilds

cat("\rEffect of protection on abundance with size classes across guilds")

tic()

## Calculate posterior distribution of counts in each deployment

posterior_guild_size <- predicted_counts_df[guild %in% c("Piscivore", "Invertivore"), .(
  mean_count = sum(count)
), by = .(draw, protection, guild, size)]

posterior_guild_size <- dcast(
  posterior_guild_size,
  draw + guild + size ~ protection,
  value.var = "mean_count",
  fill = 0
)

posterior_guild_size[, effect := (log2(Protected + 1) - log2(Unprotected + 1))]

## Summarize posterior distribution

summary_guild_size <- posterior_guild_size[, .(
  mean_effect = mean(effect),
  CI5 = quantile(effect, 0.025),
  CI25 = quantile(effect, 0.25),
  CI75 = quantile(effect, 0.75),
  CI95 = quantile(effect, 0.975)
), by = .(guild, size)]

## Plot 2: posterior distribution of counts across protection and size

p2 <- posterior_guild_size %>%
  filter(guild == "Piscivore") %>%
  ggplot(aes(x = size, y = effect)) +
  geom_boxplot(fill = "black", col = "black", width = 0.25, alpha = 0.5, outlier.shape = NA) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 1) +
  labs(x = "Size class", y = "Log fold change") +
  theme(legend.position = "none")

ggsave(here("figures/bruvs_effect_guild_size.png"), p2, width = 10, height = 6)

## Table 2: summary statistics

table_2 <- summary_guild_size %>%
  mutate(
    mean_effect = round(mean_effect, 3),
    CI5 = round(CI5, 3),
    CI25 = round(CI25, 3),
    CI75 = round(CI75, 3),
    CI95 = round(CI95, 3)
  ) %>%
  arrange(guild, protection, size) %>%
  flextable() %>%
  set_header_labels(
    protection = "Protection",
    guild = "Feeding guild",
    size = "Size class",
    mean_effect = "Effect",
    CI5 = "2.5%",
    CI25 = "25%",
    CI75 = "75%",
    CI95 = "97.5%"
  ) %>%
  bold(part = "header") %>%
  set_table_properties(width = 0.5, layout = "fixed") %>%
  theme_box()

print(table_2)

save_as_docx(table_2, path = here("figures/bruvs_effect_guild_size.docx"))

toc()

# Effect of protection on species richness of predators

cat("\rEffect of protection on species richness of predators")

tic()
