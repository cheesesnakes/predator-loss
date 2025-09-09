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

## Summarize posterior distribution

cat("\rSummarizing posterior distribution...")

summary_guild <- posterior_guild[, .(
  mean_count = mean(mean_count),
  CI5 = quantile(mean_count, 0.025),
  CI25 = quantile(mean_count, 0.25),
  CI75 = quantile(mean_count, 0.75),
  CI95 = quantile(mean_count, 0.975)
), by = .(protection, guild)]

## Plot 1: posterior distribution of counts across protection

cat("\rPlotting posterior distribution of counts across protection...")

p1 <- ggplot(posterior_guild, aes(x = protection, y = mean_count, fill = protection)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~guild, scales = "free") +
  labs(x = "Protection", y = "Predicted abundance") +
  scale_fill_brewer(palette = "Set1", name = "Protection") +
  theme(legend.position = "none")

ggsave(here("figures/bruvs_predicted_guild.png"), p1, width = 8, height = 4)

## Table 1: summary statistics

cat("\rCreating summary table...")

table_1 <- summary_guild %>%
  mutate(
    mean_count = round(mean_count, 3),
    CI5 = round(CI5, 3),
    CI25 = round(CI25, 3),
    CI75 = round(CI75, 3),
    CI95 = round(CI95, 3)
  ) %>%
  arrange(guild, protection) %>%
  flextable() %>%
  set_header_labels(
    protection = "Protection",
    guild = "Feeding guild",
    mean_count = "Mean",
    CI5 = "2.5%",
    CI25 = "25%",
    CI75 = "75%",
    CI95 = "97.5%"
  ) %>%
  bold(part = "header") %>%
  set_table_properties(width = 0.5, layout = "fixed") %>%
  theme_box()

print(table_1)

save_as_docx(table_1, path = here("figures/bruvs_predicted_guild.docx"))

toc()

# Effect of protection on abundance with size classes across guilds

cat("\rEffect of protection on abundance with size classes across guilds")

tic()

## Calculate posterior distribution of counts in each deployment

posterior_guild_size <- predicted_counts_df[guild %in% c("Piscivore", "Invertivore"), .(
  mean_count = sum(count)
), by = .(draw, protection, guild, size)]

## Summarize posterior distribution

summary_guild_size <- posterior_guild_size[, .(
  mean_count = median(mean_count),
  CI5 = quantile(mean_count, 0.10),
  CI25 = quantile(mean_count, 0.25),
  CI75 = quantile(mean_count, 0.75),
  CI95 = quantile(mean_count, 0.90)
), by = .(protection, guild, size)]

## Plot 2: posterior distribution of counts across protection and size

p2 <- summary_guild_size %>%
  filter(guild == "Piscivore") %>%
  mutate(protection = ifelse(protection == "Protected", "Inside PA", "Outside PA")) %>%
  ggplot(aes(x = size, y = mean_count, col = protection)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI5, ymax = CI95), linewidth = 1, position = position_dodge(width = 0.5), width = 0.25) +
  geom_errorbar(aes(ymin = CI25, ymax = CI75), width = 0, position = position_dodge(width = 0.5), linewidth = 2) +
  labs(x = "Size class", y = "Predicted abundance") +
  scale_color_brewer(palette = "Set1", name = "Protection") +
  theme(legend.position = "top")

ggsave(here("figures/bruvs_predicted_guild_size.png"), p2, width = 8, height = 6)

## Table 2: summary statistics

table_2 <- summary_guild_size %>%
  mutate(
    mean_count = round(mean_count, 3),
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
    mean_count = "Mean",
    CI5 = "2.5%",
    CI25 = "25%",
    CI75 = "75%",
    CI95 = "97.5%"
  ) %>%
  bold(part = "header") %>%
  set_table_properties(width = 0.5, layout = "fixed") %>%
  theme_box()

print(table_2)

save_as_docx(table_2, path = here("figures/bruvs_predicted_guild_size.docx"))

toc()
