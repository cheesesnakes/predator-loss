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
  here("outputs/bruvs-family_posteriors.fst")
)

setDT(predicted_counts_df)

toc()

# Get family guild information

guilds <- read.csv(here("outputs/standardised_data/traits_bruvs.csv"))

guilds %>%
  select(family, guild) %>%
  filter(
    guild != "Unknown",
    !(family == "Lethrinidae" & guild == "Invertivore"),
    !(family == "Labridae" & guild == "Piscivore")
  ) %>%
  distinct() -> guilds

predicted_counts_df <- predicted_counts_df[guilds, on = .(family)]

# Effect of protection on abundance across families

cat("\rEffect of protection on abundance across families")

tic()

## Calculate posterior distribution of counts in each deployment

cat("\rCalculating posterior distribution of counts in each deployment...")

posterior_family <- predicted_counts_df[, .(
  mean_count = sum(count)
), by = .(draw, protection, transect, family, guild)]

## Summarize posterior distribution

cat("\rSummarizing posterior distribution...")

summary_family <- posterior_family[, .(
  mean_count = median(mean_count),
  CI5 = quantile(mean_count, 0.10),
  CI25 = quantile(mean_count, 0.25),
  CI75 = quantile(mean_count, 0.75),
  CI95 = quantile(mean_count, 0.90)
), by = .(protection, family, guild)]

## Plot 1: posterior distribution of counts across protection

cat("\rPlotting posterior distribution of counts across protection...")

p1 <- ggplot(summary_family %>% filter(guild != "Invertivore"), aes(x = reorder(family, mean_count), y = mean_count, col = protection)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = CI5, ymax = CI95), linewidth = 1, width = 0.5) +
  geom_errorbar(aes(ymin = CI25, ymax = CI75), linewidth = 2, width = 0) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = "Family", y = "Predicted abundance") +
  theme(legend.position = "none") +
  coord_flip() +
  facet_wrap(~protection, labeller = as_labeller(c("Protected" = "Inside PA", "Unprotected" = "Outside PA")))

ggsave(here("figures/bruvs_predicted_family.png"), p1, width = 8, height = 4)

## Table 1: summary statistics

cat("\rCreating summary table...")

table_1 <- summary_family %>%
  mutate(
    mean_count = round(mean_count, 3),
    CI5 = round(CI5, 3),
    CI25 = round(CI25, 3),
    CI75 = round(CI75, 3),
    CI95 = round(CI95, 3)
  ) %>%
  arrange(family, protection) %>%
  flextable() %>%
  set_header_labels(
    protection = "Protection",
    family = "Family",
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

save_as_docx(table_1, path = here("figures/bruvs_predicted_family.docx"))

toc()

# Effect of protection on abundance with size classes across families

cat("\rEffect of protection on abundance with size classes across families")

tic()

## Calculate posterior distribution of counts in each deployment

posterior_family_size <- predicted_counts_df[, .(
  mean_count = sum(count)
), by = .(draw, protection, family, size)]

## Summarize posterior distribution

summary_family_size <- posterior_family_size[, .(
  mean_count = mean(mean_count),
  CI5 = quantile(mean_count, 0.025),
  CI25 = quantile(mean_count, 0.25),
  CI75 = quantile(mean_count, 0.75),
  CI95 = quantile(mean_count, 0.975)
), by = .(protection, family, size)]

## Plot 2: posterior distribution of counts across protection and size

p2 <- posterior_family_size %>%
  mutate(protection = ifelse(protection == "Protected", "Inside PA", "Outside PA")) %>%
  ggplot(aes(x = size, y = mean_count, fill = protection)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  labs(x = "Size class", y = "Predicted abundance") +
  scale_fill_brewer(palette = "Set1", name = "Protection") +
  scale_y_continuous(limits = c(0, 750)) +
  theme(legend.position = "top") +
  facet_wrap(~family, ncol = 4, scales = "free_y")

ggsave(here("figures/bruvs_predicted_family_size.png"), p2, width = 8, height = 6)

## Table 2: summary statistics

table_2 <- summary_family_size %>%
  mutate(
    mean_count = round(mean_count, 3),
    CI5 = round(CI5, 3),
    CI25 = round(CI25, 3),
    CI75 = round(CI75, 3),
    CI95 = round(CI95, 3)
  ) %>%
  arrange(family, protection, size) %>%
  flextable() %>%
  set_header_labels(
    protection = "Protection",
    family = "Family",
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

save_as_docx(table_2, path = here("figures/bruvs_predicted_family_size.docx"))

toc()
