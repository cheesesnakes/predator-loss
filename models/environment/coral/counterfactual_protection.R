pacman::p_load(
  here,
  dplyr,
  tidyr,
  flextable,
  ggplot2,
  data.table,
  stringr
)

theme_set(theme_bw() +
  theme(legend.position = "top", legend.text = element_text(size = 20)))

## counterfactual protection

fit <- readRDS(here("outputs", "standardised_coral_fit.rds"))

protection <- fit$draws("coral_cf", format = "matrix")

protection_df <- as.data.table(protection) %>%
  mutate(Draw = 1:nrow(protection)) %>%
  pivot_longer(
    cols = -c(Draw),
    names_to = "variable",
    values_to = "value"
  )

rm(protection, fit)
gc()

setDT(protection_df)
protection_df[, variable := str_remove(variable, "coral_cf\\[")]
protection_df[, variable := str_remove(variable, "\\]")]
protection_df[, c("protection", "plot_id") := tstrsplit(variable, ",")]
protection_df[, protection := factor(protection,
  levels = 1:2,
  labels = c("Outside", "Inside")
)]
protection_df[, variable := NULL]

dcast(protection_df,
  Draw + plot_id ~ protection,
  value.var = "value"
) -> protection_df
protection_df[, Delta := log2(Inside + 1) - log2(Outside + 1)]
protection_df[, P := Inside > Outside] # P(Inside > Outside)

summary_protection <- protection_df[, .(
  mean = mean(Delta, na.rm = TRUE),
  median = median(Delta, na.rm = TRUE),
  q025 = quantile(Delta, 0.025, na.rm = TRUE),
  q250 = quantile(Delta, 0.25, na.rm = TRUE),
  q750 = quantile(Delta, 0.75, na.rm = TRUE),
  q975 = quantile(Delta, 0.975, na.rm = TRUE),
  P = sum(P, na.rm = TRUE) / .N
), by = .(Draw)]

# plot counterfactual protection
ggplot(summary_protection, aes(x = mean)) +
  geom_density(fill = "black", alpha = 0.5) +
  geom_vline(
    aes(xintercept = 0),
    linetype = "dashed",
    color = "black"
  ) +
  geom_vline(
    aes(xintercept = mean(summary_protection$mean)),
    linetype = "dotted",
    color = "red"
  ) +
  labs(
    x = "Log-fold change",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
  ) -> p1

ggsave(
  here("figures", "counterfactual_coral.png"),
  plot = p1,
  width = 6,
  height = 4
)

# save summary
head(summary_protection)

summary_protection %>%
  summarise(
    median = median(mean),
    q025 = quantile(mean, 0.025),
    q250 = quantile(mean, 0.25),
    q750 = quantile(mean, 0.75),
    q975 = quantile(mean, 0.975),
    mean = mean(mean),
    P = mean(P, na.rm = TRUE)
  )%>%
  print(.)
