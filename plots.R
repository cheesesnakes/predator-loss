pacman::p_load(here, patchwork)

source(here("models", "behaviour-prob", "counterfactuals.R"))
source(here("models", "behaviour-time", "counterfactuals.R"))

p <- p1 + p5 + plot_layout(guides = "collect") & theme(legend.position = "top")

ggsave(here("figures", "figure3.png"), plot = p, height = 6, width = 12)

# Figure 4

p4 <- (p41 + p42) / (p43 + p44) + plot_layout(guides = "collect") & theme(legend.position = "top")

ggsave(here("figures", "figure4.png"), plot = p4, height = 12, width = 12)
