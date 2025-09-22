pacman::p_load(here, patchwork)

source(here("models", "behaviour-prob", "counterfactuals.R"))
source(here("models", "behaviour-time", "counterfactuals.R"))

p <- p1 + p5 + plot_layout(guides = "collect") & theme(legend.position = "top")

ggsave(here("figures", "figure3.png"), plot = p, height = 6, width = 12)
