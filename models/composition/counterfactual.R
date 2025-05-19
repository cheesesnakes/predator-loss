pacman::p_load(here, dplyr, tidyr, flextable, patchwork, stringr, ggplot2)

theme_set(
    theme_bw() +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 20)
        )
)

fit <- readRDS(here("outputs/validation_community_fit.rds"))

# covariates
fish <- read.csv(here("outputs/simulated_data/fish.csv")) %>%
    mutate(observation = as.character(1:nrow(.)))
transects <- read.csv(here("outputs/simulated_data/transects.csv"))
species <- read.csv(here("outputs/simulated_data/species.csv"))

# output

predicted_counts <- fit$draws("lambda", format = "matrix")

predicted_counts <- array(predicted_counts, dim = c(nrow(transects), nrow(species), n_distinct(fish$Size), 4000))

predicted_counts_df <- as.data.frame.table(predicted_counts) %>%
    rename(Transect = Var1, Species = Var2, Size = Var3, count = Freq) %>%
    mutate(
        Transect = as.numeric(Transect),
        Species = as.numeric(Species),
        Size = as.numeric(Size)
    ) %>%
    left_join(transects, by = c("Transect" = "Transect")) %>%
    left_join(species, by = c("Species" = "Species"))

predicted_counts_df <- predicted_counts_df %>%
    mutate(
        Feeding_guild = factor(as.numeric(Feeding_guild), levels = c(1, 2, 3, 4)),
        Protection = factor(as.numeric(Protection), levels = c(1, 2)),
        Size = factor(as.numeric(Size), levels = c(1, 2, 3, 4, 5, 6))
    )

ggplot(predicted_counts_df, aes(x = Size, y = count, fill = Protection)) +
    geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.5)) +
    geom_violin(alpha = 0.5, position = position_dodge(width = 0.5), width = 0.5) +
    facet_wrap(~Feeding_guild, ncol = 1)

ggsave(here("figures/counterfactual_abundance.png"), width = 6, height = 10)
