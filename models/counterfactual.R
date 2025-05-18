pacman::p_load(here, dplyr, tidyr, flextable, patchwork, stringr, ggplot2)

theme_set(
    theme_bw() +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 20)
        )
)

fit <- readRDS(here("outputs/validation_fit.rds"))

# load data 
individuals <- read.csv(here("outputs/simulated_data/individuals.csv"))
predictors <- read.csv(here("outputs/simulated_data/predictors.csv"))

covariate <- individuals%>%
left_join(predictors, by = c("plot_id"))

# extract counterfactual draws

# predictive checks

D_pred <- fit$draws("D_pred", format = "matrix")

D_pred <- array(D_pred, dim = c(500, 3, 4000))

D_pred_long <- as.data.frame.table(D_pred, responseName = "value")%>%
    rename(individual = Var1, response = Var2, draw = Var3)%>%
    mutate(
        draw = as.numeric(draw),
        response = recode(as.numeric(response),
                        "1" = "Foraging",
                        "2" = "Vigilance",
                        "3" = "Movement"),
        individual = as.numeric(individual)
    )%>%
    group_by(response, draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )


ggplot(D_pred_long, aes(x = response, y = mean, fill = response)) +
    geom_boxplot(width = 0.05, outlier.shape = NA) +
    geom_violin(aes(y = mean), alpha = 0.5, width = 0.5) +
    labs(x = "Behaviour", y = "Duration (s)")+
    scale_y_log10()+
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(
        legend.position = "none"
    ) 

ggsave(here("figures", "D_pred.png"), width = 10, height = 8)

bites_pred <- fit$draws("bites_pred", format = "matrix")
bites_pred <- array(bites_pred, dim = c(500, 4000))
bites_pred_long <- as.data.frame.table(bites_pred, responseName = "value")%>%
    rename(individual = Var1, draw = Var2)%>%
    mutate(
        draw = as.numeric(draw),
        individual = as.numeric(individual)
    )%>%
    group_by(draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )

ggplot(bites_pred_long, aes(x = mean)) +
geom_density(aes(y = ..density..), fill = "blue", alpha = 0.5) +
    labs(x = "Bites", y = "Density") +
    theme_bw()

ggsave(here("figures", "bites_pred.png"), width = 10, height = 8)

## counterfactual protection

protection <- fit$draws("D_pred_protection", format = "matrix")
protection <- array(protection, dim = c(2, 500, 3, 4000))
diff_1 <- protection[1,,,] - protection[2,,,]
protection_long <- as.data.frame.table(diff_1, responseName = "value")%>%
    rename(individual = Var1, response = Var2, draw = Var3)%>%
    mutate(
        draw = as.numeric(draw),
        response = recode(as.numeric(response),
                        "1" = "Foraging",
                        "2" = "Vigilance",
                        "3" = "Movement"),
        individual = as.numeric(individual)
    )%>%
    group_by(response, draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )

p1 <- ggplot(protection_long, aes(x = response, y = mean, fill = response)) +
    geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), outlier.shape = NA) +
    geom_violin(aes(y = mean), alpha = 0.5, width = 0.5) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    labs(x = "Behaviour", y = expression(Delta~"Duration (s)"))+
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(
        legend.position = "none"
    )

# counterfactual protection bites
protection_bites <- fit$draws("bites_pred_protection", format = "matrix")
protection_bites <- array(protection_bites, dim = c(2, 500, 4000))
diff_2 <- protection_bites[1,,] - protection_bites[2,,]
protection_bites_long <- as.data.frame.table(diff_2, responseName = "value")%>%
    rename(individual = Var1, draw = Var2)%>%
    mutate(
        draw = as.numeric(draw),
        individual = as.numeric(individual)
    )%>%
    group_by(draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )
p2 <- ggplot(protection_bites_long, aes(x = mean)) +
    geom_density(aes(y = ..density..), fill = "blue", alpha = 0.5) +
    labs(x = expression(Delta~"Bites"), y = "Density") +
    theme_bw() +
    theme(
        legend.position = "top"
    )

p1/p2

ggsave(here("figures", "counterfactual_protection.png"), width = 6, height = 8)

# counterfactual protection varying biomass

protection_biomass <- fit$draws("D_pred_protection_biomass", format = "matrix")
protection_biomass <- array(protection_biomass, dim = c(2, 500, 3, 4000))
diff_3 <- protection_biomass[1,,,] - protection_biomass[2,,,]
protection_biomass_long <- as.data.frame.table(diff_3, responseName = "value")%>%
    rename(individual = Var1, response = Var2, draw = Var3)%>%
    mutate(
        draw = as.numeric(draw),
        response = recode(as.numeric(response),
                        "1" = "Foraging",
                        "2" = "Vigilance",
                        "3" = "Movement"),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, response, value, biomass))

p1 <- protection_biomass_long%>%
  filter(draw %in% sample(1:4000, 50))%>%
  ggplot(aes(x = biomass, y = value, col = response)) +
  geom_smooth(aes(group = draw), method = "loess", se = FALSE, alpha = 0.05, color = "gray") +  # raw posterior draws
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +  # mean trend
  facet_wrap(~ response, ncol=1) +
  labs(x = "Biomass", y = expression(Delta~"Duration (s)"))+
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(
        legend.position = "none"
    )

# counterfactual protection bites
protection_bites_biomass <- fit$draws("bites_pred_protection_biomass", format = "matrix")
protection_bites_biomass <- array(protection_bites_biomass, dim = c(2, 500, 4000))
diff_4 <- protection_bites_biomass[1,,] - protection_bites_biomass[2,,]
protection_bites_biomass_long <- as.data.frame.table(diff_4, responseName = "value")%>%
    rename(individual = Var1, draw = Var2)%>%
    mutate(
        draw = as.numeric(draw),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, value, biomass))

p2 <- protection_bites_biomass_long%>%
    filter(draw %in% sample(1:4000, 50))%>%
    ggplot(aes(x = biomass, y = value)) +
    geom_smooth(aes(group = draw), method = "loess", se = FALSE, alpha = 0.05, color = "gray") +  # raw posterior draws
    geom_smooth(method = "loess", se = FALSE, size = 1.2) +  # mean trend
    labs(x = "Biomass", y = expression(Delta~"Bites")) +
    theme_bw() +
    theme(
        legend.position = "top"
    )

p1/p2 + plot_layout(heights = c(0.75, 0.25))

ggsave(here("figures", "counterfactual_protection_biomass.png"), width = 6, height = 16)

# counterfactual protection varying rugosity
protection_rugosity <- fit$draws("D_pred_protection_rugosity", format = "matrix")
protection_rugosity <- array(protection_rugosity, dim = c(2, 500, 3, 4000))
diff_5 <- protection_rugosity[1,,,] - protection_rugosity[2,,,]
protection_rugosity_long <- as.data.frame.table(diff_5, responseName = "value")%>%
    rename(individual = Var1, response = Var2, draw = Var3)%>%
    mutate(
        draw = as.numeric(draw),
        response = recode(as.numeric(response),
                        "1" = "Foraging",
                        "2" = "Vigilance",
                        "3" = "Movement"),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, response, value, mean_rugosity))
p1 <- protection_rugosity_long%>%
    filter(draw %in% sample(1:4000, 50))%>%
    ggplot(aes(x = mean_rugosity, y = value, col = response)) +
    geom_smooth(aes(group = draw), method = "loess", se = FALSE, alpha = 0.05, color = "gray") +  # raw posterior draws
    geom_smooth(method = "loess", se = FALSE, size = 1.2) +  # mean trend
    facet_wrap(~ response, ncol=1) +
    labs(x = "Rugosity", y = expression(Delta~"Duration (s)"))+
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(
        legend.position = "none"
    )
# counterfactual protection bites
protection_bites_rugosity <- fit$draws("bites_pred_protection_rugosity", format = "matrix")
protection_bites_rugosity <- array(protection_bites_rugosity, dim = c(2, 500, 4000))
diff_6 <- protection_bites_rugosity[1,,] - protection_bites_rugosity[2,,]
protection_bites_rugosity_long <- as.data.frame.table(diff_6, responseName = "value")%>%
    rename(individual = Var1, draw = Var2)%>%
    mutate(
        draw = as.numeric(draw),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, value, mean_rugosity))
p2 <- protection_bites_rugosity_long%>%
    filter(draw %in% sample(1:4000, 50))%>%
    ggplot(aes(x = mean_rugosity, y = value)) +
    geom_smooth(aes(group = draw), method = "loess", se = FALSE, alpha = 0.05, color = "gray") +  # raw posterior draws
    geom_smooth(method = "loess", se = FALSE, size = 1.2) +  # mean trend
    labs(x = "Rugosity", y = expression(Delta~"Bites")) +
    theme_bw() +
    theme(
        legend.position = "top"
    )
p1/p2 + plot_layout(heights = c(0.75, 0.25))
ggsave(here("figures", "counterfactual_protection_rugosity.png"), width = 6, height = 16)

# counterfactual protection predator presence
protection_predator <- fit$draws("D_pred_protection_predator", format = "matrix")
protection_predator <- array(protection_predator, dim = c(2, 500, 3, 4000))
diff_7 <- protection_predator[1,,,] - protection_predator[2,,,]
protection_predator_long <- as.data.frame.table(diff_7, responseName = "value")%>%
    rename(individual = Var1, response = Var2, draw = Var3)%>%
    mutate(
        draw = as.numeric(draw),
        response = recode(as.numeric(response),
                        "1" = "Foraging",
                        "2" = "Vigilance",
                        "3" = "Movement"),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, response, value, predator))%>%
    mutate(predator = recode(as.numeric(predator),
                            "1" = "Absent",
                            "2" = "Present"))%>%
    group_by(predator, response, draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )

p1 <- ggplot(protection_predator_long, aes(x = response, y = mean, fill = response)) +
    geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), outlier.shape = NA) +
    geom_violin(aes(y = mean), alpha = 0.5, width = 0.5) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    labs(x = "Behaviour", y = expression(Delta~"Duration (s)"))+
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(
        legend.position = "none"
    )
# counterfactual protection bites
protection_bites_predator <- fit$draws("bites_pred_protection_predator", format = "matrix")
protection_bites_predator <- array(protection_bites_predator, dim = c(2, 500, 4000))
diff_8 <- protection_bites_predator[1,,] - protection_bites_predator[2,,]
protection_bites_predator_long <- as.data.frame.table(diff_8, responseName = "value")%>%
    rename(individual = Var1, draw = Var2)%>%
    mutate(
        draw = as.numeric(draw),
        individual = as.numeric(individual)
    )%>%
    left_join(covariate, by = c("individual" = "ind_id"))%>%
    select(c(individual, draw, value, predator))%>%
    mutate(predator = recode(as.numeric(predator),
                            "1" = "Absent",
                            "2" = "Present"))%>%
    group_by(predator, draw)%>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        median = median(value),
        q025 = quantile(value, 0.025),
        q975 = quantile(value, 0.975)
    )
p2 <- ggplot(protection_bites_predator_long, aes(x = mean, fill = predator)) +
    geom_density(aes(y = ..density..), alpha = 0.5) +
    labs(x = expression(Delta~"Bites"), y = "Density") +
    theme_bw() +
    theme(
        legend.position = "top",
        legend.title = element_blank()
    )
p1/p2 
ggsave(here("figures", "counterfactual_protection_predator.png"), width = 6, height = 8)
