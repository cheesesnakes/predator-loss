pacman::p_load(here, data.table, dplyr, tidyr, flextable, ggplot2, stringr, ggridges, patchwork)

source(here("models", "composition", "functions.R"))
source(here("analysis/standardise_fish.R"))


# Extract summary of parameters

model_summary <- data.table()

cov_levels <- c("Overall", "Foraging guild", "Size class")

# coral cover

data = extract(p = "omega")

data[interaction == "Foraging guild", variable := factor(as.integer(variable), levels = 1:length(guild_levels), labels = guild_levels)]
data[interaction == "Size class", variable := factor(as.integer(variable), levels = 1:length(size_levels), labels = size_levels)]
data[, interaction := factor(interaction, levels = cov_levels)]

p1 <- ggplot(data[interaction == "Overall"], aes(x = value)) +
    geom_density(alpha = 0.5, fill = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Density") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16), axis.title.x = element_blank())

p2 <- ggplot(data[interaction == "Foraging guild"], aes(x = value, y = variable)) +
    geom_density_ridges(aes(fill = variable), alpha = 0.5, scale = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Foraging guild") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16),axis.title.x = element_blank())

p3 <- ggplot(data[interaction == "Size class"], aes(x = value, y = variable)) +
    geom_density_ridges(aes(fill = variable), alpha = 0.5, scale = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Size class") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

p <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1, 1)) 

ggsave(here("figures/composition/summary/omega.png"), p, width = 8, height = 9, dpi = 300)

# Save summary table
summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("model", "interaction", "variable")]
summary[, Covariate := "Coral cover"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Biomass

data = extract(p = "tau")
data[interaction == "Foraging guild", variable := factor(as.integer(variable), levels = 1:length(guild_levels), labels = guild_levels)]
data[interaction == "Size class", variable := factor(as.integer(variable), levels = 1:length(size_levels), labels = size_levels)]
data[, interaction := factor(interaction, levels = cov_levels)]

p1 <- ggplot(data[interaction == "Overall"], aes(x = value)) +
    geom_density(alpha = 0.5, fill = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Density") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16), axis.title.x = element_blank())

p2 <- ggplot(data[interaction == "Foraging guild"], aes(x = value, y = variable)) +
    geom_density_ridges(aes(fill = variable), alpha = 0.5, scale = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Foraging guild") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16),axis.title.x = element_blank())

p3 <- ggplot(data[interaction == "Size class"], aes(x = value, y = variable)) +
    geom_density_ridges(aes(fill = variable), alpha = 0.5, scale = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Size class") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

p <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1, 1)) 

ggsave(here("figures/composition/summary/tau.png"), p, width = 8, height = 9, dpi = 300)

# Save summary table

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("model", "interaction", "variable")]
summary[, Covariate := "Biomass"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()


# Full model summary

summary_table <- model_summary %>%
    select(Covariate, interaction,  variable, model, mean, q025, q975, P) %>%
    arrange(Covariate, interaction,  variable, model) %>%
    mutate_if(is.numeric, round, 3) %>%
    mutate(model = str_to_title(model), interaction = str_to_title(interaction), variable = str_to_title(variable)) %>%
    flextable() %>%
    merge_v(j = c("Covariate" , "interaction", "model", "variable")) %>%
    set_header_labels(interaction = "Interaction", variable = "Level", model = "Model", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0.95, layout = "fixed") %>%
    theme_box()

print(summary_table)

#save_as_docx(summary, path = here("outputs/composition_summary.docx"))
save_as_image(summary_table, path = here("outputs/composition_summary.png"))
