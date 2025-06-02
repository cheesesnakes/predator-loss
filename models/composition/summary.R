pacman::p_load(here, data.table, dplyr, tidyr, flextable, ggplot2, stringr)

source(here("models", "composition", "functions.R"))
source(here("analysis/standardise_fish.R"))


# Extract summary of parameters

model_summary <- data.table()

# coral cover

data = extract(p = "omega")

head(data)

p <- ggplot(data, aes(x = value)) +
    geom_density(alpha = 0.5, fill = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(omega), y = "Density") +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ model, scales = "free") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/composition/summary/omega.png"), p, width = 6, height = 3, dpi = 300)

# Save summary table
summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = model]
summary[, Variable := "Coral cover"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Biomass

data = extract(p = "tau")

p <- ggplot(data, aes(x = value)) +
    geom_density(fill = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(tau), y = "Density") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))+
    facet_wrap(~ model, scales = "free")

ggsave(here("figures/composition/summary/tau.png"), p, width = 6, height = 3, dpi = 300)

# Save summary table

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = model]
summary[, Variable := "Biomass"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()


# Full model summary

summary_table <- model_summary %>%
    select(Variable, model, mean, q025, q975, P) %>%
    arrange(Variable, model) %>%
    mutate_if(is.numeric, round, 3) %>%
    mutate(model = str_to_title(model)) %>%
    flextable() %>%
    merge_v(j = c("Variable" , "model")) %>%
    set_header_labels(model = "Model", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0.95, layout = "fixed") %>%
    theme_box()

print(summary_table)

#save_as_docx(summary, path = here("outputs/composition_summary.docx"))
save_as_image(summary_table, path = here("outputs/composition_summary.png"))
