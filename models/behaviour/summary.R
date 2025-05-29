pacman::p_load(here, data.table, dplyr, tidyr, flextable, ggplot2, stringr)

source(here("models/behaviour/functions.R"))
source(here("analysis/standardise_behaviour.R"))

guild_levels <- colnames(traits)[-1]
rm(fish, predictors, responses, traits, rugosity)

# Extract summary of parameters

model_summary <- data.table()

# Protection

data = extract(p = "beta")

data[, c("protection", "response") := tstrsplit(variable, ",")]
data[, protection := factor(protection, levels = c(1, 2), label = c("Unprotected", "Protected"))]
data[, response := factor(response, levels =c(1, 2, 3), label = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]

## calculate effect
data <- dcast(data, draw + response + model ~ protection, value.var = "value", fun.aggregate = mean)
data[, value := Protected - Unprotected]
data[, `:=`(Protected = NULL, Unprotected = NULL)]

## rename columns
data[, response := str_to_title(response)]
data[, model := str_to_title(model)]

## Plot 

p <- ggplot(data, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(beta), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/beta_D.png"), p, width = 12, height = 8, dpi = 300)

## Save summary table

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "Protection"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)

rm(data, p, summary)
gc()

# Rugosity

data = extract(p = "gamma")

head(data)

data[, response := factor(variable, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]
data[, variable := NULL]

p <- ggplot(data, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(gamma), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/gamma_D.png"), p, width = 12, height = 8, dpi = 300)

# Save summary table
summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "Rugosity"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Biomass

data = extract(p = "zeta")
data[, response := factor(variable, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]
data[, variable := NULL]

p <- ggplot(data, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(zeta), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/zeta_D.png"), p, width = 12, height = 8, dpi = 300)

# Save summary table

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "Biomass"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Predator presence

data = extract(p = "eta")

data[, c("predator", "response") := tstrsplit(variable, ",")]
data[, predator := factor(predator, levels = c(1, 2), label = c("Absent", "Present"))]
data[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]

data <- dcast(data, draw + response + model ~ predator, value.var = "value", fun.aggregate = mean)
data[, value := Present - Absent]
data[, `:=`(Present = NULL, Absent = NULL)]
data[, response := str_to_title(response)]
data[, model := str_to_title(model)]

p <- ggplot(data, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(eta), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/eta_D.png"), p, width = 12, height = 8, dpi = 300)

# Save summary table
summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "Predator presence"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Effect of grouping 

data = extract(p = "delta")

data[, c("grouping", "response") := tstrsplit(variable, ",")]
data[, grouping := factor(grouping, levels = c(1, 2), label = c("Solitary", "Grouped"))]
data[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]

data <- dcast(data, draw + response + model ~ grouping, value.var = "value", fun.aggregate = mean)
data[, value := Grouped - Solitary]
data[, `:=`(Grouped = NULL, Solitary = NULL)]
data[, response := str_to_title(response)]
data[, model := str_to_title(model)]

p <- ggplot(data, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(delta), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/delta_D.png"), p, width = 12, height = 8, dpi = 300)

# Save summary table

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "Grouping"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, p, summary)
gc()

# Size class

data = extract(p = "omega")

data[, c("size_class", "response") := tstrsplit(variable, ",")]
data[, size_class := factor(size_class, levels = 1:length(size_levels), labels = size_levels)]
data[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("size_class", "response", "model")]

summary_table <- summary%>%
select(size_class, response, model, mean, q025, q975, P) %>%
    arrange(size_class, response, desc(model)) %>%
    mutate_if(is.numeric, round, 3) %>%
    flextable() %>%
    merge_v(j = c("size_class", "response")) %>%
    set_header_labels(model = "Model", response = "Behaviour", size_class = "Size class", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0, layout = "autofit") %>%
    theme_vanilla()

print(summary_table)

save_as_docx(summary_table, path = here("outputs/behaviour_summary_size.docx"))

# Guild

data = extract(p = "epsilon")

data[, c("guild", "response") := tstrsplit(variable, ",")]
data[, guild := factor(guild, levels = 1:length(guild_levels), labels = guild_levels)]
data[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[is.na(response), response := "Bite rate"]

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("guild", "response", "model")]
summary_table <- summary %>%
    select(guild, response, model, mean, q025, q975, P) %>%
    arrange(guild, response, desc(model)) %>%
    mutate_if(is.numeric, round, 3) %>%
    flextable() %>%
    merge_v(j = c("guild", "response")) %>%
    set_header_labels(model = "Model", response = "Behaviour", guild = "Guild", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0, layout = "autofit") %>%
    theme_vanilla()

print(summary_table)
save_as_docx(summary_table, path = here("outputs/behaviour_summary_guild.docx"))

# SD duration
data = extract(p = "sigma_D")
data[, response := factor(variable, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]
data[, variable := NULL]

summary <- data[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]
summary[, Variable := "SD Duration"]
model_summary <- rbind(model_summary, summary, use.names = TRUE, fill = TRUE)
rm(data, summary)
gc()

# Full model summary

summary_table <- model_summary %>%
    select(Variable, response, model, mean, q025, q975, P) %>%
    arrange(Variable, response, desc(model)) %>%
    mutate_if(is.numeric, round, 3) %>%
    flextable() %>%
    merge_v(j = c("Variable", "response")) %>%
    set_header_labels(model = "Model", response = "Behaviour", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0, layout = "autofit") %>%
    theme_vanilla()

print(summary_table)

save_as_docx(summary, path = here("outputs/behaviour_summary.docx"))
