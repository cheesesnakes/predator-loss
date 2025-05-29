pacman::p_load(here, data.table, dplyr, tidyr, flextable, ggplot2, stringr)

# Extract summary of parameters

parameters <- c("beta", "gamma", "delta", "epsilon", "zeta", "eta", "omega")
responses <- c("pi", "D", "bites")
hyper <- c("sigma_D")

fit <- readRDS(here("outputs/standardised_behaviour_fit.rds"))

draw <- fit$draws("beta_D")

df <- as.data.table(draw)

df[, draw := as.numeric(iteration) * as.numeric(chain)]
df[, iteration := NULL]
df[, chain := NULL]
df[, variable := str_remove(variable, "beta_D\\[")]
df[, variable := str_remove(variable, "\\]")]
df[, c("protection", "response") := tstrsplit(variable, ",")]
df[, protection := as.integer(protection)]
df[, response := as.integer(response)]
df[, variable := NULL]

df <- dcast(df, draw + response ~ protection, value.var = "value", fun.aggregate = mean)

df[, value := `2` - `1`]
df[, `:=`(`1` = NULL, `2` = NULL)]

df[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]

df_bites <- fit$draws("beta_bites")
df_bites <- as.data.table(df_bites)
df_bites[, draw := as.numeric(iteration) * as.numeric(chain)]
df_bites[, iteration := NULL]
df_bites[, chain := NULL]
df_bites[, variable := str_remove(variable, "beta_bites\\[")]
df_bites[, variable := str_remove(variable, "\\]")]
df_bites[, c("protection") := tstrsplit(variable, ",")]
df_bites[, protection := as.integer(protection)]
df_bites[, response := as.integer(response)]
df_bites[, variable := NULL]
df_bites <- dcast(df_bites, draw ~ protection, value.var = "value", fun.aggregate = mean)
df_bites[, value := `2` - `1`]
df_bites[, `:=`(`1` = NULL, `2` = NULL)]
df_bites[, response := "Bites"]

df <- rbind(df, df_bites, use.names = TRUE, fill = TRUE)
rm(df_bites)
gc()

df[, model := "Behaviour"]

draw <- fit$draws("beta_pi")

df_pi <- as.data.table(draw)

df_pi[, draw := as.numeric(iteration) * as.numeric(chain)]
df_pi[, iteration := NULL]
df_pi[, chain := NULL]
df_pi[, variable := str_remove(variable, "beta_pi\\[")]
df_pi[, variable := str_remove(variable, "\\]")]
df_pi[, c("protection", "response") := tstrsplit(variable, ",")]
df_pi[, protection := as.integer(protection)]
df_pi[, response := as.integer(response)]
df_pi[, variable := NULL]
df_pi <- dcast(df_pi, draw + response ~ protection, value.var = "value", fun.aggregate = mean)
df_pi[, value := `2` - `1`]
df_pi[, `:=`(`1` = NULL, `2` = NULL)]   
df_pi[, response := factor(response, levels = c(1, 2, 3), labels = c("Foraging", "Vigilance", "Movement"))]

df_pi[, model := "Zero-Inflation"]

df <- rbind(df, df_pi, use.names = TRUE, fill = TRUE)
rm(df_pi)
gc()

head(df)

p <- ggplot(df, aes(x = value)) +
    geom_density(aes(fill = response), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(beta[D]), y = "Density") +
    facet_grid(model ~ response) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/behaviour/summary/beta_D.png"), p, width = 12, height = 8, dpi = 300)

# Save summary table

summary <- df[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975), P = mean(value > 0, na.rm = TRUE)), by = c("response", "model")]

summary <- summary %>%
    arrange(response, desc(model)) %>%
    mutate_if(is.numeric, round, 3) %>%
    flextable() %>%
    merge_v(j = "response") %>%
    set_header_labels(model = "Model",response = "Behaviour", mean = "Mean", q025 = "2.5% CrI", q975 = "97.5% CrI") %>%
    set_table_properties(width = 0, layout = "autofit") %>%
    theme_vanilla()

print(summary)

save_as_docx(summary, path = here("outputs/behaviour_summary.docx"))
