pacman::p_load(here, data.table, dplyr, tidyr, flextable, ggplot2, stringr)

# Extract summary of parameters

fit <- readRDS(here("outputs/standardised_rugosity_fit.rds"))
draw <- fit$draws("beta")
data <- as.data.table(draw)

# Create draw index
data[, draw := as.numeric(iteration) * as.numeric(chain)]
data[, c("iteration", "chain") := NULL]

# clean variable names

data[, variable := str_remove(variable, "beta\\[")]
data[, variable := str_remove(variable, "\\]")]

# Protection

data[, protection := factor(variable, levels = c(1, 2), label = c("Unprotected", "Protected"))]
data[, variable := NULL]

## calculate effect
data <- dcast(data, draw ~ protection, value.var = "value", fun.aggregate = mean)
data[, value := Protected - Unprotected]
data[, `:=`(Protected = NULL, Unprotected = NULL)]

head(data)

## Plot 

p <- ggplot(data, aes(x = value)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(x = expression(beta), y = "Density") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 16))

ggsave(here("figures/environment/summary/beta.png"), p, width = 12, height = 8, dpi = 300)

## Save summary table

summary <- data[, .(
    mean = mean(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    P = mean(value > 0, na.rm = TRUE)
)]
print(summary)
