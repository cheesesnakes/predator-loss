# Load required packages
library(tictoc)

cat("Initialising...")
tic()

pacman::p_load(
    here, dplyr, tidyr, flextable, patchwork, stringr,
    ggplot2, readr, data.table, ggridges, fst
)

# Set default ggplot theme
theme_set(
    theme_bw() +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 14),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )
)

toc()

# Load fitted model
cat("\rLoading counterfactual abundance model...")
tic()

fit <- readRDS(here("outputs/standardised_community_fit.rds"))

toc()

# Load covariates

cat("\rLoading covariates...")
tic()

source(here("analysis/standardise_fish.R"))

rm(guild)

toc()

# Extract predicted counts (lambda)
cat("\rExtracting predicted counts...")
tic()

predicted_counts_raw <- fit$draws("abundance", format = "matrix")

predicted_counts_df <- as.data.table(predicted_counts_raw)

rm(predicted_counts_raw, fit)

gc()

predicted_counts_df[, Draw := .I]

# Melt to long format using data.table's melt (faster than pivot_longer)
predicted_counts_df <- melt(predicted_counts_df,
    id.vars = "Draw",
    variable.name = "variable", value.name = "count"
)

# Remove prefix/suffix and split into multiple columns in one go
predicted_counts_df[, variable := str_remove_all(variable, "abundance\\[|\\]")]
predicted_counts_df[, c("protection", "Transect", "Guild", "Size") :=
    tstrsplit(variable, ",", type.convert = TRUE)]

# Remove the original 'variable' column
predicted_counts_df[, variable := NULL]
colnames(predicted_counts_df) <- str_to_lower(colnames(predicted_counts_df))

toc()

# set factors

tic("Setting factors...")

predicted_counts_df[, guild := factor(guild, levels = 1:length(guild_levels), labels = guild_levels)]

predicted_counts_df[, protection := factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected"))]

predicted_counts_df[, size := factor(size, levels = 1:length(size_levels), labels = size_levels)]

toc()

# sum count

cat("Summarising predicted counts...")
tic()

predicted_counts_df <- predicted_counts_df[, .(
    count = sum(count)
), by = .(draw, protection, transect, guild, size)]

toc()

# Save

print(head(predicted_counts_df))

cat("Saving predicted counts...")

write.fst(
    predicted_counts_df,
    here("outputs/transect_posteriors.fst")
)
