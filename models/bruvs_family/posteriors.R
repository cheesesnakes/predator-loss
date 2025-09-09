# Load required packages
library(tictoc)

cat("Initialising...")
tic()

pacman::p_load(
    here, dplyr, tidyr, flextable, patchwork, stringr,
    ggplot2, readr, data.table, ggridges, fst
)

toc()

# Load fitted model
cat("\rLoading counterfactual abundance model...")

tic()

fit <- readRDS(here("outputs/standardised_bruvs-family_fit.rds"))

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

cat("\nMelting data...\n")

predicted_counts_df <- melt(predicted_counts_df,
    id.vars = "Draw",
    variable.name = "variable", value.name = "count"
)

# Filter count > 0

predicted_counts_df <- predicted_counts_df[count > 0]

print(paste("Number of rows:", nrow(predicted_counts_df)))

# Remove prefix/suffix and split into multiple columns in one go
cat("\nFixing variable names...")

predicted_counts_df[, variable := gsub("abundance\\[|\\]", "", variable, perl = TRUE)]
predicted_counts_df[, c("protection", "transect", "family", "size") := tstrsplit(variable, ",", type.convert = TRUE)]
predicted_counts_df[, variable := NULL]

gc()

toc()

# Load covariates

cat("\rLoading covariates...")

tic()

source(here("analysis/standardise_bruvs.R"))

rm(fish, guild, individuals, metadata, no_guild, predictors, traits)

gc()

toc()

# Join metadata

colnames(predicted_counts_df) <- str_to_lower(colnames(predicted_counts_df))


# set factors

cat("\rSetting factors...")

tic()

predicted_counts_df[, family := factor(family, levels = 1:length(family_levels), labels = family_levels)]

predicted_counts_df[, protection := factor(protection, levels = c(1, 2), labels = c("Unprotected", "Protected"))]

predicted_counts_df[, size := factor(size, levels = 1:length(size_levels), labels = size_levels)]

toc()

# summary of abundance by protection status

print(head(predicted_counts_df))

cat("Saving predicted counts...")

write.fst(
    predicted_counts_df,
    here("outputs/bruvs-family_posteriors.fst")
)
