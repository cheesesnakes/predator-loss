pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_fish.R"))

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

dir.create(here("outputs", "habitat"), showWarnings = FALSE)
dir.create(here("figures", "habitat"), showWarnings = FALSE)

# make data frame

data <- benthic %>%
    left_join(transects, by = c("transect_id" = "transect"))

# load models

model <- readRDS(here("outputs", "habitat", paste0("model.rds")))

# Counterfactual effect of protection on habitat composition

posterior_protection <- data %>%
ungroup() %>%
select(location, site, transect_id, quadrat, aspect)%>%
distinct() %>%
expand_grid(protection = c(1, 2)) %>%
add_predicted_draws(model, value = "prop") %>%
ungroup()

effect <- posterior_protection %>%
    pivot_wider(names_from = protection, values_from = prop, id_cols = c(transect_id, quadrat, .draw, .category)) %>%
    rename("Unprotected" = `1`, "Protected" = `2`) %>%
    mutate(
        effect = Protected - Unprotected
    )

P <- effect %>%
    group_by(.category) %>%
    summarise(prob = mean(effect > 0))

# Summary table

effect_summary <- effect %>%
group_by(.category) %>%
median_qi(effect, .width = c(0.5, 0.9)) %>%
pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
left_join(P, by = ".category") %>%
select(.category, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, prob) %>%
mutate(across(where(is.numeric), ~ round(.x, 3)))

table3 <- effect_summary %>%
    mutate(.category = str_to_title(.category),
           effect= effect * 100,
           .lower_0.9 = .lower_0.9 * 100,
           .lower_0.5 = .lower_0.5 * 100,
           .upper_0.5 = .upper_0.5 * 100,
           .upper_0.9 = .upper_0.9 * 100) %>%
    flextable() %>%
    set_header_labels(
        .category = "Substrate Type",
        effect = "Median (%)",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI",
        prob = "P(effect > 0)"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "habitat", paste0("effect_protection.docx")))


# plot effect of protection on habitat composition

effect %>%
mutate(.category = str_to_title(.category),
effect= effect * 100) %>%
group_by(.category) %>%
median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = .category, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper, color = .category), position = position_dodge(width = 0.5)) +
    labs(
        x = "Substrate Type",
        y = expression(Delta ~ "Cover (%)")
    ) +
    theme(legend.position = "none")

ggsave(here("figures", "habitat", "effect_protection.png"), width = 8, height = 6)

# Counterfactual effect of aspect on habitat composition

posterior_aspect <- data %>%
    ungroup() %>%
    select(location, site, transect_id, quadrat, protection) %>%
    distinct() %>%
    expand_grid(aspect = c(1,2)) %>%
    add_predicted_draws(model, value = "prop") %>%
    ungroup()

effect <- posterior_aspect %>%
    pivot_wider(names_from = aspect, values_from = prop, id_cols = c(transect_id, quadrat, .draw, .category)) %>%
    rename("Flat" = `1`, "Slope" = `2`) %>%
    mutate(
        effect = Slope - Flat
    )

P <- effect %>%
    group_by(.category) %>%
    summarise(prob = mean(effect > 0))

# Summary table

effect_summary <- effect %>%
group_by(.category) %>%
median_qi(effect, .width = c(0.5, 0.9)) %>%
pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
left_join(P, by = ".category") %>%
select(.category, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9, prob) %>%
mutate(across(where(is.numeric), ~ round(.x, 3)))

table3 <- effect_summary %>%
    mutate(.category = str_to_title(.category),
           effect= effect * 100,
           .lower_0.9 = .lower_0.9 * 100,
           .lower_0.5 = .lower_0.5 * 100,
           .upper_0.5 = .upper_0.5 * 100,
           .upper_0.9 = .upper_0.9 * 100) %>%
    flextable() %>%
    set_header_labels(
        .category = "Substrate Type",
        effect = "Median (%)",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI",
        prob = "P(effect > 0)"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "habitat", paste0("effect_aspect.docx")))

# Plot effect of aspect on habitat composition

effect %>%
    mutate(.category = str_to_title(.category),
           effect= effect * 100) %>%
    group_by(.category) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ggplot(aes(x = .category, y = effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper, color = .category), position = position_dodge(width = 0.5)) +
    labs(
        x = "Substrate Type",
        y = expression(Delta ~ "Cover (%)")
    ) +
    theme(legend.position = "none")

ggsave(here("figures", "habitat", "effect_aspect.png"), width = 8, height = 6)
