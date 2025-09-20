pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_bruvs.R"))
source(here("functions", "summaries.R"))

dir.create(here("outputs", "bruvs"), showWarnings = FALSE)
dir.create(here("figures", "bruvs"), showWarnings = FALSE)

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

# standardise fish data

data <- counts %>%
uncount(abundance) %>%
left_join(predictors, by = c("deployment_id"))%>%
mutate(size = as.numeric(size))

# define formula

formula <- bf(size ~ protection + aspect + (1|site))

# priors

priors <- c(set_prior("normal(0,0.5)", class = "Intercept"),
          set_prior("normal(0,0.5)", class = "b"),
          set_prior("exponential(2)", class = "sd"))

# fit model

model <- brm(formula, data = data,
              family = sratio(), chains = 4, cores =4, iter = 2000, file = here("outputs", "bruvs", "model_2"))

# generate stan code from formula, not fitted model
stan_code <- make_stancode(formula, data = data)
writeLines(stan_code, con = here("models", "bruvs", "model_2.stan"))

# save model
saveRDS(model, here("outputs", "bruvs", "model_2.rds"))

# summarise model

print(summary(model))

# posterior predictive check

p2 <- pp_check(model)

ggsave(here("figures", "bruvs", paste0("pp_check_model_2.png")), plot = p2)

# Random effects

random_effects <- model %>%
    spread_draws(`r_site`[l_s, ]) %>%
    summarise_draws() %>%
    ungroup() %>%
    mutate(
        Site = l_s
    ) %>%
    select(-c(l_s, variable)) %>%
    select(Site,  mean, median, sd, q5, q95)%>%
    mutate(Site = str_to_title(site_levels[Site]))

# make table for random intercepts

table_re(random_effects, var_names = "Site", path = here("outputs", "bruvs", "model_2_intercepts.docx"))

# Parameters
get_variables(model)

summarise_params(model, type = "Parameter", vars = c("b_protection", "b_aspect"), path = here("outputs", "bruvs", "model_2_parameters.docx"))

# Goodness of fit

summarise_r2(model, path = here("outputs", "bruvs", "model_2_R2.docx"))

# plot posterior proportions by protection

posterior_prediction <- data %>%
    ungroup() %>%
    select(site, aspect) %>%
    distinct() %>%
    expand_grid(protection = c(0, 1)) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

effect <- posterior_prediction %>%
        pivot_wider(names_from = protection, values_from = prop, id_cols = c(site, aspect, .category, .draw)) %>%
        rename("Unprotected" = `0`, "Protected" = `1`) %>%
        mutate(
            effect = Protected - Unprotected
        )

P <- effect %>%
group_by(.category) %>%
rename(size = .category) %>%
summarise(P = mean(effect > 0))
    
## Across guilds and size classes

posterior_prediction%>%
rename(size = .category) %>%
group_by(protection, size)%>%
median_qi(prop, .width = c(0.5, 0.9)) %>%
mutate(size = size_levels[as.numeric(size)],
protection = case_when(protection == 0 ~ "Outside PA",
TRUE ~ "Inside PA")) %>%
ggplot(aes(x = size, y = prop, col = protection)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), width = 0.2, position = position_dodge(width = 0.5)) +
    labs(
        x = "Size",
        y = "Expected Proportion",
        col = ""
    )+
    scale_color_brewer(palette = "Set1")

ggsave(here("figures", "bruvs", "effect_model_2.png"), plot = last_plot(), height = 6, width = 8)

## Table
effect_summary <- effect %>%
    rename(size = .category) %>%
    group_by(size) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    mutate(size = size_levels[as.numeric(size)],) %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(size, effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    left_join(P, by = "size") 

table3 <- effect_summary %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        size = "Size",
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "bruvs", paste0("effect_protection_size.docx")))
