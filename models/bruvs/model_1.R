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

data <- fish %>%
mutate(deployment_id = as.numeric(factor(deployment_id, levels = deployment_levels))) %>%
left_join(predictors, by = c("deployment_id"))%>%
group_by(protection, location, site, deployment_id, aspect) %>%
summarise(abundance = sum(max_n))%>%
ungroup()%>%
mutate(protection = protection - 1)

# define formula

formula <- bf(abundance ~ protection + aspect + (1|site))

# priors

priors <- c(set_prior("normal(0,0.5)", class = "Intercept"),
          set_prior("normal(0,0.5)", class = "b"),
          set_prior("exponential(2)", class = "sd"))

# fit model

model <- brm(formula, data = data, family = "poisson", chains = 4, cores =4, iter = 2000, file = here("outputs", "bruvs", "model_1"))

# generate stan code from formula, not fitted model
stan_code <- make_stancode(formula, data = data)
writeLines(stan_code, con = here("models", "bruvs", "model_1.stan"))

# save model
saveRDS(model, here("outputs", "bruvs", "model_1.rds"))

# summarise model

print(summary(model))

# posterior predictive check

p2 <- pp_check(model)

ggsave(here("figures", "bruvs", paste0("pp_check_model_1.png")), plot = p2)

# Random effects

random_effects <- model %>%
    spread_draws(`r_site`[l_s, ]) %>%
    summarise_draws() %>%
    ungroup() %>%
    mutate(
        Site = l_s
    ) %>%
    select(-c(l_s, variable)) %>%
    select(Site,  mean, median, sd, q5, q95) %>%
    mutate(Site = str_to_title(site_levels[Site]))

# make table for random intercepts

table_re(random_effects, var_names = "Site", path = here("outputs", "bruvs", "model_1_intercepts.docx"))

# Parameters

summarise_params(model, type = "Parameter", vars = c("b_Intercept", "b_protection", "b_aspect"), path = here("outputs", "bruvs", "model_1_parameters.docx"))

# Goodness of fit

summarise_r2(model, path = here("outputs", "bruvs", "model_1_R2.docx"))

# plot posterior proportions by protection

posterior_protection <- data %>%
    ungroup() %>%
    select(site, aspect) %>%
    distinct() %>%
    expand_grid(protection = c(0, 1)) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

effect <- posterior_protection %>%
        pivot_wider(names_from = protection, values_from = prop, id_cols = c(site, .draw)) %>%
        rename("Unprotected" = `0`, "Protected" = `1`) %>%
        mutate(
            effect = log2(Protected) - log2(Unprotected)
        )

P <- effect %>%
summarise(P = mean(effect > 0))

## Across guilds and size classes

effect %>%
ggplot(aes(x = effect)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = median(effect$effect), linetype = "dashed", color = "red", linewidth = 1) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    labs(
        x = "Effect of Protection",
        y = "Density"
    )

ggsave(here("figures", "bruvs", "effect_model_1.png"), plot = last_plot(), height = 6, width = 8)

## Table

effect_summary <- effect %>%
    cross_join(P) %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    cross_join(P)

table3 <- effect_summary %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI",
        P = "P(effect > 0)"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "bruvs", paste0("effect_protection.docx")))
