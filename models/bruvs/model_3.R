pacman::p_load(brms, here, ordbetareg, tidybayes, ggplot2, marginaleffects, flextable)

source(here("analysis", "standardise_bruvs.R"))

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

data <- family_counts %>%
left_join(predictors, by = c("deployment_id"))%>%
group_by(protection, location, site, deployment_id, family) %>%
summarise(abundance = sum(abundance))%>%
ungroup()%>%
mutate(protection = protection - 1,
family = as.factor(family))


# define formula

formula <- bf(abundance ~ family +(1|site))

# priors

priors <- c(set_prior("normal(0,0.5)", class = "Intercept"),
          set_prior("normal(0,0.5)", class = "b"),
          set_prior("exponential(2)", class = "sd"))

# fit model

model <- brm(formula, data = data, family = "poisson", chains = 4, cores =4, iter = 2000)

# generate stan code from formula, not fitted model
stan_code <- make_stancode(formula, data = data)
writeLines(stan_code, con = here("models", "bruvs", "model_3.stan"))

# save model
saveRDS(model, here("outputs", "bruvs", "model_3.rds"))

# summarise model

print(summary(model))

# posterior predictive check

p2 <- pp_check(model)

ggsave(here("figures", "bruvs", paste0("pp_check_model_3.png")), plot = p2)

# Random effects

random_effects <- model %>%
    spread_draws(`r_site`[l_s, ]) %>%
    summarise_draws() %>%
    ungroup() %>%
    mutate(
        Site = l_s
    ) %>%
    select(-c(l_s, variable)) %>%
    select(Site,  mean, median, sd, q5, q95)

# make table for random intercepts

table1 <- random_effects %>%
select(-c(sd)) %>%
mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        Site = "Site",
        mean = "Mean",
        median = "Median",
        sd = "SD",
        q5 = "Lower CI",
        q95 = "Upper CI"
    ) %>%
    merge_v(~ Site) %>%
    set_table_properties(width = 0.5, layout = "fixed") %>%
    theme_box()

print(table1)

save_as_docx(table1, path = here("outputs", "bruvs", "model_3_intercepts.docx"))


# Parameters

parameters <- model %>%
    spread_draws(b_Intercept, regex = TRUE) %>%
    summarise_draws() %>%
    mutate(
        type = "Parameter",
        variable = str_replace(variable, "b_", "")
    )

# make table

table2 <- parameters %>%
    select(c(variable, mean, median, q5, q95)) %>%
    mutate(variable = str_to_title(variable)) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        variable = "Variable",
        mean = "Mean",
        median = "Median",
        sd = "SD",
        q5 = "Lower CI",
        q95 = "Upper CI"
    ) %>%
    set_table_properties(width = 0.5, layout = "fixed") %>%
    theme_box()

print(table2)

save_as_docx(table2, path = here("outputs", "bruvs", "model_3_parameters.docx"))

# plot posterior proportions by protection

posterior_prediction <- data %>%
    ungroup() %>%
    select(site) %>%
    distinct() %>%
    expand_grid(family = c(1:12)) %>%
    add_epred_draws(model, value = "prop") %>%
    ungroup()

## Across guilds and size classes

posterior_prediction%>%
group_by(family)%>%
median_qi(prop, .width = c(0.5, 0.9)) %>%
mutate(family = family_levels[as.numeric(family)]) %>%
ggplot(aes(x = reorder(family, prop), y = prop)) +
    geom_pointinterval(aes(ymin = .lower, ymax = .upper), width = 0.2) +
    labs(
        x = "Family",
        y = "Expected Abundance"
    )+
    coord_flip()

ggsave(here("figures", "bruvs", "effect_model_3.png"), plot = last_plot(), height = 6, width = 8)

## Table

table3 <- effect %>%
    median_qi(effect, .width = c(0.5, 0.9)) %>%
    ungroup() %>%
    pivot_wider(names_from = .width, values_from = c(.lower, .upper)) %>%
    select(effect, .lower_0.9, .lower_0.5, .upper_0.5, .upper_0.9) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    flextable() %>%
    set_header_labels(
        effect = "Median",
        .lower_0.9 = "90% Lower CI",
        .lower_0.5 = "50% Lower CI",
        .upper_0.5 = "50% Upper CI",
        .upper_0.9 = "90% Upper CI"
    ) %>%
    theme_box() %>%
    set_table_properties(width = 0.5, layout = "fixed")

print(table3)

save_as_docx(table3, path = here("outputs", "bruvs", paste0("effect_protection_guild_size.docx")))
