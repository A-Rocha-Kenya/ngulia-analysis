library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(splines)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "04_age_sex_demography")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

cli_h1("Assess demographic-model sensitivity")

age_data <- read_csv(
  file.path(model_data_dir, "age_daily_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
sex_data <- read_csv(
  file.path(model_data_dir, "adult_sex_daily_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
age_models <- readRDS(file.path(model_dir, "age_composition_models.rds"))
sex_models <- readRDS(file.path(model_dir, "adult_sex_composition_models.rds"))

fit_sensitivities <- function(data, models, response) {
  groups <- data |> group_split(avibase_id, common_name)
  map_dfr(groups, function(species_data) {
    species_id <- species_data$avibase_id[[1]]
    selected_rhs <- models[[species_id]]$formula_rhs
    linear_rhs <- sub(
      "splines::ns\\(season_centered, df = 2\\)",
      "I(season_centered / 10)", selected_rhs
    )
    top_threshold <- quantile(species_data$n_total, 0.99)
    specifications <- list(
      capped_primary = list(data = species_data, weight = pmin(species_data$n_total, 50)),
      equal_date = list(data = species_data, weight = rep(1, nrow(species_data))),
      full_individual = list(data = species_data, weight = species_data$n_total),
      remove_top_1_percent = list(
        data = species_data |> filter(n_total <= top_threshold),
        weight = pmin(species_data$n_total[species_data$n_total <= top_threshold], 50)
      )
    )

    imap_dfr(specifications, function(specification, sensitivity) {
      model_weights <- specification$weight
      model <- glm(
        as.formula(paste(response, "~", linear_rhs)),
        data = specification$data,
        weights = model_weights,
        family = quasibinomial()
      )
      term <- "I(season_centered/10)"
      estimate <- coef(model)[term]
      standard_error <- sqrt(vcov(model)[term, term])
      tibble(
        avibase_id = species_id,
        common_name = species_data$common_name[[1]],
        sensitivity,
        n_dates = nrow(specification$data),
        effective_sample_size = sum(model_weights),
        odds_ratio_per_decade = exp(estimate),
        lower = exp(estimate - 1.96 * standard_error),
        upper = exp(estimate + 1.96 * standard_error),
        direction = case_when(
          lower > 1 ~ "increase",
          upper < 1 ~ "decrease",
          TRUE ~ "uncertain"
        )
      )
    })
  })
}

age_sensitivity <- fit_sensitivities(age_data, age_models, "first_year_fraction")

# The sex primary cap is 30 rather than 50; replace its primary rows.
sex_sensitivity <- fit_sensitivities(sex_data, sex_models, "male_fraction") |>
  mutate(sensitivity = if_else(sensitivity == "capped_primary", "capped_50", sensitivity))

sex_primary <- sex_data |>
  group_split(avibase_id, common_name) |>
  map_dfr(function(species_data) {
    species_id <- species_data$avibase_id[[1]]
    linear_rhs <- sub(
      "splines::ns\\(season_centered, df = 2\\)",
      "I(season_centered / 10)", sex_models[[species_id]]$formula_rhs
    )
    model_weights <- pmin(species_data$n_total, 30)
    model <- glm(
      as.formula(paste("male_fraction ~", linear_rhs)),
      data = species_data,
      weights = model_weights,
      family = quasibinomial()
    )
    term <- "I(season_centered/10)"
    estimate <- coef(model)[term]
    standard_error <- sqrt(vcov(model)[term, term])
    tibble(
      avibase_id = species_id,
      common_name = species_data$common_name[[1]],
      sensitivity = "capped_primary",
      n_dates = nrow(species_data),
      effective_sample_size = sum(model_weights),
      odds_ratio_per_decade = exp(estimate),
      lower = exp(estimate - 1.96 * standard_error),
      upper = exp(estimate + 1.96 * standard_error),
      direction = case_when(
        lower > 1 ~ "increase",
        upper < 1 ~ "decrease",
        TRUE ~ "uncertain"
      )
    )
  })

sex_sensitivity <- bind_rows(sex_primary, sex_sensitivity) |>
  arrange(common_name, match(
    sensitivity,
    c("capped_primary", "equal_date", "full_individual", "remove_top_1_percent", "capped_50")
  ))

age_robustness <- age_sensitivity |>
  group_by(avibase_id, common_name) |>
  summarise(
    primary_direction = direction[sensitivity == "capped_primary"],
    n_sensitivity_directions = n_distinct(direction),
    minimum_odds_ratio = min(odds_ratio_per_decade),
    maximum_odds_ratio = max(odds_ratio_per_decade),
    conclusion_stable = n_sensitivity_directions == 1,
    .groups = "drop"
  )

write_csv(age_sensitivity, file.path(table_dir, "age_trend_weighting_sensitivity.csv"))
write_csv(sex_sensitivity, file.path(table_dir, "adult_sex_trend_weighting_sensitivity.csv"))
write_csv(age_robustness, file.path(table_dir, "age_trend_robustness_summary.csv"))

cli_alert_success("Compared capped, equal-date, full-individual and high-catch-deletion trends.")
