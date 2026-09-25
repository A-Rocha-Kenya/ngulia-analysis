library(dplyr)
library(tidyr)
library(readr)
library(splines)
library(purrr)
library(mgcv)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
source(file.path(project_dir, "scripts", "helpers", "weighted_quantile.R"))

analysis_dir <- file.path(paths$analysis_output_dir, "04_age_sex_demography")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cli_h1("Fit age and adult sex composition models")

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
day_support <- read_csv(file.path(table_dir, "day_of_season_support.csv"), show_col_types = FALSE)

# Shared model definitions ------------------------------------------------

condition_terms <- paste(
  "moon_distance_from_new_moon_z + mist_probability_light_patchy +",
  "mist_probability_good + era5_rain_log_z + wind_speed_10m_mean_ms_z +",
  "temperature_2m_mean_c_z + surface_pressure_mean_hpa_z"
)

model_specifications <- tribble(
  ~model_id, ~model, ~formula_rhs,
  "M0", "Constant composition", "1",
  "M1", "+ within-season timing", "splines::ns(day_of_season, df = 3)",
  "M2", "+ broad long-term change", paste(
    "splines::ns(day_of_season, df = 3)",
    "splines::ns(season_centered, df = 2)", sep = " + "
  ),
  "M3", "+ daily grounding conditions", paste(
    "splines::ns(day_of_season, df = 3)",
    "splines::ns(season_centered, df = 2)", condition_terms, sep = " + "
  )
)

set.seed(402)
season_folds <- tibble(season = sort(unique(age_data$season))) |>
  mutate(fold = sample(rep(seq_len(5), length.out = n())))
write_csv(season_folds, file.path(table_dir, "season_validation_folds.csv"))

fit_quasibinomial <- function(data, response, weight, formula_rhs) {
  model_formula <- as.formula(paste(response, "~", formula_rhs))
  model_weights <- data[[weight]]
  glm(model_formula, data = data, weights = model_weights, family = quasibinomial())
}

log_loss <- function(observed, predicted, weight) {
  predicted <- pmin(pmax(predicted, 1e-8), 1 - 1e-8)
  -weighted.mean(observed * log(predicted) + (1 - observed) * log(1 - predicted), weight)
}

validate_species <- function(data, response, weight) {
  data <- data |> left_join(season_folds, by = "season")
  map_dfr(seq_len(nrow(model_specifications)), function(i) {
    heldout <- rep(NA_real_, nrow(data))
    for (fold in sort(unique(data$fold))) {
      fit <- fit_quasibinomial(
        data |> filter(.data$fold != .env$fold), response, weight,
        model_specifications$formula_rhs[[i]]
      )
      heldout[data$fold == fold] <- predict(
        fit, newdata = data[data$fold == fold, ], type = "response"
      )
    }
    tibble(
      model_id = model_specifications$model_id[[i]],
      model = model_specifications$model[[i]],
      heldout_log_loss = log_loss(data[[response]], heldout, data[[weight]])
    )
  })
}

fit_analysis <- function(data, analysis, response, count_success, count_failure, weight) {
  cli_h2("{analysis} models")

  validation <- data |>
    group_split(avibase_id, common_name) |>
    map_dfr(function(species_data) {
      validate_species(species_data, response, weight) |>
        mutate(
          avibase_id = species_data$avibase_id[[1]],
          common_name = species_data$common_name[[1]],
          .before = 1
        )
    })

  selections <- validation |>
    filter(model_id %in% c("M2", "M3")) |>
    group_by(avibase_id, common_name) |>
    arrange(heldout_log_loss, .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    select(avibase_id, common_name, selected_model_id = model_id,
           selected_model = model, selected_heldout_log_loss = heldout_log_loss)

  species_groups <- data |> group_split(avibase_id, common_name)
  names(species_groups) <- map_chr(species_groups, ~ .x$avibase_id[[1]])
  fitted <- map(species_groups, function(species_data) {
      species_id <- species_data$avibase_id[[1]]
      selected_id <- selections$selected_model_id[selections$avibase_id == species_id]
      selected_rhs <- model_specifications$formula_rhs[model_specifications$model_id == selected_id]
      annual_rhs <- sub(
        "splines::ns\\(season_centered, df = 2\\)",
        "s(season_centered, k = 5) + s(season_factor, bs = 're')",
        selected_rhs
      )
      linear_rhs <- sub(
        "splines::ns\\(season_centered, df = 2\\)",
        "I(season_centered / 10)", selected_rhs
      )
      annual_data <- species_data |>
        mutate(season_factor = factor(season))
      annual_weights <- annual_data[[weight]]
      list(
        smooth = fit_quasibinomial(species_data, response, weight, selected_rhs),
        annual = mgcv::gam(
          as.formula(paste(response, "~", annual_rhs)),
          data = annual_data,
          weights = annual_weights,
          family = quasibinomial(),
          method = "REML"
        ),
        linear = fit_quasibinomial(species_data, response, weight, linear_rhs),
        model_id = selected_id,
        formula_rhs = selected_rhs
      )
    })

  diagnostics <- imap_dfr(fitted, function(bundle, species_id) {
    model <- bundle$smooth
    species_data <- data |> filter(avibase_id == species_id)
    tibble(
      avibase_id = species_id,
      common_name = species_data$common_name[[1]],
      selected_model_id = bundle$model_id,
      n_dates = nrow(species_data),
      n_seasons = n_distinct(species_data$season),
      n_individuals = sum(species_data$n_total),
      effective_sample_size = sum(species_data[[weight]]),
      residual_deviance = deviance(model),
      residual_df = df.residual(model),
      dispersion = summary(model)$dispersion,
      converged = model$converged
    )
  })

  trend_summary <- imap_dfr(fitted, function(bundle, species_id) {
    model <- bundle$linear
    term <- "I(season_centered/10)"
    estimate <- coef(model)[term]
    standard_error <- sqrt(vcov(model)[term, term])
    species_data <- data |> filter(avibase_id == species_id)
    tibble(
      avibase_id = species_id,
      common_name = species_data$common_name[[1]],
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

  reference_days <- day_support |>
    filter(.data$analysis == !!analysis, in_reference_window) |>
    pull(day_of_season)

  annual_ratios <- imap_dfr(fitted, function(bundle, species_id) {
    species_data <- data |> filter(avibase_id == species_id)
    model <- bundle$annual
    seasons <- sort(unique(species_data$season))
    grid <- expand_grid(season = seasons, day_of_season = reference_days) |>
      mutate(
        season_factor = factor(season, levels = levels(model$model$season_factor)),
        season_centered = season - median(data$season),
        moon_distance_from_new_moon_z = 0,
        mist_probability_light_patchy = mean(species_data$mist_probability_light_patchy),
        mist_probability_good = mean(species_data$mist_probability_good),
        era5_rain_log_z = 0,
        wind_speed_10m_mean_ms_z = 0,
        temperature_2m_mean_c_z = 0,
        surface_pressure_mean_hpa_z = 0
      )
    design <- predict(model, newdata = grid, type = "lpmatrix")
    set.seed(404)
    coefficient_draws <- MASS::mvrnorm(400, coef(model), vcov(model))
    probability_draws <- plogis(design %*% t(coefficient_draws))
    point_probability <- plogis(drop(design %*% coef(model)))
    annual_estimates <- grid |>
      mutate(point_probability = point_probability, row_id = row_number()) |>
      group_by(season) |>
      summarise(
        estimate = mean(point_probability),
        row_ids = list(row_id),
        .groups = "drop"
      ) |>
      rowwise() |>
      mutate(
        draws = list(colMeans(probability_draws[unlist(row_ids), , drop = FALSE])),
        lower = quantile(unlist(draws), 0.025),
        upper = quantile(unlist(draws), 0.975)
      ) |>
      ungroup() |>
      transmute(
        avibase_id = species_id,
        common_name = species_data$common_name[[1]],
        season, estimate, lower, upper
      )
    annual_sample <- species_data |>
      group_by(season) |>
      summarise(n_individuals = sum(n_total), .groups = "drop")
    left_join(annual_estimates, annual_sample, by = "season")
  })

  # Predict composition across the central observed passage window.
  phenology <- imap_dfr(fitted, function(bundle, species_id) {
    species_data <- data |> filter(avibase_id == species_id)
    model <- bundle$smooth
    day_values <- seq(
      floor(quantile(species_data$day_of_season, 0.02)),
      ceiling(quantile(species_data$day_of_season, 0.98)),
      length.out = 100
    )
    grid <- tibble(
      day_of_season = day_values,
      season_centered = 0,
      moon_distance_from_new_moon_z = 0,
      mist_probability_light_patchy = mean(species_data$mist_probability_light_patchy),
      mist_probability_good = mean(species_data$mist_probability_good),
      era5_rain_log_z = 0,
      wind_speed_10m_mean_ms_z = 0,
      temperature_2m_mean_c_z = 0,
      surface_pressure_mean_hpa_z = 0
    )
    prediction <- predict(model, newdata = grid, type = "link", se.fit = TRUE)
    tibble(
      avibase_id = species_id,
      common_name = species_data$common_name[[1]],
      day_of_season = day_values,
      estimate = plogis(prediction$fit),
      lower = plogis(prediction$fit - 1.96 * prediction$se.fit),
      upper = plogis(prediction$fit + 1.96 * prediction$se.fit)
    )
  })

  timing_summary <- data |>
    group_by(avibase_id, common_name) |>
    summarise(
      first_group_median_day = weighted_quantile(day_of_season, .data[[count_success]], 0.5),
      second_group_median_day = weighted_quantile(day_of_season, .data[[count_failure]], 0.5),
      median_day_difference = first_group_median_day - second_group_median_day,
      .groups = "drop"
    )

  list(
    models = fitted,
    validation = validation,
    selections = selections,
    diagnostics = diagnostics,
    trends = trend_summary,
    annual = annual_ratios,
    phenology = phenology,
    timing = timing_summary
  )
}

age_results <- fit_analysis(
  age_data, "age", "first_year_fraction", "n_first_year", "n_adult", "effective_n"
)
sex_results <- fit_analysis(
  sex_data, "sex", "male_fraction", "n_male", "n_female", "effective_n"
)

# Write outputs ----------------------------------------------------------

saveRDS(age_results$models, file.path(model_dir, "age_composition_models.rds"))
saveRDS(sex_results$models, file.path(model_dir, "adult_sex_composition_models.rds"))

write_csv(model_specifications, file.path(table_dir, "demographic_model_specifications.csv"))
write_csv(age_results$validation, file.path(table_dir, "age_model_validation.csv"))
write_csv(sex_results$validation, file.path(table_dir, "adult_sex_model_validation.csv"))
write_csv(age_results$selections, file.path(table_dir, "age_model_selections.csv"))
write_csv(sex_results$selections, file.path(table_dir, "adult_sex_model_selections.csv"))
write_csv(age_results$diagnostics, file.path(table_dir, "age_model_diagnostics.csv"))
write_csv(sex_results$diagnostics, file.path(table_dir, "adult_sex_model_diagnostics.csv"))
write_csv(age_results$trends, file.path(table_dir, "age_decadal_trends.csv"))
write_csv(sex_results$trends, file.path(table_dir, "adult_sex_decadal_trends.csv"))
write_csv(age_results$annual, file.path(table_dir, "standardized_annual_first_year_proportions.csv"))
write_csv(sex_results$annual, file.path(table_dir, "standardized_annual_adult_male_proportions.csv"))
write_csv(age_results$phenology, file.path(table_dir, "age_phenology_predictions.csv"))
write_csv(sex_results$phenology, file.path(table_dir, "adult_sex_phenology_predictions.csv"))
write_csv(age_results$timing, file.path(table_dir, "raw_age_timing_summary.csv"))
write_csv(sex_results$timing, file.path(table_dir, "raw_adult_sex_timing_summary.csv"))

cli_alert_success("Fitted standardized demographic models and annual composition indices.")
