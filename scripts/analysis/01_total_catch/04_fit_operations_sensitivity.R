library(dplyr)
library(readr)
library(stringr)
library(mgcv)
library(nnet)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
source(file.path(project_dir, "scripts", "helpers", "mist_state.R"))

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
model_data_path <- file.path(analysis_dir, "model_data", "documented_operation_model_data.csv")
calibration_data_path <- file.path(analysis_dir, "model_data", "positive_catch_model_data.csv")
daily_coverage_path <- file.path(paths$curated_dir, "daily_coverage.csv")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read positive-catch dates in stable post-transition configurations ------

cli_h1("Compare configuration and playback models")

model_data <- read_csv(
  model_data_path,
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(
    bush_period = factor(bush_net_configuration, levels = c("back_bush", "front_bush")),
    playback_used = playback_nocturnal_observed
  )

daily_coverage <- read_csv(
  daily_coverage_path,
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)

calibration_data <- read_csv(
  calibration_data_path,
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(
    mist_state = factor(mist_observation, levels = c("none", "light_patchy", "good")),
    cloud_base_height_km = cloud_base_height_mean_m / 1000
  ) |>
  filter(
    mist_observation %in% c("none", "light_patchy", "good"),
    complete.cases(total_cloud_cover_mean, cloud_base_height_km, relative_humidity_mean_pct, wind_u_10m_mean_ms)
  )

operation_coverage <- tribble(
  ~covariate, ~available,
  "Minimum recorded team size", !is.na(daily_coverage$djp_team_size_minimum),
  "Nocturnal playback", !is.na(daily_coverage$playback_nocturnal_observed),
  "Operated net sites", !is.na(daily_coverage$net_sites_observed)
) |>
  rowwise() |>
  summarise(
    covariate,
    n_positive_dates_available = sum(available & daily_coverage$ringing_happened),
    n_positive_dates = sum(daily_coverage$ringing_happened),
    coverage = n_positive_dates_available / n_positive_dates,
    n_seasons_available = n_distinct(daily_coverage$season[available & daily_coverage$ringing_happened]),
    .groups = "drop"
  )

# Draw the unified mist state --------------------------------------------

mist_levels <- c("none", "light_patchy", "good")

mist_calibration_formula <- mist_state ~ total_cloud_cover_mean + cloud_base_height_km +
  relative_humidity_mean_pct + wind_u_10m_mean_ms

# Specify M4-S, M5 and M6 on identical dates -----------------------------

m4_formula <- total_birds_ringed ~
  s(season_day, k = 12) +
  s(moon_distance_from_new_moon, k = 6) +
  mist_state +
  s(era5_rain_log, k = 6) +
  s(wind_speed_10m_mean_ms, k = 6) +
  s(temperature_2m_mean_c, k = 6) +
  s(surface_pressure_mean_hpa, k = 6)

m4_trend_formula <- update(m4_formula, . ~ . + s(season, k = 10))
m5_formula <- update(m4_trend_formula, . ~ . + bush_period)
m6_formula <- update(m5_formula, . ~ . + playback_used)

validation_formulas <- list(
  `M4-S Smooth-year baseline` = m4_trend_formula,
  `M5 + bush configuration` = m5_formula,
  `M6 + playback` = m6_formula
)

# Hold out complete seasons on an identical subset -----------------------

set.seed(73)
season_folds <- tibble(season = sample(unique(model_data$season))) |>
  mutate(fold = rep(1:5, length.out = n()))
model_data <- left_join(model_data, season_folds, by = "season")

poisson_deviance <- function(observed, expected) {
  2 * sum(if_else(observed == 0, expected, observed * log(observed / expected) - observed + expected))
}

score_model <- function(formula) {
  set.seed(731)
  prediction <- matrix(NA_real_, nrow(model_data), 5)
  baseline <- rep(NA_real_, nrow(model_data))

  for (fold in 1:5) {
    training <- filter(model_data, fold != .env$fold)
    testing <- filter(model_data, fold == .env$fold)
    baseline[model_data$fold == fold] <- mean(training$total_birds_ringed)
    mist_model <- multinom(
      mist_calibration_formula,
      data = filter(calibration_data, !season %in% testing$season),
      trace = FALSE
    )
    training_probability <- predict_mist_probabilities(mist_model, training)
    testing_probability <- predict_mist_probabilities(mist_model, testing)

    for (imputation in 1:5) {
      training_imputed <- draw_mist_state(training, training_probability)
      testing_imputed <- draw_mist_state(testing, testing_probability)
      fit <- bam(formula, data = training_imputed, family = nb(), method = "fREML", discrete = TRUE)
      prediction[model_data$fold == fold, imputation] <- predict(fit, newdata = testing_imputed, type = "response")
    }
  }

  prediction <- rowMeans(prediction)
  tibble(
    cv_deviance_reduction = 1 - poisson_deviance(model_data$total_birds_ringed, prediction) /
      poisson_deviance(model_data$total_birds_ringed, baseline),
    log_rmse = sqrt(mean((log1p(model_data$total_birds_ringed) - log1p(prediction))^2)),
    log_mae = mean(abs(log1p(model_data$total_birds_ringed) - log1p(prediction)))
  )
}

validation <- bind_rows(lapply(validation_formulas, score_model), .id = "model")

operation_contribution <- validation |>
  mutate(
    covariate_block = c("Weather and year", "Bush configuration", "Nocturnal playback"),
    heldout_log_rmse_change_from_previous = c(NA_real_, diff(log_rmse)),
    heldout_deviance_change_from_previous = c(NA_real_, diff(cv_deviance_reduction))
  )

# Fit the smooth-year models ----------------------------------------------

set.seed(315)
imputed_data <- lapply(1:20, function(imputation) draw_mist_state(model_data))
fitted_formulas <- list(
  `M4-S Smooth-year baseline` = m4_trend_formula,
  `M5 + bush configuration` = m5_formula,
  `M6 + playback` = m6_formula
)

fitted_models <- lapply(fitted_formulas, function(formula) {
  lapply(imputed_data, function(data) bam(formula, data = data, family = nb(), method = "fREML", discrete = TRUE))
})

model_comparison <- bind_rows(lapply(names(fitted_models), function(model_name) {
  models <- fitted_models[[model_name]]
  tibble(
    model = model_name,
    n_dates = nrow(model_data),
    n_seasons = n_distinct(model_data$season),
    n_parameters = mean(vapply(models, function(model) length(coef(model)), numeric(1))),
    aic = mean(vapply(models, AIC, numeric(1))),
    deviance_explained = mean(vapply(models, function(model) summary(model)$dev.expl, numeric(1)))
  )
})) |>
  left_join(validation, by = "model") |>
  mutate(delta_aic = aic - min(aic))

# Write model outputs -----------------------------------------------------

saveRDS(fitted_models, file.path(model_dir, "configuration_playback_gam_mi.rds"))
write_csv(model_comparison, file.path(table_dir, "documented_operations_model_comparison.csv"))
write_csv(validation, file.path(table_dir, "documented_operations_season_block_validation.csv"))
write_csv(operation_contribution, file.path(table_dir, "documented_operations_covariate_contribution.csv"))
write_csv(operation_coverage, file.path(table_dir, "documented_operations_covariate_coverage.csv"))
unlink(file.path(model_dir, c("documented_operations_gam_mi.rds",
  "documented_operations_by_bush_period_gam_mi.rds")))
unlink(file.path(table_dir, "documented_operations_by_bush_period_comparison.csv"))

cli_alert_success("Wrote documented-operations sensitivity models to {analysis_dir}")
