library(dplyr)
library(tidyr)
library(readr)
library(mgcv)
library(nnet)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
source(file.path(project_dir, "scripts", "helpers", "mist_state.R"))

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
model_data_path <- file.path(analysis_dir, "model_data", "positive_catch_model_data.csv")
daily_coverage_path <- file.path(paths$curated_dir, "daily_coverage.csv")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read model data ---------------------------------------------------------

cli_h1("Fit adjusted annual positive-catch index")

model_data <- read_csv(
  model_data_path,
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(season_factor = factor(season), cloud_base_height_km = cloud_base_height_mean_m / 1000)

mist_levels <- c("none", "light_patchy", "good")
n_imputations <- 20L
n_validation_imputations <- 5L

calibration_data <- read_csv(
  daily_coverage_path,
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(
    cloud_base_height_km = cloud_base_height_mean_m / 1000,
    mist_state = factor(mist_observation, levels = mist_levels)
  ) |>
  filter(
    mist_observation %in% mist_levels,
    complete.cases(total_cloud_cover_mean, cloud_base_height_km, relative_humidity_mean_pct, wind_u_10m_mean_ms)
  )

mist_calibration_formula <- mist_state ~ total_cloud_cover_mean + cloud_base_height_km +
  relative_humidity_mean_pct + wind_u_10m_mean_ms

# Define nested candidate models -----------------------------------------

model_labels <- c(
  "M1 Seasonal timing",
  "M2 + Mist and rain",
  "M3 + Moon",
  "M4 + Local weather"
)

model_additions <- c(
  "Within-season timing",
  "Latent mist state + rainfall",
  "Moon distance",
  "Wind speed + temperature + pressure"
)

validation_formulas <- list(
  total_birds_ringed ~ s(season_day, k = 12),
  total_birds_ringed ~ s(season_day, k = 12) +
    mist_state +
    s(era5_rain_log, k = 6),
  total_birds_ringed ~ s(season_day, k = 12) +
    s(moon_distance_from_new_moon, k = 6) +
    mist_state +
    s(era5_rain_log, k = 6),
  total_birds_ringed ~ s(season_day, k = 12) +
    s(moon_distance_from_new_moon, k = 6) +
    mist_state +
    s(era5_rain_log, k = 6) +
    s(wind_speed_10m_mean_ms, k = 6) +
    s(temperature_2m_mean_c, k = 6) +
    s(surface_pressure_mean_hpa, k = 6)
)

fitted_formulas <- lapply(validation_formulas, update, . ~ season_factor + .)

# Compare covariate sets with seasons held out ----------------------------

score_season_block_cv <- function(data, formula, model_name, model_number, folds) {
  prediction <- rep(NA_real_, nrow(data))
  baseline_prediction <- rep(NA_real_, nrow(data))

  for (fold in sort(unique(folds))) {
    training <- data[folds != fold, ]
    testing <- data[folds == fold, ]

    if (model_number < 2) {
      fit <- bam(formula, data = training, family = nb(), method = "fREML", discrete = TRUE)
      prediction[folds == fold] <- predict(fit, newdata = testing, type = "response")
    } else {
      mist_model <- multinom(
        mist_calibration_formula,
        data = filter(calibration_data, !season %in% testing$season),
        trace = FALSE
      )
      training_probability <- predict_mist_probabilities(mist_model, training)
      testing_probability <- predict_mist_probabilities(mist_model, testing)
      imputed_prediction <- matrix(NA_real_, nrow(testing), n_validation_imputations)

      for (imputation in seq_len(n_validation_imputations)) {
        training_imputed <- draw_mist_state(training, training_probability)
        testing_imputed <- draw_mist_state(testing, testing_probability)
        fit <- bam(formula, data = training_imputed, family = nb(), method = "fREML", discrete = TRUE)
        imputed_prediction[, imputation] <- predict(fit, newdata = testing_imputed, type = "response")
      }
      prediction[folds == fold] <- rowMeans(imputed_prediction)
    }
    baseline_prediction[folds == fold] <- mean(data$total_birds_ringed[folds != fold])
  }

  poisson_deviance <- function(observed, expected) {
    2 * sum(if_else(observed == 0, expected, observed * log(observed / expected) - observed + expected))
  }

  tibble(
    model = model_name,
    cv_deviance_reduction = 1 - poisson_deviance(data$total_birds_ringed, prediction) /
      poisson_deviance(data$total_birds_ringed, baseline_prediction),
    log_rmse = sqrt(mean((log1p(data$total_birds_ringed) - log1p(prediction))^2)),
    log_mae = mean(abs(log1p(data$total_birds_ringed) - log1p(prediction)))
  )
}

set.seed(73)
season_folds <- tibble(season = sort(unique(model_data$season))) |>
  mutate(fold = sample(rep(1:5, length.out = n())))

model_data <- model_data |>
  left_join(season_folds, by = "season")

season_block_validation <- bind_rows(lapply(seq_along(model_labels), function(i) {
  score_season_block_cv(model_data, validation_formulas[[i]], model_labels[[i]], i, model_data$fold)
}))

# Fit candidate models across mist-state imputations ---------------------

set.seed(314)
imputed_data <- lapply(seq_len(n_imputations), function(i) draw_mist_state(model_data))

candidate_models <- lapply(imputed_data, function(data) {
  lapply(fitted_formulas, function(formula) gam(formula, data = data, family = nb(), method = "ML"))
})

full_weather_models <- lapply(imputed_data, function(data) {
  bam(fitted_formulas[[4]], data = data, family = nb(), method = "fREML", discrete = TRUE)
})

candidate_metrics <- bind_rows(lapply(seq_len(n_imputations), function(imputation) {
  bind_rows(lapply(seq_along(model_labels), function(model_number) {
    fit <- candidate_models[[imputation]][[model_number]]
    tibble(
      imputation,
      model = model_labels[[model_number]],
      n_parameters = length(coef(fit)),
      aic = AIC(fit),
      deviance_explained = summary(fit)$dev.expl
    )
  }))
}))

model_comparison <- candidate_metrics |>
  group_by(model) |>
  summarise(
    n_mist_imputations = n(),
    n_parameters = mean(n_parameters),
    aic = mean(aic),
    imputation_deviance_sd = sd(deviance_explained, na.rm = TRUE),
    deviance_explained = mean(deviance_explained),
    .groups = "drop"
  ) |>
  mutate(model = factor(model, levels = model_labels)) |>
  arrange(model) |>
  mutate(model = as.character(model), added_terms = model_additions) |>
  left_join(season_block_validation, by = "model") |>
  mutate(
    deviance_gain_pp = 100 * (deviance_explained - lag(deviance_explained, default = 0)),
    cv_deviance_gain_pp = 100 * (cv_deviance_reduction - lag(cv_deviance_reduction, default = 0)),
    delta_aic = aic - min(aic),
    primary_index_model = model == model_labels[[4]]
  )

# Write model outputs -----------------------------------------------------

saveRDS(full_weather_models, file.path(model_dir, "adjusted_positive_catch_gam_mi.rds"))
write_csv(model_comparison, file.path(table_dir, "annual_index_model_comparison.csv"))
write_csv(season_block_validation, file.path(table_dir, "season_block_model_validation.csv"))

cli_alert_success("Wrote primary daily-count models to {analysis_dir}")
