library(dplyr)
library(readr)
library(mgcv)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
model_data_path <- file.path(analysis_dir, "model_data", "positive_catch_model_data.csv")
primary_model_path <- file.path(analysis_dir, "models", "adjusted_positive_catch_gam_mi.rds")
table_dir <- file.path(analysis_dir, "tables")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read model data and imputed primary models -----------------------------

cli_h1("Validate adjusted positive-catch models")

model_data <- read_csv(
  model_data_path,
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)

primary_models <- readRDS(primary_model_path)

# Pool residual and smooth-basis diagnostics -----------------------------

model_diagnostics_by_imputation <- bind_rows(lapply(seq_along(primary_models), function(imputation) {
  model <- primary_models[[imputation]]
  pearson_residuals <- residuals(model, type = "pearson")
  deviance_residuals <- residuals(model, type = "deviance")
  tibble(
    imputation,
    pearson_dispersion = sum(pearson_residuals^2) / model$df.residual,
    deviance_residual_median = median(deviance_residuals),
    deviance_residual_q025 = quantile(deviance_residuals, 0.025),
    deviance_residual_q975 = quantile(deviance_residuals, 0.975),
    correlation_fitted_observed_log = cor(log1p(fitted(model)), log1p(model$y))
  )
}))

model_diagnostics <- model_diagnostics_by_imputation |>
  summarise(
    n_dates = nrow(model_data),
    n_seasons = n_distinct(model_data$season),
    n_mist_imputations = n(),
    across(-imputation, mean)
  )

smooth_basis_check <- bind_rows(lapply(seq_along(primary_models), function(imputation) {
  as.data.frame(k.check(primary_models[[imputation]])) |>
    tibble::rownames_to_column("smooth") |>
    as_tibble() |>
    mutate(imputation)
}))

smooth_basis_summary <- smooth_basis_check |>
  group_by(smooth) |>
  summarise(
    k_prime = median(`k'`),
    edf = mean(edf),
    k_index = mean(`k-index`),
    minimum_p_value = min(`p-value`),
    .groups = "drop"
  )

# Write outputs -----------------------------------------------------------

write_csv(model_diagnostics, file.path(table_dir, "primary_model_diagnostics.csv"))
write_csv(model_diagnostics_by_imputation, file.path(table_dir, "primary_model_diagnostics_by_imputation.csv"))
write_csv(smooth_basis_summary, file.path(table_dir, "primary_model_smooth_basis_check.csv"))

cli_alert_success("Wrote total-catch model validation to {table_dir}")
