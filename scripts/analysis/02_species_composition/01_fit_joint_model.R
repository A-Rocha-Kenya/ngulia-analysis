library(dplyr)
library(tidyr)
library(readr)
library(nnet)
library(splines)
library(purrr)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "02_species_composition")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read and reshape model data --------------------------------------------

cli_h1("Fit joint species-composition model")

model_data <- read_csv(
  file.path(model_data_dir, "composition_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)

categories <- read_csv(
  file.path(model_data_dir, "composition_categories.csv"),
  show_col_types = FALSE
) |>
  arrange(category_order)

predictor_columns <- c(
  "ringing_date", "season", "season_centered", "day_of_season",
  "moon_distance_from_new_moon_z",
  "mist_probability_light_patchy", "mist_probability_good",
  "era5_rain_log_z", "wind_speed_10m_mean_ms_z",
  "temperature_2m_mean_c_z", "surface_pressure_mean_hpa_z"
)

daily_predictors <- model_data |>
  distinct(across(all_of(predictor_columns))) |>
  arrange(ringing_date)

count_data <- model_data |>
  select(ringing_date, common_name, count) |>
  pivot_wider(names_from = common_name, values_from = count) |>
  arrange(ringing_date)

stopifnot(identical(daily_predictors$ringing_date, count_data$ringing_date))

count_matrix <- as.matrix(count_data[, categories$common_name])
storage.mode(count_matrix) <- "double"
stopifnot(all(rowSums(count_matrix) > 0))

# Define staged candidate models -----------------------------------------

# Restrict flexibility to broad biological patterns rather than allowing
# small within-season or multi-year fluctuations to drive interpretation.
season_spline_df <- 2L
phenology_spline_df <- 3L

model_specifications <- tribble(
  ~model_id, ~model, ~added_information, ~formula_rhs,
  "M0", "M0 Constant composition", "No temporal adjustment", "1",
  "M1", "M1 Broad long-term season", "Broad long-term season pattern", paste0(
    "splines::ns(season_centered, df = ", season_spline_df, ")"
  ),
  "M2", "M2 + Within-season timing", "Species-specific day-of-season phenology", paste(
    paste0("splines::ns(season_centered, df = ", season_spline_df, ")"),
    paste0("splines::ns(day_of_season, df = ", phenology_spline_df, ")"),
    sep = " + "
  ),
  "M3", "M3 + Daily conditions", "Moon, mist, rain, wind, temperature and pressure", paste(
    paste0("splines::ns(season_centered, df = ", season_spline_df, ")"),
    paste0("splines::ns(day_of_season, df = ", phenology_spline_df, ")"),
    "splines::ns(moon_distance_from_new_moon_z, df = 3)",
    "mist_probability_light_patchy + mist_probability_good",
    "era5_rain_log_z + wind_speed_10m_mean_ms_z",
    "temperature_2m_mean_c_z + surface_pressure_mean_hpa_z",
    sep = " + "
  )
)

fit_multinomial <- function(counts, data, formula_rhs, hessian = FALSE) {
  model_formula <- as.formula(paste("counts ~", formula_rhs), env = environment())
  multinom(
    model_formula,
    data = data,
    trace = FALSE,
    Hess = hessian,
    maxit = 1500,
    MaxNWts = 20000
  )
}

multinomial_scores <- function(observed, probability) {
  probability <- pmax(probability[, colnames(observed), drop = FALSE], 1e-12)
  probability <- probability / rowSums(probability)
  daily_log_loss <- -rowSums(observed * log(probability)) / rowSums(observed)
  tibble(
    bird_weighted_log_loss = -sum(observed * log(probability)) / sum(observed),
    mean_daily_log_loss = mean(daily_log_loss),
    median_daily_log_loss = median(daily_log_loss)
  )
}

# Compare candidates with complete seasons held out ----------------------

set.seed(2209)
season_folds <- tibble(season = sort(unique(daily_predictors$season))) |>
  mutate(fold = sample(rep(seq_len(5), length.out = n())))

daily_predictors <- daily_predictors |>
  left_join(season_folds, by = "season")

score_candidate <- function(formula_rhs) {
  heldout_probability <- matrix(
    NA_real_,
    nrow = nrow(count_matrix),
    ncol = ncol(count_matrix),
    dimnames = list(NULL, colnames(count_matrix))
  )

  for (fold in sort(unique(daily_predictors$fold))) {
    training_rows <- daily_predictors$fold != fold
    testing_rows <- daily_predictors$fold == fold
    fit <- fit_multinomial(
      count_matrix[training_rows, , drop = FALSE],
      daily_predictors[training_rows, ],
      formula_rhs
    )
    heldout_probability[testing_rows, ] <- predict(
      fit,
      newdata = daily_predictors[testing_rows, ],
      type = "probs"
    )[, colnames(count_matrix), drop = FALSE]
  }

  multinomial_scores(count_matrix, heldout_probability)
}

cli_h2("Season-block model comparison")

season_block_validation <- model_specifications |>
  mutate(scores = map(formula_rhs, score_candidate)) |>
  select(-formula_rhs) |>
  unnest(scores)

# Fit full-data candidates and select one primary model ------------------

cli_h2("Full-data candidate fits")

candidate_models <- map(
  model_specifications$formula_rhs,
  ~ fit_multinomial(count_matrix, daily_predictors, .x)
)
names(candidate_models) <- model_specifications$model_id

candidate_metrics <- model_specifications |>
  mutate(
    n_parameters = map_int(candidate_models, ~ length(coef(.x))),
    aic = map_dbl(candidate_models, AIC),
    convergence = map_int(candidate_models, "convergence")
  ) |>
  left_join(season_block_validation, by = c("model_id", "model", "added_information"))

# M2 is the minimum acceptable model because timing adjustment is required.
# M3 becomes primary only when it improves held-out bird-weighted log loss.
eligible_primary <- candidate_metrics |>
  filter(model_id %in% c("M2", "M3")) |>
  arrange(bird_weighted_log_loss, mean_daily_log_loss)

primary_model_id <- eligible_primary$model_id[[1]]
primary_specification <- model_specifications |>
  filter(model_id == primary_model_id)

primary_model <- fit_multinomial(
  count_matrix,
  daily_predictors,
  primary_specification$formula_rhs,
  hessian = TRUE
)

primary_probability <- predict(primary_model, newdata = daily_predictors, type = "probs")
primary_probability <- primary_probability[, colnames(count_matrix), drop = FALSE]
primary_scores <- multinomial_scores(count_matrix, primary_probability)

expected_counts <- rowSums(count_matrix) * primary_probability
n_residual_cells <- nrow(count_matrix) * (ncol(count_matrix) - 1)
pearson_dispersion <- sum((count_matrix - expected_counts)^2 / pmax(expected_counts, 1e-8)) /
  (n_residual_cells - length(coef(primary_model)))

untempered_diagnostics <- tibble(
  primary_model_id,
  primary_model_label = primary_specification$model,
  n_dates = nrow(count_matrix),
  n_seasons = n_distinct(daily_predictors$season),
  n_focal_species = sum(categories$is_focal_species),
  n_model_categories = nrow(categories),
  n_birds = sum(count_matrix),
  n_parameters = length(coef(primary_model)),
  convergence = primary_model$convergence,
  aic = AIC(primary_model),
  pearson_dispersion,
  fitted_bird_weighted_log_loss = primary_scores$bird_weighted_log_loss,
  fitted_mean_daily_log_loss = primary_scores$mean_daily_log_loss,
  max_probability_sum_error = max(abs(rowSums(primary_probability) - 1))
)

untempered_bundle <- list(
  model = primary_model,
  model_id = primary_model_id,
  specification = primary_specification,
  model_specifications = model_specifications,
  categories = categories,
  predictor_columns = predictor_columns,
  category_order = colnames(count_matrix),
  season_spline_df = season_spline_df,
  phenology_spline_df = phenology_spline_df
)

# Write model outputs -----------------------------------------------------

saveRDS(untempered_bundle, file.path(model_dir, "joint_multinomial_untempered.rds"))
saveRDS(candidate_models, file.path(model_dir, "joint_multinomial_candidates.rds"))
write_csv(candidate_metrics, file.path(table_dir, "composition_model_comparison.csv"))
write_csv(season_folds, file.path(table_dir, "season_validation_folds.csv"))
write_csv(untempered_diagnostics, file.path(table_dir, "untempered_model_diagnostics.csv"))

cli_alert_success("Selected {primary_specification$model} as the composition mean structure.")
cli_alert_info(
  "Held-out bird-weighted log loss: {signif(eligible_primary$bird_weighted_log_loss[[1]], 5)}."
)
