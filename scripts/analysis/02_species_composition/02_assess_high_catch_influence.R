library(dplyr)
library(tidyr)
library(readr)
library(nnet)
library(splines)
library(cli)

analysis_dir <- "outputs/analysis/02_species_composition"
model_data <- read_csv(
  file.path(analysis_dir, "model_data", "composition_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
categories <- read_csv(
  file.path(analysis_dir, "model_data", "composition_categories.csv"),
  show_col_types = FALSE
) |>
  arrange(category_order)
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

cli_h1("Assess high-catch influence and calibrate daily information")

bundle <- readRDS(file.path(model_dir, "joint_multinomial_untempered.rds"))
formula_rhs <- bundle$specification$formula_rhs

predictor_columns <- c(
  "ringing_date", "season", "season_centered", "day_of_season",
  "moon_distance_from_new_moon_z", "mist_probability_light_patchy",
  "mist_probability_good", "era5_rain_log_z", "wind_speed_10m_mean_ms_z",
  "temperature_2m_mean_c_z", "surface_pressure_mean_hpa_z"
)
predictors <- model_data |>
  distinct(across(all_of(predictor_columns))) |>
  arrange(ringing_date) |>
  left_join(
    read_csv(file.path(analysis_dir, "tables", "season_validation_folds.csv"), show_col_types = FALSE),
    by = "season"
  )
count_data <- model_data |>
  select(ringing_date, common_name, count) |>
  pivot_wider(names_from = common_name, values_from = count) |>
  arrange(ringing_date)
counts <- as.matrix(count_data[, categories$common_name])
storage.mode(counts) <- "double"
daily_total <- rowSums(counts)
focal_names <- categories$common_name[categories$is_focal_species]

fit_multinomial <- function(response, data, hessian = FALSE, formula_rhs_override = formula_rhs) {
  model_formula <- as.formula(paste("response ~", formula_rhs_override), env = environment())
  multinom(
    model_formula,
    data = data,
    trace = FALSE,
    Hess = hessian,
    maxit = 1500,
    MaxNWts = 20000
  )
}

probability_scores <- function(observed, probability) {
  probability <- pmax(probability[, colnames(observed), drop = FALSE], 1e-12)
  daily_log_loss <- -rowSums(observed * log(probability)) / rowSums(observed)
  tibble(
    bird_weighted_log_loss = sum(daily_log_loss * rowSums(observed)) / sum(observed),
    mean_daily_log_loss = mean(daily_log_loss),
    high_catch_log_loss = weighted.mean(
      daily_log_loss,
      pmin(rowSums(observed), quantile(rowSums(observed), 0.95))
    )
  )
}

softmax <- function(linear_predictor) {
  row_maximum <- apply(linear_predictor, 1, max)
  exponentiated <- exp(linear_predictor - row_maximum)
  exponentiated / rowSums(exponentiated)
}

predict_from_coefficients <- function(model_matrix, coefficients, category_names) {
  linear_predictor <- cbind(0, model_matrix %*% t(coefficients))
  colnames(linear_predictor) <- category_names
  softmax(linear_predictor)
}

dm_objective <- function(parameters, model_matrix, response, include_constant = FALSE) {
  n_categories <- ncol(response)
  n_coefficients <- ncol(model_matrix)
  coefficient_matrix <- matrix(
    parameters[-length(parameters)],
    nrow = n_categories - 1,
    byrow = TRUE
  )
  concentration <- exp(parameters[[length(parameters)]])
  probability <- predict_from_coefficients(
    model_matrix,
    coefficient_matrix,
    colnames(response)
  )
  probability <- pmax(probability, 1e-15)
  probability <- probability / rowSums(probability)
  alpha <- concentration * probability
  total <- rowSums(response)
  log_likelihood <- sum(
    lgamma(concentration) - lgamma(concentration + total) +
      rowSums(lgamma(alpha + response) - lgamma(alpha))
  )
  if (include_constant) {
    log_likelihood <- log_likelihood + sum(lgamma(total + 1) - rowSums(lgamma(response + 1)))
  }

  digamma_difference <- digamma(alpha + response) - digamma(alpha)
  average_digamma_difference <- rowSums(probability * digamma_difference)
  eta_gradient <- concentration * probability *
    (digamma_difference - average_digamma_difference)
  coefficient_gradient <- t(eta_gradient[, -1, drop = FALSE]) %*% model_matrix
  concentration_gradient <- concentration * sum(
    digamma(concentration) - digamma(concentration + total) +
      rowSums(probability * digamma_difference)
  )

  list(
    value = -log_likelihood,
    gradient = -c(as.vector(t(coefficient_gradient)), concentration_gradient),
    probability = probability,
    concentration = concentration
  )
}

fit_dirichlet_multinomial <- function(response, data) {
  multinomial_fit <- fit_multinomial(response, data)
  model_matrix <- model.matrix(
    delete.response(terms(multinomial_fit)),
    data = data,
    contrasts.arg = multinomial_fit$contrasts,
    xlev = multinomial_fit$xlevels
  )
  initial_parameters <- c(as.vector(t(coef(multinomial_fit))), log(50))
  optimization <- optim(
    initial_parameters,
    fn = function(parameters) dm_objective(parameters, model_matrix, response)$value,
    gr = function(parameters) dm_objective(parameters, model_matrix, response)$gradient,
    method = "L-BFGS-B",
    lower = c(rep(-30, length(initial_parameters) - 1), log(0.1)),
    upper = c(rep(30, length(initial_parameters) - 1), log(1e5)),
    control = list(maxit = 1000, factr = 1e7, pgtol = 1e-5)
  )
  coefficient_matrix <- matrix(
    optimization$par[-length(optimization$par)],
    nrow = ncol(response) - 1,
    byrow = TRUE,
    dimnames = dimnames(coef(multinomial_fit))
  )
  list(
    coefficients = coefficient_matrix,
    concentration = exp(tail(optimization$par, 1)),
    convergence = optimization$convergence,
    value = optimization$value,
    terms = terms(multinomial_fit),
    contrasts = multinomial_fit$contrasts,
    xlevels = multinomial_fit$xlevels
  )
}

predict_dm <- function(model, newdata) {
  model_matrix <- model.matrix(
    delete.response(model$terms),
    data = newdata,
    contrasts.arg = model$contrasts,
    xlev = model$xlevels
  )
  predict_from_coefficients(model_matrix, model$coefficients, categories$common_name)
}

dm_log_score <- function(observed, probability, concentration) {
  alpha <- concentration * probability
  total <- rowSums(observed)
  log_probability <-
    lgamma(total + 1) - rowSums(lgamma(observed + 1)) +
    lgamma(concentration) - lgamma(concentration + total) +
    rowSums(lgamma(alpha + observed) - lgamma(alpha))
  tibble(
    mean_daily_dm_log_score = -mean(log_probability),
    median_daily_dm_log_score = -median(log_probability)
  )
}

cap_95 <- unname(quantile(daily_total, 0.95))
median_total <- median(daily_total)
dm_full_model <- fit_dirichlet_multinomial(counts, predictors)
dm_concentration <- dm_full_model$concentration
effective_counts <- list(
  multinomial = counts,
  capped_95 = counts * pmin(1, cap_95 / daily_total),
  square_root = counts * (sqrt(daily_total * median_total) / daily_total),
  variance_tempered = counts * ((dm_concentration + 1) / (daily_total + dm_concentration)),
  equal_date = counts * (median_total / daily_total)
)

# Held-out whole-season validation ---------------------------------------

heldout_probabilities <- lapply(effective_counts, function(x) {
  matrix(
    NA_real_, nrow = nrow(counts), ncol = ncol(counts),
    dimnames = list(NULL, colnames(counts))
  )
})
dm_heldout_probability <- heldout_probabilities$multinomial
dm_heldout_log_probability <- rep(NA_real_, nrow(counts))
dm_fold_concentration <- numeric()

for (fold in sort(unique(predictors$fold))) {
  training <- predictors$fold != fold
  testing <- predictors$fold == fold
  dm_model <- fit_dirichlet_multinomial(
    counts[training, , drop = FALSE],
    predictors[training, ]
  )
  dm_fold_concentration[as.character(fold)] <- dm_model$concentration

  for (model_id in names(effective_counts)) {
    training_response <- effective_counts[[model_id]][training, , drop = FALSE]
    if (model_id == "variance_tempered") {
      training_total <- daily_total[training]
      training_response <- counts[training, , drop = FALSE] *
        ((dm_model$concentration + 1) / (training_total + dm_model$concentration))
    }
    model <- fit_multinomial(
      training_response,
      predictors[training, ]
    )
    heldout_probabilities[[model_id]][testing, ] <- predict(
      model,
      newdata = predictors[testing, ],
      type = "probs"
    )[, colnames(counts), drop = FALSE]
  }

  dm_probability <- predict_dm(dm_model, predictors[testing, ])
  dm_heldout_probability[testing, ] <- dm_probability
  alpha <- dm_model$concentration * dm_probability
  test_counts <- counts[testing, , drop = FALSE]
  test_total <- rowSums(test_counts)
  dm_heldout_log_probability[testing] <-
    lgamma(test_total + 1) - rowSums(lgamma(test_counts + 1)) +
    lgamma(dm_model$concentration) - lgamma(dm_model$concentration + test_total) +
    rowSums(lgamma(alpha + test_counts) - lgamma(alpha))
}

validation_scores <- bind_rows(lapply(names(heldout_probabilities), function(model_id) {
  probability_scores(counts, heldout_probabilities[[model_id]]) |>
    mutate(model_id = model_id, .before = 1)
})) |>
  bind_rows(
    probability_scores(counts, dm_heldout_probability) |>
      mutate(model_id = "dirichlet_multinomial", .before = 1)
  )

# Confirm the mean structure under the final variance-tempered weighting.
tempered_mean_probabilities <- list(M3 = heldout_probabilities$variance_tempered)
for (model_id in c("M0", "M1", "M2")) {
  heldout_probability <- matrix(
    NA_real_, nrow = nrow(counts), ncol = ncol(counts),
    dimnames = list(NULL, colnames(counts))
  )
  candidate_formula_rhs <- bundle$model_specifications$formula_rhs[
    bundle$model_specifications$model_id == model_id
  ]
  for (fold in sort(unique(predictors$fold))) {
    training <- predictors$fold != fold
    testing <- predictors$fold == fold
    fold_concentration <- dm_fold_concentration[as.character(fold)]
    training_total <- daily_total[training]
    training_response <- counts[training, , drop = FALSE] *
      ((fold_concentration + 1) / (training_total + fold_concentration))
    model <- fit_multinomial(
      training_response,
      predictors[training, ],
      formula_rhs_override = candidate_formula_rhs
    )
    heldout_probability[testing, ] <- predict(
      model,
      newdata = predictors[testing, ],
      type = "probs"
    )[, colnames(counts), drop = FALSE]
  }
  tempered_mean_probabilities[[model_id]] <- heldout_probability
}

tempered_mean_structure_comparison <- bind_rows(lapply(
  bundle$model_specifications$model_id,
  function(model_id) {
    probability_scores(counts, tempered_mean_probabilities[[model_id]]) |>
      mutate(model_id = model_id, .before = 1)
  }
)) |>
  left_join(
    bundle$model_specifications |> select(model_id, model, added_information),
    by = "model_id"
  ) |>
  arrange(match(model_id, bundle$model_specifications$model_id))

selected_tempered_mean <- tempered_mean_structure_comparison |>
  filter(model_id %in% c("M2", "M3")) |>
  arrange(mean_daily_log_loss, bird_weighted_log_loss) |>
  slice(1)
stopifnot(selected_tempered_mean$model_id == bundle$specification$model_id)

multinomial_dm_log_probability <-
  lgamma(daily_total + 1) - rowSums(lgamma(counts + 1)) +
  rowSums(counts * log(pmax(heldout_probabilities$multinomial, 1e-12)))
count_distribution_scores <- tibble(
  model_id = c("multinomial", "dirichlet_multinomial"),
  mean_daily_count_log_score = c(
    -mean(multinomial_dm_log_probability),
    -mean(dm_heldout_log_probability)
  ),
  median_daily_count_log_score = c(
    -median(multinomial_dm_log_probability),
    -median(dm_heldout_log_probability)
  )
)

# Full fits and trend sensitivity ----------------------------------------

reference_days <- read_csv(
  file.path(analysis_dir, "tables", "day_of_season_support.csv"),
  show_col_types = FALSE
) |>
  filter(in_reference_window) |>
  pull(day_of_season)
seasons <- sort(unique(predictors$season))
prediction_grid <- expand_grid(season = seasons, day_of_season = reference_days) |>
  mutate(
    season_centered = season - median(model_data$season),
    moon_distance_from_new_moon_z = 0,
    mist_probability_light_patchy = mean(predictors$mist_probability_light_patchy),
    mist_probability_good = mean(predictors$mist_probability_good),
    era5_rain_log_z = 0,
    wind_speed_10m_mean_ms_z = 0,
    temperature_2m_mean_c_z = 0,
    surface_pressure_mean_hpa_z = 0
  )

community_index <- function(probability) {
  annual_probability <- probability |>
    as.data.frame() |>
    mutate(season = prediction_grid$season, .before = 1) |>
    group_by(season) |>
    summarise(across(everything(), mean), .groups = "drop")
  focal_probability <- as.matrix(annual_probability[, focal_names])
  centered_log_ratio <- log(focal_probability) - rowMeans(log(focal_probability))
  exp(sweep(centered_log_ratio, 2, colMeans(centered_log_ratio), FUN = "-"))
}

full_models <- lapply(effective_counts, fit_multinomial, data = predictors)
full_probability <- lapply(full_models, function(model) {
  predict(model, newdata = prediction_grid, type = "probs")[, colnames(counts), drop = FALSE]
})

full_probability$dirichlet_multinomial <- predict_dm(dm_full_model, prediction_grid)
indices <- lapply(full_probability, community_index)

high_catch_thresholds <- quantile(daily_total, c(0.99, 0.95))
for (threshold_name in names(high_catch_thresholds)) {
  keep <- daily_total <= high_catch_thresholds[[threshold_name]]
  model_id <- paste0("exclude_top_", 100 - as.numeric(sub("%", "", threshold_name)), "pct")
  exclusion_model <- fit_multinomial(counts[keep, , drop = FALSE], predictors[keep, ])
  probability <- predict(
    exclusion_model,
    newdata = prediction_grid,
    type = "probs"
  )[, colnames(counts), drop = FALSE]
  indices[[model_id]] <- community_index(probability)
}

baseline_index <- indices$multinomial
trend_sensitivity <- bind_rows(lapply(names(indices), function(model_id) {
  alternative_index <- indices[[model_id]]
  log_difference <- log(alternative_index) - log(baseline_index)
  endpoint_ratio <- alternative_index[match(2022, seasons), ] /
    alternative_index[match(1977, seasons), ]
  baseline_endpoint_ratio <- baseline_index[match(2022, seasons), ] /
    baseline_index[match(1977, seasons), ]
  tibble(
    model_id,
    rms_log_index_difference = sqrt(mean(log_difference^2)),
    maximum_fold_index_difference = max(exp(abs(log_difference))),
    median_endpoint_fold_difference = median(exp(abs(log(endpoint_ratio / baseline_endpoint_ratio)))),
    maximum_endpoint_fold_difference = max(exp(abs(log(endpoint_ratio / baseline_endpoint_ratio)))),
    curve_log_correlation = cor(as.vector(log(alternative_index)), as.vector(log(baseline_index)))
  )
}))

endpoint_comparison <- bind_rows(lapply(names(indices), function(model_id) {
  alternative_index <- indices[[model_id]]
  tibble(
    model_id,
    common_name = focal_names,
    relative_change_1977_2022 =
      alternative_index[match(2022, seasons), ] /
      alternative_index[match(1977, seasons), ]
  )
})) |>
  left_join(categories |> select(common_name, abundance_rank), by = "common_name") |>
  arrange(abundance_rank, model_id)

# One-at-a-time deletion of the ten largest dates ------------------------

largest_dates <- order(daily_total, decreasing = TRUE)[1:10]
top_date_influence <- bind_rows(lapply(largest_dates, function(row_index) {
  keep <- seq_len(nrow(counts)) != row_index
  model <- fit_multinomial(counts[keep, , drop = FALSE], predictors[keep, ])
  probability <- predict(model, newdata = prediction_grid, type = "probs")
  index <- community_index(probability[, colnames(counts), drop = FALSE])
  log_difference <- log(index) - log(baseline_index)
  tibble(
    ringing_date = predictors$ringing_date[[row_index]],
    season = predictors$season[[row_index]],
    total_comparable_count = daily_total[[row_index]],
    rms_log_index_difference = sqrt(mean(log_difference^2)),
    maximum_fold_index_difference = max(exp(abs(log_difference)))
  )
})) |>
  arrange(desc(rms_log_index_difference))

catch_concentration <- tibble(
  threshold = c("top_1_percent_dates", "top_5_percent_dates", "top_10_percent_dates", "top_20_percent_dates"),
  n_dates = ceiling(nrow(counts) * c(0.01, 0.05, 0.10, 0.20))
) |>
  rowwise() |>
  mutate(
    share_of_birds = sum(sort(daily_total, decreasing = TRUE)[seq_len(n_dates)]) / sum(daily_total),
    smallest_daily_total = sort(daily_total, decreasing = TRUE)[[n_dates]]
  ) |>
  ungroup()

# Fit the variance-tempered primary model --------------------------------

primary_model <- fit_multinomial(
  effective_counts$variance_tempered,
  predictors,
  hessian = TRUE
)
primary_probability <- predict(primary_model, newdata = predictors, type = "probs")
primary_probability <- primary_probability[, colnames(counts), drop = FALSE]
effective_total <- rowSums(effective_counts$variance_tempered)
expected_effective_counts <- effective_total * primary_probability
n_residual_cells <- nrow(counts) * (ncol(counts) - 1)
pearson_dispersion <- sum(
  (effective_counts$variance_tempered - expected_effective_counts)^2 /
    pmax(expected_effective_counts, 1e-8)
) / (n_residual_cells - length(coef(primary_model)))
primary_scores <- probability_scores(counts, primary_probability)

primary_diagnostics <- tibble(
  primary_model_id = "M3_variance_tempered",
  primary_model_label = "M3 + Daily conditions, variance-tempered dates",
  n_dates = nrow(counts),
  n_seasons = n_distinct(predictors$season),
  n_focal_species = sum(categories$is_focal_species),
  n_model_categories = nrow(categories),
  n_birds = sum(counts),
  effective_n_birds = sum(effective_total),
  median_daily_count = median(daily_total),
  median_effective_daily_count = median(effective_total),
  maximum_daily_count = max(daily_total),
  maximum_effective_daily_count = max(effective_total),
  dm_concentration,
  n_parameters = length(coef(primary_model)),
  convergence = primary_model$convergence,
  aic = AIC(primary_model),
  pearson_dispersion,
  uncertainty_dispersion = max(1, pearson_dispersion),
  fitted_bird_weighted_log_loss = primary_scores$bird_weighted_log_loss,
  fitted_mean_daily_log_loss = primary_scores$mean_daily_log_loss,
  max_probability_sum_error = max(abs(rowSums(primary_probability) - 1))
)

primary_bundle <- list(
  model = primary_model,
  model_id = "M3_variance_tempered",
  specification = bundle$specification,
  categories = categories,
  predictor_columns = predictor_columns,
  category_order = colnames(counts),
  season_spline_df = bundle$season_spline_df,
  phenology_spline_df = bundle$phenology_spline_df,
  dm_concentration = dm_concentration,
  weighting = "Dirichlet-multinomial variance-matched effective daily totals"
)

weighting_descriptions <- tibble(
  model_id = c(
    "multinomial", "capped_95", "square_root", "variance_tempered",
    "equal_date", "dirichlet_multinomial"
  ),
  approach = c(
    "Ordinary multinomial: each bird has equal weight",
    "Daily effective total capped at the 95th percentile",
    "Daily influence grows with the square root of total catch",
    "Daily effective total matched to Dirichlet-multinomial variance",
    "Every positive-catch date has equal total weight",
    "Full Dirichlet-multinomial likelihood"
  ),
  role = c("comparison", "sensitivity", "sensitivity", "primary", "sensitivity", "distributional sensitivity")
)
validation_scores <- validation_scores |>
  left_join(weighting_descriptions, by = "model_id") |>
  select(model_id, approach, role, everything())

model_decision <- tibble(
  primary_model = "variance_tempered",
  decision = "Use variance-tempered multinomial mean model",
  rationale = paste(
    "Exceptional dates do not materially determine the curves, but daily species vectors are strongly overdispersed.",
    "The Dirichlet-multinomial concentration provides a principled effective daily sample size.",
    "The resulting mean model has the best held-out mean-daily log loss while retaining nearly unchanged bird-weighted prediction.",
    "The full Dirichlet-multinomial is retained as a distributional sensitivity because it predicts count-vector dispersion well but predicts held-out mean composition less accurately."
  )
)

saveRDS(primary_bundle, file.path(model_dir, "joint_composition_primary.rds"))
saveRDS(dm_full_model, file.path(model_dir, "dirichlet_multinomial_sensitivity.rds"))
write_csv(primary_diagnostics, file.path(table_dir, "primary_model_diagnostics.csv"))
write_csv(validation_scores, file.path(table_dir, "high_catch_model_comparison.csv"))
write_csv(
  tempered_mean_structure_comparison,
  file.path(table_dir, "tempered_mean_structure_comparison.csv")
)
write_csv(count_distribution_scores, file.path(table_dir, "count_distribution_scores.csv"))
write_csv(trend_sensitivity, file.path(table_dir, "high_catch_trend_sensitivity.csv"))
write_csv(endpoint_comparison, file.path(table_dir, "high_catch_endpoint_comparison.csv"))
write_csv(top_date_influence, file.path(table_dir, "top_date_influence.csv"))
write_csv(catch_concentration, file.path(table_dir, "catch_concentration.csv"))
write_csv(
  tibble(fold = names(dm_fold_concentration), concentration = dm_fold_concentration),
  file.path(table_dir, "dm_fold_concentration.csv")
)
write_csv(model_decision, file.path(table_dir, "high_catch_model_decision.csv"))

cli_alert_success("Selected the variance-tempered M3 model as primary.")
cli_alert_info(
  "Dirichlet-multinomial concentration: {signif(dm_concentration, 4)}; maximum effective daily total: {signif(max(effective_total), 4)}."
)
