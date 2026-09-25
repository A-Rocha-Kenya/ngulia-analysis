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
source(file.path(project_dir, "scripts", "helpers", "weighted_quantile.R"))

analysis_dir <- file.path(paths$analysis_output_dir, "03_migration_phenology")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")

n_draws <- 1000L
primary_start_season <- 1977L
full_record_start_season <- 1969L
fit_first_day <- 35L
fit_last_day <- 90L
reference_first_day <- 42L
reference_last_day <- 81L
quantile_probabilities <- c(q25 = 0.25, q50 = 0.50, q75 = 0.75)

# Read species-level daily counts ----------------------------------------

cli_h1("Fit within-species passage-date models")

daily_data <- read_csv(
  file.path(model_data_dir, "species_daily_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
) |>
  filter(between(day_of_season, fit_first_day, fit_last_day))

date_covariates <- daily_data |>
  distinct(
    ringing_date, season, day_of_season, moon_distance_from_new_moon,
    mist_probability_light_patchy, mist_probability_good
  ) |>
  filter(season >= primary_start_season)

moon_center <- mean(date_covariates$moon_distance_from_new_moon)
moon_scale <- sd(date_covariates$moon_distance_from_new_moon)
season_center <- median(date_covariates$season)

daily_data <- daily_data |>
  mutate(
    season_decades = (season - season_center) / 10,
    moon_distance_z = (moon_distance_from_new_moon - moon_center) / moon_scale,
    pre_night_net = as.integer(capture_era == "dawn-net era"),
    transition = as.integer(capture_era == "transition")
  )

primary_data <- daily_data |>
  filter(season >= primary_start_season)

daily_caps <- primary_data |>
  filter(count > 0) |>
  group_by(avibase_id, common_name) |>
  summarise(daily_count_cap = quantile(count, 0.95, type = 1), .groups = "drop")

daily_data <- daily_data |>
  left_join(daily_caps, by = c("avibase_id", "common_name")) |>
  mutate(model_count = pmin(count, daily_count_cap))

primary_data <- daily_data |>
  filter(season >= primary_start_season)

day_basis <- ns(primary_data$day_of_season, df = 4)
moon_basis <- ns(primary_data$moon_distance_z, df = 3)

make_design <- function(data, changing = TRUE, adjusted = TRUE, era_adjusted = FALSE) {
  day <- predict(day_basis, data$day_of_season)
  colnames(day) <- paste0("day", seq_len(ncol(day)))
  design <- day

  if (changing) {
    trend <- day * data$season_decades
    colnames(trend) <- paste0("trend_day", seq_len(ncol(trend)))
    design <- cbind(design, trend)
  }

  if (adjusted) {
    moon <- predict(moon_basis, data$moon_distance_z)
    colnames(moon) <- paste0("moon", seq_len(ncol(moon)))
    design <- cbind(
      design, moon,
      mist_light_patchy = data$mist_probability_light_patchy,
      mist_good = data$mist_probability_good
    )
  }

  if (era_adjusted) {
    pre_era <- day * data$pre_night_net
    transition_era <- day * data$transition
    colnames(pre_era) <- paste0("pre_era_day", seq_len(ncol(pre_era)))
    colnames(transition_era) <- paste0("transition_day", seq_len(ncol(transition_era)))
    design <- cbind(design, pre_era, transition_era)
  }

  unname(as.matrix(design))
}

# The likelihood conditions on each species' annual total. This removes annual
# abundance and any effort multiplier that is constant within season. Daily
# counts are capped at the species-specific 95th percentile so that mass falls
# do not act like thousands of independent observations of passage date.
conditional_objective <- function(parameters, design, count, season) {
  linear_predictor <- as.vector(design %*% parameters)
  negative_log_likelihood <- 0

  for (season_value in unique(season)) {
    rows <- season == season_value
    season_count <- count[rows]
    season_total <- sum(season_count)
    if (season_total == 0) next

    season_predictor <- linear_predictor[rows]
    log_normalizer <- max(season_predictor) + log(sum(exp(season_predictor - max(season_predictor))))
    negative_log_likelihood <- negative_log_likelihood -
      sum(season_count * season_predictor) + season_total * log_normalizer
  }

  negative_log_likelihood
}

conditional_gradient <- function(parameters, design, count, season) {
  linear_predictor <- as.vector(design %*% parameters)
  gradient <- rep(0, ncol(design))

  for (season_value in unique(season)) {
    rows <- season == season_value
    season_count <- count[rows]
    season_total <- sum(season_count)
    if (season_total == 0) next
    probability <- softmax(linear_predictor[rows])
    gradient <- gradient + crossprod(design[rows, , drop = FALSE], season_total * probability - season_count)
  }

  as.vector(gradient)
}

conditional_hessian <- function(parameters, design, count, season) {
  linear_predictor <- as.vector(design %*% parameters)
  hessian <- matrix(0, ncol(design), ncol(design))

  for (season_value in unique(season)) {
    rows <- season == season_value
    season_total <- sum(count[rows])
    if (season_total == 0) next
    season_design <- design[rows, , drop = FALSE]
    probability <- softmax(linear_predictor[rows])
    weighted_crossproduct <- crossprod(season_design, season_design * probability)
    mean_design <- colSums(season_design * probability)
    hessian <- hessian + season_total * (weighted_crossproduct - tcrossprod(mean_design))
  }

  hessian
}

clustered_covariance <- function(parameters, design, count, season) {
  hessian <- conditional_hessian(parameters, design, count, season)
  bread <- tryCatch(solve(hessian), error = function(error) MASS::ginv(hessian))
  linear_predictor <- as.vector(design %*% parameters)
  scores <- lapply(unique(season), function(season_value) {
    rows <- season == season_value
    season_count <- count[rows]
    season_total <- sum(season_count)
    if (season_total == 0) return(rep(0, ncol(design)))
    probability <- softmax(linear_predictor[rows])
    as.vector(crossprod(design[rows, , drop = FALSE], season_count - season_total * probability))
  }) |>
    (\(x) do.call(rbind, x))()
  meat <- crossprod(scores)
  n_cluster <- nrow(scores)
  bread %*% meat %*% bread * n_cluster / (n_cluster - 1)
}

fit_species_model <- function(data, changing = TRUE, adjusted = TRUE, era_adjusted = FALSE) {
  design <- make_design(data, changing, adjusted, era_adjusted)
  initial <- rep(0, ncol(design))
  fit <- optim(
    initial, conditional_objective, conditional_gradient,
    design = design, count = data$model_count, season = data$season,
    method = "BFGS",
    control = list(maxit = 3000, reltol = 1e-9)
  )
  covariance <- clustered_covariance(fit$par, design, data$model_count, data$season)

  list(
    coefficients = fit$par,
    covariance = covariance,
    negative_log_likelihood = fit$value,
    aic = 2 * fit$value + 2 * length(fit$par),
    convergence = fit$convergence,
    changing = changing,
    adjusted = adjusted,
    era_adjusted = era_adjusted
  )
}

softmax <- function(linear_predictor) {
  probability <- exp(linear_predictor - max(linear_predictor))
  probability / sum(probability)
}

curve_quantiles <- function(probability, day) {
  cumulative <- cumsum(probability)
  vapply(quantile_probabilities, function(tau) {
    approx(c(0, cumulative), c(day[[1]] - 1, day), xout = tau, ties = "ordered")$y
  }, numeric(1))
}

prediction_grid <- expand_grid(
  season = sort(unique(primary_data$season)),
  day_of_season = seq(reference_first_day, reference_last_day)
) |>
  mutate(
    season_decades = (season - season_center) / 10,
    moon_distance_z = 0,
    mist_probability_light_patchy = mean(primary_data$mist_probability_light_patchy),
    mist_probability_good = mean(primary_data$mist_probability_good),
    pre_night_net = 0L,
    transition = 0L
  )

model_quantiles <- function(model, grid, coefficient = model$coefficients) {
  design <- make_design(grid, model$changing, model$adjusted, model$era_adjusted)
  linear_predictor <- as.vector(design %*% coefficient)

  split(seq_len(nrow(grid)), grid$season) |>
    map_dfr(function(rows) {
      values <- curve_quantiles(softmax(linear_predictor[rows]), grid$day_of_season[rows])
      tibble(
        season = grid$season[rows[[1]]],
        quantile = names(values),
        passage_day = unname(values)
      )
    })
}

quantile_slopes <- function(quantiles) {
  quantiles |>
    group_by(quantile) |>
    summarise(
      slope_days_per_decade = 10 * cov(season, passage_day) / var(season),
      .groups = "drop"
    )
}

draw_coefficients <- function(model, n) {
  covariance <- (model$covariance + t(model$covariance)) / 2
  eigen_decomposition <- eigen(covariance, symmetric = TRUE)
  eigen_decomposition$values <- pmax(eigen_decomposition$values, 1e-10)
  root <- eigen_decomposition$vectors %*%
    diag(sqrt(eigen_decomposition$values), nrow = length(eigen_decomposition$values))
  draws <- matrix(rnorm(n * length(model$coefficients)), nrow = n) %*% t(root)
  sweep(draws, 2, model$coefficients, "+")
}

# Fit stable, changing and condition-adjusted models ---------------------

species_metadata <- primary_data |>
  distinct(avibase_id, common_name, abundance_rank) |>
  arrange(abundance_rank)

species_models <- vector("list", nrow(species_metadata))
names(species_models) <- species_metadata$common_name
model_comparison <- vector("list", nrow(species_metadata))
annual_standardized_quantiles <- vector("list", nrow(species_metadata))
trend_results <- vector("list", nrow(species_metadata))
condition_sensitivity <- vector("list", nrow(species_metadata))

set.seed(1977)
for (species_row in seq_len(nrow(species_metadata))) {
  species_name <- species_metadata$common_name[[species_row]]
  cli_alert_info("Fitting {species_name}")
  species_data <- primary_data |> filter(common_name == species_name)

  stable_adjusted <- fit_species_model(species_data, changing = FALSE, adjusted = TRUE)
  changing_unadjusted <- fit_species_model(species_data, changing = TRUE, adjusted = FALSE)
  changing_adjusted <- fit_species_model(species_data, changing = TRUE, adjusted = TRUE)
  species_models[[species_name]] <- changing_adjusted

  model_comparison[[species_row]] <- tibble(
    avibase_id = species_metadata$avibase_id[[species_row]],
    common_name = species_name,
    abundance_rank = species_metadata$abundance_rank[[species_row]],
    model = c("Stable phenology + moon and mist", "Changing phenology, unadjusted", "Changing phenology + moon and mist"),
    aic = c(stable_adjusted$aic, changing_unadjusted$aic, changing_adjusted$aic),
    n_parameters = c(length(stable_adjusted$coefficients), length(changing_unadjusted$coefficients), length(changing_adjusted$coefficients)),
    convergence = c(stable_adjusted$convergence, changing_unadjusted$convergence, changing_adjusted$convergence)
  )

  point_quantiles <- model_quantiles(changing_adjusted, prediction_grid) |>
    mutate(
      avibase_id = species_metadata$avibase_id[[species_row]],
      common_name = species_name,
      abundance_rank = species_metadata$abundance_rank[[species_row]],
      .before = 1
    )
  annual_standardized_quantiles[[species_row]] <- point_quantiles
  point_slopes <- quantile_slopes(point_quantiles)

  condition_sensitivity[[species_row]] <- bind_rows(
    quantile_slopes(model_quantiles(changing_unadjusted, prediction_grid)) |>
      mutate(model = "Unadjusted for moon and mist", .before = 1),
    point_slopes |>
      mutate(model = "Moon- and mist-standardized", .before = 1)
  ) |>
    mutate(
      avibase_id = species_metadata$avibase_id[[species_row]],
      common_name = species_name,
      abundance_rank = species_metadata$abundance_rank[[species_row]],
      .before = 1
    )

  coefficient_draws <- draw_coefficients(changing_adjusted, n_draws)
  slope_draws <- matrix(NA_real_, nrow = n_draws, ncol = length(quantile_probabilities), dimnames = list(NULL, names(quantile_probabilities)))
  for (draw in seq_len(n_draws)) {
    draw_quantiles <- model_quantiles(changing_adjusted, prediction_grid, coefficient_draws[draw, ])
    draw_slopes <- quantile_slopes(draw_quantiles)
    slope_draws[draw, draw_slopes$quantile] <- draw_slopes$slope_days_per_decade
  }

  trend_results[[species_row]] <- point_slopes |>
    rowwise() |>
    mutate(
      avibase_id = species_metadata$avibase_id[[species_row]],
      common_name = species_name,
      abundance_rank = species_metadata$abundance_rank[[species_row]],
      lower = quantile(slope_draws[, quantile], 0.025),
      upper = quantile(slope_draws[, quantile], 0.975),
      coefficient_p = min(1, 2 * min(
        (sum(slope_draws[, quantile] <= 0) + 1) / (n_draws + 1),
        (sum(slope_draws[, quantile] >= 0) + 1) / (n_draws + 1)
      )),
      .before = 1
    ) |>
    ungroup()
}

model_comparison <- bind_rows(model_comparison)
annual_standardized_quantiles <- bind_rows(annual_standardized_quantiles)
trend_results <- bind_rows(trend_results)
condition_sensitivity <- bind_rows(condition_sensitivity)

# Define supported quantiles and multiplicity families -------------------

quantile_support <- primary_data |>
  group_by(avibase_id, common_name, abundance_rank) |>
  summarise(
    pooled_q25 = weighted_quantile(day_of_season, count, 0.25),
    pooled_q50 = weighted_quantile(day_of_season, count, 0.50),
    pooled_q75 = weighted_quantile(day_of_season, count, 0.75),
    .groups = "drop"
  ) |>
  pivot_longer(starts_with("pooled_q"), names_prefix = "pooled_", names_to = "quantile", values_to = "pooled_passage_day") |>
  mutate(
    reference_first_day,
    reference_last_day,
    supported = pooled_passage_day >= reference_first_day + 2 & pooled_passage_day <= reference_last_day - 2
  )

trend_results <- trend_results |>
  left_join(
    quantile_support |> select(avibase_id, quantile, pooled_passage_day, supported),
    by = c("avibase_id", "quantile")
  ) |>
  group_by(test_family = if_else(quantile == "q50", "primary median", "secondary quartiles")) |>
  mutate(adjusted_p = p.adjust(if_else(supported, coefficient_p, NA_real_), method = "BH")) |>
  ungroup() |>
  mutate(
    direction = case_when(
      !supported ~ "outside reference-window support",
      adjusted_p < 0.05 & slope_days_per_decade < 0 ~ "earlier",
      adjusted_p < 0.05 & slope_days_per_decade > 0 ~ "later",
      TRUE ~ "no clear change"
    )
  ) |>
  arrange(abundance_rank, match(quantile, names(quantile_probabilities)))

# Full-record sensitivity with capture-regime adjustment ----------------

full_grid <- expand_grid(
  season = seq(full_record_start_season, max(primary_data$season)),
  day_of_season = seq(reference_first_day, reference_last_day)
) |>
  mutate(
    season_decades = (season - season_center) / 10,
    moon_distance_z = 0,
    mist_probability_light_patchy = mean(primary_data$mist_probability_light_patchy),
    mist_probability_good = mean(primary_data$mist_probability_good),
    pre_night_net = 0L,
    transition = 0L
  )

full_record_sensitivity <- vector("list", nrow(species_metadata))
full_record_models <- vector("list", nrow(species_metadata))
names(full_record_models) <- species_metadata$common_name

for (species_row in seq_len(nrow(species_metadata))) {
  species_name <- species_metadata$common_name[[species_row]]
  species_data <- daily_data |> filter(common_name == species_name)
  full_model <- fit_species_model(species_data, changing = TRUE, adjusted = TRUE, era_adjusted = TRUE)
  full_record_models[[species_name]] <- full_model
  full_slopes <- model_quantiles(full_model, full_grid) |> quantile_slopes()
  full_record_sensitivity[[species_row]] <- full_slopes |>
    mutate(
      avibase_id = species_metadata$avibase_id[[species_row]],
      common_name = species_name,
      abundance_rank = species_metadata$abundance_rank[[species_row]],
      model = "1969–2023, capture-regime adjusted",
      .before = 1
    )
}

full_record_sensitivity <- bind_rows(full_record_sensitivity) |>
  bind_rows(
    trend_results |>
      select(avibase_id, common_name, abundance_rank, quantile, slope_days_per_decade) |>
      mutate(model = "1977–2023 primary", .before = 1)
  ) |>
  arrange(abundance_rank, quantile, model)

# Save outputs ------------------------------------------------------------

model_diagnostics <- tibble(
  n_species = nrow(species_metadata),
  primary_first_season = min(primary_data$season),
  primary_last_season = max(primary_data$season),
  primary_n_seasons = n_distinct(primary_data$season),
  primary_n_sampling_dates = n_distinct(primary_data$ringing_date),
  full_record_first_season = min(daily_data$season),
  full_record_n_sampling_dates = n_distinct(daily_data$ringing_date),
  fit_first_day, fit_last_day, reference_first_day, reference_last_day,
  n_draws
)

saveRDS(
  list(
    models = species_models,
    full_record_models = full_record_models,
    day_basis = day_basis,
    moon_basis = moon_basis,
    moon_center = moon_center,
    moon_scale = moon_scale,
    season_center = season_center,
    reference_conditions = c(
      moon_distance_z = 0,
      mist_probability_light_patchy = mean(primary_data$mist_probability_light_patchy),
      mist_probability_good = mean(primary_data$mist_probability_good)
    )
  ),
  file.path(model_dir, "within_species_phenology_models.rds")
)

write_csv(annual_standardized_quantiles, file.path(model_data_dir, "species_standardized_annual_quantiles.csv"))
write_csv(trend_results, file.path(table_dir, "species_phenology_trends.csv"))
write_csv(model_comparison, file.path(table_dir, "species_phenology_model_comparison.csv"))
write_csv(quantile_support, file.path(table_dir, "species_phenology_quantile_support.csv"))
write_csv(full_record_sensitivity, file.path(table_dir, "full_record_phenology_sensitivity.csv"))
write_csv(condition_sensitivity, file.path(table_dir, "moon_mist_phenology_sensitivity.csv"))
write_csv(model_diagnostics, file.path(table_dir, "species_phenology_model_diagnostics.csv"))
write_csv(daily_caps, file.path(table_dir, "species_daily_count_caps.csv"))

cli_alert_success("Fitted within-species passage models for {nrow(species_metadata)} species.")
cli_alert_info("Primary inference uses q50; q25 and q75 form a separate secondary family.")
