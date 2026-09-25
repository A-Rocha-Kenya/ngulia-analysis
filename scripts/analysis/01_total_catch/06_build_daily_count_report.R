library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(mgcv)
library(MASS)
library(scales)
library(patchwork)
library(htmltools)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- ngulia_figure_dir(file.path(analysis_dir, "figures"))

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read fitted models and validation summaries ----------------------------

cli_h1("Build daily-count story report")

primary_models <- readRDS(file.path(model_dir, "adjusted_positive_catch_gam_mi.rds"))
operation_models <- readRDS(file.path(model_dir, "configuration_playback_gam_mi.rds"))
m4_subset_models <- operation_models[["M4-S Smooth-year baseline"]]
m5_models <- operation_models[["M5 + bush configuration"]]
m6_models <- operation_models[["M6 + playback"]]
model_progression <- read_csv(file.path(table_dir, "annual_index_model_comparison.csv"), show_col_types = FALSE)
operation_comparison <- read_csv(file.path(table_dir, "documented_operations_model_comparison.csv"), show_col_types = FALSE)
operation_contribution <- read_csv(file.path(table_dir, "documented_operations_covariate_contribution.csv"), show_col_types = FALSE)
operation_coverage <- read_csv(file.path(table_dir, "documented_operations_covariate_coverage.csv"), show_col_types = FALSE)
diagnostics <- read_csv(file.path(table_dir, "primary_model_diagnostics.csv"), show_col_types = FALSE)
daily_coverage <- read_csv(file.path(paths$curated_dir, "daily_coverage.csv"), show_col_types = FALSE)
operation_source <- read_csv(file.path(analysis_dir, "model_data", "documented_operation_model_data.csv"), show_col_types = FALSE)

primary_data <- lapply(primary_models, function(model) {
  model$model |>
    mutate(season = as.integer(as.character(season_factor)))
})
operation_data <- lapply(m6_models, function(model) model$model)

# Pool adjusted covariate effects ----------------------------------------

pool_effect_curve <- function(models, variable, values, reference) {
  estimates <- bind_rows(lapply(seq_along(models), function(imputation) {
    model <- models[[imputation]]
    prediction_data <- model$model[rep(1, length(values)), ]
    prediction_data[[variable]] <- values
    reference_data <- prediction_data[1, , drop = FALSE]
    reference_data[[variable]] <- reference
    contrast <- sweep(
      predict(model, newdata = prediction_data, type = "lpmatrix"),
      2,
      drop(predict(model, newdata = reference_data, type = "lpmatrix"))
    )
    tibble(
      imputation,
      value = values,
      estimate = drop(contrast %*% coef(model)),
      variance = rowSums((contrast %*% vcov(model)) * contrast)
    )
  }))

  estimates |>
    group_by(value) |>
    summarise(
      pooled_estimate = mean(estimate),
      within_variance = mean(variance),
      between_variance = var(estimate),
      n_imputations = n(),
      .groups = "drop"
    ) |>
    mutate(
      total_variance = within_variance + (1 + 1 / n_imputations) * between_variance,
      relative_effect = exp(pooled_estimate),
      lower = exp(pooled_estimate - 1.96 * sqrt(total_variance)),
      upper = exp(pooled_estimate + 1.96 * sqrt(total_variance))
    )
}

core_variables <- tribble(
  ~variable, ~label,
  "season_day", "Date within season",
  "moon_distance_from_new_moon", "Distance from new moon (days)",
  "era5_rain_log", "Rainfall: log(1 + mm)",
  "wind_speed_10m_mean_ms", "Wind speed (m/s)",
  "temperature_2m_mean_c", "Temperature (°C)",
  "surface_pressure_mean_hpa", "Surface pressure (hPa)"
)

core_effects <- bind_rows(lapply(seq_len(nrow(core_variables)), function(i) {
  variable <- core_variables$variable[[i]]
  observed <- primary_data[[1]][[variable]]
  values <- seq(quantile(observed, 0.02), quantile(observed, 0.98), length.out = 100)
  pool_effect_curve(primary_models, variable, values, median(observed)) |>
    mutate(variable = variable, label = core_variables$label[[i]])
})) |>
  mutate(label = factor(label, levels = core_variables$label))

mist_levels <- c("none", "light_patchy", "good")
mist_effects_raw <- bind_rows(lapply(seq_along(primary_models), function(imputation) {
  model <- primary_models[[imputation]]
  reference_data <- model$model[1, , drop = FALSE]
  reference_data$mist_state <- factor("none", levels = mist_levels)
  reference_matrix <- predict(model, newdata = reference_data, type = "lpmatrix")

  bind_rows(lapply(mist_levels, function(state) {
    prediction_data <- reference_data
    prediction_data$mist_state <- factor(state, levels = mist_levels)
    contrast <- predict(model, newdata = prediction_data, type = "lpmatrix") - reference_matrix
    tibble(
      imputation,
      mist_state = state,
      estimate = drop(contrast %*% coef(model)),
      variance = drop(contrast %*% vcov(model) %*% t(contrast))
    )
  }))
}))

mist_effects <- mist_effects_raw |>
  group_by(mist_state) |>
  summarise(
    pooled_estimate = mean(estimate),
    within_variance = mean(variance),
    between_variance = var(estimate),
    n_imputations = n(),
    .groups = "drop"
  ) |>
  mutate(
    total_variance = within_variance + (1 + 1 / n_imputations) * between_variance,
    relative_effect = exp(pooled_estimate),
    lower = exp(pooled_estimate - 1.96 * sqrt(total_variance)),
    upper = exp(pooled_estimate + 1.96 * sqrt(total_variance)),
    label = factor(mist_state, levels = mist_levels, labels = c("No mist", "Light/patchy", "Good mist"))
  )

operation_binary_effects <- bind_rows(
  pool_effect_curve(m5_models, "bush_period", "front_bush", "back_bush") |>
    mutate(variable = "bush_period", label = "Front versus back bush (M5)") |>
    dplyr::select(-value),
  pool_effect_curve(m6_models, "bush_period", "front_bush", "back_bush") |>
    mutate(variable = "bush_period_m6", label = "Front versus back bush (M6)") |>
    dplyr::select(-value),
  pool_effect_curve(m6_models, "playback_used", 1, 0) |>
    mutate(variable = "playback_used", label = "Nocturnal playback (M6)") |>
    dplyr::select(-value)
)

# Estimate annual fluctuations from season effects -----------------------

draw_annual_index <- function(models, model_name, draws_per_imputation = 200) {
  index_draws <- vector("list", length(models))
  point_estimates <- vector("list", length(models))

  for (imputation in seq_along(models)) {
    model <- models[[imputation]]
    seasons <- levels(model$model$season_factor)
    coefficient_names <- c("(Intercept)", grep("^season_factor", names(coef(model)), value = TRUE))
    design <- matrix(0, nrow = length(seasons), ncol = length(coefficient_names))
    colnames(design) <- coefficient_names
    design[, "(Intercept)"] <- 1
    for (i in seq_along(seasons)[-1]) {
      design[i, paste0("season_factor", seasons[[i]])] <- 1
    }

    coefficient_draws <- mvrnorm(
      draws_per_imputation,
      coef(model)[coefficient_names],
      vcov(model)[coefficient_names, coefficient_names, drop = FALSE]
    )
    expected_draws <- exp(design %*% t(coefficient_draws))
    relative_draws <- sweep(expected_draws, 2, colMeans(expected_draws), "/")
    point <- exp(drop(design %*% coef(model)[coefficient_names]))

    colnames(relative_draws) <- paste0("draw_", seq_len(ncol(relative_draws)))
    index_draws[[imputation]] <- as_tibble(relative_draws) |>
      mutate(season = as.integer(seasons), .before = 1) |>
      pivot_longer(-season, names_to = "draw", values_to = "relative_index") |>
      mutate(draw = paste(imputation, draw, sep = "_"))
    point_estimates[[imputation]] <- tibble(
      imputation,
      season = as.integer(seasons),
      relative_index = point / mean(point)
    )
  }

  draws <- bind_rows(index_draws)
  summary <- draws |>
    group_by(season) |>
    summarise(
      relative_index = median(relative_index),
      lower = quantile(relative_index, 0.025),
      upper = quantile(relative_index, 0.975),
      .groups = "drop"
    ) |>
    mutate(model = model_name)

  list(summary = summary, draws = draws, point = bind_rows(point_estimates))
}

set.seed(910)
primary_annual <- draw_annual_index(primary_models, "M4 post-transition")
annual_sample_size <- primary_data[[1]] |>
  count(season, name = "n_positive_catch_dates") |>
  mutate(model = "M4 post-transition")
annual_index <- primary_annual$summary |>
  left_join(annual_sample_size, by = c("season", "model"))

# Fit smooth and linear long-term trends ---------------------------------

m4_trend_formula <- total_birds_ringed ~
  s(season, k = 10) + s(season_day, k = 12) +
  s(moon_distance_from_new_moon, k = 6) + mist_state +
  s(era5_rain_log, k = 6) + s(wind_speed_10m_mean_ms, k = 6) +
  s(temperature_2m_mean_c, k = 6) + s(surface_pressure_mean_hpa, k = 6)

m4_linear_formula <- update(m4_trend_formula, . ~ . - s(season, k = 10) + season_centered)

fit_trend_models <- function(data, smooth_formula, linear_formula) {
  smooth_models <- vector("list", length(data))
  linear_models <- vector("list", length(data))
  for (imputation in seq_along(data)) {
    model_data <- data[[imputation]] |>
      mutate(season_centered = season - mean(season))
    smooth_models[[imputation]] <- bam(smooth_formula, data = model_data, family = nb(), method = "fREML", discrete = TRUE)
    linear_models[[imputation]] <- gam(linear_formula, data = model_data, family = nb(), method = "REML")
  }
  list(smooth = smooth_models, linear = linear_models)
}

primary_trend_models <- fit_trend_models(primary_data, m4_trend_formula, m4_linear_formula)
fit_linear_models <- function(formula) {
  lapply(operation_data, function(data) {
    gam(formula, data = mutate(data, season_centered = season - mean(season)),
      family = nb(), method = "REML")
  })
}
subset_weather_trend_models <- list(smooth = m4_subset_models, linear = fit_linear_models(m4_linear_formula))
configuration_trend_models <- list(smooth = m5_models,
  linear = fit_linear_models(update(m4_linear_formula, . ~ . + bush_period)))
playback_trend_models <- list(smooth = m6_models,
  linear = fit_linear_models(update(m4_linear_formula, . ~ . + bush_period + playback_used)))

pool_trend_curve <- function(models, model_name) {
  seasons <- seq(
    min(models[[1]]$model$season),
    max(models[[1]]$model$season)
  )
  estimates <- bind_rows(lapply(seq_along(models), function(imputation) {
    model <- models[[imputation]]
    prediction_data <- model$model[rep(1, length(seasons)), ]
    prediction_data$season <- seasons
    prediction <- predict(model, newdata = prediction_data, type = "terms", terms = "s(season)", se.fit = TRUE)
    tibble(
      imputation,
      season = seasons,
      estimate = drop(prediction$fit),
      variance = drop(prediction$se.fit)^2
    )
  }))

  pooled <- estimates |>
    group_by(season) |>
    summarise(
      pooled_estimate = mean(estimate),
      within_variance = mean(variance),
      between_variance = var(estimate),
      n_imputations = n(),
      .groups = "drop"
    ) |>
    mutate(total_variance = within_variance + (1 + 1 / n_imputations) * between_variance)

  normalizer <- mean(exp(pooled$pooled_estimate))
  pooled |>
    mutate(
      relative_trend = exp(pooled_estimate) / normalizer,
      lower = exp(pooled_estimate - 1.96 * sqrt(total_variance)) / normalizer,
      upper = exp(pooled_estimate + 1.96 * sqrt(total_variance)) / normalizer,
      model = model_name
    )
}

pool_linear_trend <- function(models, model_name, first_season, last_season) {
  estimates <- bind_rows(lapply(seq_along(models), function(imputation) {
    coefficients <- summary(models[[imputation]])$p.table
    tibble(
      imputation,
      estimate = coefficients["season_centered", "Estimate"],
      variance = coefficients["season_centered", "Std. Error"]^2
    )
  }))

  pooled <- estimates |>
    summarise(
      pooled_estimate = mean(estimate),
      within_variance = mean(variance),
      between_variance = var(estimate),
      n_imputations = n()
    ) |>
    mutate(total_variance = within_variance + (1 + 1 / n_imputations) * between_variance)

  tibble(
    model = model_name,
    first_season = first_season,
    last_season = last_season,
    average_annual_change = exp(pooled$pooled_estimate) - 1,
    lower = exp(pooled$pooled_estimate - 1.96 * sqrt(pooled$total_variance)) - 1,
    upper = exp(pooled$pooled_estimate + 1.96 * sqrt(pooled$total_variance)) - 1
  )
}

trend_curves <- bind_rows(
  pool_trend_curve(primary_trend_models$smooth, "M4 post-transition"),
  pool_trend_curve(subset_weather_trend_models$smooth, "M4-S common dates"),
  pool_trend_curve(configuration_trend_models$smooth, "M5 + bush configuration"),
  pool_trend_curve(playback_trend_models$smooth, "M6 + playback")
)
trend_summary <- bind_rows(
  pool_linear_trend(primary_trend_models$linear, "M4 post-transition", 1977, 2023),
  pool_linear_trend(subset_weather_trend_models$linear, "M4-S common dates", 1977, 2014),
  pool_linear_trend(configuration_trend_models$linear, "M5 + bush configuration", 1977, 2014),
  pool_linear_trend(playback_trend_models$linear, "M6 + playback", 1977, 2014)
)

# Build separate story figures -------------------------------------------

daily_counts <- primary_data[[1]] |>
  transmute(total_birds_ringed) |>
  arrange(total_birds_ringed) |>
  mutate(
    night_percentile = 100 * row_number() / n(),
    cumulative_bird_percent = 100 * cumsum(total_birds_ringed) / sum(total_birds_ringed)
  )

median_catch <- median(daily_counts$total_birds_ringed)
maximum_catch <- max(daily_counts$total_birds_ringed)
top_five_cutoff <- quantile(daily_counts$total_birds_ringed, 0.95)
top_five_share <- daily_counts |>
  slice_max(total_birds_ringed, n = ceiling(0.05 * nrow(daily_counts))) |>
  summarise(share = sum(total_birds_ringed) / sum(daily_counts$total_birds_ringed)) |>
  pull(share)

dark_theme <- ngulia_theme(dark = TRUE)

distribution_plot <- ggplot(daily_counts, aes(total_birds_ringed)) +
  geom_histogram(binwidth = 100, boundary = 0, fill = "#35546D", colour = "#0B1320", linewidth = 0.25) +
  geom_histogram(
    data = filter(daily_counts, total_birds_ringed >= top_five_cutoff),
    binwidth = 100, boundary = 0, fill = "#FFB84D", colour = "#0B1320", linewidth = 0.25
  ) +
  geom_vline(xintercept = median_catch, colour = "#58D6E7", linewidth = 0.9) +
  annotate(
    "label", x = median_catch + 55, y = Inf, label = paste0("Median\n", comma(median_catch), " birds"),
    hjust = 0, vjust = 1.1, colour = "#0B1320", fill = "#58D6E7", linewidth = 0, size = 3.5
  ) +
  annotate(
    "text", x = 2450, y = Inf, label = "BUSIEST 5% OF NIGHTS", hjust = 0.5,
    vjust = 1.5, colour = "#FFB84D", fontface = "bold", size = 3.4
  ) +
  scale_x_continuous(labels = label_comma(), breaks = seq(0, 3000, 500), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Most nights have modest catches", x = "Birds ringed in one night", y = "Number of nights") +
  dark_theme +
  theme(panel.grid.major.x = element_blank())

stat_cards <- ggplot() +
  annotate("text", x = 0, y = 3.0, label = comma(nrow(daily_counts)), hjust = 0, colour = "#58D6E7", fontface = "bold", size = 10) +
  annotate("text", x = 0, y = 2.55, label = "positive-catch nights", hjust = 0, colour = "#B9C8D4", size = 4.2) +
  annotate("segment", x = 0, xend = 1, y = 2.15, yend = 2.15, colour = "#26384A", linewidth = 0.6) +
  annotate("text", x = 0, y = 1.55, label = comma(maximum_catch), hjust = 0, colour = "#FFB84D", fontface = "bold", size = 10) +
  annotate("text", x = 0, y = 1.1, label = "birds on the busiest night", hjust = 0, colour = "#B9C8D4", size = 4.2) +
  annotate("text", x = 0, y = 0.35, label = paste0(percent(top_five_share, accuracy = 1), " of all birds"), hjust = 0, colour = "#FFB84D", fontface = "bold", size = 5.5) +
  annotate("text", x = 0, y = 0.0, label = "came from the busiest 5% of nights", hjust = 0, colour = "#B9C8D4", size = 3.8) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.2, 3.35), clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA))

concentration_plot <- ggplot(daily_counts, aes(night_percentile, cumulative_bird_percent)) +
  annotate("rect", xmin = 95, xmax = 100, ymin = -Inf, ymax = Inf, fill = "#FFB84D", alpha = 0.15) +
  geom_area(fill = "#58D6E7", alpha = 0.16) +
  geom_line(colour = "#58D6E7", linewidth = 1.05) +
  annotate(
    "segment", x = 95, xend = 95, y = 0, yend = 100 * (1 - top_five_share),
    colour = "#FFB84D", linetype = "dashed", linewidth = 0.7
  ) +
  annotate(
    "segment", x = 95, xend = 100, y = 100 * (1 - top_five_share), yend = 100 * (1 - top_five_share),
    colour = "#FFB84D", linetype = "dashed", linewidth = 0.7
  ) +
  annotate(
    "label", x = 94, y = 47,
    label = paste0("The final 5% of nights\nadd ", percent(top_five_share, accuracy = 1), " of all birds"),
    hjust = 1, colour = "#0B1320", fill = "#FFB84D", linewidth = 0, size = 4, fontface = "bold"
  ) +
  scale_x_continuous(labels = label_percent(scale = 1), breaks = c(0, 25, 50, 75, 95, 100)) +
  scale_y_continuous(labels = label_percent(scale = 1), breaks = seq(0, 100, 25)) +
  labs(
    title = "A small number of exceptional nights carry a large share of the season",
    subtitle = "Nights are ordered from the smallest catch to the largest",
    x = "Recorded nights, from quietest to busiest",
    y = "Cumulative share of all birds"
  ) +
  dark_theme

count_distribution_plot <- (distribution_plot | stat_cards) / concentration_plot +
  plot_layout(heights = c(1, 0.9), widths = c(2.25, 1)) +
  plot_annotation(
    title = "Daily catch is strongly right-skewed",
    subtitle = "Most ringing nights record hundreds of birds; a few exceptional nights record thousands.",
    caption = paste0("Ngulia non-swallow catch · ", nrow(primary_data[[1]]),
      " positive-catch dates · 1976 and 1994–1995 transitions excluded from modeling"),
    theme = theme(
      plot.background = element_rect(fill = "#0B1320", colour = NA),
      plot.title = element_text(colour = "#F5F8FA", face = "bold", size = 23),
      plot.subtitle = element_text(colour = "#B9C8D4", size = 13),
      plot.caption = element_text(colour = "#8295A5", size = 9, hjust = 0)
    )
  )

progression_plot_data <- model_progression |>
  transmute(
    model = factor(model, levels = model),
    `Held-out deviance reduction` = 100 * cv_deviance_reduction,
    `Held-out log RMSE (lower is better)` = log_rmse
  ) |>
  pivot_longer(-model, names_to = "metric", values_to = "value")

progression_plot <- ggplot(progression_plot_data, aes(model, value, group = 1)) +
  geom_line(colour = "#65777D", linewidth = 0.7) +
  geom_point(colour = ngulia_colours[["teal"]], size = 3) +
  geom_text(aes(label = if_else(value > 10, sprintf("%.1f%%", value), sprintf("%.2f", value))), vjust = -0.7, size = 3.3) +
  facet_wrap(vars(metric), scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = c("M1\nTiming", "M2\n+ Mist/rain", "M3\n+ Moon", "M4\n+ Weather")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) +
  labs(
    title = "Each covariate block improves prediction of unseen seasons",
    subtitle = "Complete seasons are omitted from model fitting and then predicted",
    x = NULL,
    y = NULL
  ) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), strip.text = element_text(face = "bold"))

plot_core_effect <- function(variable, title, x_label) {
  ggplot(filter(core_effects, .data$variable == .env$variable), aes(value, relative_effect)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "#71838B") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ngulia_colours[["pale_teal"]], alpha = 0.3) +
    geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.9) +
    geom_rug(data = tibble(value = primary_data[[1]][[variable]]), aes(x = value),
      inherit.aes = FALSE, sides = "b", colour = "#415963", alpha = 0.22, linewidth = 0.2) +
    labs(
      title = title,
      subtitle = "Adjusted association; rug marks show observed dates (overlap appears darker)",
      x = x_label,
      y = "Relative expected positive catch"
    ) +
    ngulia_theme(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

season_timing_plot <- plot_core_effect("season_day", "M1: catch changes strongly through the season", "Day since 20 October")
rain_plot <- plot_core_effect("era5_rain_log", "M2: rainfall has a nonlinear association with catch", "log(1 + rainfall from 00:00–08:00, mm)")
moon_plot <- plot_core_effect("moon_distance_from_new_moon", "M3: catch decreases away from new moon", "Days from new moon")
wind_plot <- plot_core_effect("wind_speed_10m_mean_ms", "M4: wind adds a comparatively small association", "Mean wind speed from 00:00–08:00 (m/s)")
temperature_plot <- plot_core_effect("temperature_2m_mean_c", "M4: warmer mornings are associated with higher catch", "Mean temperature from 00:00–08:00 (°C)")
pressure_plot <- plot_core_effect("surface_pressure_mean_hpa", "M4: pressure has a shallow nonlinear association", "Mean surface pressure from 00:00–08:00 (hPa)")

mist_plot <- ggplot(mist_effects, aes(label, relative_effect)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#71838B") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, colour = ngulia_colours[["teal"]]) +
  geom_point(colour = ngulia_colours[["teal"]], size = 3) +
  geom_text(data = primary_data[[1]] |> count(mist_state) |>
      mutate(label = factor(mist_state, levels = mist_levels,
        labels = c("No mist", "Light/patchy", "Good mist"))),
    aes(x = label, y = min(mist_effects$lower) * 0.72, label = paste0("n=", n)),
    inherit.aes = FALSE, colour = "#415963", size = 3.2) +
  labs(
    title = "M2: mist is the strongest measured daily association",
    subtitle = "No mist is the reference; n shows dates in one mist-state imputation",
    x = NULL,
    y = "Relative expected positive catch"
  ) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

team_by_season <- daily_coverage |>
  filter(ringing_happened, season >= 1977, season != 1994, !is.na(djp_team_size_minimum)) |>
  group_by(season) |>
  summarise(
    q25 = quantile(djp_team_size_minimum, 0.25),
    median = median(djp_team_size_minimum),
    q75 = quantile(djp_team_size_minimum, 0.75),
    n = n(), .groups = "drop"
  )
team_plot <- ggplot(team_by_season, aes(season, median)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = ngulia_colours[["pale_teal"]], alpha = 0.4) +
  geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.8) +
  geom_point(aes(size = n), colour = ngulia_colours[["teal"]], alpha = 0.8) +
  scale_size_continuous(name = "Recorded dates", range = c(1.5, 4)) +
  labs(title = "Recorded team size increased over time",
    subtitle = "Median and middle half of daily minimum counts; descriptive only",
    x = "Season", y = "Minimum recorded team size") +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank())

net_periods <- daily_coverage |>
  filter(season >= 1977, season <= 2014) |>
  distinct(season, bush_net_configuration) |>
  mutate(bush_net_configuration = factor(bush_net_configuration,
    levels = c("back_bush", "transition", "front_bush"),
    labels = c("Back bush", "1994–1995 transition", "Front bush")))
net_period_plot <- ggplot(net_periods, aes(season, 1, fill = bush_net_configuration)) +
  geom_tile(width = 0.95, height = 0.7) +
  scale_fill_manual(values = c("Back bush" = ngulia_colours[["blue"]],
    "1994–1995 transition" = ngulia_colours[["gold"]],
    "Front bush" = ngulia_colours[["teal"]]), name = "Bush position") +
  scale_x_continuous(limits = c(1976.5, 2014.5), breaks = seq(1980, 2010, 5)) +
  labs(title = "Bush-net position changed once", x = NULL, y = NULL) +
  ngulia_theme(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "bottom")

net_use_by_season <- daily_coverage |>
  filter(season >= 1977, season <= 2014, ringing_happened,
    !is.na(net_sites_observed), net_sites_observed != "none") |>
  transmute(
    season,
    back_bush = grepl("back_bush", net_sites_observed, fixed = TRUE),
    front_bush = grepl("front_bush", net_sites_observed, fixed = TRUE),
    outside_night_nets = grepl("outside_night_nets", net_sites_observed, fixed = TRUE)
  ) |>
  group_by(season) |>
  summarise(
    n_dates = n(),
    across(c(back_bush, front_bush, outside_night_nets), sum),
    .groups = "drop"
  ) |>
  pivot_longer(c(back_bush, front_bush, outside_night_nets), names_to = "site", values_to = "n_used") |>
  mutate(
    share = n_used / n_dates,
    site = recode(site, back_bush = "Back bush", front_bush = "Front bush",
      outside_night_nets = "Night nets")
  )
net_date_range <- range(distinct(net_use_by_season, season, n_dates)$n_dates)
net_use_plot <- ggplot(net_use_by_season, aes(season, share, colour = site)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_dates), alpha = 0.85) +
  scale_colour_manual(values = c("Back bush" = ngulia_colours[["blue"]],
    "Front bush" = ngulia_colours[["teal"]], "Night nets" = ngulia_colours[["purple"]]),
    name = "Recorded site") +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
  scale_x_continuous(limits = c(1976.5, 2014.5), breaks = seq(1980, 2010, 5)) +
  scale_size_continuous(name = "Dates in denominator", range = c(1.4, 4)) +
  guides(size = guide_legend(override.aes = list(colour = ngulia_palette()[["muted"]]))) +
  labs(title = "Daily operation still varied within each period",
    subtitle = paste0("Share among positive-catch dates with a recorded operating site (",
      net_date_range[[1]], "–", net_date_range[[2]], " dates per season); categories can overlap"),
    x = "Season", y = "Share of recorded-site dates") +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
net_position_plot <- net_period_plot / net_use_plot + plot_layout(heights = c(0.45, 1))

playback_by_season <- operation_source |>
  group_by(season, bush_net_configuration) |>
  summarise(n_dates = n(), n_playback = sum(playback_nocturnal_observed == 1),
    share = n_playback / n_dates, .groups = "drop") |>
  mutate(bush_net_configuration = recode(bush_net_configuration,
    back_bush = "Back bush", front_bush = "Front bush"))
playback_overlap <- operation_source |>
  group_by(bush_net_configuration) |>
  summarise(n_dates = n(), n_playback = sum(playback_nocturnal_observed == 1),
    share = n_playback / n_dates, .groups = "drop")
playback_overlap_plot <- ggplot(playback_by_season,
    aes(season, share, fill = bush_net_configuration)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c("Back bush" = ngulia_colours[["blue"]],
    "Front bush" = ngulia_colours[["teal"]]), name = "Bush position") +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
  labs(title = "Playback is concentrated in later seasons",
    subtitle = "Share of comparable positive-catch dates with recorded nocturnal playback",
    x = "Season", y = "Dates with playback") +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

bush_effect_plot_data <- bind_rows(
  tibble(label = "Back bush (reference)", relative_effect = 1, lower = 1, upper = 1),
  filter(operation_binary_effects, variable == "bush_period") |>
    transmute(label = "Front bush", relative_effect, lower, upper)
)
bush_effect_plot <- ggplot(bush_effect_plot_data,
    aes(relative_effect, factor(label, levels = c("Front bush", "Back bush (reference)")))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#71838B") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.13, colour = "#A86509") +
  geom_point(colour = "#A86509", size = 2.8) +
  labs(
    title = "M5: front versus back bush",
    subtitle = "A period contrast estimated under the smooth-year assumption",
    x = "Relative expected positive catch",
    y = NULL
  ) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank())

m6_effect_plot_data <- bind_rows(
  tibble(label = c("Back bush (reference)", "No playback (reference)"),
    relative_effect = 1, lower = 1, upper = 1),
  filter(operation_binary_effects, variable == "bush_period_m6") |>
    transmute(label = "Front bush", relative_effect, lower, upper),
  filter(operation_binary_effects, variable == "playback_used") |>
    transmute(label = "Playback used", relative_effect, lower, upper)
)
m6_effect_plot <- ggplot(m6_effect_plot_data,
    aes(relative_effect, factor(label, levels = c("Playback used", "No playback (reference)",
      "Front bush", "Back bush (reference)")))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#71838B") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.13,
    colour = ngulia_colours[["purple"]]) +
  geom_point(colour = ngulia_colours[["purple"]], size = 2.8) +
  labs(title = "M6: playback and the remaining bush-period contrast",
    subtitle = "Both adjusted associations are conditional and share historical variation",
    x = "Relative expected positive catch", y = NULL) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank())

operation_validation_plot <- operation_comparison |>
  transmute(model = factor(model, levels = model),
    `Held-out log RMSE` = log_rmse,
    `Held-out deviance reduction` = cv_deviance_reduction) |>
  pivot_longer(-model, names_to = "metric", values_to = "value") |>
  ggplot(aes(model, value, group = 1)) +
  geom_line(colour = "#65777D", linewidth = 0.7) +
  geom_point(colour = ngulia_colours[["teal"]], size = 3) +
  facet_wrap(vars(metric), scales = "free_y") +
  labs(title = "M4-S, M5 and M6 prediction on the same dates",
    subtitle = "Entire seasons held out; lower log RMSE and higher deviance reduction are better",
    x = NULL, y = NULL) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

trend_plot <- ggplot() +
  geom_ribbon(
    data = trend_curves,
    aes(season, ymin = lower, ymax = upper, fill = model),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(data = trend_curves, aes(season, relative_trend, colour = model), linewidth = 1) +
  geom_errorbar(
    data = filter(annual_index, model == "M4 post-transition"),
    aes(season, ymin = lower, ymax = upper),
    colour = ngulia_colours[["teal"]],
    alpha = 0.25,
    width = 0
  ) +
  geom_point(
    data = filter(annual_index, model == "M4 post-transition"),
    aes(season, relative_index, size = n_positive_catch_dates),
    colour = ngulia_colours[["teal"]],
    alpha = 0.65
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#71838B") +
  scale_colour_manual(values = c(
    "M4 post-transition" = ngulia_colours[["teal"]],
    "M4-S common dates" = "#65777D",
    "M5 + bush configuration" = "#A86509",
    "M6 + playback" = ngulia_colours[["purple"]]
  )) +
  scale_fill_manual(values = c(
    "M4 post-transition" = ngulia_colours[["pale_teal"]],
    "M4-S common dates" = "#AAB7BC",
    "M5 + bush configuration" = "#E7A84B",
    "M6 + playback" = "#D6C5E8"
  )) +
  scale_size_continuous(range = c(1, 3.3)) +
  labs(
    title = "Annual pattern under successive assumptions",
    subtitle = "Points are M4 annual effects; smooth curves are normalized within each fitted sample",
    x = "Season",
    y = "Adjusted positive-catch index",
    colour = NULL,
    fill = NULL,
    size = "Positive-catch dates"
  ) +
  ngulia_theme(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# Create compact result tables -------------------------------------------

model_summary <- bind_rows(
  model_progression |>
    transmute(
      model,
      scope = paste0("Post-transition: ", nrow(primary_data[[1]]),
        " dates, ", n_distinct(primary_data[[1]]$season), " seasons"),
      added_terms,
      heldout_deviance_reduction = cv_deviance_reduction,
      heldout_log_rmse = log_rmse
    ),
  operation_comparison |>
    transmute(
      model,
      scope = paste0("Comparable post-1976 operations: ", nrow(operation_data[[1]]),
        " dates, ", n_distinct(operation_data[[1]]$season), " seasons"),
      added_terms = case_when(
        model == "M4-S Smooth-year baseline" ~ "M4 covariates + smooth year on common dates",
        model == "M5 + bush configuration" ~ "M4 + back/front bush period",
        TRUE ~ "M5 + nocturnal playback"
      ),
      heldout_deviance_reduction = cv_deviance_reduction,
      heldout_log_rmse = log_rmse
    )
)

effect_summary <- core_effects |>
  group_by(variable, label) |>
  summarise(
    central_observed_range = paste0(signif(min(value), 3), " to ", signif(max(value), 3)),
    fitted_relative_effect_range = paste0(signif(min(relative_effect), 3), " to ", signif(max(relative_effect), 3)),
    .groups = "drop"
  ) |>
  mutate(label = as.character(label)) |>
  bind_rows(
    mist_effects |>
      transmute(
        variable = "mist_state",
        label = as.character(label),
        central_observed_range = NA_character_,
        fitted_relative_effect_range = paste0(signif(relative_effect, 3), " (", signif(lower, 3), "–", signif(upper, 3), ")")
      ),
    operation_binary_effects |>
      transmute(
        variable = "operation_binary",
        label,
        central_observed_range = if_else(grepl("bush_period", variable), "back vs front", "absent vs used"),
        fitted_relative_effect_range = paste0(signif(relative_effect, 3), " (", signif(lower, 3), "–", signif(upper, 3), ")")
      )
  )

# Build one explanatory HTML report --------------------------------------

html_table <- function(data) {
  tags$table(
    tags$thead(tags$tr(lapply(names(data), tags$th))),
    tags$tbody(lapply(seq_len(nrow(data)), function(i) {
      tags$tr(lapply(data, function(column) tags$td(as.character(column[[i]]))))
    }))
  )
}

report_model_table <- model_summary |>
  mutate(
    `Held-out deviance reduction` = percent(heldout_deviance_reduction, accuracy = 0.1),
    `Held-out log RMSE` = sprintf("%.3f", heldout_log_rmse)
  ) |>
  dplyr::select(Model = model, Scope = scope, `Added information` = added_terms, `Held-out deviance reduction`, `Held-out log RMSE`)

report_trend_table <- trend_summary |>
  transmute(
    Model = model,
    Period = paste(first_season, last_season, sep = "–"),
    `Linear-year annual change` = paste0(
      percent(average_annual_change, accuracy = 0.1),
      " (95% interval ", percent(lower, accuracy = 0.1), " to ", percent(upper, accuracy = 0.1), ")"
    )
  )

report_playback_table <- playback_overlap |>
  transmute(
    `Bush period` = recode(bush_net_configuration,
      back_bush = "Back bush, 1977–1993", front_bush = "Front bush, 1996–2014"),
    `Comparable dates` = n_dates,
    `Playback dates` = n_playback,
    `Playback share` = percent(share, accuracy = 1)
  )

m4_result <- filter(model_progression, primary_index_model)
m4_subset <- filter(operation_comparison, model == "M4-S Smooth-year baseline")
m5_result <- filter(operation_comparison, model == "M5 + bush configuration")
m6_result <- filter(operation_comparison, model == "M6 + playback")
primary_change <- filter(trend_summary, model == "M4 post-transition")
subset_weather_change <- filter(trend_summary, model == "M4-S common dates")
m5_change <- filter(trend_summary, model == "M5 + bush configuration")
m6_change <- filter(trend_summary, model == "M6 + playback")
m5_effect <- filter(operation_binary_effects, variable == "bush_period")
m6_bush_effect <- filter(operation_binary_effects, variable == "bush_period_m6")
m6_playback_effect <- filter(operation_binary_effects, variable == "playback_used")
operation_year_correlations <- read_csv(
  file.path(analysis_dir, "model_data", "documented_operation_model_data.csv"),
  show_col_types = FALSE
) |>
  summarise(
    team_size = cor(season, djp_team_size_minimum, method = "spearman", use = "complete.obs"),
    playback = cor(season, playback_nocturnal_observed, method = "spearman")
  )

model_report <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$title("Ngulia daily count analysis"),
    tags$style(HTML(paste(readLines(file.path(project_dir, "assets", "report.css")), collapse = "\n")))
  ),
  tags$body(
    tags$h1("What explains daily catch at Ngulia?"),
    tags$p(class = "lede", "A step-by-step model of daily non-swallow catch, followed by the adjusted annual pattern."),
    tags$div(class = "note", tags$b("Response: "), "non-swallow catch on dates with at least one recorded capture. Incomplete operated zeroes and net-hours mean that this is positive-catch intensity, not abundance."),

    tags$h2("The daily counts are strongly skewed"),
    tags$p("Most positive-catch nights contain hundreds of birds, while a small number contain thousands. The count model must preserve those exceptional nights while allowing variability to increase with expected catch."),
    tags$img(class = "figure wide", src = "figures/00_daily_count_distribution.png", alt = "Dark infographic showing the right-skewed distribution of positive daily bird catches"),

    tags$h2("1. Models M1–M4 and their adjusted associations"),
    tags$div(class = "equation", HTML("Y<sub>d,y</sub> &sim; Negative binomial(&mu;<sub>d,y</sub>, &theta;); fitted on positive-catch dates")),
    tags$p(HTML("Let &eta;<sub>d,y</sub> = log(&mu;<sub>d,y</sub>). Here d indexes a date and y its season. The negative-binomial dispersion parameter &theta; allows the variance to exceed the mean. Every s<sub>k</sub>(x) below is a penalized smooth; k is its basis dimension, which limits flexibility but is not the number of fitted bends.")),
    tags$p("Each effect figure in M1–M4 is calculated from the final M4 fit with the other covariates held fixed, and is placed beside the model stage where that variable first enters. Rug marks and category counts show the observed support; bands pool coefficient and mist-imputation uncertainty."),

    tags$h3("M1 — seasonal timing"),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M1</sup> = &alpha;<sub>y</sub> + s<sub>12</sub>(season day<sub>d</sub>)")),
    tags$p(HTML("Season day is the number of days since 20 October and is fitted as a nonlinear smooth. The unrestricted season intercept &alpha;<sub>y</sub> gives every year its own catch level in the final annual-index fit.")),
    tags$p("Expected catch peaks near the middle of the ringing season and is much lower near either end. Rug marks show the observed season days."),
    tags$img(class = "figure", src = "figures/02_season_timing_effect.png", alt = "Adjusted seasonal timing association"),

    tags$h3("M2 — seasonal timing, mist and rain"),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M2</sup> = &eta;<sub>d,y</sub><sup>M1</sup> + &gamma;<sub>light</sub>I(light/patchy mist) + &gamma;<sub>good</sub>I(good mist) + s<sub>6</sub>[log(1 + rain<sub>d</sub>)]")),
    tags$p("Mist is a three-level categorical effect with no mist as the reference. Observed mist fixes the category. Where mist is unobserved, the category is drawn from ERA5-calibrated probabilities based on cloud cover, cloud-base height, humidity and zonal wind; estimates pool 20 completed datasets. ERA5 rainfall from 00:00–08:00 local time is transformed as log(1 + mm) and then fitted nonlinearly."),
    tags$p(paste0(
      "Relative to no mist, expected catch is ",
      sprintf("%.1f", mist_effects$relative_effect[mist_effects$mist_state == "light_patchy"]),
      " times higher in light or patchy mist and ",
      sprintf("%.1f", mist_effects$relative_effect[mist_effects$mist_state == "good"]),
      " times higher in good mist. Counts beneath the categories show their representation."
    )),
    tags$img(class = "figure", src = "figures/03_mist_effect.png", alt = "Adjusted mist association"),
    tags$p("The rainfall association rises from the driest conditions and then changes little across the wetter part of the observed range. Rug marks show the rainfall distribution."),
    tags$img(class = "figure", src = "figures/04_rain_effect.png", alt = "Adjusted rainfall association"),

    tags$h3("M3 — add moon"),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M3</sup> = &eta;<sub>d,y</sub><sup>M2</sup> + s<sub>6</sub>(distance from new moon<sub>d</sub>)")),
    tags$p("Moon is the absolute number of days from new moon, so new moon is 0 and dates approaching full moon have larger values. The smooth permits a nonlinear association and treats waxing and waning dates at the same distance equally."),
    tags$p("Expected catch is highest near new moon and decreases as the date approaches full moon."),
    tags$img(class = "figure", src = "figures/05_moon_effect.png", alt = "Adjusted lunar association"),

    tags$h3("M4 — add local weather"),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M4</sup> = &eta;<sub>d,y</sub><sup>M3</sup> + s<sub>6</sub>(wind speed<sub>d</sub>) + s<sub>6</sub>(temperature<sub>d</sub>) + s<sub>6</sub>(surface pressure<sub>d</sub>)")),
    tags$p("Wind speed, 2-m air temperature and surface pressure are ERA5 means for 00:00–08:00 local time. Each is modeled with its own penalized smooth. Their effects are additive on the log scale; no interactions are included."),
    tags$p("Wind has a comparatively small association over the central observed range."),
    tags$img(class = "figure", src = "figures/06_wind_effect.png", alt = "Adjusted wind association"),
    tags$p("Warmer early-morning conditions are associated with higher positive catch, with wider uncertainty at the warm end."),
    tags$img(class = "figure", src = "figures/07_temperature_effect.png", alt = "Adjusted temperature association"),
    tags$p("Pressure has a shallow U-shaped association and contributes as part of the local-weather block."),
    tags$img(class = "figure", src = "figures/08_pressure_effect.png", alt = "Adjusted surface-pressure association"),

    tags$div(class = "note", HTML("<b>Validation versus final estimation:</b> complete seasons are held out when comparing M1–M4, so the held-out predictions omit &alpha;<sub>y</sub>; an unseen year has no fitted annual coefficient. The final M1–M4 fits include &alpha;<sub>y</sub> to estimate the adjusted annual index. M4-S, M5 and M6 instead use s<sub>10</sub>(year) and are validated on identical dates and season folds.")),
    tags$img(class = "figure wide", src = "figures/01_model_progression.png", alt = "Held-out performance as covariate blocks are added"),
    html_table(report_model_table),

    tags$h2("2. Team size: an exploratory record"),
    tags$p("Recorded team size rose as the project grew. The figure shows its distribution by season, including minimum counts when Earthwatch numbers were unspecified. Team size is not a covariate in M4–M6 because its long-term increase could absorb the pattern under study."),
    tags$img(class = "figure wide", src = "figures/09_team_size_descriptive.png", alt = "Recorded team-size distribution by season"),

    tags$h2("3. M4-S and M5: back versus front bush"),
    tags$p("The top band shows the historical bush-net position. The lower panel shows recorded daily site use, including night nets. Point size gives the number of positive-catch dates in the yearly denominator. Daily opening varies with mist, rain and workload, so these percentages are descriptive rather than a measure of the physical position."),
    tags$img(class = "figure wide", src = "figures/10_net_position_by_year.png", alt = "Historical bush-net position and yearly recorded net-site use"),
    tags$h3("M4-S — smooth-year sensitivity baseline"),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M4-S</sup> = &beta;<sub>0</sub> + s<sub>10</sub>(year<sub>y</sub>) + all M4 daily terms")),
    tags$p(HTML("M4-S uses the common documented-operation subset and replaces the unrestricted &alpha;<sub>y</sub> values with one smooth function of year. This restriction is necessary because a back/front period indicator cannot be estimated separately from a complete set of annual fixed effects.")),
    tags$p("M4-S introduces no new daily covariate effect to plot; its smooth annual pattern is shown in the final model-comparison section."),
    tags$h3("M5 — add bush-net configuration"),
    tags$p("M5 adds one fixed category for the stable back-bush (1977–1993) and front-bush (1996–2014) configurations. Both 1994 and 1995 are excluded as transition years."),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M5</sup> = &eta;<sub>d,y</sub><sup>M4-S</sup> + &beta;<sub>front</sub>I(front bush<sub>y</sub>)")),
    tags$p(HTML("Back bush is the reference, with I(front bush) = 0 and relative effect 1. Front bush has I(front bush) = 1 and a front/back ratio of exp(&beta;<sub>front</sub>). The model estimates one contrast rather than two unrelated coefficients.")),
    tags$img(class = "figure", src = "figures/11_bush_effect.png", alt = "M5 back versus front bush contrast"),
    tags$p(paste0("A fixed bush category cannot be separated from unrestricted annual fixed effects. M4-S, M5 and M6 therefore use a smooth year term on the same ", nrow(operation_data[[1]]), " dates. The size of the bush-period contrast depends on that smoothness assumption and also contains any other change coincident with the move.")),
    tags$p(class = "result", paste0("Held-out log RMSE: M4-S ", sprintf("%.3f", m4_subset$log_rmse),
      "; M5 ", sprintf("%.3f", m5_result$log_rmse), ".")),
    tags$p(paste0("The fitted front/back contrast is ", sprintf("%.2f", m5_effect$relative_effect),
      " (95% interval ", sprintf("%.2f", m5_effect$lower), "–", sprintf("%.2f", m5_effect$upper),
      "). Its broad interval and the slightly worse held-out scores do not establish a configuration effect.")),

    tags$h2("4. M6: playback and its overlap with net position"),
    tags$p("Playback is not spread evenly across the two bush periods. The yearly bars show when it was recorded on the comparable dates."),
    tags$img(class = "figure wide", src = "figures/12_playback_by_year.png", alt = "Yearly playback frequency in each bush-net period"),
    html_table(report_playback_table),
    tags$p("All 38 workbook playback rows in the back-bush subset occur in 1993 and carry the same behind-lodge tape code. The historical synthesis independently confirms experiments began that year in the southern bush. This is one season-wide exposure pattern, not 38 independent seasons. Playback is more common in the stable front-bush subset and varies within later seasons, so separation remains limited."),
    tags$p("M6 adds documented nocturnal playback to M5. The workbook distinguishes tapes at the front from tapes behind the lodge, but does not identify night-net versus bush-net speakers. Their locations remain unresolved and are not separate model levels."),
    tags$div(class = "equation", HTML("&eta;<sub>d,y</sub><sup>M6</sup> = &eta;<sub>d,y</sub><sup>M5</sup> + &beta;<sub>play</sub>I(playback<sub>d</sub>)")),
    tags$p(HTML("No documented playback is the reference, with relative effect 1. A playback date has I(playback) = 1 and a playback/no-playback ratio of exp(&beta;<sub>play</sub>). The figure also repeats the conditional front/back contrast from M6.")),
    tags$img(class = "figure", src = "figures/13_m6_effects.png", alt = "M6 playback and bush-period contrasts"),
    tags$p(class = "result", paste0("Held-out log RMSE: M5 ", sprintf("%.3f", m5_result$log_rmse),
      "; M6 ", sprintf("%.3f", m6_result$log_rmse), ".")),
    tags$p(paste0("The M6 front/back contrast is ", sprintf("%.2f", m6_bush_effect$relative_effect),
      " after accounting for playback (M5: ", sprintf("%.2f", m5_effect$relative_effect),
      "). This change is a diagnostic of overlap, not evidence that one variable caused the other's association.")),

    tags$h2("5. Model comparison and annual pattern"),
    tags$p("M4-S, M5 and M6 use identical dates and season folds. Log RMSE improves with playback, while held-out deviance reduction worsens with each added term. The scores do not select a unique adjustment."),
    tags$img(class = "figure wide", src = "figures/14_operations_validation.png", alt = "Held-out comparison of M4-S, M5 and M6"),
    tags$p("The plot compares M4 on all post-transition dates with M4-S, M5 and M6 on the common documented-operation subset. The common-date curves isolate the added assumptions; the wider M4 curve also reflects a longer observation period. Each curve is normalized within its own sample, so slopes and shape are more comparable than absolute vertical position."),
    tags$img(class = "figure wide", src = "figures/15_adjusted_trend.png", alt = "Adjusted annual fluctuations and trends"),
    tags$p("The curves use flexible smooth-year models. The table below is a separate linear-year sensitivity, fitted to express one annual percentage. A step effect and a straight-line trend compete strongly for the same historical change; the table's slopes need not match the local shape of the smooth curves."),
    html_table(report_trend_table),
    tags$p(class = "result", paste0("On common dates, estimated annual change is ",
      percent(subset_weather_change$average_annual_change, accuracy = 0.1), " for M4, ",
      percent(m5_change$average_annual_change, accuracy = 0.1), " for M5 and ",
      percent(m6_change$average_annual_change, accuracy = 0.1), " for M6.")),
    tags$p(class = "note", HTML(paste0(
      "<b>These are sensitivity trends, not identified abundance trends.</b> M4's positive linear slope weakens and its interval crosses zero once a bush-period step is allowed, and remains near zero after playback is added. Yet M5 does not improve held-out prediction, so the weaker slope cannot be credited specifically to net position. The M5 step can also capture unrelated changes near 1994–1995, and M6 playback may be selected in response to conditions. Correlations with year are ",
      sprintf("%.2f", operation_year_correlations$team_size), " for team size and ",
      sprintf("%.2f", operation_year_correlations$playback), " for playback. No model can uniquely allocate the historical change between birds and catchability."
    ))),

    tags$h2("6. Conclusion"),
    tags$p("M4 describes adjusted positive-catch intensity. M5 shows sensitivity to a back/front bush period contrast under a smooth-year assumption; M6 adds recorded playback. Team size remains descriptive. Agreement or disagreement among the trend curves describes model sensitivity, not proof of migrant abundance change."),
    tags$p(paste0(
      "Checks: pooled Pearson dispersion ", sprintf("%.2f", diagnostics$pearson_dispersion),
      "; fitted-versus-observed log-count correlation ", sprintf("%.2f", diagnostics$correlation_fitted_observed_log), "."
    )),
    tags$p("Tables: ",
      tags$a(href = "tables/daily_count_model_summary.csv", "models"), " · ",
      tags$a(href = "tables/daily_count_covariate_dictionary.csv", "variables"), " · ",
      tags$a(href = "tables/daily_count_covariate_coverage.csv", "coverage"), " · ",
      tags$a(href = "tables/daily_count_net_site_use_by_season.csv", "net-site evidence"), " · ",
      tags$a(href = "tables/daily_count_covariate_effect_summary.csv", "effects"), " · ",
      tags$a(href = "tables/daily_count_adjusted_annual_index.csv", "annual index"), " · ",
      tags$a(href = "tables/daily_count_trend_summary.csv", "trends")
    )
  )
)

# Write outputs and remove superseded report products --------------------

story_figures <- list(
  "00_daily_count_distribution.png" = list(count_distribution_plot, 12, 9),
  "01_model_progression.png" = list(progression_plot, 10, 5.6),
  "02_season_timing_effect.png" = list(season_timing_plot, 8, 5),
  "03_mist_effect.png" = list(mist_plot, 8, 5),
  "04_rain_effect.png" = list(rain_plot, 8, 5),
  "05_moon_effect.png" = list(moon_plot, 8, 5),
  "06_wind_effect.png" = list(wind_plot, 8, 5),
  "07_temperature_effect.png" = list(temperature_plot, 8, 5),
  "08_pressure_effect.png" = list(pressure_plot, 8, 5),
  "09_team_size_descriptive.png" = list(team_plot, 10, 5.5),
  "10_net_position_by_year.png" = list(net_position_plot, 10, 7),
  "11_bush_effect.png" = list(bush_effect_plot, 8, 4.5),
  "12_playback_by_year.png" = list(playback_overlap_plot, 10, 5.5),
  "13_m6_effects.png" = list(m6_effect_plot, 8, 5),
  "14_operations_validation.png" = list(operation_validation_plot, 10, 5.2),
  "15_adjusted_trend.png" = list(trend_plot, 11, 6.5)
)
for (file_name in names(story_figures)) {
  figure <- story_figures[[file_name]]
  ngulia_save(file.path(figure_dir, file_name), figure[[1]], width = figure[[2]], height = figure[[3]], dpi = 220)
}
write_csv(model_summary, file.path(table_dir, "daily_count_model_summary.csv"))
write_csv(net_use_by_season, file.path(table_dir, "daily_count_net_site_use_by_season.csv"))
write_csv(effect_summary, file.path(table_dir, "daily_count_covariate_effect_summary.csv"))
write_csv(annual_index, file.path(table_dir, "daily_count_adjusted_annual_index.csv"))
write_csv(trend_summary, file.path(table_dir, "daily_count_trend_summary.csv"))
rendered_report <- renderTags(model_report)
report_html <- sub(
  "<html>",
  paste0("<!doctype html>\n<html>\n<head>\n", rendered_report$head, "\n</head>"),
  rendered_report$html,
  fixed = TRUE
)
writeLines(report_html, file.path(analysis_dir, "daily_count_analysis.html"))

unlink(file.path(analysis_dir, "daily_count_model_report.html"))
unlink(file.path(analysis_dir, "covariate_exploration"), recursive = TRUE)
unlink(file.path(figure_dir, c(
  "daily_count_model_story.png",
  "adjusted_annual_positive_catch_index.png",
  "daily_count_covariate_effects.png",
  "daily_count_mist_effect.png",
  "daily_count_model_progression.png",
  "documented_operations_adjusted_effects.png",
  "documented_operations_coverage.png",
  "documented_operations_model_comparison.png",
  "10_playback_net_effects.png",
  "11_operations_evidence.png",
  "10_configuration_playback_effects.png",
  "11_operations_validation.png",
  "12_adjusted_trend.png"
)))
unlink(file.path(model_dir, c(
  "core_weather_gam_mi.rds",
  "daily_count_trend_m4_gam_mi.rds",
  "daily_count_trend_m5_gam_mi.rds"
)))
unlink(file.path(figure_dir, "legacy_2014"), recursive = TRUE)
unlink(file.path(table_dir, c(
  "adjusted_annual_positive_catch_index.csv",
  "adjusted_positive_catch_mist_effect.csv",
  "adjusted_positive_catch_parametric_terms.csv",
  "adjusted_positive_catch_smooth_terms.csv",
  "annual_index_sample_size.csv",
  "documented_operations_parametric_terms.csv",
  "documented_operations_smooth_terms.csv"
)))

cli_alert_success("Wrote daily-count story report to {file.path(analysis_dir, 'daily_count_analysis.html')}")
