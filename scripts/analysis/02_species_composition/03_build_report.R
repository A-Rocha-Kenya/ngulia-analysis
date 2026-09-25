library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(nnet)
library(splines)
library(scales)
library(htmltools)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "02_species_composition")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- ngulia_figure_dir(file.path(analysis_dir, "figures"))

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read model products ----------------------------------------------------

cli_h1("Build species-composition analysis report")

primary_bundle <- readRDS(file.path(model_dir, "joint_composition_primary.rds"))
primary_model <- primary_bundle$model
categories <- primary_bundle$categories |>
  arrange(category_order)

model_data <- read_csv(
  file.path(model_data_dir, "composition_model_data.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
model_comparison <- read_csv(
  file.path(table_dir, "composition_model_comparison.csv"),
  show_col_types = FALSE
)
tempered_mean_structure_comparison <- read_csv(
  file.path(table_dir, "tempered_mean_structure_comparison.csv"),
  show_col_types = FALSE
)
diagnostics <- read_csv(
  file.path(table_dir, "primary_model_diagnostics.csv"),
  show_col_types = FALSE
)
coverage_by_season <- read_csv(
  file.path(table_dir, "composition_coverage_by_season.csv"),
  show_col_types = FALSE
)
day_support <- read_csv(
  file.path(table_dir, "day_of_season_support.csv"),
  show_col_types = FALSE
)
species_support <- read_csv(
  file.path(table_dir, "species_model_support.csv"),
  show_col_types = FALSE
)
untempered_diagnostics <- read_csv(
  file.path(table_dir, "untempered_model_diagnostics.csv"),
  show_col_types = FALSE
)
high_catch_model_comparison <- read_csv(
  file.path(table_dir, "high_catch_model_comparison.csv"),
  show_col_types = FALSE
)
count_distribution_scores <- read_csv(
  file.path(table_dir, "count_distribution_scores.csv"),
  show_col_types = FALSE
)
high_catch_trend_sensitivity <- read_csv(
  file.path(table_dir, "high_catch_trend_sensitivity.csv"),
  show_col_types = FALSE
)
high_catch_endpoint_comparison <- read_csv(
  file.path(table_dir, "high_catch_endpoint_comparison.csv"),
  show_col_types = FALSE
)
top_date_influence <- read_csv(
  file.path(table_dir, "top_date_influence.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
catch_concentration <- read_csv(
  file.path(table_dir, "catch_concentration.csv"),
  show_col_types = FALSE
)
analysis_period_decision <- read_csv(
  file.path(table_dir, "analysis_period_decision.csv"),
  show_col_types = FALSE
)
excluded_early_seasons <- read_csv(
  file.path(table_dir, "excluded_early_seasons.csv"),
  show_col_types = FALSE
)

focal_species <- categories |>
  filter(is_focal_species) |>
  arrange(abundance_rank)
focal_names <- focal_species$common_name
all_category_names <- categories$common_name
seasons <- sort(unique(model_data$season))
trend_start <- max(species_support$first_season[species_support$selected])
trend_end <- min(species_support$last_season[species_support$selected])
trend_start_index <- match(trend_start, seasons)
trend_end_index <- match(trend_end, seasons)
reference_days <- day_support |>
  filter(in_reference_window) |>
  pull(day_of_season)
daily_predictors <- model_data |>
  distinct(ringing_date, .keep_all = TRUE)

# Prediction helpers -----------------------------------------------------

softmax <- function(linear_predictor) {
  row_maximum <- apply(linear_predictor, 1, max)
  exponentiated <- exp(linear_predictor - row_maximum)
  exponentiated / rowSums(exponentiated)
}

predict_from_coefficient_matrix <- function(model_matrix, coefficient_matrix, category_names) {
  nonreference_names <- rownames(coefficient_matrix)
  reference_name <- setdiff(category_names, nonreference_names)
  stopifnot(length(reference_name) == 1)
  linear_predictor <- cbind(0, model_matrix %*% t(coefficient_matrix))
  colnames(linear_predictor) <- c(reference_name, nonreference_names)
  probability <- softmax(linear_predictor)
  probability[, category_names, drop = FALSE]
}

average_probability_by_season <- function(probability, prediction_grid) {
  probability |>
    as.data.frame() |>
    mutate(season = prediction_grid$season, .before = 1) |>
    group_by(season) |>
    summarise(across(everything(), mean), .groups = "drop")
}

community_center_index <- function(focal_probability) {
  centered_log_ratio <- log(focal_probability) - rowMeans(log(focal_probability))
  exp(sweep(centered_log_ratio, 2, colMeans(centered_log_ratio), FUN = "-"))
}

# Standardize every season over the same well-supported dates ------------

prediction_grid <- expand_grid(
  season = seasons,
  day_of_season = reference_days
) |>
  mutate(
    season_centered = season - median(model_data$season),
    moon_distance_from_new_moon_z = 0,
    mist_probability_light_patchy = mean(daily_predictors$mist_probability_light_patchy),
    mist_probability_good = mean(daily_predictors$mist_probability_good),
    era5_rain_log_z = 0,
    wind_speed_10m_mean_ms_z = 0,
    temperature_2m_mean_c_z = 0,
    surface_pressure_mean_hpa_z = 0
  )

prediction_probability <- predict(primary_model, newdata = prediction_grid, type = "probs")
prediction_probability <- prediction_probability[, all_category_names, drop = FALSE]
annual_probability <- average_probability_by_season(prediction_probability, prediction_grid)
focal_probability <- as.matrix(annual_probability[, focal_names])
point_index <- community_center_index(focal_probability)

# Propagate coefficient, normalization and overdispersion uncertainty ----

prediction_model_matrix <- model.matrix(
  delete.response(terms(primary_model)),
  data = prediction_grid,
  contrasts.arg = primary_model$contrasts,
  xlev = primary_model$xlevels
)

coefficient_matrix <- coef(primary_model)
coefficient_mean <- as.vector(t(coefficient_matrix))
adjusted_covariance <- vcov(primary_model) * diagnostics$uncertainty_dispersion
n_draws <- 500L

set.seed(1905)
coefficient_draws <- MASS::mvrnorm(
  n_draws,
  mu = coefficient_mean,
  Sigma = adjusted_covariance
)

index_draws <- array(
  NA_real_,
  dim = c(n_draws, length(seasons), length(focal_names)),
  dimnames = list(NULL, seasons, focal_names)
)
endpoint_probability_draws <- array(
  NA_real_,
  dim = c(n_draws, 2, length(focal_names)),
  dimnames = list(NULL, c("start", "end"), focal_names)
)

for (draw in seq_len(n_draws)) {
  draw_coefficient_matrix <- matrix(
    coefficient_draws[draw, ],
    nrow = nrow(coefficient_matrix),
    byrow = TRUE,
    dimnames = dimnames(coefficient_matrix)
  )
  draw_probability <- predict_from_coefficient_matrix(
    prediction_model_matrix,
    draw_coefficient_matrix,
    all_category_names
  )
  draw_annual_probability <- average_probability_by_season(draw_probability, prediction_grid)
  draw_focal_probability <- as.matrix(draw_annual_probability[, focal_names])
  index_draws[draw, , ] <- community_center_index(draw_focal_probability)
  endpoint_probability_draws[draw, 1, ] <- draw_focal_probability[trend_start_index, ]
  endpoint_probability_draws[draw, 2, ] <- draw_focal_probability[trend_end_index, ]
}

index_lower <- apply(index_draws, c(2, 3), quantile, probs = 0.025)
index_upper <- apply(index_draws, c(2, 3), quantile, probs = 0.975)

matrix_to_long <- function(matrix, value_name) {
  matrix |>
    as.data.frame() |>
    mutate(season = seasons, .before = 1) |>
    pivot_longer(-season, names_to = "common_name", values_to = value_name)
}

annual_index <- matrix_to_long(point_index, "relative_index") |>
  left_join(matrix_to_long(index_lower, "lower"), by = c("season", "common_name")) |>
  left_join(matrix_to_long(index_upper, "upper"), by = c("season", "common_name")) |>
  left_join(
    focal_species |> select(avibase_id, common_name, total_count, abundance_rank),
    by = "common_name"
  ) |>
  select(avibase_id, common_name, total_count, abundance_rank, season, relative_index, lower, upper) |>
  arrange(abundance_rank, season)

# Summarize long-term and pairwise relative change -----------------------

trend_summary <- lapply(seq_along(focal_names), function(species_index) {
  species_name <- focal_names[[species_index]]
  point_ratio <- point_index[trend_end_index, species_index] / point_index[trend_start_index, species_index]
  draw_ratio <- index_draws[, trend_end_index, species_index] / index_draws[, trend_start_index, species_index]
  tibble(
    common_name = species_name,
    first_season = trend_start,
    last_season = trend_end,
    relative_change_ratio = point_ratio,
    lower = quantile(draw_ratio, 0.025),
    upper = quantile(draw_ratio, 0.975),
    relative_change_percent = 100 * (point_ratio - 1),
    direction = case_when(
      quantile(draw_ratio, 0.025) > 1 ~ "increase relative to community",
      quantile(draw_ratio, 0.975) < 1 ~ "decrease relative to community",
      TRUE ~ "uncertain relative change"
    )
  )
}) |>
  bind_rows() |>
  left_join(
    focal_species |> select(avibase_id, common_name, total_count, abundance_rank),
    by = "common_name"
  ) |>
  select(avibase_id, common_name, total_count, abundance_rank, everything()) |>
  arrange(abundance_rank)

# Re-center the same fitted endpoint changes on abundant-species subsets.
reference_sets <- tibble(
  reference = c("All 19 focal species", "Top 10 by catch", "Top 5 by catch"),
  n_species = c(19L, 10L, 5L)
)

reference_sensitivity <- bind_rows(lapply(seq_len(nrow(reference_sets)), function(i) {
  reference_shift <- exp(mean(log(trend_summary$relative_change_ratio[
    trend_summary$abundance_rank <= reference_sets$n_species[[i]]
  ])))
  trend_summary |>
    transmute(
      common_name, abundance_rank,
      reference = reference_sets$reference[[i]],
      reference_species = reference_sets$n_species[[i]],
      relative_change_ratio = relative_change_ratio / reference_shift
    )
}))

species_pairs <- combn(focal_names, 2, simplify = FALSE)
pairwise_trend_contrasts <- lapply(species_pairs, function(pair) {
  first_index <- match(pair[[1]], focal_names)
  second_index <- match(pair[[2]], focal_names)
  point_log_change <-
    log(focal_probability[trend_end_index, first_index] / focal_probability[trend_end_index, second_index]) -
    log(focal_probability[trend_start_index, first_index] / focal_probability[trend_start_index, second_index])
  draw_log_change <-
    log(endpoint_probability_draws[, 2, first_index] / endpoint_probability_draws[, 2, second_index]) -
    log(endpoint_probability_draws[, 1, first_index] / endpoint_probability_draws[, 1, second_index])
  tibble(
    species_1 = pair[[1]],
    species_2 = pair[[2]],
    relative_change_ratio = exp(point_log_change),
    lower = exp(quantile(draw_log_change, 0.025)),
    upper = exp(quantile(draw_log_change, 0.975))
  )
}) |>
  bind_rows() |>
  left_join(
    focal_species |> select(species_1 = common_name, species_1_rank = abundance_rank),
    by = "species_1"
  ) |>
  left_join(
    focal_species |> select(species_2 = common_name, species_2_rank = abundance_rank),
    by = "species_2"
  ) |>
  arrange(species_1_rank, species_2_rank)

marsh_thrush_contrast <- pairwise_trend_contrasts |>
  filter(species_1 == "Marsh Warbler", species_2 == "Thrush Nightingale")

# Build descriptive and phenology summaries -----------------------------

raw_annual_composition <- model_data |>
  group_by(season, common_name) |>
  summarise(annual_count = sum(count), .groups = "drop") |>
  group_by(season) |>
  mutate(
    annual_comparable_count = sum(annual_count),
    raw_annual_proportion = annual_count / annual_comparable_count
  ) |>
  ungroup() |>
  semi_join(focal_species, by = "common_name") |>
  left_join(
    focal_species |> select(avibase_id, common_name, total_count, abundance_rank),
    by = "common_name"
  ) |>
  select(
    avibase_id, common_name, total_count, abundance_rank, season,
    annual_count, annual_comparable_count, raw_annual_proportion
  ) |>
  arrange(abundance_rank, season)

phenology_days <- reference_days

phenology_grid <- tibble(
  season = median(seasons),
  season_centered = 0,
  day_of_season = phenology_days,
  moon_distance_from_new_moon_z = 0,
  mist_probability_light_patchy = mean(daily_predictors$mist_probability_light_patchy),
  mist_probability_good = mean(daily_predictors$mist_probability_good),
  era5_rain_log_z = 0,
  wind_speed_10m_mean_ms_z = 0,
  temperature_2m_mean_c_z = 0,
  surface_pressure_mean_hpa_z = 0
)

phenology_probability <- predict(primary_model, newdata = phenology_grid, type = "probs")
phenology_predictions <- phenology_probability[, focal_names, drop = FALSE] |>
  as.data.frame() |>
  mutate(day_of_season = phenology_days, .before = 1) |>
  pivot_longer(-day_of_season, names_to = "common_name", values_to = "predicted_proportion") |>
  group_by(common_name) |>
  mutate(relative_phenology = predicted_proportion / max(predicted_proportion)) |>
  ungroup() |>
  left_join(
    focal_species |> select(common_name, total_count, abundance_rank),
    by = "common_name"
  ) |>
  arrange(abundance_rank, day_of_season)

# Figures ----------------------------------------------------------------

plot_theme <- ngulia_theme()

sampling_window_plot <- ggplot(coverage_by_season, aes(season)) +
  geom_ribbon(
    aes(ymin = first_day_of_season, ymax = last_day_of_season),
    fill = "#9ecae1",
    alpha = 0.45
  ) +
  geom_line(aes(y = median_day_of_season), colour = ngulia_colours[["blue"]], linewidth = 0.8) +
  geom_hline(yintercept = range(reference_days), linetype = "dashed", colour = "#238b45") +
  labs(
    title = "The observed ringing window changes among seasons",
    subtitle = "Ribbon: first to last positive-catch date; line: median date; dashed: common prediction window",
    x = "Season",
    y = "Days since 1 October"
  ) +
  plot_theme

model_progression_plot <- tempered_mean_structure_comparison |>
  mutate(model = factor(model_id, levels = model_id)) |>
  ggplot(aes(model, mean_daily_log_loss, group = 1)) +
  geom_line(colour = "#65777D", linewidth = 0.7) +
  geom_point(aes(colour = model_id == primary_bundle$specification$model_id), size = 3) +
  geom_text(aes(label = sprintf("%.3f", mean_daily_log_loss)), vjust = -0.8, size = 3.4) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = ngulia_colours[["teal"]]), guide = "none") +
  scale_x_discrete(labels = c("M0\nconstant", "M1\nyear", "M2\n+ timing", "M3\n+ conditions")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    title = "Phenology and daily conditions improve unseen-season prediction",
    subtitle = "Five-fold validation holds out complete seasons; lower log loss is better",
    x = NULL,
    y = "Held-out mean-daily log loss"
  ) +
  plot_theme

phenology_plot <- ggplot(
  phenology_predictions |>
    mutate(common_name = factor(common_name, levels = focal_names)),
  aes(day_of_season, relative_phenology, colour = common_name)
) +
  geom_line(linewidth = 0.75, show.legend = FALSE) +
  facet_wrap(vars(common_name), ncol = 4) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "Broad seasonal timing differs among species",
    subtitle = paste0(
      "Common day ", min(reference_days), "–", max(reference_days),
      " window; panels rank total retained catch from largest to smallest"
    ),
    x = "Days since 1 October",
    y = "Relative fitted representation"
  ) +
  plot_theme

trend_plot <- annual_index |>
  filter(between(season, trend_start, trend_end)) |>
  mutate(common_name = factor(common_name, levels = focal_names)) |>
  ggplot(aes(season, relative_index)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9ecae1", alpha = 0.38) +
  geom_line(colour = ngulia_colours[["blue"]], linewidth = 0.75) +
  facet_wrap(vars(common_name), scales = "free_y", ncol = 4) +
  scale_y_log10() +
  labs(
    title = "Broad species-specific change relative to the equal-weight community center",
    subtitle = paste0(
      trend_start, "–", trend_end,
      "; panels rank total retained catch from largest to smallest"
    ),
    x = "Season",
    y = "Community-centered relative index (log scale)"
  ) +
  plot_theme

reference_sensitivity_plot <- reference_sensitivity |>
  mutate(
    common_name = factor(common_name, levels = rev(focal_names)),
    reference = factor(reference, levels = reference_sets$reference)
  ) |>
  ggplot(aes(relative_change_ratio, common_name, colour = reference)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(aes(group = common_name), colour = "grey75", linewidth = 0.5) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  scale_colour_manual(values = c("#2166ac", "#b2182b", "#1b7837")) +
  labs(
    title = paste0("Endpoint changes depend on the reference species, ", trend_start, "–", trend_end),
    subtitle = "Same fitted model for all three series; only the equal-weight reference changes",
    x = "Fitted change relative to reference (fold, log scale)",
    y = NULL,
    colour = "Reference"
  ) +
  plot_theme

pairwise_plot_data <- bind_rows(
  pairwise_trend_contrasts |>
    transmute(species_1, species_2, log2_change = log2(relative_change_ratio)),
  pairwise_trend_contrasts |>
    transmute(species_1 = species_2, species_2 = species_1, log2_change = -log2(relative_change_ratio)),
  tibble(species_1 = focal_names, species_2 = focal_names, log2_change = 0)
) |>
  mutate(
    species_1 = factor(species_1, levels = rev(focal_names)),
    species_2 = factor(species_2, levels = focal_names)
  )

pairwise_limit <- max(abs(pairwise_plot_data$log2_change), na.rm = TRUE)
pairwise_heatmap <- ggplot(pairwise_plot_data, aes(species_2, species_1, fill = log2_change)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-pairwise_limit, pairwise_limit),
    name = "log2 relative\nchange"
  ) +
  labs(
    title = paste0("Pairwise change from ", trend_start, " to ", trend_end),
    subtitle = "Each row is the change of that species relative to the column species",
    x = "Comparison species",
    y = "Focal species"
  ) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

# Report tables ----------------------------------------------------------

html_table <- function(data) {
  tags$table(
    tags$thead(tags$tr(lapply(names(data), tags$th))),
    tags$tbody(lapply(seq_len(nrow(data)), function(i) {
      tags$tr(lapply(data, function(column) tags$td(as.character(column[[i]]))))
    }))
  )
}

report_model_table <- tempered_mean_structure_comparison |>
  transmute(
    Model = model,
    `Added information` = added_information,
    `Mean-daily log loss` = sprintf("%.3f", mean_daily_log_loss),
    `Bird-weighted log loss` = sprintf("%.3f", bird_weighted_log_loss),
    `Selected mean structure` = if_else(model_id == primary_bundle$specification$model_id, "yes", "")
  )

report_weighting_table <- high_catch_model_comparison |>
  transmute(
    Approach = approach,
    Role = role,
    `Bird-weighted log loss` = sprintf("%.3f", bird_weighted_log_loss),
    `Mean daily log loss` = sprintf("%.3f", mean_daily_log_loss)
  )

timing_gain <- tempered_mean_structure_comparison$mean_daily_log_loss[
  tempered_mean_structure_comparison$model_id == "M1"
] - tempered_mean_structure_comparison$mean_daily_log_loss[
  tempered_mean_structure_comparison$model_id == "M2"
]
condition_gain <- tempered_mean_structure_comparison$mean_daily_log_loss[
  tempered_mean_structure_comparison$model_id == "M2"
] - tempered_mean_structure_comparison$mean_daily_log_loss[
  tempered_mean_structure_comparison$model_id == "M3"
]

top_five_endpoint_sensitivity <- high_catch_endpoint_comparison |>
  filter(model_id %in% c("multinomial", "exclude_top_5pct")) |>
  select(model_id, common_name, relative_change_1977_2022) |>
  pivot_wider(names_from = model_id, values_from = relative_change_1977_2022) |>
  mutate(fold_difference = pmax(
    exclude_top_5pct / multinomial,
    multinomial / exclude_top_5pct
  )) |>
  arrange(desc(fold_difference))

excluded_period_table <- excluded_early_seasons |>
  transmute(
    Season = season,
    `Positive-catch dates` = n_positive_dates,
    `Comparable birds` = comma(total_comparable_count),
    Reason = exclusion_reason
  )

# Build one explanatory HTML report --------------------------------------

model_report <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$title("Ngulia species-composition analysis"),
    tags$style(HTML(paste(readLines(file.path(project_dir, "assets", "report.css")), collapse = "\n")))
  ),
    tags$body(
    tags$h1("Which species changed relative to the Ngulia migrant community?"),
    tags$p(class = "lede", paste0(
      "A joint composition model of ", diagnostics$n_focal_species,
      " comparably caught species across ", diagnostics$n_seasons,
      " post-transition ringing seasons."
    )),
    tags$div(
      class = "note",
      tags$b("Interpretation: "),
      "These are trends in relative representation within the Ngulia catch, not absolute population trends. A shared increase or decrease across every species is not identifiable from composition alone."
    ),

    tags$h2("1. Analysis period and comparable catches"),
    tags$p(
      "The primary analysis starts in 1977. Historical operations evidence shows that catches through 1975 came mainly from six to ten dawn nets south of the lodge. Night nets below the floodlights were introduced in 1976, and the contemporary 1976/77 report describes increased night effort while some dates were still caught almost entirely at dawn. The 1976 season is therefore treated as the transition, and 1969–1976 is excluded before species selection or model fitting."
    ),
    html_table(excluded_period_table),
    tags$p(
      "Barn Swallows, Bank Swallows, Western House-Martins, Red-rumped Swallows and the recorded hybrid are excluded because targeted daytime catching is a different observation process. Species enter the model with at least 500 birds, 10 seasons and 100 positive dates. This support threshold controls precision; it is not a percentage-of-catch cutoff."
    ),
    tags$p(class = "result", paste0(
      comma(diagnostics$n_birds), " birds on ", comma(diagnostics$n_dates),
      " positive-catch dates remain in 19 focal species plus one pooled other-taxa category."
    )),

    tags$h2("2. Day of season must be standardized"),
    tags$p(paste0(
      "The dates sampled differ markedly among seasons. All annual predictions are therefore averaged over the same day ",
      min(reference_days), "–", max(reference_days), " window, where every day is represented in at least ",
      min(day_support$n_seasons[day_support$in_reference_window]),
      " seasons (40% of the retained record)."
    )),
    tags$img(class = "figure wide", src = "figures/00_sampling_window.png", alt = "Ringing date coverage by season"),

    tags$h2("3. Joint model structure"),
    tags$div(class = "equation", HTML("N<sup>eff</sup><sub>d</sub> = N<sub>d</sub>(&phi; + 1)/(N<sub>d</sub> + &phi;)")),
    tags$div(class = "equation", HTML("Y<sup>*</sup><sub>d</sub> &sim; working Multinomial(N<sup>eff</sup><sub>d</sub>, p<sub>1d</sub>, &hellip;, p<sub>Kd</sub>)")),
    tags$div(class = "equation", HTML("log(p<sub>sd</sub>/p<sub>rd</sub>) = &alpha;<sub>s</sub> + f<sub>s</sub>(season) + g<sub>s</sub>(day) + &beta;<sub>s</sub>X<sub>d</sub>")),
    tags$p(
      paste0(
        "The model conditions on the total comparable catch N on each date and estimates how that total is divided among species. Counts are scaled without changing their proportions so that daily information matches the variance of a fitted Dirichlet-multinomial model. Its concentration φ is ",
        sprintf("%.1f", diagnostics$dm_concentration),
        ": the median catch of ", comma(diagnostics$median_daily_count),
        " birds contributes ", sprintf("%.1f", diagnostics$median_effective_daily_count),
        " effective observations, while the largest catch of ", comma(diagnostics$maximum_daily_count),
        " contributes ", sprintf("%.1f", diagnostics$maximum_effective_daily_count),
        ". This is a variance-matched working likelihood rather than a claim that fractional birds were observed. The complete species vector remains joint, so probabilities sum to one."
      )
    ),
    tags$ul(
      tags$li(paste0(
        "Each species has its own broad long-term curve (", primary_bundle$season_spline_df,
        " spline degrees of freedom) and seasonal-timing curve (",
        primary_bundle$phenology_spline_df, " degrees of freedom)."
      )),
      tags$li("Moon, mist, rain, wind, temperature and pressure may affect species differently and are included as daily covariates."),
      tags$li("Annual indices average predictions over the same supported dates and hold daily conditions constant, so changing sampling dates do not create a trend."),
      tags$li("Dates are treated as conditionally independent. Within-date clustering is represented by the reduced effective daily total; residual Pearson dispersion is retained for uncertainty but is never used to shrink intervals."),
      tags$li("Five-fold validation holds out complete seasons and re-estimates the tempering concentration within each training fold; M2 is the minimum candidate because adjustment for species phenology is required by design.")
    ),
    tags$img(class = "figure", src = "figures/01_model_progression.png", alt = "Held-out model comparison"),
    html_table(report_model_table),
    tags$p(class = "result", paste0(
      "Under the final tempered weighting, adding phenology reduces held-out mean-daily log loss by ", sprintf("%.3f", timing_gain),
      "; daily conditions reduce it by another ", sprintf("%.3f", condition_gain),
      ". M3 therefore remains the selected mean structure after calibrating daily information."
    )),

    tags$h2("4. Are exceptional high-catch dates driving the trends?"),
    tags$p(paste0(
      "The largest 1% of dates contain ",
      percent(catch_concentration$share_of_birds[catch_concentration$threshold == "top_1_percent_dates"], accuracy = 0.1),
      " of birds, and the largest 5% contain ",
      percent(catch_concentration$share_of_birds[catch_concentration$threshold == "top_5_percent_dates"], accuracy = 0.1),
      ". Removing the largest date one at a time changes any fitted index by at most ",
      percent(max(top_date_influence$maximum_fold_index_difference - 1), accuracy = 0.1),
      ". Removing the largest 5% of dates changes the median species endpoint ratio by ",
      percent(
        high_catch_trend_sensitivity$median_endpoint_fold_difference[
          high_catch_trend_sensitivity$model_id == "exclude_top_5pct"
        ] - 1,
        accuracy = 0.1
      ),
      ". The largest species-specific endpoint change is ",
      percent(top_five_endpoint_sensitivity$fold_difference[[1]] - 1, accuracy = 0.1),
      " for ", top_five_endpoint_sensitivity$common_name[[1]],
      ". Exceptional dates therefore do not determine the general curves, although individual-species sensitivity remains visible."
    )),
    html_table(report_weighting_table),
    tags$p(class = "result", paste0(
      "Daily vectors are nevertheless strongly overdispersed: the ordinary multinomial count-vector log score is ",
      sprintf("%.1f", count_distribution_scores$mean_daily_count_log_score[count_distribution_scores$model_id == "multinomial"]),
      " versus ",
      sprintf("%.1f", count_distribution_scores$mean_daily_count_log_score[count_distribution_scores$model_id == "dirichlet_multinomial"]),
      " for the Dirichlet-multinomial (lower is better). Variance-tempered weighting has the best held-out mean-daily composition score and nearly unchanged bird-weighted performance, so it is the primary estimator. The full Dirichlet-multinomial remains a sensitivity because it represents dispersion coherently but predicts held-out mean composition less accurately."
    )),

    tags$h2("5. Species have different seasonal timing"),
    tags$p(
      paste0(
        "The phenology curves explain why raw annual proportions cannot be compared safely when the ringing window changes. ",
        primary_bundle$phenology_spline_df,
        " spline degrees of freedom retain broad passage timing while suppressing short wiggles that have no strong biological interpretation. Species are ordered by total catch in the retained 1977–2023 data."
      )
    ),
    tags$img(class = "figure wide", src = "figures/02_species_phenology.png", alt = "Species-specific within-season timing"),

    tags$h2("6. Community-centered relative trends"),
    tags$p(
      paste0(
        "For each season, fitted species probabilities are converted to log ratios against the geometric mean of the 19 focal species. Every species therefore contributes equally to the reference center: Marsh Warbler is one of 19 contributors, not 41% of the reference. Each species curve is then scaled to a long-term geometric mean of one. ",
        primary_bundle$season_spline_df,
        " spline degrees of freedom emphasize the general multi-decadal pattern rather than short fluctuations."
      )
    ),
    tags$img(class = "figure wide", src = "figures/03_species_relative_trends.png", alt = "Community-centered relative species trends"),
    tags$p(paste0(
      "The equal-weight center is sensitive to which species enter it. Re-centering the same fitted endpoint predictions on the ten most-caught species changes Marsh Warbler from ",
      sprintf("%.2f", reference_sensitivity$relative_change_ratio[reference_sensitivity$common_name == "Marsh Warbler" & reference_sensitivity$reference == "All 19 focal species"]),
      " to ",
      sprintf("%.2f", reference_sensitivity$relative_change_ratio[reference_sensitivity$common_name == "Marsh Warbler" & reference_sensitivity$reference == "Top 10 by catch"]),
      " times its reference. This is a sensitivity of the display reference, not a refit with a different focal-species list. Pairwise species contrasts are unchanged by this re-centering."
    )),
    tags$img(class = "figure wide", src = "figures/03_reference_species_sensitivity.png", alt = "Species endpoint changes under three reference species sets"),

    tags$h2("7. Direct pairwise comparisons"),
    tags$p(
      paste0(
        "The pairwise table and heatmap remove the community center entirely. For example, the Marsh Warbler–Thrush Nightingale contrast answers whether their fitted relative representation changed against one another. A value of two means that the row species' odds doubled relative to the column species between ",
        trend_start, " and ", trend_end, ". This common supported period avoids endpoint extrapolation for species absent from the first or last project seasons."
      )
    ),
    tags$p(class = "result", paste0(
      "From ", trend_start, " to ", trend_end,
      ", the Marsh Warbler-to-Thrush Nightingale relative odds changed to ",
      sprintf("%.2fx", marsh_thrush_contrast$relative_change_ratio),
      " (95% interval ", sprintf("%.2f", marsh_thrush_contrast$lower),
      "–", sprintf("%.2f", marsh_thrush_contrast$upper),
      "). This comparison is independent of how either species relates to the other 17 focal species."
    )),
    tags$img(class = "figure wide", src = "figures/04_pairwise_change_heatmap.png", alt = "Pairwise long-term relative species changes"),

    tags$h2("8. Diagnostics and limits"),
    tags$p(paste0(
      "The model converged with probability-sum error below ",
      format(diagnostics$max_probability_sum_error, scientific = TRUE, digits = 2),
      ". The untempered multinomial Pearson dispersion was ",
      sprintf("%.2f", untempered_diagnostics$pearson_dispersion),
      "; after variance-tempering it is ", sprintf("%.2f", diagnostics$pearson_dispersion),
      ". The effective-count covariance is not shrunk when dispersion is below one, and dates retain the requested conditional-independence assumption."
    )),
    tags$p(
      "Intervals condition on the fitted Dirichlet-multinomial concentration rather than propagating uncertainty in that parameter. Its estimates across held-out training folds are stable, but this remains a small unrepresented source of uncertainty."
    ),
    tags$p(
      "The model resolves changes among species, but composition alone cannot distinguish a focal increase from a collective decline of its comparison species. Excluding 1969–1976 removes the largest early change in capture regime, but it does not make later operations constant. Weaker lamps in 1984 and the 1994–95 move to more productive northern dawn nets may still alter species composition if taxa differ in night-versus-dawn catchability. Daily weather adjustment cannot remove that protocol bias; pairwise outputs only make the biological comparison explicit."
    ),
    tags$p(
      "Tables: ",
      tags$a(href = "tables/species_community_centered_index.csv", "annual index"), " · ",
      tags$a(href = "tables/species_relative_trend_summary.csv", "trend summary"), " · ",
      tags$a(href = "tables/species_reference_sensitivity.csv", "reference sensitivity"), " · ",
      tags$a(href = "tables/species_pairwise_trend_contrasts.csv", "pairwise contrasts"), " · ",
      tags$a(href = "tables/tempered_mean_structure_comparison.csv", "mean-structure comparison"), " · ",
      tags$a(href = "tables/species_model_support.csv", "species support"), " · ",
      tags$a(href = "tables/analysis_period_decision.csv", "period decision"), " · ",
      tags$a(href = "tables/high_catch_model_comparison.csv", "high-catch model comparison"), " · ",
      tags$a(href = "tables/high_catch_trend_sensitivity.csv", "trend sensitivity")
    )
  )
)

# Write outputs ----------------------------------------------------------

story_figures <- list(
  "00_sampling_window.png" = list(sampling_window_plot, 10, 5.5),
  "01_model_progression.png" = list(model_progression_plot, 8, 5),
  "02_species_phenology.png" = list(phenology_plot, 12, 14),
  "03_species_relative_trends.png" = list(trend_plot, 12, 14),
  "04_pairwise_change_heatmap.png" = list(pairwise_heatmap, 12, 10),
  "03_reference_species_sensitivity.png" = list(reference_sensitivity_plot, 11, 8)
)

for (file_name in names(story_figures)) {
  figure <- story_figures[[file_name]]
  ngulia_save(
    file.path(figure_dir, file_name),
    figure[[1]],
    width = figure[[2]],
    height = figure[[3]],
    dpi = 220
  )
}

write_csv(raw_annual_composition, file.path(table_dir, "raw_annual_species_composition.csv"))
write_csv(phenology_predictions, file.path(table_dir, "species_phenology_predictions.csv"))
write_csv(annual_index, file.path(table_dir, "species_community_centered_index.csv"))
write_csv(trend_summary, file.path(table_dir, "species_relative_trend_summary.csv"))
write_csv(pairwise_trend_contrasts, file.path(table_dir, "species_pairwise_trend_contrasts.csv"))
write_csv(reference_sensitivity, file.path(table_dir, "species_reference_sensitivity.csv"))

rendered_report <- renderTags(model_report)
report_html <- sub(
  "<html>",
  paste0("<!doctype html>\n<html>\n<head>\n", rendered_report$head, "\n</head>"),
  rendered_report$html,
  fixed = TRUE
)
writeLines(report_html, file.path(analysis_dir, "species_composition_analysis.html"))

cli_alert_success(
  "Wrote species-composition report to {file.path(analysis_dir, 'species_composition_analysis.html')}"
)
