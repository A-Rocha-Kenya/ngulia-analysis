library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
model_path <- file.path(analysis_dir, "models", "adjusted_positive_catch_gam_mi.rds")
figure_path <- file.path(ngulia_figure_dir(file.path(analysis_dir, "figures")), "13_predicted_best_capture_windows.png")

cli_h1("Predict best future capture windows")

# Build future season calendar -------------------------------------------

synodic_month_days <- 29.530588853
reference_new_moon <- ymd_hms("2000-01-06 18:14:00", tz = "UTC")
future_seasons <- 2026:2035

future_dates <- bind_rows(lapply(future_seasons, function(season) {
  tibble(ringing_date = seq(ymd(paste0(season, "-10-20")), ymd(paste0(season + 1, "-01-12")), by = "day")) |>
    mutate(season = season)
})) |>
  mutate(
    season_day = as.integer(ringing_date - make_date(season, 10L, 20L)) + 1L,
    moon_age_days = as.numeric(difftime(as_datetime(ringing_date, tz = "UTC") + hours(12), reference_new_moon, units = "days")) %%
      synodic_month_days,
    moon_days_from_new_moon_exact = if_else(
      moon_age_days <= synodic_month_days / 2,
      moon_age_days,
      moon_age_days - synodic_month_days
    ),
    moon_distance_from_new_moon = abs(as.integer(round(moon_days_from_new_moon_exact))),
    season_label = paste0(season, "–", substr(season + 1, 3, 4))
  )

# Predict opportunity from seasonal timing and moon ----------------------

models <- readRDS(model_path)

daily_predictions <- bind_rows(lapply(seq_along(models), function(imputation) {
  model <- models[[imputation]]
  prediction_data <- model$model[rep(1, nrow(future_dates)), ]
  prediction_data$season_day <- future_dates$season_day
  prediction_data$moon_distance_from_new_moon <- future_dates$moon_distance_from_new_moon

  timing_and_moon <- predict(
    model,
    newdata = prediction_data,
    type = "terms",
    terms = c("s(season_day)", "s(moon_distance_from_new_moon)")
  )

  future_dates |>
    transmute(imputation, ringing_date, relative_opportunity = exp(rowSums(timing_and_moon)))
}))

future_scores <- future_dates |>
  left_join(
    daily_predictions |>
      group_by(ringing_date) |>
      summarise(relative_opportunity = mean(relative_opportunity), .groups = "drop"),
    by = "ringing_date"
  ) |>
  mutate(relative_opportunity = relative_opportunity / mean(relative_opportunity))

# Select the best continuous 14-day window -------------------------------

window_days <- 14L

window_candidates <- future_scores |>
  group_by(season) |>
  arrange(ringing_date, .by_group = TRUE) |>
  mutate(
    window_score = vapply(
      seq_len(n()),
      function(i) {
        if (i + window_days - 1L > n()) return(NA_real_)
        sum(relative_opportunity[i:(i + window_days - 1L)])
      },
      numeric(1)
    ),
    window_days = window_days,
    end_date = ringing_date + window_days - 1L
  ) |>
  filter(!is.na(window_score)) |>
  ungroup()

selected_windows <- window_candidates |>
  group_by(season, window_days) |>
  slice_max(window_score, n = 1, with_ties = FALSE) |>
  ungroup() |>
  rename(start_date = ringing_date) |>
  mutate(
    midpoint_date = start_date + floor((window_days - 1L) / 2),
    mean_daily_relative_opportunity = window_score / window_days
  )

nearest_new_moons <- selected_windows |>
  select(season, window_days, start_date, end_date) |>
  inner_join(
    select(future_scores, season, ringing_date, moon_days_from_new_moon_exact),
    by = "season",
    relationship = "many-to-many"
  ) |>
  filter(ringing_date >= start_date, ringing_date <= end_date) |>
  group_by(season, window_days) |>
  slice_min(abs(moon_days_from_new_moon_exact), n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(season, window_days, nearest_new_moon = ringing_date)

season_average <- future_scores |>
  group_by(season) |>
  summarise(season_average_opportunity = mean(relative_opportunity), .groups = "drop")

best_windows <- selected_windows |>
  left_join(nearest_new_moons, by = c("season", "window_days")) |>
  left_join(season_average, by = "season") |>
  mutate(advantage_over_average_season_date = mean_daily_relative_opportunity / season_average_opportunity - 1) |>
  select(
    season, window_days, start_date, end_date, midpoint_date, nearest_new_moon,
    mean_daily_relative_opportunity, advantage_over_average_season_date
  )

# Plot future planning calendar ------------------------------------------

new_moons <- future_scores |>
  group_by(season) |>
  arrange(ringing_date, .by_group = TRUE) |>
  filter(
    abs(moon_days_from_new_moon_exact) == min(abs(moon_days_from_new_moon_exact)) |
      abs(moon_days_from_new_moon_exact) < lag(abs(moon_days_from_new_moon_exact), default = Inf) &
      abs(moon_days_from_new_moon_exact) < lead(abs(moon_days_from_new_moon_exact), default = Inf)
  ) |>
  ungroup()

plot_dates <- future_scores

plot_windows <- best_windows |>
  mutate(
    start_day = as.integer(start_date - make_date(season, 10L, 20L)) + 1L,
    end_day = as.integer(end_date - make_date(season, 10L, 20L)) + 1L,
    midpoint_day = (start_day + end_day) / 2,
    date_label = paste(format(start_date, "%d %b"), format(end_date, "%d %b"), sep = " – ")
  )

plot_new_moons <- new_moons

date_breaks <- sort(unique(c(seq(1, 85, by = 5), 85)))
axis_dates <- as.Date("2000-10-20") + date_breaks - 1L
date_labels <- if_else(
  seq_along(axis_dates) == 1 | month(axis_dates) != lag(month(axis_dates), default = month(axis_dates[[1]])),
  format(axis_dates, "%d %b"),
  format(axis_dates, "%d")
)

planning_plot <- ggplot(plot_dates, aes(season_day, factor(season))) +
  geom_tile(aes(fill = relative_opportunity), height = 0.8) +
  geom_vline(xintercept = date_breaks, colour = "#D5DFE6", alpha = 0.22, linewidth = 0.3) +
  geom_rect(
    data = plot_windows,
    aes(xmin = start_day - 0.5, xmax = end_day + 0.5, ymin = as.numeric(factor(season)) - 0.42,
        ymax = as.numeric(factor(season)) + 0.42),
    fill = NA, colour = "#58D6E7", linewidth = 1.15, inherit.aes = FALSE
  ) +
  geom_label(
    data = plot_windows,
    aes(midpoint_day, factor(season), label = date_label),
    fill = "#0B1320", colour = "#F5F8FA", alpha = 0.82,
    linewidth = 0, size = 3.1, fontface = "bold", inherit.aes = FALSE,
    position = position_nudge(y = 0.24)
  ) +
  geom_point(
    data = plot_new_moons,
    aes(season_day, factor(season)),
    shape = 23, fill = "#F5F8FA", colour = "#0B1320", stroke = 0.4, size = 2.5,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = as.integer(as.Date("2000-12-01") - as.Date("2000-10-20")) + 1L,
             linetype = "dotted", colour = "#D5DFE6", linewidth = 0.55) +
  scale_fill_gradientn(
    colours = c("#172538", "#24506A", "#35A7B6", "#F1C75B", "#F0714F"),
    values = rescale(c(0.1, 0.6, 1, 1.8, 3.3)),
    name = "Relative daily\nopportunity"
  ) +
  scale_x_continuous(breaks = date_breaks, labels = date_labels, expand = c(0, 0)) +
  labs(
    title = "Best continuous 14-day ringing windows",
    subtitle = "Labels give the recommended start and end dates · vertical guides and ticks are every 5 days · diamonds mark new moon · dotted line marks 1 December",
    x = NULL,
    y = "Ringing season",
    caption = "Scheduling guide from timing and moon effects only. Mist, weather, team size and net operation are not included."
  ) +
  ngulia_theme(dark = TRUE, base_size = 12) +
  theme(
    text = element_text(colour = "#EAF2F7"),
    plot.background = element_rect(fill = "#0B1320", colour = NA),
    panel.background = element_rect(fill = "#0B1320", colour = NA),
    panel.grid = element_blank(),
    axis.text = element_text(colour = "#B9C8D4"),
    axis.title = element_text(colour = "#D9E5EC"),
    legend.position = "right",
    legend.title = element_text(colour = "#D9E5EC"),
    legend.text = element_text(colour = "#B9C8D4"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(colour = "#B9C8D4", size = 11),
    plot.caption = element_text(colour = "#8295A5", hjust = 0, size = 9)
  )

# Write outputs -----------------------------------------------------------

ngulia_save(figure_path, planning_plot, width = 12, height = 6.75, dpi = 220)
unlink(file.path(analysis_dir, "predict_best_dates.html"))
unlink(file.path(analysis_dir, "tables", "predicted_best_capture_windows.csv"))

cli_alert_success("Wrote future-window figure to {figure_path}")
