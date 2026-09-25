library(dplyr)
library(readr)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

daily_coverage_path <- file.path(paths$curated_dir, "daily_coverage.csv")
model_data_dir <- file.path(paths$analysis_output_dir, "01_total_catch", "model_data")

dir.create(model_data_dir, recursive = TRUE, showWarnings = FALSE)

# Read and prepare data ---------------------------------------------------

cli_h1("Prepare total-catch model data")

daily_coverage <- read_csv(
  daily_coverage_path,
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(
    era5_rain_log = log1p(total_precipitation_00_08_mm),
    documented_active_net_site = effort_status == "documented_operation" & !is.na(total_birds_ringed),
    eligible_positive_catch_model = ringing_happened & season >= 1977L & !season %in% 1994:1995,
    eligible_documented_operations_model = ringing_happened & documented_active_net_site &
      bush_net_configuration %in% c("back_bush", "front_bush"),
    positive_catch_sampling_basis = if_else(documented_active_net_site, "documented_net_operation", "positive_catch_record_only")
  )

positive_model_columns <- c(
  "ringing_date", "season", "season_day", "total_birds_ringed",
  "moon_distance_from_new_moon",
  "mist_observation", "mist_probability_none", "mist_probability_light_patchy", "mist_probability_good", "era5_rain_log",
  "total_cloud_cover_mean", "cloud_base_height_mean_m", "relative_humidity_mean_pct", "wind_u_10m_mean_ms",
  "wind_speed_10m_mean_ms", "temperature_2m_mean_c",
  "surface_pressure_mean_hpa"
)

operations_model_columns <- c(
  "ringing_date", "season", "season_day", "total_birds_ringed",
  "moon_distance_from_new_moon", "mist_observation", "mist_probability_none", "mist_probability_light_patchy", "mist_probability_good", "era5_rain_log",
  "total_cloud_cover_mean", "cloud_base_height_mean_m", "relative_humidity_mean_pct", "wind_u_10m_mean_ms",
  "wind_speed_10m_mean_ms", "temperature_2m_mean_c", "surface_pressure_mean_hpa",
  "djp_team_size_minimum", "djp_team_size_interpretation",
  "daily_count_status", "net_sites_observed", "djp_tape", "playback_nocturnal_observed",
  "bush_net_configuration", "night_net_configuration", "rain_observed",
  "night_net_operation", "dawn_net_operation"
)

positive_catch_model_data <- daily_coverage |>
  filter(
    eligible_positive_catch_model,
    if_all(
      c(
        season_day,
        moon_distance_from_new_moon,
        mist_probability_none,
        mist_probability_light_patchy,
        mist_probability_good,
        era5_rain_log,
        wind_speed_10m_mean_ms,
        temperature_2m_mean_c,
        surface_pressure_mean_hpa
      ),
      ~ !is.na(.x)
    )
  ) |>
  select(all_of(positive_model_columns))

documented_operation_model_data <- daily_coverage |>
  filter(
    eligible_documented_operations_model,
    if_all(
      c(season_day, moon_distance_from_new_moon, mist_probability_none, mist_probability_light_patchy, mist_probability_good,
        era5_rain_log, wind_speed_10m_mean_ms, temperature_2m_mean_c, surface_pressure_mean_hpa,
        net_sites_observed, playback_nocturnal_observed),
      ~ !is.na(.x)
    )
  ) |>
  select(all_of(operations_model_columns))

model_coverage <- daily_coverage |>
  group_by(season) |>
  summarise(
    n_calendar_dates = n(),
    n_positive_catch_dates = sum(ringing_happened),
    n_documented_net_operation_dates = sum(documented_active_net_site),
    n_documented_operation_zero_catch_dates = sum(documented_active_net_site & !ringing_happened),
    n_positive_catch_model_dates = sum(ringing_date %in% positive_catch_model_data$ringing_date),
    n_documented_operation_model_dates = sum(ringing_date %in% documented_operation_model_data$ringing_date),
    .groups = "drop"
  )

# Write outputs -----------------------------------------------------------

write_csv(positive_catch_model_data, file.path(model_data_dir, "positive_catch_model_data.csv"))
write_csv(documented_operation_model_data, file.path(model_data_dir, "documented_operation_model_data.csv"))
write_csv(model_coverage, file.path(model_data_dir, "model_coverage_by_season.csv"))

cli_alert_success("Wrote total-catch model inputs to {model_data_dir}")
