library(dplyr)
library(tidyr)
library(readr)
library(cli)

# Set paths ---------------------------------------------------------------

source(here::here("scripts/helpers/data_paths.R"))
paths <- get_data_paths()
analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
table_dir <- file.path(analysis_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read canonical daily data ----------------------------------------------

daily <- read_csv(
  file.path(paths$curated_dir, "daily_coverage.csv"),
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
)

# Document model-facing covariates ---------------------------------------

registry <- tribble(
  ~covariate, ~label, ~type, ~model_use, ~interpretation,
  "season", "Season", "numeric", "Annual effect and trend", "Captures residual differences among seasons after daily adjustment.",
  "season_day", "Date within season", "numeric", "M1-M6", "Flexible seasonal timing effect.",
  "mist_state", "Unified mist state", "categorical", "M2-M6", "Observed no/light/good mist where known; otherwise drawn from the ERA5-calibrated state probabilities.",
  "total_precipitation_00_08_mm", "Rainfall", "numeric", "M2-M6", "ERA5 rainfall from 00:00 to 08:00 local time; modeled as log(1 + mm).",
  "moon_distance_from_new_moon", "Distance from new moon", "numeric", "M3-M6", "Days from new moon, the selected lunar representation.",
  "wind_speed_10m_mean_ms", "Wind speed", "numeric", "M4-M6", "ERA5 mean from 00:00 to 08:00 local time.",
  "temperature_2m_mean_c", "Temperature", "numeric", "M4-M6", "ERA5 mean from 00:00 to 08:00 local time.",
  "surface_pressure_mean_hpa", "Surface pressure", "numeric", "M4-M6", "ERA5 mean from 00:00 to 08:00 local time.",
  "djp_team_size_minimum", "Minimum recorded team size", "numeric", "Exploratory only", "Exact count or a lower bound where Earthwatch participant numbers are unspecified; not used to adjust the trend.",
  "playback_nocturnal_observed", "Nocturnal playback", "categorical", "M6 sensitivity", "Missing is unknown; blank cells inside a documented source block mean no playback.",
  "net_sites_observed", "Daily operated net sites", "categorical", "Descriptive only", "Daily opening reflects weather and operations, not a fixed configuration effect.",
  "bush_net_configuration", "Bush-net configuration", "categorical", "M5 fixed period", "Back bush in 1977–1993, transition in 1994–1995, front bush from 1996; identified under a smooth-year assumption.",
  "night_net_configuration", "Night-net configuration", "categorical", "Period restriction", "Established from 1977; pre-1977 omitted from count models."
)

coverage <- bind_rows(lapply(seq_len(nrow(registry)), function(i) {
  variable <- registry$covariate[[i]]
  if (variable == "mist_state") {
    available <- complete.cases(
      daily$mist_probability_none,
      daily$mist_probability_light_patchy,
      daily$mist_probability_good
    )
  } else {
    value <- daily[[variable]]
    available <- !is.na(value) & !as.character(value) %in% c("unknown", "not_recorded", "recorded_unknown")
  }
  tibble(
    covariate = variable,
    n_positive_dates_available = sum(daily$ringing_happened & available),
    n_positive_dates = sum(daily$ringing_happened),
    positive_date_coverage = n_positive_dates_available / n_positive_dates,
    n_seasons_available = n_distinct(daily$season[daily$ringing_happened & available])
  )
})) |>
  left_join(registry, by = "covariate")

# Audit model inputs ------------------------------------------------------

range_checks <- daily |>
  transmute(
    ringing_date,
    invalid_mist_probability = !between(mist_probability_none, 0, 1) |
      !between(mist_probability_light_patchy, 0, 1) |
      !between(mist_probability_good, 0, 1),
    invalid_mist_probability_sum = abs(mist_probability_none + mist_probability_light_patchy + mist_probability_good - 1) > 1e-8,
    negative_rain = total_precipitation_00_08_mm < 0,
    negative_wind_speed = wind_speed_10m_mean_ms < 0,
    negative_team_size = djp_team_size_minimum < 0,
    invalid_count = total_birds_ringed < 0 | total_birds_ringed != round(total_birds_ringed),
    invalid_swallow_subtraction = total_birds_ringed + swallow_birds_ringed != all_birds_ringed
  ) |>
  pivot_longer(-ringing_date, names_to = "check", values_to = "failed")

stopifnot(!any(range_checks$failed, na.rm = TRUE))

# Write compact reference tables -----------------------------------------

write_csv(registry, file.path(table_dir, "daily_count_covariate_dictionary.csv"))
write_csv(coverage, file.path(table_dir, "daily_count_covariate_coverage.csv"))

cli_alert_success("Validated and documented daily-count covariates in {table_dir}")
