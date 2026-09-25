library(dplyr)
library(tidyr)
library(readr)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "04_age_sex_demography")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")

dir.create(model_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_season <- 1991L

# Read data ---------------------------------------------------------------

cli_h1("Prepare age and sex demography model data")

ring_events <- read_csv(
  file.path(paths$curated_dir, "ring_events.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date()),
  col_select = c(
    season, ringing_date, retrap, avibase_id, common_name, species_code,
    age, sex
  )
) |>
  filter(season >= analysis_start_season, retrap %in% FALSE, !is.na(avibase_id), !is.na(common_name)) |>
  mutate(
    age_class = case_when(
      age %in% c(3, 5) ~ "first_year",
      age %in% c(4, 6) ~ "adult",
      TRUE ~ NA_character_
    ),
    sex_certain = if_else(sex %in% c("M", "F"), sex, NA_character_),
    day_of_season = as.integer(ringing_date - as.Date(sprintf("%s-10-01", season)))
  )

daily_coverage <- read_csv(
  file.path(paths$curated_dir, "daily_coverage.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
) |>
  transmute(
    ringing_date,
    moon_distance_from_new_moon,
    mist_probability_light_patchy,
    mist_probability_good,
    era5_rain_log = log1p(total_precipitation_00_08_mm),
    wind_speed_10m_mean_ms,
    temperature_2m_mean_c,
    surface_pressure_mean_hpa
  )

excluded_capture_groups <- read_csv(
  file.path(project_dir, "config", "analysis", "excluded_capture_groups.csv"),
  show_col_types = FALSE
)

comparable_events <- ring_events |>
  anti_join(excluded_capture_groups, by = "avibase_id")

# Audit age and sex recording --------------------------------------------

species_recording_coverage <- comparable_events |>
  group_by(avibase_id, common_name) |>
  summarise(
    n_initial_captures = n(),
    n_ageable = sum(!is.na(age_class)),
    age_recording_fraction = mean(!is.na(age_class)),
    n_adults = sum(age_class == "adult", na.rm = TRUE),
    n_adults_sexed = sum(age_class == "adult" & !is.na(sex_certain), na.rm = TRUE),
    adult_sex_recording_fraction = n_adults_sexed / n_adults,
    n_seasons = n_distinct(season),
    n_dates = n_distinct(ringing_date),
    .groups = "drop"
  )

age_support <- species_recording_coverage |>
  mutate(
    selected_age = n_ageable >= 1000 & n_seasons >= 20 & n_dates >= 100 &
      age_recording_fraction >= 0.80,
    age_selection_rule = paste(
      "ageable >= 1,000; seasons >= 20; dates >= 100;",
      "age recording fraction >= 0.80"
    )
  ) |>
  arrange(desc(n_ageable), common_name)

sex_support <- species_recording_coverage |>
  mutate(
    selected_sex = n_adults_sexed >= 500 & n_seasons >= 15 &
      adult_sex_recording_fraction >= 0.80,
    sex_selection_rule = paste(
      "certainly sexed adults >= 500; seasons >= 15;",
      "adult sex recording fraction >= 0.80"
    )
  ) |>
  arrange(desc(n_adults_sexed), common_name)

focal_age_species <- age_support |>
  filter(selected_age) |>
  arrange(desc(n_ageable)) |>
  mutate(species_order = row_number())

focal_sex_species <- sex_support |>
  filter(selected_sex) |>
  arrange(desc(n_adults_sexed)) |>
  mutate(species_order = row_number())

annual_recording_coverage <- comparable_events |>
  semi_join(focal_age_species, by = "avibase_id") |>
  group_by(season) |>
  summarise(
    n_initial_captures = n(),
    age_recording_fraction = mean(!is.na(age_class)),
    adult_sex_recording_fraction = sum(
      age_class == "adult" & !is.na(sex_certain), na.rm = TRUE
    ) / sum(age_class == "adult", na.rm = TRUE),
    .groups = "drop"
  )

# Aggregate first captures by species and ringing date -------------------

age_daily <- comparable_events |>
  filter(!is.na(age_class)) |>
  semi_join(focal_age_species, by = "avibase_id") |>
  count(
    avibase_id, common_name, species_code, season, ringing_date,
    day_of_season, age_class, name = "n"
  ) |>
  complete(
    nesting(avibase_id, common_name, species_code, season, ringing_date, day_of_season),
    age_class = c("first_year", "adult"),
    fill = list(n = 0)
  ) |>
  pivot_wider(names_from = age_class, values_from = n, names_prefix = "n_") |>
  mutate(
    n_total = n_first_year + n_adult,
    first_year_fraction = n_first_year / n_total,
    effective_n = pmin(n_total, 50),
    season_centered = season - median(ring_events$season)
  ) |>
  left_join(daily_coverage, by = "ringing_date") |>
  arrange(avibase_id, ringing_date)

sex_daily <- comparable_events |>
  filter(age_class == "adult", !is.na(sex_certain)) |>
  semi_join(focal_sex_species, by = "avibase_id") |>
  count(
    avibase_id, common_name, species_code, season, ringing_date,
    day_of_season, sex_certain, name = "n"
  ) |>
  complete(
    nesting(avibase_id, common_name, species_code, season, ringing_date, day_of_season),
    sex_certain = c("M", "F"),
    fill = list(n = 0)
  ) |>
  pivot_wider(names_from = sex_certain, values_from = n, names_prefix = "n_") |>
  transmute(
    avibase_id, common_name, species_code, season, ringing_date, day_of_season,
    n_male = n_M,
    n_female = n_F,
    n_total = n_M + n_F,
    male_fraction = n_M / n_total,
    effective_n = pmin(n_total, 30),
    season_centered = season - median(ring_events$season)
  ) |>
  left_join(daily_coverage, by = "ringing_date") |>
  arrange(avibase_id, ringing_date)

# Scale daily conditions once for fitting and prediction -----------------

condition_variables <- c(
  "moon_distance_from_new_moon", "era5_rain_log", "wind_speed_10m_mean_ms",
  "temperature_2m_mean_c", "surface_pressure_mean_hpa"
)

condition_scaling <- bind_rows(
  age = age_daily,
  sex = sex_daily,
  .id = "analysis"
) |>
  distinct(analysis, ringing_date, across(all_of(condition_variables))) |>
  group_by(analysis) |>
  summarise(across(
    all_of(condition_variables),
    list(center = ~ mean(.x, na.rm = TRUE), scale = ~ sd(.x, na.rm = TRUE)),
    .names = "{.col}_{.fn}"
  ), .groups = "drop") |>
  pivot_longer(
    -analysis,
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(center|scale)$"
  )

scale_conditions <- function(data, analysis_name) {
  scaling <- condition_scaling |> filter(analysis == analysis_name)
  for (variable in condition_variables) {
    center <- scaling$center[scaling$variable == variable]
    scale <- scaling$scale[scaling$variable == variable]
    data[[paste0(variable, "_z")]] <- (data[[variable]] - center) / scale
  }
  data
}

age_daily <- scale_conditions(age_daily, "age")
sex_daily <- scale_conditions(sex_daily, "sex")

# Define well-supported common prediction windows -----------------------

day_support <- bind_rows(
  age = age_daily |> distinct(season, day_of_season),
  sex = sex_daily |> distinct(season, day_of_season),
  .id = "analysis"
) |>
  count(analysis, day_of_season, name = "n_seasons") |>
  group_by(analysis) |>
  mutate(in_reference_window = n_seasons >= ceiling(0.40 * n_distinct(ring_events$season))) |>
  ungroup()

# QA checks and outputs --------------------------------------------------

stopifnot(all(age_daily$n_total > 0), all(sex_daily$n_total > 0))
stopifnot(all(age_daily$n_first_year + age_daily$n_adult == age_daily$n_total))
stopifnot(all(sex_daily$n_male + sex_daily$n_female == sex_daily$n_total))
stopifnot(all(age_daily$first_year_fraction >= 0 & age_daily$first_year_fraction <= 1))
stopifnot(all(sex_daily$male_fraction >= 0 & sex_daily$male_fraction <= 1))

write_csv(age_daily, file.path(model_data_dir, "age_daily_model_data.csv"))
write_csv(sex_daily, file.path(model_data_dir, "adult_sex_daily_model_data.csv"))
write_csv(condition_scaling, file.path(model_data_dir, "condition_scaling.csv"))
write_csv(age_support, file.path(table_dir, "age_species_support.csv"))
write_csv(sex_support, file.path(table_dir, "adult_sex_species_support.csv"))
write_csv(annual_recording_coverage, file.path(table_dir, "annual_recording_coverage.csv"))
write_csv(day_support, file.path(table_dir, "day_of_season_support.csv"))
write_csv(excluded_capture_groups, file.path(table_dir, "excluded_capture_groups.csv"))

cli_alert_success(
  "Prepared age models for {nrow(focal_age_species)} species and adult sex models for {nrow(focal_sex_species)} species."
)
cli_alert_info(
  "Daily information is capped at 50 ageable birds and 30 sexed adults per species-date to limit mass-fall pseudo-replication."
)
