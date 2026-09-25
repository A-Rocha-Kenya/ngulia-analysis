library(dplyr)
library(tidyr)
library(readr)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "02_species_composition")
model_data_dir <- file.path(analysis_dir, "model_data")
model_dir <- file.path(analysis_dir, "models")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")

analysis_start_season <- 1977L

dir.create(model_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read data ---------------------------------------------------------------

cli_h1("Prepare joint species-composition model data")

daily_counts <- read_csv(
  file.path(paths$curated_dir, "daily_counts.csv"),
  show_col_types = FALSE,
  col_types = cols(
    ringing_date = col_date(),
    season = col_integer(),
    avibase_id = col_character(),
    common_name = col_character(),
    n_records = col_double()
  )
)

daily_coverage_all <- read_csv(
  file.path(paths$curated_dir, "daily_coverage.csv"),
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
) |>
  mutate(era5_rain_log = log1p(total_precipitation_00_08_mm))

daily_coverage <- daily_coverage_all |>
  filter(ringing_happened, season >= analysis_start_season)

# Define comparable capture process --------------------------------------

# Swallows and martins were sometimes caught with a separate targeted
# daytime process. They are excluded from both focal species and the daily
# comparison total.
excluded_capture_groups <- read_csv(
  file.path(paths$dataset_dir, "config", "daily_counts", "targeted_capture_groups.csv"),
  show_col_types = FALSE
)

comparable_counts_all <- daily_counts |>
  filter(is.na(avibase_id) | !avibase_id %in% excluded_capture_groups$avibase_id)

comparable_counts <- comparable_counts_all |>
  filter(season >= analysis_start_season)

analysis_period_decision <- tibble(
  primary_start_season = analysis_start_season,
  excluded_seasons = "1969–1976",
  decision = "Begin with the first full season after the 1976 introduction of intensive night netting",
  evidence = paste(
    "1969–1975 catches were mainly from southern dawn nets (operations evidence OH001).",
    "The historical synthesis states that most birds were caught at night from 1976 (OH004).",
    "The contemporary 1976/77 report describes increased night effort during that first season,",
    "with some dates still virtually confined to dawn catching; 1976 is therefore treated as transitional."
  )
)

excluded_early_seasons <- comparable_counts_all |>
  filter(season < analysis_start_season) |>
  group_by(season) |>
  summarise(
    n_positive_dates = n_distinct(ringing_date),
    total_comparable_count = sum(n_records),
    .groups = "drop"
  ) |>
  mutate(exclusion_reason = case_when(
    season <= 1975 ~ "pre-night-net capture regime",
    season == 1976 ~ "first night-net season; mixed transition in contemporary report"
  ))

minimum_total_count <- 500
minimum_seasons <- 10
minimum_positive_days <- 100

species_support <- comparable_counts |>
  filter(!is.na(avibase_id), !is.na(common_name)) |>
  group_by(avibase_id, common_name) |>
  summarise(
    total_count = sum(n_records),
    n_seasons_observed = n_distinct(season),
    n_positive_days = n_distinct(ringing_date),
    first_season = min(season),
    last_season = max(season),
    .groups = "drop"
  ) |>
  mutate(
    selected = total_count >= minimum_total_count &
      n_seasons_observed >= minimum_seasons &
      n_positive_days >= minimum_positive_days,
    selection_rule = paste0(
      "total_count >= ", minimum_total_count,
      "; seasons >= ", minimum_seasons,
      "; positive_days >= ", minimum_positive_days
    )
  ) |>
  arrange(desc(total_count), common_name)

focal_species <- species_support |>
  filter(selected) |>
  arrange(desc(total_count), common_name) |>
  mutate(
    abundance_rank = row_number(),
    category_order = abundance_rank
  )

other_category <- tibble(
  avibase_id = "other_comparable_taxa",
  common_name = "Other comparable taxa",
  total_count = NA_real_,
  abundance_rank = NA_integer_,
  category_order = nrow(focal_species) + 1L
)

category_metadata <- bind_rows(
  focal_species |>
    transmute(
      avibase_id, common_name, total_count, abundance_rank,
      category_order, is_focal_species = TRUE
    ),
  other_category |>
    transmute(
      avibase_id, common_name, total_count, abundance_rank,
      category_order, is_focal_species = FALSE
    )
)

# Complete the daily count matrix ----------------------------------------

model_days <- daily_coverage |>
  transmute(
    ringing_date,
    season,
    season_centered = season - median(season),
    day_of_season = as.integer(ringing_date - as.Date(sprintf("%s-10-01", season))),
    moon_distance_from_new_moon,
    mist_probability_light_patchy,
    mist_probability_good,
    era5_rain_log,
    wind_speed_10m_mean_ms,
    temperature_2m_mean_c,
    surface_pressure_mean_hpa
  )

focal_daily_counts <- expand_grid(
  ringing_date = model_days$ringing_date,
  avibase_id = focal_species$avibase_id
) |>
  left_join(
    comparable_counts |>
      semi_join(focal_species, by = "avibase_id") |>
      select(ringing_date, avibase_id, count = n_records),
    by = c("ringing_date", "avibase_id")
  ) |>
  mutate(count = replace_na(count, 0)) |>
  left_join(focal_species |> select(avibase_id, common_name, category_order), by = "avibase_id")

other_daily_counts <- comparable_counts |>
  filter(ringing_date %in% model_days$ringing_date) |>
  anti_join(focal_species, by = "avibase_id") |>
  group_by(ringing_date) |>
  summarise(count = sum(n_records), .groups = "drop") |>
  right_join(model_days |> select(ringing_date), by = "ringing_date") |>
  mutate(
    count = replace_na(count, 0),
    avibase_id = other_category$avibase_id,
    common_name = other_category$common_name,
    category_order = other_category$category_order
  )

model_data <- bind_rows(focal_daily_counts, other_daily_counts) |>
  left_join(model_days, by = "ringing_date") |>
  group_by(ringing_date) |>
  mutate(total_comparable_count = sum(count)) |>
  ungroup() |>
  filter(total_comparable_count > 0) |>
  arrange(ringing_date, category_order)

# Standardize continuous daily covariates once for fitting and prediction.
standardized_covariates <- c(
  "moon_distance_from_new_moon",
  "era5_rain_log",
  "wind_speed_10m_mean_ms",
  "temperature_2m_mean_c",
  "surface_pressure_mean_hpa"
)

daily_covariate_scaling <- model_data |>
  distinct(ringing_date, across(all_of(standardized_covariates))) |>
  summarise(across(
    all_of(standardized_covariates),
    list(center = mean, scale = sd),
    .names = "{.col}_{.fn}"
  )) |>
  pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(center|scale)$"
  )

for (variable in standardized_covariates) {
  center <- daily_covariate_scaling$center[daily_covariate_scaling$variable == variable]
  scale <- daily_covariate_scaling$scale[daily_covariate_scaling$variable == variable]
  model_data[[paste0(variable, "_z")]] <- (model_data[[variable]] - center) / scale
}

# Define the common within-season prediction window from empirical support.
minimum_reference_seasons <- ceiling(0.40 * n_distinct(model_data$season))

day_support <- model_data |>
  distinct(season, day_of_season) |>
  count(day_of_season, name = "n_seasons") |>
  mutate(in_reference_window = n_seasons >= minimum_reference_seasons)

reference_days <- day_support |>
  filter(in_reference_window) |>
  pull(day_of_season)

coverage_by_season <- model_data |>
  distinct(ringing_date, season, day_of_season, total_comparable_count) |>
  group_by(season) |>
  summarise(
    n_positive_dates = n(),
    first_day_of_season = min(day_of_season),
    median_day_of_season = median(day_of_season),
    last_day_of_season = max(day_of_season),
    total_comparable_count = sum(total_comparable_count),
    .groups = "drop"
  )

# QA checks ---------------------------------------------------------------

expected_rows <- n_distinct(model_data$ringing_date) * nrow(category_metadata)
stopifnot(nrow(model_data) == expected_rows)
stopifnot(all(model_data$total_comparable_count > 0))
stopifnot(all(reference_days == seq(min(reference_days), max(reference_days))))
stopifnot(all(model_data$count >= 0))

# Write model inputs ------------------------------------------------------

write_csv(model_data, file.path(model_data_dir, "composition_model_data.csv"))
write_csv(category_metadata, file.path(model_data_dir, "composition_categories.csv"))
write_csv(daily_covariate_scaling, file.path(model_data_dir, "daily_covariate_scaling.csv"))
write_csv(day_support, file.path(table_dir, "day_of_season_support.csv"))
write_csv(coverage_by_season, file.path(table_dir, "composition_coverage_by_season.csv"))
write_csv(species_support, file.path(table_dir, "species_model_support.csv"))
write_csv(excluded_capture_groups, file.path(table_dir, "excluded_capture_groups.csv"))
write_csv(analysis_period_decision, file.path(table_dir, "analysis_period_decision.csv"))
write_csv(excluded_early_seasons, file.path(table_dir, "excluded_early_seasons.csv"))

cli_alert_success(
  "Prepared {nrow(focal_species)} focal species across {n_distinct(model_data$ringing_date)} positive-catch dates."
)
cli_alert_info(
  "The common prediction window is day {min(reference_days)} to {max(reference_days)} after 1 October."
)
cli_alert_info("Primary inference begins in season {analysis_start_season}; seasons 1969–1976 are excluded.")
