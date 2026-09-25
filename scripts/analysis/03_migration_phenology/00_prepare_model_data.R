library(dplyr)
library(tidyr)
library(readr)
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
figure_dir <- file.path(analysis_dir, "figures")

dir.create(model_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

analysis_start_season <- 1977L
minimum_total_count <- 500
minimum_seasons <- 10
minimum_positive_days <- 100
minimum_annual_count <- 20
minimum_annual_positive_days <- 5
minimum_boundary_buffer_days <- 3
full_record_start_season <- 1969L

# Read data ---------------------------------------------------------------

cli_h1("Prepare species passage-date data")

daily_counts <- read_csv(
  file.path(paths$curated_dir, "daily_counts.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)

daily_coverage <- read_csv(
  file.path(paths$curated_dir, "daily_coverage.csv"),
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
)

# Define the comparable capture process ----------------------------------

excluded_capture_groups <- read_csv(
  file.path(paths$dataset_dir, "config", "analysis", "excluded_capture_groups.csv"),
  show_col_types = FALSE
)

comparable_counts_all <- daily_counts |>
  filter(
    is.na(avibase_id) | !avibase_id %in% excluded_capture_groups$avibase_id
  ) |>
  mutate(day_of_season = as.integer(ringing_date - as.Date(sprintf("%s-10-01", season))))

comparable_counts <- comparable_counts_all |>
  filter(season >= analysis_start_season)

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
    selected_for_annual_summary = total_count >= minimum_total_count &
      n_seasons_observed >= minimum_seasons &
      n_positive_days >= minimum_positive_days,
    species_selection_rule = paste0(
      "total count >= ", minimum_total_count,
      "; observed seasons >= ", minimum_seasons,
      "; positive dates >= ", minimum_positive_days
    )
  ) |>
  arrange(desc(total_count), common_name) |>
  mutate(abundance_rank = row_number())

focal_species <- species_support |>
  filter(selected_for_annual_summary)

phenology_daily_counts <- comparable_counts |>
  semi_join(focal_species, by = "avibase_id") |>
  select(ringing_date, season, day_of_season, avibase_id, common_name, n_records) |>
  arrange(avibase_id, season, ringing_date)

season_coverage <- daily_coverage |>
  filter(season >= analysis_start_season, ringing_happened) |>
  mutate(day_of_season = as.integer(ringing_date - as.Date(sprintf("%s-10-01", season)))) |>
  group_by(season) |>
  summarise(
    first_positive_catch_day = min(day_of_season),
    last_positive_catch_day = max(day_of_season),
    n_positive_catch_dates = n(),
    .groups = "drop"
  )

# Estimate annual passage quantiles --------------------------------------

annual_quantiles <- phenology_daily_counts |>
  group_by(avibase_id, common_name, season) |>
  summarise(
    n_captures = sum(n_records),
    n_species_positive_dates = n(),
    passage_q25 = weighted_quantile(day_of_season, n_records, 0.25),
    passage_q50 = weighted_quantile(day_of_season, n_records, 0.50),
    passage_q75 = weighted_quantile(day_of_season, n_records, 0.75),
    .groups = "drop"
  ) |>
  left_join(season_coverage, by = "season") |>
  mutate(
    lower_boundary_buffer_days = passage_q25 - first_positive_catch_day,
    upper_boundary_buffer_days = last_positive_catch_day - passage_q75,
    minimum_sample_supported = n_captures >= minimum_annual_count &
      n_species_positive_dates >= minimum_annual_positive_days &
      !is.na(first_positive_catch_day) & !is.na(last_positive_catch_day),
    eligible_q25 = minimum_sample_supported &
      passage_q25 - first_positive_catch_day >= minimum_boundary_buffer_days,
    eligible_q50 = minimum_sample_supported &
      passage_q50 - first_positive_catch_day >= minimum_boundary_buffer_days &
      last_positive_catch_day - passage_q50 >= minimum_boundary_buffer_days,
    eligible_q75 = minimum_sample_supported &
      last_positive_catch_day - passage_q75 >= minimum_boundary_buffer_days,
    exclusion_reason = case_when(
      n_captures < minimum_annual_count ~ "fewer than 20 captures",
      n_species_positive_dates < minimum_annual_positive_days ~ "fewer than 5 positive species dates",
      !eligible_q25 & !eligible_q50 & !eligible_q75 ~ "all quantiles close to an observed season boundary",
      TRUE ~ NA_character_
    )
  ) |>
  left_join(focal_species |> select(avibase_id, abundance_rank), by = "avibase_id") |>
  arrange(abundance_rank, season)

trend_support <- annual_quantiles |>
  select(
    avibase_id, common_name, abundance_rank, season, n_captures,
    starts_with("eligible_q")
  ) |>
  pivot_longer(starts_with("eligible_q"), names_prefix = "eligible_", names_to = "quantile", values_to = "eligible") |>
  group_by(avibase_id, common_name, abundance_rank, quantile) |>
  summarise(
    total_count = sum(n_captures),
    n_seasons_with_captures = n(),
    n_eligible_seasons = sum(eligible),
    first_eligible_season = if_else(any(eligible), min(season[eligible]), NA_integer_),
    last_eligible_season = if_else(any(eligible), max(season[eligible]), NA_integer_),
    selected_for_trend = n_eligible_seasons >= minimum_seasons,
    .groups = "drop"
  ) |>
  arrange(abundance_rank)

# Complete species counts on known ringing dates -------------------------

model_dates <- comparable_counts_all |>
  filter(season >= full_record_start_season) |>
  group_by(ringing_date, season, day_of_season) |>
  summarise(total_comparable_count = sum(n_records), .groups = "drop") |>
  filter(total_comparable_count > 0) |>
  left_join(
    daily_coverage |>
      select(
        ringing_date, moon_distance_from_new_moon,
        mist_probability_light_patchy, mist_probability_good,
        total_precipitation_00_08_mm, wind_speed_10m_mean_ms,
        temperature_2m_mean_c, surface_pressure_mean_hpa,
        djp_team_size_minimum, playback_nocturnal_observed, net_sites_observed
      ),
    by = "ringing_date"
  ) |>
  mutate(
    capture_era = case_when(
      season <= 1975 ~ "dawn-net era",
      season == 1976 ~ "transition",
      TRUE ~ "night-net era"
    )
  )

species_daily_model_data <- expand_grid(
  ringing_date = model_dates$ringing_date,
  avibase_id = focal_species$avibase_id
) |>
  left_join(
    comparable_counts_all |>
      semi_join(focal_species, by = "avibase_id") |>
      select(ringing_date, avibase_id, count = n_records),
    by = c("ringing_date", "avibase_id")
  ) |>
  mutate(count = replace_na(count, 0)) |>
  left_join(focal_species |> select(avibase_id, common_name, abundance_rank), by = "avibase_id") |>
  left_join(model_dates, by = "ringing_date") |>
  arrange(abundance_rank, season, ringing_date)

analysis_decisions <- tribble(
  ~decision, ~value, ~reason,
  "First season", as.character(analysis_start_season), "First full season after the 1976 transition to intensive night netting",
  "Passage summaries", "25th, 50th and 75th percentiles", "Estimate early, central and late passage without relying on sparsely sampled distribution tails",
  "Annual minimum", paste(minimum_annual_count, "captures on", minimum_annual_positive_days, "dates"), "Avoid unstable annual quantiles",
  "Boundary buffer", paste(minimum_boundary_buffer_days, "days"), "Exclude quantiles pressed against the observed ringing window",
  "Trend form", "linear days per decade", "Estimate an interpretable broad multi-decadal shift",
  "Uncertainty", "season bootstrap", "Treat seasons, rather than individual birds, as independent replicates"
)

# QA checks ---------------------------------------------------------------

stopifnot(all(phenology_daily_counts$n_records > 0))
stopifnot(all(annual_quantiles$passage_q25 <= annual_quantiles$passage_q50))
stopifnot(all(annual_quantiles$passage_q50 <= annual_quantiles$passage_q75))
stopifnot(all(trend_support$n_eligible_seasons <= trend_support$n_seasons_with_captures))
stopifnot(all(species_daily_model_data$count >= 0))
stopifnot(nrow(species_daily_model_data) == nrow(model_dates) * nrow(focal_species))

# Write model inputs ------------------------------------------------------

write_csv(phenology_daily_counts, file.path(model_data_dir, "phenology_daily_counts.csv"))
write_csv(species_daily_model_data, file.path(model_data_dir, "species_daily_model_data.csv"))
write_csv(annual_quantiles, file.path(model_data_dir, "annual_passage_quantiles.csv"))
write_csv(species_support, file.path(table_dir, "species_support.csv"))
write_csv(trend_support, file.path(table_dir, "trend_support.csv"))
write_csv(season_coverage, file.path(table_dir, "season_coverage.csv"))
write_csv(excluded_capture_groups, file.path(table_dir, "excluded_capture_groups.csv"))
write_csv(analysis_decisions, file.path(table_dir, "analysis_decisions.csv"))

cli_alert_success(
  "Prepared annual passage quantiles for {nrow(focal_species)} focal species; {sum(trend_support$selected_for_trend)} species-quantile trends have sufficient support."
)
