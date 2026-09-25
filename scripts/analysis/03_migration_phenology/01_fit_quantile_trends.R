library(dplyr)
library(tidyr)
library(readr)
library(purrr)
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

n_bootstrap <- 2000L
quantile_columns <- c("passage_q25", "passage_q50", "passage_q75")
quantile_labels <- c(
  q25 = "Early passage (25th percentile)",
  q50 = "Median passage (50th percentile)",
  q75 = "Late passage (75th percentile)"
)

# Read prepared data ------------------------------------------------------

cli_h1("Fit species passage-date trends")

annual_quantiles <- read_csv(
  file.path(model_data_dir, "annual_passage_quantiles.csv"),
  show_col_types = FALSE
)
daily_counts <- read_csv(
  file.path(model_data_dir, "phenology_daily_counts.csv"),
  show_col_types = FALSE,
  col_types = cols(ringing_date = col_date())
)
trend_support <- read_csv(
  file.path(table_dir, "trend_support.csv"),
  show_col_types = FALSE
)

model_quantiles <- annual_quantiles |>
  select(
    avibase_id, common_name, abundance_rank, season,
    all_of(quantile_columns), starts_with("eligible_q")
  ) |>
  pivot_longer(
    matches("^(passage|eligible)_q"),
    names_to = c(".value", "quantile"),
    names_pattern = "^(passage|eligible)_(q[0-9]+)$"
  ) |>
  filter(eligible) |>
  semi_join(trend_support |> filter(selected_for_trend), by = c("avibase_id", "quantile")) |>
  rename(passage_day = passage) |>
  mutate(quantile_label = unname(quantile_labels[quantile]))

# Fit broad linear changes and bootstrap seasons -------------------------

bootstrap_slopes <- function(season, passage_day, n_bootstrap) {
  n <- length(season)
  sampled_rows <- matrix(sample.int(n, n * n_bootstrap, replace = TRUE), nrow = n)
  sampled_season <- matrix(season[sampled_rows], nrow = n)
  sampled_day <- matrix(passage_day[sampled_rows], nrow = n)
  centered_season <- sweep(sampled_season, 2, colMeans(sampled_season), "-")
  centered_day <- sweep(sampled_day, 2, colMeans(sampled_day), "-")
  10 * colSums(centered_season * centered_day) / colSums(centered_season^2)
}

fit_trend <- function(data, n_bootstrap) {
  model <- lm(passage_day ~ I((season - mean(season)) / 10), data = data)
  slope <- unname(coef(model)[2])
  draws <- bootstrap_slopes(data$season, data$passage_day, n_bootstrap)
  draws <- draws[is.finite(draws)]
  bootstrap_tail <- min(
    (sum(draws <= 0) + 1) / (length(draws) + 1),
    (sum(draws >= 0) + 1) / (length(draws) + 1)
  )

  tibble(
    n_seasons = nrow(data),
    first_season = min(data$season),
    last_season = max(data$season),
    slope_days_per_decade = slope,
    lower = unname(quantile(draws, 0.025)),
    upper = unname(quantile(draws, 0.975)),
    bootstrap_p = min(1, 2 * bootstrap_tail),
    total_fitted_change_days = slope * (max(data$season) - min(data$season)) / 10,
    r_squared = summary(model)$r.squared,
    model = list(model),
    bootstrap_draws = list(draws)
  )
}

set.seed(1977)
trend_fits <- model_quantiles |>
  group_by(avibase_id, common_name, abundance_rank, quantile, quantile_label) |>
  group_modify(~ fit_trend(.x, n_bootstrap)) |>
  ungroup() |>
  mutate(
    adjusted_p = p.adjust(bootstrap_p, method = "BH"),
    direction = case_when(
      adjusted_p < 0.05 & slope_days_per_decade < 0 ~ "earlier",
      adjusted_p < 0.05 & slope_days_per_decade > 0 ~ "later",
      TRUE ~ "no clear linear change"
    )
  ) |>
  arrange(abundance_rank, match(quantile, names(quantile_labels)))

# Test sensitivity to exceptional mass-catch dates -----------------------

daily_caps <- daily_counts |>
  group_by(avibase_id) |>
  summarise(daily_count_cap = quantile(n_records, 0.95, type = 1), .groups = "drop")

capped_annual_quantiles <- daily_counts |>
  inner_join(daily_caps, by = "avibase_id") |>
  mutate(capped_count = pmin(n_records, daily_count_cap)) |>
  group_by(avibase_id, common_name, season) |>
  summarise(
    passage_q25 = weighted_quantile(day_of_season, capped_count, 0.25),
    passage_q50 = weighted_quantile(day_of_season, capped_count, 0.50),
    passage_q75 = weighted_quantile(day_of_season, capped_count, 0.75),
    .groups = "drop"
  ) |>
  pivot_longer(all_of(quantile_columns), names_prefix = "passage_", names_to = "quantile", values_to = "passage_day") |>
  semi_join(model_quantiles, by = c("avibase_id", "season", "quantile"))

capped_trends <- capped_annual_quantiles |>
  group_by(avibase_id, common_name, quantile) |>
  summarise(
    capped_slope_days_per_decade = 10 * cov(season, passage_day) / var(season),
    .groups = "drop"
  )

trend_summary <- trend_fits |>
  select(-model, -bootstrap_draws) |>
  left_join(capped_trends, by = c("avibase_id", "common_name", "quantile")) |>
  mutate(
    capped_same_direction = sign(slope_days_per_decade) == sign(capped_slope_days_per_decade)
  )

median_summary <- trend_summary |>
  filter(quantile == "q50") |>
  arrange(slope_days_per_decade) |>
  select(
    avibase_id, common_name, abundance_rank, n_seasons, first_season, last_season,
    slope_days_per_decade, lower, upper, adjusted_p, direction,
    total_fitted_change_days, capped_slope_days_per_decade, capped_same_direction
  )

# Describe change in the central 50% passage duration --------------------

duration_data <- annual_quantiles |>
  filter(eligible_q25, eligible_q75) |>
  mutate(passage_day = passage_q75 - passage_q25) |>
  group_by(avibase_id, common_name, abundance_rank) |>
  filter(n() >= 10) |>
  ungroup()

set.seed(1978)
duration_trends <- duration_data |>
  group_by(avibase_id, common_name, abundance_rank) |>
  group_modify(~ fit_trend(.x, n_bootstrap)) |>
  ungroup() |>
  mutate(
    adjusted_p = p.adjust(bootstrap_p, method = "BH"),
    direction = case_when(
      adjusted_p < 0.05 & slope_days_per_decade < 0 ~ "shorter passage window",
      adjusted_p < 0.05 & slope_days_per_decade > 0 ~ "longer passage window",
      TRUE ~ "no clear linear change"
    )
  ) |>
  select(-model, -bootstrap_draws) |>
  rename(
    slope_duration_days_per_decade = slope_days_per_decade,
    total_fitted_duration_change_days = total_fitted_change_days
  ) |>
  arrange(abundance_rank)

# Save model objects and tables ------------------------------------------

saveRDS(
  list(
    trend_fits = trend_fits,
    duration_trends = duration_trends,
    quantile_labels = quantile_labels,
    n_bootstrap = n_bootstrap,
    daily_caps = daily_caps
  ),
  file.path(model_dir, "passage_date_quantile_trends.rds")
)

write_csv(model_quantiles, file.path(model_data_dir, "eligible_annual_passage_quantiles_long.csv"))
write_csv(capped_annual_quantiles, file.path(model_data_dir, "capped_annual_passage_quantiles_long.csv"))
write_csv(trend_summary, file.path(table_dir, "passage_quantile_trends.csv"))
write_csv(median_summary, file.path(table_dir, "median_passage_trends.csv"))
write_csv(duration_trends, file.path(table_dir, "passage_duration_trends.csv"))
write_csv(daily_caps, file.path(table_dir, "species_daily_count_caps.csv"))

cli_alert_success(
  "Fitted early, median and late passage trends for {n_distinct(trend_summary$avibase_id)} species using {n_bootstrap} season bootstrap replicates."
)
