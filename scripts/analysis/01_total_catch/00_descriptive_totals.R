library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "01_total_catch")
figure_dir <- ngulia_figure_dir(file.path(analysis_dir, "figures"))
table_dir <- file.path(analysis_dir, "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read data ---------------------------------------------------------------

cli_h1("Analyse total catches")

daily_counts <- read_csv(
  file.path(paths$curated_dir, "daily_counts.csv"),
  show_col_types = FALSE,
  guess_max = Inf
)

daily_coverage <- read_csv(
  file.path(paths$curated_dir, "daily_coverage.csv"),
  show_col_types = FALSE,
  guess_max = Inf,
  col_types = cols(ringing_date = col_date())
)

# Summarize annual catches ------------------------------------------------

excluded_capture_groups <- read_csv(
  file.path(project_dir, "config", "analysis", "excluded_capture_groups.csv"),
  show_col_types = FALSE
)

annual_totals <- daily_counts |>
  filter(!avibase_id %in% excluded_capture_groups$avibase_id) |>
  group_by(season) |>
  summarise(total_birds = sum(n_records), .groups = "drop") |>
  left_join(
    daily_coverage |>
      filter(ringing_happened) |>
      count(season, name = "positive_catch_nights"),
    by = "season"
  ) |>
  mutate(birds_per_positive_catch_night = total_birds / positive_catch_nights)

write_csv(annual_totals, file.path(table_dir, "annual_catch_totals.csv"))

# Plotting ----------------------------------------------------------------

annual_total_plot <- ggplot(annual_totals, aes(season, total_birds)) +
  geom_col(fill = ngulia_colours[["blue"]], width = 0.82) +
  scale_x_continuous(breaks = seq(1970, 2025, by = 5)) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.07))) +
  labs(
    title = "Annual nocturnal-migrant catch totals at Ngulia",
    subtitle = "Swallow and martin catches are excluded because they come from a separate targeted process",
    x = "Ringing season",
    y = "Birds recorded"
  ) +
  ngulia_theme()

positive_night_plot <- ggplot(annual_totals, aes(season, birds_per_positive_catch_night)) +
  geom_col(fill = ngulia_colours[["pale_teal"]], width = 0.82) +
  scale_x_continuous(breaks = seq(1970, 2025, by = 5)) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.07))) +
  labs(
    title = "Annual catch per positive-catch night",
    subtitle = "This is not effort-standardised because quantitative effort and operated zero-catch coverage are incomplete",
    x = "Ringing season",
    y = "Birds per positive-catch night"
  ) +
  ngulia_theme()

ngulia_save(file.path(figure_dir, "annual_catch_totals.png"), annual_total_plot, width = 11, height = 5.8, dpi = 220)
ngulia_save(file.path(figure_dir, "annual_catch_per_positive_night.png"), positive_night_plot, width = 11, height = 5.8, dpi = 220)

cli_alert_warning("Catch per positive-catch night is descriptive, not an effort-corrected abundance index.")
cli_alert_success("Wrote catch-total results to {.file {analysis_dir}}")
