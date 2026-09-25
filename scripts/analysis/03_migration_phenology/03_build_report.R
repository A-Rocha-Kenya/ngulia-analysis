library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(htmltools)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "03_migration_phenology")
model_data_dir <- file.path(analysis_dir, "model_data")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- ngulia_figure_dir(file.path(analysis_dir, "figures"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read results ------------------------------------------------------------

cli_h1("Build within-species migration-phenology report")

daily_data <- read_csv(file.path(model_data_dir, "species_daily_model_data.csv"), show_col_types = FALSE, col_types = cols(ringing_date = col_date()))
raw_annual <- read_csv(file.path(model_data_dir, "annual_passage_quantiles.csv"), show_col_types = FALSE)
standardized_annual <- read_csv(file.path(model_data_dir, "species_standardized_annual_quantiles.csv"), show_col_types = FALSE)
trends <- read_csv(file.path(table_dir, "species_phenology_trends.csv"), show_col_types = FALSE)
model_comparison <- read_csv(file.path(table_dir, "species_phenology_model_comparison.csv"), show_col_types = FALSE)
quantile_support <- read_csv(file.path(table_dir, "species_phenology_quantile_support.csv"), show_col_types = FALSE)
full_record_sensitivity <- read_csv(file.path(table_dir, "full_record_phenology_sensitivity.csv"), show_col_types = FALSE)
condition_sensitivity <- read_csv(file.path(table_dir, "moon_mist_phenology_sensitivity.csv"), show_col_types = FALSE)
model_diagnostics <- read_csv(file.path(table_dir, "species_phenology_model_diagnostics.csv"), show_col_types = FALSE)
species_support <- read_csv(file.path(table_dir, "species_support.csv"), show_col_types = FALSE) |>
  filter(selected_for_annual_summary)

species_order <- species_support |>
  arrange(desc(abundance_rank)) |>
  pull(common_name)
quantile_labels <- c(q25 = "Early passage (q25)", q50 = "Median passage (q50)", q75 = "Late passage (q75)")
date_labels <- function(day) format(as.Date("2000-10-01") + day, "%d %b")

# Figures -----------------------------------------------------------------

pooled_timing <- quantile_support |>
  select(common_name, quantile, pooled_passage_day) |>
  pivot_wider(names_from = quantile, values_from = pooled_passage_day) |>
  mutate(common_name = factor(common_name, levels = species_order))

timing_plot <- ggplot(pooled_timing, aes(y = common_name)) +
  annotate("rect", xmin = model_diagnostics$reference_first_day, xmax = model_diagnostics$reference_last_day, ymin = -Inf, ymax = Inf, fill = "#E8F1FA") +
  geom_segment(aes(x = q25, xend = q75, yend = common_name), linewidth = 1.1, colour = "#5B7083") +
  geom_point(aes(x = q50), size = 2.5, colour = ngulia_colours[["red"]]) +
  geom_vline(xintercept = c(model_diagnostics$reference_first_day, model_diagnostics$reference_last_day), linetype = 2, colour = ngulia_colours[["blue"]]) +
  scale_x_continuous(labels = date_labels) +
  labs(title = "Seasonal timing", subtitle = "Pooled q25–q75 and median; shading is the prediction window", x = "Date", y = NULL) +
  ngulia_theme(base_size = 10) +
  theme(panel.grid.major.y = element_blank())

abundance_plot <- species_support |>
  mutate(common_name = factor(common_name, levels = species_order)) |>
  ggplot(aes(total_count, common_name)) +
  geom_segment(aes(x = 1, xend = total_count, yend = common_name), colour = "grey78") +
  geom_point(size = 2.4, colour = ngulia_colours[["blue"]]) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  labs(title = "Capture volume", subtitle = "Comparable captures, 1977–2023", x = "Total captures (log scale)", y = NULL) +
  ngulia_theme(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

descriptive_plot <- timing_plot + abundance_plot +
  plot_layout(widths = c(2.2, 1)) +
  plot_annotation(title = "Observed central passage differs among the 19 focal species")

coverage_plot <- daily_data |>
  filter(season >= 1977) |>
  distinct(ringing_date, season, day_of_season) |>
  ggplot(aes(day_of_season, season)) +
  annotate("rect", xmin = model_diagnostics$reference_first_day, xmax = model_diagnostics$reference_last_day, ymin = -Inf, ymax = Inf, fill = "#E8F1FA") +
  geom_point(size = 1.1, alpha = 0.8, colour = ngulia_colours[["blue"]]) +
  scale_x_continuous(labels = date_labels) +
  labs(
    title = "Ringing dates vary among seasons",
    subtitle = "Each point is a date with a positive comparable catch; the model conditions on the dates available within each season",
    x = "Date", y = "Season"
  ) +
  ngulia_theme()

trend_plot <- trends |>
  filter(supported) |>
  mutate(
    common_name = factor(common_name, levels = species_order),
    quantile_label = factor(unname(quantile_labels[quantile]), levels = unname(quantile_labels)),
    result = case_when(
      adjusted_p < 0.05 & slope_days_per_decade < 0 ~ "FDR-supported earlier",
      adjusted_p < 0.05 & slope_days_per_decade > 0 ~ "FDR-supported later",
      lower > 0 | upper < 0 ~ "nominal interval excludes zero",
      TRUE ~ "interval includes zero"
    )
  ) |>
  ggplot(aes(slope_days_per_decade, common_name)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, colour = "grey45") +
  geom_point(aes(fill = result), shape = 21, size = 2.7, colour = "black", stroke = 0.3) +
  facet_wrap(vars(quantile_label), nrow = 1) +
  scale_fill_manual(
    values = c(
      "FDR-supported earlier" = "#2C7BB6", "FDR-supported later" = "#D7191C",
      "nominal interval excludes zero" = "#E69F00", "interval includes zero" = "#D9D9D9"
    ),
    name = NULL
  ) +
  labs(
    title = "Estimated change in passage date",
    subtitle = "Within-species date model; moon and mist standardized; season-clustered 95% intervals",
    x = "Change in passage date (days per decade; negative is earlier)", y = NULL
  ) +
  ngulia_theme(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

annual_plot <- ggplot() +
  geom_point(
    data = raw_annual |> filter(eligible_q50) |> mutate(common_name = factor(common_name, levels = rev(species_order))),
    aes(season, passage_q50), colour = "grey60", alpha = 0.6, size = 1
  ) +
  geom_line(
    data = standardized_annual |> filter(quantile == "q50") |> mutate(common_name = factor(common_name, levels = rev(species_order))),
    aes(season, passage_day), colour = ngulia_colours[["red"]], linewidth = 0.75
  ) +
  facet_wrap(vars(common_name), ncol = 4) +
  scale_y_continuous(labels = date_labels) +
  labs(
    title = "Annual median passage dates",
    subtitle = "Grey: empirical median; red: fitted median under common moon and mist conditions",
    x = "Season", y = "Median passage date"
  ) +
  ngulia_theme(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 8))

condition_plot_data <- condition_sensitivity |>
  filter(quantile == "q50") |>
  select(common_name, model, slope_days_per_decade) |>
  pivot_wider(names_from = model, values_from = slope_days_per_decade) |>
  mutate(
    difference = `Moon- and mist-standardized` - `Unadjusted for moon and mist`,
    label = if_else(abs(difference) >= 0.25, common_name, ""),
    label_hjust = if_else(`Unadjusted for moon and mist` > 2.5, 1.05, -0.05)
  )

condition_plot <- condition_plot_data |>
  ggplot(aes(`Unadjusted for moon and mist`, `Moon- and mist-standardized`)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 2.6, colour = ngulia_colours[["blue"]]) +
  geom_text(aes(label = label, hjust = label_hjust), check_overlap = TRUE, nudge_y = 0.08, size = 3) +
  labs(
    title = "Moon and mist adjustment changes some median trends",
    subtitle = "Dashed line indicates no change after standardization",
    x = "Unadjusted median trend (days per decade)",
    y = "Standardized median trend (days per decade)"
  ) +
  ngulia_theme() +
  theme(plot.margin = margin(8, 28, 8, 8))

full_record_plot <- full_record_sensitivity |>
  filter(quantile == "q50") |>
  mutate(common_name = factor(common_name, levels = species_order)) |>
  ggplot(aes(slope_days_per_decade, common_name, colour = model)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(
    title = "Including 1969–1976 does not materially change most point estimates",
    subtitle = "The full-record model adjusts the seasonal curve for the dawn-net and 1976 transition regimes",
    x = "Median trend (days per decade)", y = NULL
  ) +
  ngulia_theme(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

figures <- list(
  "00_descriptive_phenology.png" = list(descriptive_plot, 13, 8),
  "01_sampling_dates.png" = list(coverage_plot, 11, 7),
  "02_species_phenology_trends.png" = list(trend_plot, 15, 9),
  "03_annual_median_passage.png" = list(annual_plot, 12, 14),
  "04_moon_mist_sensitivity.png" = list(condition_plot, 9, 7),
  "05_full_record_sensitivity.png" = list(full_record_plot, 10, 8)
)

for (file_name in names(figures)) {
  figure <- figures[[file_name]]
  ngulia_save(file.path(figure_dir, file_name), figure[[1]], width = figure[[2]], height = figure[[3]], dpi = 220)
}

# Report ------------------------------------------------------------------

median_table <- function(data) {
  data <- data |> filter(quantile == "q50") |> arrange(abundance_rank)
  tags$table(
    tags$thead(tags$tr(
      tags$th("Species"), tags$th("Days/decade"), tags$th("95% interval"), tags$th("FDR-adjusted p")
    )),
    tags$tbody(lapply(seq_len(nrow(data)), function(row) {
      tags$tr(
        tags$td(data$common_name[row]),
        tags$td(sprintf("%+.2f", data$slope_days_per_decade[row])),
        tags$td(sprintf("%+.2f to %+.2f", data$lower[row], data$upper[row])),
        tags$td(sprintf("%.3f", data$adjusted_p[row]))
      )
    }))
  )
}

n_median_nominal <- trends |> filter(quantile == "q50", supported, lower > 0 | upper < 0) |> nrow()
n_median_fdr <- trends |> filter(quantile == "q50", supported, adjusted_p < 0.05) |> nrow()
n_quartile_nominal <- trends |> filter(quantile != "q50", supported, lower > 0 | upper < 0) |> nrow()
n_quartile_fdr <- trends |> filter(quantile != "q50", supported, adjusted_p < 0.05) |> nrow()
n_change_aic <- model_comparison |>
  select(common_name, model, aic) |>
  pivot_wider(names_from = model, values_from = aic) |>
  summarise(n = sum(`Changing phenology + moon and mist` + 2 < `Stable phenology + moon and mist`)) |>
  pull(n)
median_full_difference <- full_record_sensitivity |>
  filter(quantile == "q50") |>
  select(common_name, model, slope_days_per_decade) |>
  pivot_wider(names_from = model, values_from = slope_days_per_decade) |>
  summarise(value = median(abs(`1969–2023, capture-regime adjusted` - `1977–2023 primary`))) |>
  pull(value)

report <- tags$html(
  tags$head(
    tags$title("Migration phenology at Ngulia, 1969–2023"),
    tags$style(HTML(paste(readLines(file.path(project_dir, "assets", "report.css")), collapse = "\n")))
  ),
  tags$body(
    tags$h1("Migration phenology at Ngulia"),
    tags$p(class = "lead", "Do individual migrant species now pass Ngulia earlier or later than they did in the past? This analysis models the distribution of capture dates separately for each species. The median passage date is the primary endpoint; the 25th and 75th percentiles describe changes in earlier and later passage."),
    tags$div(class = "result",
      tags$p(paste0(
        "Main result: from 1977–2023, ", n_median_nominal,
        " species had median-trend intervals that excluded zero before multiple-testing correction, but ", n_median_fdr,
        " retained support after controlling the false discovery rate across the 19 primary species tests."
      )),
      tags$p("The data therefore show substantial year-to-year movement in capture dates, but no species has a statistically resolved long-term median shift after accounting for moon, mist, mass-fall nights, uncertainty among seasons and the number of species examined.")),

    tags$h2("1. Observed passage"),
    tags$p("The 19 focal species were selected using abundance and temporal-support rules defined before trend fitting. The central 50% of passage is summarized by q25–q75, with q50 marking the date by which half of that season's captures had occurred."),
    tags$img(src = "figures/00_descriptive_phenology.png", alt = "Observed central passage timing and capture volume"),
    tags$p("Ringing did not occur on the same dates every year. The likelihood therefore compares the distribution of a species' captures only among the dates available in that season; it does not turn unknown dates into zero catches."),
    tags$img(src = "figures/01_sampling_dates.png", alt = "Ringing dates by season"),

    tags$h2("2. Why the median is primary"),
    tags$p("Both a mean and a median are unchanged by multiplying every daily count in a season by the same effort factor. The difference is sensitivity to the shape of the catch distribution. A few exceptionally large falls, or a small number of very early or late captures, can move the mean strongly. The median depends on the cumulative halfway point and is more stable for skewed or multi-peaked migration."),
    tags$p("The median is not universally superior: the mean uses all observations and can be more precise for a symmetric, consistently sampled distribution. Ngulia's mass-fall catches and incomplete tails favour the median. The q25 and q75 analyses retain information about changes in the first and second halves of passage without relying on the poorly observed q10 and q90 tails."),

    tags$h2("3. Species-level date model"),
    tags$h3("Conditioning within a species and season"),
    tags$p("For each species and season, the response is its vector of daily catches over known ringing dates from 5 November–30 December. The model conditions on that species' total within this fitting window. Consequently, annual abundance and any effort multiplier that is constant through a season affect the total but not the estimated passage-date distribution."),
    tags$p(tags$code("daily species catches | annual species total ~ conditional date distribution")),
    tags$p("The expected date pattern contains a smooth baseline passage curve and a year × date interaction. That interaction allows the passage distribution to move and change shape over the decades. Moon distance and the probabilities of light/patchy and good mist are fitted as observation-condition effects."),

    tags$h3("Mass falls and uncertainty"),
    tags$p("For each species, daily counts are capped at its 95th percentile before model fitting. This preserves ordinary variation in catch size but prevents a single exceptional fall from acting like hundreds or thousands of independent measurements of passage date. Confidence intervals use a season-clustered covariance estimate, treating seasons rather than individual birds as the principal independent replicates."),

    tags$h3("Standardized annual curves"),
    tags$p(paste0(
      "Each season is predicted over the same ", date_labels(model_diagnostics$reference_first_day), "–",
      date_labels(model_diagnostics$reference_last_day),
      " window with moon held at its average position and mist probabilities held at their long-term means. The predicted daily curve is normalized within that species and season, then q25, q50 and q75 are calculated."
    )),
    tags$p("The q50 tests are the pre-defined primary family and are corrected across 19 species. The q25 and q75 tests form a separate secondary family. This avoids allowing the less reliable tails to determine the headline conclusion."),

    tags$h2("4. Results for 1977–2023"),
    tags$p(paste0(
      "A changing date curve improved AIC by more than two units for ", n_change_aic,
      " of 19 species, showing that interannual shape variation is real and worth modelling. A flexible curve fitting better does not, however, imply a consistent linear shift in its median."
    )),
    tags$p(paste0(
      "For median passage, ", n_median_nominal, " intervals excluded zero before multiplicity correction and ", n_median_fdr,
      " remained supported afterwards. For the secondary q25/q75 family, the corresponding numbers were ",
      n_quartile_nominal, " and ", n_quartile_fdr, "."
    )),
    tags$img(src = "figures/02_species_phenology_trends.png", alt = "Species-level q25, q50 and q75 passage trends"),
    median_table(trends),
    tags$p("The annual series make the main pattern easier to see: year-to-year differences are often much larger than the fitted multi-decadal change."),
    tags$img(src = "figures/03_annual_median_passage.png", alt = "Annual empirical and standardized median passage dates"),

    tags$h2("5. Sensitivity analyses"),
    tags$h3("Moon and mist"),
    tags$p("Moon and mist are not assumed to average out perfectly over 47 seasons. The standardized model estimates their effects from the daily data and predicts every year under the same reference conditions. The comparison below shows how much this changes each median slope."),
    tags$img(src = "figures/04_moon_mist_sensitivity.png", alt = "Median slopes before and after moon and mist standardization"),

    tags$h3("Including 1969–1976"),
    tags$p("The primary trend starts in 1977 because catches through 1975 came mainly from southern dawn nets and 1976 was a transition to intensive night netting. Constant effort within an individual season does not remove a date-dependent difference between dawn and night catching."),
    tags$p(paste0(
      "A full-record sensitivity includes 1969–2023, fits separate date-curve adjustments for the dawn-net and transition regimes, and then predicts all years under the night-net regime. The median absolute difference from the primary q50 slopes is ",
      sprintf("%.2f", median_full_difference), " days per decade. The early years therefore add useful corroboration, but they are not allowed to define the primary result."
    )),
    tags$img(src = "figures/05_full_record_sensitivity.png", alt = "Primary and full-record median trend sensitivity"),

    tags$h2("6. Interpretation and limits"),
    tags$div(class = "caution",
      tags$p("The analysis estimates capture phenology. Mist-net catches depend on migration traffic, the probability that birds descend at Ngulia, and the ringing process. They are not a direct census of all birds passing overhead."),
      tags$p("Conditioning removes abundance and effort differences that multiply all dates in a season equally. It cannot remove an unrecorded effort change that preferentially affects early or late dates."),
      tags$p("Only 11 operated zero-catch dates are documented from 1977 onward. Dates with unknown operation status are therefore omitted rather than assumed to be zero."),
      tags$p("Moon and mist are standardized because they strongly affect capture opportunity. Broader weather variables are not removed from the primary biological trend because weather may genuinely alter migration timing rather than simply observation.")),
    tags$p("The most defensible conclusion is that the record contains strong interannual variation but does not resolve a long-term species-level passage shift after multiplicity correction. The full-record sensitivity supports rather than overturns that conclusion, while the 1977 start remains methodologically cleaner."),

    tags$h2("Data products"),
    tags$p(
      tags$a(href = "tables/species_phenology_trends.csv", "species phenology trends"), " · ",
      tags$a(href = "model_data/species_standardized_annual_quantiles.csv", "annual standardized quantiles"), " · ",
      tags$a(href = "tables/species_phenology_model_comparison.csv", "model comparison"), " · ",
      tags$a(href = "tables/moon_mist_phenology_sensitivity.csv", "moon and mist sensitivity"), " · ",
      tags$a(href = "tables/full_record_phenology_sensitivity.csv", "full-record sensitivity"), " · ",
      tags$a(href = "tables/passage_quantile_trends.csv", "empirical quantile trends"))
  )
)

rendered <- renderTags(report)
report_html <- sub(
  "<html>",
  paste0("<!doctype html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n", rendered$head, "\n</head>"),
  rendered$html,
  fixed = TRUE
)
writeLines(report_html, file.path(analysis_dir, "migration_phenology_analysis.html"))

cli_alert_success("Wrote report to {file.path(analysis_dir, 'migration_phenology_analysis.html')}")
