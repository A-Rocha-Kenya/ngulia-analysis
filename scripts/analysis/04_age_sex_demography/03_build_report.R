library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(htmltools)
library(cli)

# Set paths ---------------------------------------------------------------

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)

analysis_dir <- file.path(paths$analysis_output_dir, "04_age_sex_demography")
model_data_dir <- file.path(analysis_dir, "model_data")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- ngulia_figure_dir(file.path(analysis_dir, "figures"))

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cli_h1("Build age and sex demography report")

# Read results ------------------------------------------------------------

age_data <- read_csv(file.path(model_data_dir, "age_daily_model_data.csv"), show_col_types = FALSE)
sex_data <- read_csv(file.path(model_data_dir, "adult_sex_daily_model_data.csv"), show_col_types = FALSE)
age_support <- read_csv(file.path(table_dir, "age_species_support.csv"), show_col_types = FALSE) |>
  filter(selected_age)
sex_support <- read_csv(file.path(table_dir, "adult_sex_species_support.csv"), show_col_types = FALSE)
annual_coverage <- read_csv(file.path(table_dir, "annual_recording_coverage.csv"), show_col_types = FALSE)
age_validation <- read_csv(file.path(table_dir, "age_model_validation.csv"), show_col_types = FALSE)
age_selections <- read_csv(file.path(table_dir, "age_model_selections.csv"), show_col_types = FALSE)
sex_selections <- read_csv(file.path(table_dir, "adult_sex_model_selections.csv"), show_col_types = FALSE)
age_diagnostics <- read_csv(file.path(table_dir, "age_model_diagnostics.csv"), show_col_types = FALSE)
sex_diagnostics <- read_csv(file.path(table_dir, "adult_sex_model_diagnostics.csv"), show_col_types = FALSE)
age_trends <- read_csv(file.path(table_dir, "age_decadal_trends.csv"), show_col_types = FALSE)
sex_trends <- read_csv(file.path(table_dir, "adult_sex_decadal_trends.csv"), show_col_types = FALSE)
age_annual <- read_csv(file.path(table_dir, "standardized_annual_first_year_proportions.csv"), show_col_types = FALSE)
sex_annual <- read_csv(file.path(table_dir, "standardized_annual_adult_male_proportions.csv"), show_col_types = FALSE)
age_phenology <- read_csv(file.path(table_dir, "age_phenology_predictions.csv"), show_col_types = FALSE)
sex_phenology <- read_csv(file.path(table_dir, "adult_sex_phenology_predictions.csv"), show_col_types = FALSE)
age_timing <- read_csv(file.path(table_dir, "raw_age_timing_summary.csv"), show_col_types = FALSE)
sex_timing <- read_csv(file.path(table_dir, "raw_adult_sex_timing_summary.csv"), show_col_types = FALSE)
age_sensitivity <- read_csv(file.path(table_dir, "age_trend_weighting_sensitivity.csv"), show_col_types = FALSE)
age_robustness <- read_csv(file.path(table_dir, "age_trend_robustness_summary.csv"), show_col_types = FALSE)

# Publication-ready summaries -------------------------------------------

species_order <- age_support |>
  arrange(desc(n_ageable)) |>
  pull(common_name)
sex_species_order <- sex_support |>
  filter(selected_sex) |>
  arrange(desc(n_adults_sexed)) |>
  pull(common_name)

pooled_age_summary <- age_data |>
  group_by(avibase_id, common_name) |>
  summarise(
    n_ageable = sum(n_total),
    first_year_fraction = sum(n_first_year) / sum(n_total),
    .groups = "drop"
  ) |>
  left_join(age_support |> select(avibase_id, age_recording_fraction), by = "avibase_id") |>
  arrange(desc(n_ageable))

pooled_sex_summary <- sex_data |>
  group_by(avibase_id, common_name) |>
  summarise(
    n_sexed_adults = sum(n_total),
    adult_male_fraction = sum(n_male) / sum(n_total),
    .groups = "drop"
  ) |>
  left_join(
    sex_support |> select(avibase_id, adult_sex_recording_fraction),
    by = "avibase_id"
  )

timing_result <- age_timing |>
  summarise(
    n_later = sum(median_day_difference > 0),
    n_equal = sum(median_day_difference == 0),
    maximum_delay = max(median_day_difference)
  )

robust_signals <- age_robustness |>
  filter(conclusion_stable, primary_direction != "uncertain") |>
  arrange(primary_direction, common_name)

legacy_age_summary <- tibble::tribble(
  ~common_name, ~legacy_first_year_percent,
  "Red-backed Shrike", 82,
  "Isabelline Shrike", 54,
  "River Warbler", 70,
  "Basra Reed Warbler", 49,
  "Marsh Warbler", 64,
  "Olive-tree Warbler", 55,
  "Willow Warbler", 71,
  "Garden Warbler", 62,
  "Barred Warbler", 52,
  "Greater Whitethroat", 58,
  "Thrush Nightingale", 71,
  "Common Nightingale", 74,
  "White-throated Robin", 44,
  "Spotted Flycatcher", 67
)

legacy_comparison <- legacy_age_summary |>
  left_join(
    pooled_age_summary |>
      transmute(common_name, current_first_year_percent = 100 * first_year_fraction),
    by = "common_name"
  ) |>
  mutate(change_percentage_points = current_first_year_percent - legacy_first_year_percent)

legacy_median_difference <- median(abs(legacy_comparison$change_percentage_points))
legacy_largest_changes <- legacy_comparison |>
  slice_max(abs(change_percentage_points), n = 2) |>
  arrange(desc(change_percentage_points))

weighting_report_table <- tibble::tribble(
  ~`Weighting choice`, ~`What one species-date contributes`, ~`Question it answers`,
  "Capped daily (primary)", "Up to 50 aged birds", "Uses more information from larger catches while limiting mass-fall leverage.",
  "Equal date", "One equal unit", "Would the trend remain if every capture date counted equally?",
  "Full individual", "All aged birds", "What happens if each bird is treated as an independent contribution?",
  "Remove largest 1%", "Delete the largest 1% of dates, then cap at 50", "Is the trend being driven by exceptional mass-fall dates?"
)

headline_trends <- age_trends |>
  filter(common_name %in% c("Marsh Warbler", "Willow Warbler", "Garden Warbler"))
marsh_trend <- headline_trends |> filter(common_name == "Marsh Warbler")
willow_trend <- headline_trends |> filter(common_name == "Willow Warbler")
garden_trend <- headline_trends |> filter(common_name == "Garden Warbler")

age_report_table <- pooled_age_summary |>
  transmute(
    Species = common_name,
    `Ageable birds` = comma(n_ageable),
    `Age recorded` = percent(age_recording_fraction, accuracy = 0.1),
    `First-year share` = percent(first_year_fraction, accuracy = 0.1)
  )

sex_report_table <- pooled_sex_summary |>
  transmute(
    Species = common_name,
    `Sexed adults` = comma(n_sexed_adults),
    `Adults sexed` = percent(adult_sex_recording_fraction, accuracy = 0.1),
    `Male among sexed adults` = percent(adult_male_fraction, accuracy = 0.1)
  )

# Figures ----------------------------------------------------------------

plot_theme <- ngulia_theme() + theme(legend.position = "top")

recording_plot <- annual_coverage |>
  pivot_longer(
    c(age_recording_fraction, adult_sex_recording_fraction),
    names_to = "measure", values_to = "fraction"
  ) |>
  mutate(measure = recode(
    measure,
    age_recording_fraction = "Age class recorded",
    adult_sex_recording_fraction = "Adult sex recorded"
  )) |>
  ggplot(aes(season, fraction, colour = measure)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_colour_manual(values = c(ngulia_colours[["teal"]], ngulia_colours[["gold"]])) +
  labs(
    title = "Age recording is stable; sex recording is selective",
    subtitle = "Initial captures of the 14 focal age species, 1991–2023",
    x = "Season", y = "Recording fraction", colour = NULL
  ) +
  plot_theme

age_phenology_plot <- age_phenology |>
  mutate(common_name = factor(common_name, levels = species_order)) |>
  ggplot(aes(day_of_season, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = ngulia_colours[["gold"]], alpha = 0.25) +
  geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.8) +
  facet_wrap(~ common_name, ncol = 3) +
  scale_x_continuous(
    breaks = c(31, 61, 92), labels = c("1 Nov", "1 Dec", "1 Jan")
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "Age composition changes through the passage season",
    subtitle = "Predicted first-year share at the middle season and average grounding conditions",
    x = NULL, y = "First-year share"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

timing_plot <- age_timing |>
  mutate(common_name = factor(common_name, levels = rev(species_order))) |>
  ggplot(aes(median_day_difference, common_name)) +
  geom_vline(xintercept = 0, colour = "#98a7ad") +
  geom_segment(aes(x = 0, xend = median_day_difference, yend = common_name),
               colour = "#9cc4c9", linewidth = 1.2) +
  geom_point(colour = ngulia_colours[["teal"]], size = 2.6) +
  scale_x_continuous(breaks = -2:6) +
  labs(
    title = "First-year passage is usually later than adult passage",
    subtitle = "Raw pooled median day; positive values mean first-year birds pass later",
    x = "First-year minus adult median passage day", y = NULL
  ) +
  plot_theme

age_annual_plot <- age_annual |>
  mutate(common_name = factor(common_name, levels = species_order)) |>
  ggplot(aes(season, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = ngulia_colours[["pale_teal"]], alpha = 0.28) +
  geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.65) +
  geom_point(aes(size = n_individuals), colour = ngulia_colours[["gold"]], alpha = 0.65) +
  facet_wrap(~ common_name, ncol = 3) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_size_continuous(range = c(0.4, 2.8), trans = "sqrt") +
  labs(
    title = "Standardized annual first-year proportions",
    subtitle = "Same dates and grounding conditions each season; point size shows aged sample",
    x = "Season", y = "First-year share", size = "Aged birds"
  ) +
  plot_theme

sensitivity_ranges <- age_sensitivity |>
  group_by(common_name) |>
  summarise(
    minimum = min(odds_ratio_per_decade),
    maximum = max(odds_ratio_per_decade),
    .groups = "drop"
  )

trend_order <- age_trends |>
  arrange(odds_ratio_per_decade) |>
  pull(common_name)

trend_plot <- age_trends |>
  left_join(sensitivity_ranges, by = "common_name") |>
  mutate(common_name = factor(common_name, levels = trend_order)) |>
  ggplot(aes(odds_ratio_per_decade, common_name)) +
  geom_vline(xintercept = 1, colour = "#98a7ad") +
  geom_segment(aes(x = minimum, xend = maximum, yend = common_name),
               linewidth = 2.6, colour = ngulia_colours[["gold"]], alpha = 0.55) +
  geom_segment(aes(x = lower, xend = upper, yend = common_name),
               linewidth = 0.7, colour = ngulia_colours[["teal"]]) +
  geom_point(aes(fill = direction), shape = 21, size = 2.8, colour = "white") +
  scale_x_log10(breaks = c(0.75, 0.9, 1, 1.1, 1.25, 1.5)) +
  scale_fill_manual(values = c(decrease = ngulia_colours[["red"]], uncertain = ngulia_colours[["muted"]], increase = ngulia_colours[["green"]])) +
  labs(
    title = "Long-term change in first-year odds",
    subtitle = "Points and 95% intervals use capped daily information; gold bars span four weighting sensitivities",
    x = "Odds ratio per decade", y = NULL, fill = "Primary result"
  ) +
  plot_theme

sex_phenology_plot <- sex_phenology |>
  mutate(common_name = factor(common_name, levels = sex_species_order)) |>
  ggplot(aes(day_of_season, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = ngulia_colours[["gold"]], alpha = 0.25) +
  geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.8) +
  facet_wrap(~ common_name, ncol = 2) +
  scale_x_continuous(
    breaks = c(31, 61, 92), labels = c("1 Nov", "1 Dec", "1 Jan")
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(
    title = "Adult sex composition through the passage season",
    subtitle = "Restricted to species with at least 80% of adults sexed",
    x = NULL, y = "Male share among sexed adults"
  ) +
  plot_theme

sex_annual_plot <- sex_annual |>
  mutate(common_name = factor(common_name, levels = sex_species_order)) |>
  ggplot(aes(season, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = ngulia_colours[["pale_teal"]], alpha = 0.28) +
  geom_line(colour = ngulia_colours[["teal"]], linewidth = 0.7) +
  geom_point(aes(size = n_individuals), colour = ngulia_colours[["gold"]], alpha = 0.7) +
  facet_wrap(~ common_name, ncol = 2) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_size_continuous(range = c(0.5, 3), trans = "sqrt") +
  labs(
    title = "Standardized annual adult male proportions",
    subtitle = "No species has a clear directional multi-decadal trend",
    x = "Season", y = "Male share among sexed adults", size = "Sexed adults"
  ) +
  plot_theme

figures <- list(
  "00_recording_coverage.png" = list(recording_plot, 9, 5),
  "01_age_phenology.png" = list(age_phenology_plot, 12, 14),
  "02_age_median_timing.png" = list(timing_plot, 8, 6.5),
  "03_annual_first_year_proportions.png" = list(age_annual_plot, 12, 14),
  "04_age_decadal_trends.png" = list(trend_plot, 9, 7),
  "05_adult_sex_phenology.png" = list(sex_phenology_plot, 9, 7),
  "06_annual_adult_male_proportions.png" = list(sex_annual_plot, 9, 7)
)

for (file_name in names(figures)) {
  figure <- figures[[file_name]]
  ngulia_save(
    file.path(figure_dir, file_name), figure[[1]],
    width = figure[[2]], height = figure[[3]], dpi = 220
  )
}

# Explanatory HTML report -------------------------------------------------

html_table <- function(data) {
  tags$table(
    tags$thead(tags$tr(lapply(names(data), tags$th))),
    tags$tbody(lapply(seq_len(nrow(data)), function(i) {
      tags$tr(lapply(data[i, ], function(value) tags$td(as.character(value))))
    }))
  )
}

model_report <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$title("Ngulia age and sex demography"),
    tags$style(HTML(paste(readLines(file.path(project_dir, "assets", "report.css")), collapse = "\n")))
  ),
  tags$body(
    tags$h1("Who passes Ngulia, and when?"),
    tags$p(class = "lede", paste0(
      "Age and adult sex composition of initial captures from ",
      min(age_data$season), "–", max(age_data$season),
      ", standardized for passage date and tested against mass-fall sensitivity."
    )),
    tags$div(
      class = "note",
      tags$b("Estimand: "),
      "The analysis estimates demographic composition among birds caught at Ngulia. First-year share can reflect breeding output, age-specific survival, routes, timing and capture susceptibility; it is not a direct productivity or population-abundance estimate."
    ),

    tags$h2("1. Why the 2014 summaries needed a model"),
    tags$p(
      "The main Ngulia synthesis used 1974–2010 catches to describe a roughly two-to-one excess of first-year birds, broad seasonal age bins and species-specific timing (",
      tags$a(href = "../../../data/02_reference/publications/pdfs/149917_394216_1_SM.pdf", "Pearson, Backhurst & Jackson 2014, pp. 23–25"),
      "). It reported adults featuring more strongly before mid-November, young Red-backed Shrikes peaking about two weeks after adults, and young Marsh and Olive-tree Warblers passing about five days later. River Warbler, Thrush Nightingale, White-throated Robin (then called Irania) and Spotted Flycatcher showed little or no separation."
    ),
    tags$p(
      "That summary pooled EURING ages 3/5 as first-year and 4/6 as adult, then compared broad date bins. It was biologically informative, but annual catches cover different portions of the passage season and a few mist-driven falls can contain thousands of birds. The new workflow retains the same transparent age definition while estimating passage timing continuously, standardizing every annual ratio over the same dates, and bounding the leverage of any one species-date."
    ),
    tags$p(
      "The descriptive species pattern has changed surprisingly little. Across the 14 matched species, the current 1991–2023 pooled first-year shares differ from the 2014 table by a median of only ",
      number(legacy_median_difference, accuracy = 0.1),
      " percentage points. The largest differences are Barred Warbler (",
      number(legacy_largest_changes$change_percentage_points[legacy_largest_changes$common_name == "Barred Warbler"], accuracy = 0.1, prefix = "+"),
      " points) and Willow Warbler (",
      number(legacy_largest_changes$change_percentage_points[legacy_largest_changes$common_name == "Willow Warbler"], accuracy = 0.1),
      " points). The periods overlap substantially, so this is continuity of the descriptive pattern, not an independent validation or a trend test."
    ),
    tags$p(
      "Later work strengthens the biological interpretation but also argues against a single expected direction. Ringing studies of western Palaearctic migrants found that autumn age separation depends on moult strategy and breeding latitude; a 2021 multi-species study found age separation in 11 of 25 autumn migrants, usually with adults earlier, while a 2022 analysis again identified moult, sex and food ecology as important. Tracking also shows that first-time migrants may depart later and use less direct routes. Sex differences are generally weaker or can reverse in autumn. These expectations are therefore tested species by species rather than imposed (",
      tags$a(href = "https://doi.org/10.1371/journal.pone.0147471", "Kiat & Izhaki 2016"), "; ",
      tags$a(href = "https://doi.org/10.1007/s00265-020-02957-3", "Wobker, Heim & Schmaljohann 2021"), "; ",
      tags$a(href = "https://doi.org/10.1007/s43388-022-00108-y", "Bozó et al. 2022"), "; ",
      tags$a(href = "https://doi.org/10.1371/journal.pone.0273686", "Patchett et al. 2022"), "; ",
      tags$a(href = "https://doi.org/10.1098/rspb.2018.2821", "Briedis et al. 2019"), ")."
    ),

    tags$h2("2. Recording support determines what can be tested"),
    tags$p(paste0(
      "Fourteen non-swallow species meet the preregistered support rule: at least 1,000 ageable birds, 20 seasons, 100 capture dates and 80% age recording. Their age-recording fractions range from ",
      percent(min(age_support$age_recording_fraction), accuracy = 0.1), " to ",
      percent(max(age_support$age_recording_fraction), accuracy = 0.1),
      ". Swallows and martins are excluded because targeted daytime catching is a different observation process. Retrapped or history-conflicted events are excluded."
    )),
    tags$p(
      "Sex is not missing at random: in several species adults can be sexed from plumage much more often than first-year birds. Sex models therefore use adults only and require at least 500 certainly sexed adults, 15 seasons and 80% adult sex recording. Only four species pass."
    ),
    tags$img(class = "figure", src = "figures/00_recording_coverage.png", alt = "Annual age and adult sex recording coverage"),
    html_table(age_report_table),

    tags$h2("3. Model and validation"),
    tags$div(class = "equation", HTML("Y<sub>sd</sub> / N<sub>sd</sub> &sim; quasi-Binomial(N<sup>eff</sup><sub>sd</sub>, p<sub>sd</sub>)")),
    tags$div(class = "equation", HTML("logit(p<sub>sd</sub>) = &alpha;<sub>s</sub> + f<sub>s</sub>(date) + g<sub>s</sub>(season) + &beta;<sub>s</sub>X<sub>d</sub>")),
    tags$p(
      "For each species-date, Y is the first-year count (or adult male count), N is the corresponding aged (or sexed-adult) total, and p is its composition. Information rises with sample size but is capped at 50 aged birds or 30 sexed adults per species-date. This prevents one mass fall from being treated as hundreds of independent demographic replicates. Quasi-binomial dispersion expands uncertainty for remaining within-date clustering."
    ),
    tags$ul(
      tags$li("Natural splines use three degrees of freedom for passage date and two for broad multi-decadal change."),
      tags$li("Five-fold validation holds out complete seasons. Timing and broad season change are the minimum model; moon, mist, rain, wind, temperature and pressure are retained only when they improve held-out log loss."),
      tags$li(paste0(
        "Daily grounding conditions were selected for ",
        sum(age_selections$selected_model_id == "M3"), " of ", nrow(age_selections),
        " age models and ", sum(sex_selections$selected_model_id == "M3"),
        " of ", nrow(sex_selections), " adult sex models."
      )),
      tags$li("Annual values combine a broad season smooth with partially pooled season deviations, average across the common supported date window, and hold grounding conditions fixed.")
    ),
    tags$p(
      "The capped-date, equal-date, full-individual and top-1%-date-deletion fits form the trend sensitivity. This responds directly to the well-known dependence of migrant age ratios on sampling design rather than assuming every caught bird is an independent draw (",
      tags$a(href = "https://doi.org/10.1093/condor/102.3.699", "Kelly & Finch 2000"), ")."
    ),

    tags$h2("4. First-year birds generally pass later"),
    tags$p(class = "result", paste0(
      "First-year median passage is later in ", timing_result$n_later, " of 14 focal species and equal in ",
      timing_result$n_equal, ". The largest raw median delay is ", timing_result$maximum_delay,
      " days. This repeated pattern is the clearest cross-species demographic result."
    )),
    tags$img(class = "figure", src = "figures/02_age_median_timing.png", alt = "Median timing difference between first-year and adult birds"),
    tags$img(class = "figure", src = "figures/01_age_phenology.png", alt = "Predicted within-season first-year proportions"),
    tags$p(
      "This broadly confirms the 2014 account, including little separation in River Warbler, Thrush Nightingale and White-throated Robin. The current raw median gap is two days for Marsh Warbler and three for Olive-tree Warbler, smaller than the approximately five days reported from the earlier period. Red-backed Shrike remains the clearest separation, although the six-day median gap here is not directly comparable with the earlier two-week difference between peak dates."
    ),
    tags$p(
      "The fitted curves are composition curves, not migration-intensity curves. An increasing first-year share can arise because first-year passage increases, adult passage declines, or both. They nevertheless show why a fixed late-November sample cannot be compared directly with a season extending into January."
    ),

    tags$h2("5. Annual age composition and long-term signals"),
    tags$img(class = "figure", src = "figures/03_annual_first_year_proportions.png", alt = "Standardized annual first-year proportions"),
    tags$h3("What are the four weighting choices?"),
    tags$p(
      "They are four ways of deciding how much influence each species-date has on the long-term trend. They use the same response, dates and covariates; only the leverage of small versus mass-fall catches changes. The gold bar in the next figure spans the four resulting trend estimates, not four confidence intervals."
    ),
    html_table(weighting_report_table),
    tags$img(class = "figure", src = "figures/04_age_decadal_trends.png", alt = "Decadal trends and weighting sensitivity"),
    tags$p(class = "result", paste0(
      "Three directional results persist under all four weighting choices: first-year odds decrease in Marsh Warbler and Willow Warbler and increase in Garden Warbler. The primary Basra Reed Warbler decrease is not stable to weighting choice and is therefore treated as sensitive, not a headline result."
    )),
    tags$p(
      "Under the primary capped-date analysis, the odds ratios per decade are ",
      number(marsh_trend$odds_ratio_per_decade, accuracy = 0.001), " for Marsh Warbler (",
      number(marsh_trend$lower, accuracy = 0.001), "–", number(marsh_trend$upper, accuracy = 0.001), "), ",
      number(willow_trend$odds_ratio_per_decade, accuracy = 0.001), " for Willow Warbler (",
      number(willow_trend$lower, accuracy = 0.001), "–", number(willow_trend$upper, accuracy = 0.001), "), and ",
      number(garden_trend$odds_ratio_per_decade, accuracy = 0.001), " for Garden Warbler (",
      number(garden_trend$lower, accuracy = 0.001), "–", number(garden_trend$upper, accuracy = 0.001), "). These correspond to about ",
      percent(1 - marsh_trend$odds_ratio_per_decade, accuracy = 0.1), " and ",
      percent(1 - willow_trend$odds_ratio_per_decade, accuracy = 0.1),
      " lower first-year odds per decade, versus ",
      percent(garden_trend$odds_ratio_per_decade - 1, accuracy = 0.1),
      " higher odds; they are changes in odds, not percentage-point changes in first-year share."
    ),
    tags$p(
      "These trends are best read as hypotheses about changing age composition along this migration corridor. A lower first-year share could reflect weaker production, differential juvenile mortality before Ngulia, a route shift, or a changing age-specific probability of attraction and capture. Linking the annual indices to breeding-ground climate or population monitoring should be a separate analysis with explicit source regions and lag structure."
    ),

    tags$h2("6. Adult sex composition is deliberately narrower"),
    html_table(sex_report_table),
    tags$img(class = "figure", src = "figures/05_adult_sex_phenology.png", alt = "Predicted adult male proportions through the season"),
    tags$img(class = "figure", src = "figures/06_annual_adult_male_proportions.png", alt = "Standardized annual adult male proportions"),
    tags$p(class = "result",
      "Adult male share shows no clear multi-decadal directional trend in Barred Warbler, Isabelline Shrike, Red-backed Shrike or White-throated Robin. Raw male-versus-female median passage differences range from -2 to +1 days. The scientifically useful result is the absence of a strong general sex signal in the defensible subset, not an all-catch sex ratio."
    ),

    tags$h2("7. Conclusions"),
    tags$div(
      class = "result",
      tags$ol(
        tags$li(tags$b("Timing is the strongest general result. "), "First-year birds pass later in 12 of 14 species. The 2014 pattern survives a larger, cleaned and explicitly standardized analysis, although current Marsh and Olive-tree median gaps are smaller."),
        tags$li(tags$b("The overall species composition is unexpectedly persistent. "), "Pooled first-year shares are typically within two percentage points of the 1974–2010 synthesis. That stability coexists with meaningful species-specific change."),
        tags$li(tags$b("The long-term signals are divergent, not a common drift. "), "Marsh and Willow Warblers show robust declines in first-year odds, while Garden Warbler increases under every weighting choice. Their opposite directions make a single site-wide recording change an incomplete explanation, but do not by themselves identify breeding productivity, survival or route change as the cause."),
        tags$li(tags$b("Sensitivity changes the scientific conclusion for Basra Reed Warbler. "), "Its primary decrease should remain a hypothesis because significance depends on how large fall dates are weighted."),
        tags$li(tags$b("The sex result is narrower but still useful. "), "Among the four species with consistently sexed adults, neither within-season separation nor a multi-decadal shift is strong. Extending this conclusion to all species would not be justified.")
      )
    ),
    tags$p(
      "Taken together, the data support a migration-process interpretation before a productivity interpretation: Ngulia is exceptionally informative about age-structured timing and corridor composition, while attribution to breeding output requires independent breeding-ground or survival data. A particularly useful next analysis would test whether the three robust species trends covary with source-region climate, population indices or migration-route change using pre-specified spatial origins and lags."
    ),

    tags$h2("8. Limits and reusable outputs"),
    tags$ul(
      tags$li("The record begins in 1991, so this analysis does not infer demography for the earlier Ngulia decades."),
      tags$li("Age is a plumage-based EURING observation. Codes 3/5 are first-year and 4/6 adult; unknown, uncertain and unusual older codes are excluded."),
      tags$li("Annual intervals condition on the selected model and do not include species-selection or age-classification uncertainty."),
      tags$li("Daily caps reduce pseudo-replication but do not identify family groups or shared flocks."),
      tags$li("Later lighting, vegetation and dawn-net changes can still affect age or sex groups differently; the outputs are standardized catch-composition indices, not census ratios.")
    ),
    tags$p(
      "Tables: ",
      tags$a(href = "tables/standardized_annual_first_year_proportions.csv", "annual age indices"), " · ",
      tags$a(href = "tables/age_decadal_trends.csv", "age trends"), " · ",
      tags$a(href = "tables/age_trend_weighting_sensitivity.csv", "weighting sensitivity"), " · ",
      tags$a(href = "tables/age_phenology_predictions.csv", "age phenology"), " · ",
      tags$a(href = "tables/standardized_annual_adult_male_proportions.csv", "annual adult sex indices"), " · ",
      tags$a(href = "tables/age_species_support.csv", "species support"), " · ",
      tags$a(href = "tables/age_model_validation.csv", "held-out validation")
    )
  )
)

rendered_report <- renderTags(model_report)
report_html <- sub(
  "<html>",
  paste0("<!doctype html>\n<html>\n<head>\n", rendered_report$head, "\n</head>"),
  rendered_report$html,
  fixed = TRUE
)
writeLines(report_html, file.path(analysis_dir, "age_sex_demography_analysis.html"))

cli_alert_success(
  "Wrote age and sex demography report to {file.path(analysis_dir, 'age_sex_demography_analysis.html')}"
)
