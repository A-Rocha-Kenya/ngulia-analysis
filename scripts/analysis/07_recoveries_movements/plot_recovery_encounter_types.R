# Summarize recovery condition and mortality classification ------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
recoveries_path <- file.path(paths$curated_dir, "recoveries.csv")
output_dir <- ngulia_figure_dir(file.path(project_dir, "outputs", "analysis", "07_recoveries_movements", "figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data ------------------------------------------------------------------

recoveries <- read_csv(recoveries_path, show_col_types = FALSE)

dark <- c(background = "#0B1320", text = "#EAF2F7", grid = "#26384A")

condition_labels <- c(
  alive = "Alive",
  dead = "Dead",
  unknown = "Condition unknown"
)

mortality_labels <- c(
  intentional_human = "Intentional human cause",
  unspecified_killing = "Killing; cause unspecified",
  domestic_animal = "Domestic animal",
  unintentional_human = "Unintentional human cause",
  wild_predation = "Wild predation",
  unknown = "Cause unknown"
)

condition_summary <- recoveries |>
  count(encounter_condition, name = "n") |>
  mutate(
    label = condition_labels[encounter_condition],
    label = factor(label, levels = rev(condition_labels)),
    percentage = n / sum(n),
    annotation = sprintf("%s (%.1f%%)", n, 100 * percentage)
  )

mortality_summary <- recoveries |>
  filter(encounter_condition == "dead") |>
  count(mortality_cause_class, name = "n") |>
  mutate(
    label = mortality_labels[mortality_cause_class],
    label = factor(label, levels = rev(mortality_labels)),
    percentage = n / sum(n),
    annotation = sprintf("%s (%.1f%%)", n, 100 * percentage)
  )

# Plot -----------------------------------------------------------------------

condition_plot <- ggplot(condition_summary, aes(n, label, fill = label)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = annotation), hjust = -0.08, size = 3.4, colour = dark[["text"]]) +
  scale_fill_manual(values = c(
    "Alive" = "#009E73",
    "Dead" = "#D55E00",
    "Condition unknown" = "#999999"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Condition at encounter",
    subtitle = sprintf("All %s recovery and control records", nrow(recoveries)),
    x = "Number of records", y = NULL
  ) +
  ngulia_theme(dark = TRUE, base_size = 11) +
  theme(
    plot.background = element_rect(fill = dark[["background"]], colour = NA),
    panel.background = element_rect(fill = dark[["background"]], colour = NA),
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = dark[["grid"]]),
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", colour = dark[["text"]]),
    plot.subtitle = element_text(colour = dark[["text"]]), axis.title = element_text(colour = dark[["text"]]),
    axis.text = element_text(colour = dark[["text"]])
  )

mortality_plot <- ggplot(mortality_summary, aes(n, label, fill = label)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = annotation), hjust = -0.08, size = 3.4, colour = dark[["text"]]) +
  scale_fill_manual(values = c(
    "Intentional human cause" = "#D55E00",
    "Killing; cause unspecified" = "#E69F00",
    "Domestic animal" = "#CC79A7",
    "Unintentional human cause" = "#F0E442",
    "Wild predation" = "#009E73",
    "Cause unknown" = "#999999"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Reported mortality causes",
    subtitle = sprintf("Among %s encounters recorded as dead", sum(mortality_summary$n)),
    x = "Number of records", y = NULL
  ) +
  ngulia_theme(dark = TRUE, base_size = 11) +
  theme(
    plot.background = element_rect(fill = dark[["background"]], colour = NA),
    panel.background = element_rect(fill = dark[["background"]], colour = NA),
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = dark[["grid"]]),
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", colour = dark[["text"]]),
    plot.subtitle = element_text(colour = dark[["text"]]), axis.title = element_text(colour = dark[["text"]]),
    axis.text = element_text(colour = dark[["text"]])
  )

figure <- condition_plot / mortality_plot +
  plot_annotation(
    caption = "Classifications are conservative; source wording is retained in encounter_method.",
    theme = theme(
      plot.background = element_rect(fill = dark[["background"]], colour = NA),
      plot.caption = element_text(colour = dark[["text"]], hjust = 1)
    )
  )

ngulia_save(file.path(output_dir, "recovery_encounter_types.png"), figure, width = 8.5, height = 7, dpi = 300)
