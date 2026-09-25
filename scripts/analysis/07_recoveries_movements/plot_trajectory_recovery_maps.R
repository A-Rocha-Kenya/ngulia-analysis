# Plot species trajectories, BirdLife ranges, and Ngulia recoveries --------

library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(rnaturalearth)
library(readr)
library(sf)
library(stringr)

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
trajectory_path <- file.path(paths$geolocator_dir, "ngulia_most_likely_paths.csv")
recoveries_path <- file.path(paths$curated_dir, "recoveries.csv")
crosswalk_path <- file.path(paths$dataset_dir, "config", "website", "ngulia_taxonomy_crosswalk.csv")
range_dir <- file.path(paths$reference_dir, "birdlife_ranges", "website")
output_dir <- ngulia_figure_dir(file.path(project_dir, "outputs", "analysis", "07_recoveries_movements", "figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trajectory_paths <- read_csv(trajectory_path, show_col_types = FALSE)
recoveries <- read_csv(recoveries_path, show_col_types = FALSE)
crosswalk <- read_csv(crosswalk_path, show_col_types = FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf")

species_lookup <- crosswalk |>
  filter(include_processing, include_recovery) |>
  transmute(
    scientific_name = avilist_scientific_name,
    avibase_id,
    ngulia_name = ngulia_english
  ) |>
  distinct(scientific_name, .keep_all = TRUE)

mapped_species <- trajectory_paths |>
  distinct(scientific_name) |>
  inner_join(species_lookup, by = "scientific_name") |>
  arrange(ngulia_name)

recoveries <- recoveries |>
  filter(curation_status != "needs_review", !is.na(other_latitude), !is.na(other_longitude)) |>
  transmute(
    avibase_id,
    other_longitude,
    other_latitude,
    direction_label = recode(direction, from_ngulia = "Ringed at Ngulia", to_ngulia = "Ringed elsewhere")
  )

ngulia_longitude <- 38.21098
ngulia_latitude <- -3.01361

dark <- c(
  ink = "#EAF2F7", land = "#20353D", outline = "#48646B", grid = "#26384A",
  spring = "#56B4E9", autumn = "#D1495B", unclear = "#A0A7AD",
  ringed = "#009E73", elsewhere = "#D55E00", range = "#E9C46A"
)

trajectory_points <- trajectory_paths |>
  mutate(
    point_time = as.POSIXct(start, tz = "UTC") +
      as.numeric(difftime(as.POSIXct(end, tz = "UTC"), as.POSIXct(start, tz = "UTC"), units = "secs")) / 2
  ) |>
  group_by(datapackage_id, tag_id) |>
  arrange(stap_id, .by_group = TRUE) |>
  mutate(
    next_lat = lead(lat), next_lon = lead(lon), next_time = lead(point_time),
    phase = case_when(
      !is.na(next_lat) & next_lat - lat > 0.05 ~ "Spring passage",
      !is.na(next_lat) & next_lat - lat < -0.05 ~ "Autumn passage",
      TRUE ~ "Stationary / unclear"
    )
  ) |>
  ungroup()

winter_anchors <- trajectory_points |>
  filter(stopover_days > 30) |>
  group_by(datapackage_id, tag_id) |>
  arrange(lat, desc(stopover_days), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  transmute(datapackage_id, tag_id, winter_anchor_time = point_time)

trajectory_points <- trajectory_points |>
  left_join(winter_anchors, by = c("datapackage_id", "tag_id")) |>
  mutate(
    segment_time = point_time + (next_time - point_time) / 2,
    phase = case_when(
      !is.na(next_lat) & next_lat - lat > 0.05 & segment_time >= winter_anchor_time ~ "Spring passage",
      !is.na(next_lat) & next_lat - lat < -0.05 & segment_time < winter_anchor_time ~ "Autumn passage",
      TRUE ~ "Stationary / unclear"
    )
  )

all_crossings <- trajectory_points |>
  filter(
    !is.na(next_lat), !is.na(next_lon), !is.na(next_time),
    abs(next_lat - lat) > 0.05,
    (lat - ngulia_latitude) * (next_lat - ngulia_latitude) <= 0
  ) |>
  mutate(
    fraction = (ngulia_latitude - lat) / (next_lat - lat),
    crossing_time = point_time + (next_time - point_time) * fraction,
    crossing_lon = lon + (next_lon - lon) * fraction,
    passage_season = if_else(crossing_time < winter_anchor_time, "Autumn passage", "Spring passage"),
    seasonal_date = as.Date(if_else(
      format(crossing_time, "%m-%d") == "02-29", "2001-02-28",
      paste0(if_else(month(crossing_time) >= 8, "2000-", "2001-"), format(crossing_time, "%m-%d"))
    ))
  ) |>
  filter(
    fraction >= 0, fraction <= 1,
    (crossing_time < winter_anchor_time & next_lat < lat) |
      (crossing_time >= winter_anchor_time & next_lat > lat)
  )

seasonal_limits <- all_crossings |>
  group_by(passage_season) |>
  summarise(date_min = min(seasonal_date) - 7, date_max = max(seasonal_date) + 7, .groups = "drop")

for (i in seq_len(nrow(mapped_species))) {
  species <- mapped_species[i, ]
  points <- trajectory_points |>
    filter(scientific_name == species$scientific_name)

  stops <- points |>
    filter(!is.na(stopover_days)) |>
    mutate(stopover_days_plot = pmin(stopover_days, 90))

  crossings <- all_crossings |>
    filter(scientific_name == species$scientific_name) |>
    mutate(tag_row = as.integer(factor(tag_id)))

  range <- st_read(file.path(range_dir, paste0(species$avibase_id, ".geojson")), quiet = TRUE)
  species_recoveries <- recoveries |>
    filter(avibase_id == species$avibase_id)

  p <- ggplot() +
    geom_sf(data = world, fill = dark[["land"]], colour = dark[["outline"]], linewidth = 0.14) +
    geom_sf(data = range, fill = scales::alpha(dark[["range"]], 0.34), colour = NA) +
    geom_path(
      data = points,
      aes(lon, lat, group = interaction(datapackage_id, tag_id), colour = phase),
      linewidth = 0.5, alpha = 0.72
    ) +
    geom_point(
      data = stops,
      aes(lon, lat, size = stopover_days_plot, fill = phase),
      shape = 21, colour = dark[["ink"]], stroke = 0.2, alpha = 0.86
    ) +
    geom_point(
      data = crossings,
      aes(x = crossing_lon, y = ngulia_latitude, colour = passage_season),
      shape = 4, size = 2.4, stroke = 0.8, inherit.aes = FALSE
    ) +
    geom_segment(
      data = species_recoveries,
      aes(ngulia_longitude, ngulia_latitude, xend = other_longitude, yend = other_latitude, colour = direction_label),
      linewidth = 0.42, alpha = 0.7
    ) +
    geom_point(
      data = species_recoveries,
      aes(other_longitude, other_latitude, colour = direction_label),
      size = 1.5, alpha = 0.84
    ) +
    geom_point(
      data = tibble(longitude = ngulia_longitude, latitude = ngulia_latitude),
      aes(longitude, latitude), shape = 8, size = 2.7, colour = dark[["ink"]]
    ) +
    annotate("text", x = -22, y = 71.5, label = species$ngulia_name, hjust = 0, vjust = 1,
      colour = dark[["ink"]], fontface = "bold", size = 5.2) +
    coord_sf(xlim = c(-25, 100), ylim = c(-40, 75), expand = FALSE) +
    scale_colour_manual(
      values = c(
        "Spring passage" = dark[["spring"]], "Autumn passage" = dark[["autumn"]],
        "Stationary / unclear" = dark[["unclear"]],
        "Ringed at Ngulia" = dark[["ringed"]], "Ringed elsewhere" = dark[["elsewhere"]]
      ), breaks = c("Spring passage", "Autumn passage", "Stationary / unclear", "Ringed at Ngulia", "Ringed elsewhere"), name = NULL, drop = FALSE
    ) +
    scale_fill_manual(
      values = c("Spring passage" = dark[["spring"]], "Autumn passage" = dark[["autumn"]], "Stationary / unclear" = dark[["unclear"]]),
      name = "Trajectory period"
    ) +
    scale_size_continuous(name = "Stationary period (days)", range = c(0.5, 6), trans = "sqrt", breaks = c(1, 7, 30, 90), labels = c("1", "7", "30", "90+")) +
    guides(
      colour = guide_legend(ncol = 1, byrow = TRUE),
      fill = guide_legend(order = 1, ncol = 1),
      size = guide_legend(order = 2, ncol = 1)
    ) +
    ngulia_theme(dark = TRUE, base_size = 11) +
    theme(
      plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA),
      axis.title = element_blank(), axis.text = element_text(colour = dark[["ink"]]),
      panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.35), legend.position = "right", legend.box = "vertical",
      legend.spacing.y = grid::unit(0.08, "cm"), legend.key.height = grid::unit(0.35, "cm"), legend.key.width = grid::unit(0.35, "cm"),
      legend.text = element_text(colour = dark[["ink"]], size = 8), legend.title = element_text(colour = dark[["ink"]], face = "bold", size = 8)
    )

  timeline <- crossings |>
    arrange(seasonal_date, tag_id) |>
    mutate(tag_row = as.integer(factor(tag_id, levels = unique(tag_id))))

  make_timeline <- function(data, season_label, season_title) {
    season_limits <- seasonal_limits |>
      filter(passage_season == season_label)
    if (nrow(data) == 0) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste("No", str_to_lower(season_title), "passage"), colour = dark[["ink"]], size = 3.2) +
        xlim(0, 1) + ylim(0, 1) + theme_void() +
        theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA)))
    }
    data <- data |>
      arrange(seasonal_date, tag_id) |>
      mutate(tag_id = factor(tag_id, levels = rev(unique(tag_id))), tag_row = as.integer(tag_id))
    ggplot(data, aes(seasonal_date, tag_row)) +
      geom_hline(yintercept = seq_along(levels(data$tag_id)), colour = dark[["grid"]], linewidth = 0.2) +
      geom_point(aes(colour = passage_season), size = 2.5) +
      scale_colour_manual(values = c("Spring passage" = dark[["spring"]], "Autumn passage" = dark[["autumn"]]), guide = "none") +
      scale_x_date(
        limits = c(season_limits$date_min, season_limits$date_max),
        date_breaks = "1 month", date_labels = "%b", expand = expansion(mult = c(0.01, 0.03))
      ) +
      scale_y_continuous(breaks = seq_along(levels(data$tag_id)), labels = levels(data$tag_id), expand = expansion(mult = c(0.08, 0.18))) +
      labs(title = season_title, x = NULL, y = NULL) +
      ngulia_theme(dark = TRUE, base_size = 8.5) +
      theme(
        plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA),
        plot.title = element_text(colour = dark[["ink"]], face = "bold", size = 9),
        panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25),
        axis.text = element_text(colour = dark[["ink"]]), axis.title = element_text(colour = dark[["ink"]])
      )
  }

  spring_timeline <- make_timeline(filter(timeline, passage_season == "Spring passage"), "Spring passage", "Spring")
  autumn_timeline <- make_timeline(filter(timeline, passage_season == "Autumn passage"), "Autumn passage", "Autumn")
  timeline_plot <- wrap_elements(
    full = (spring_timeline | autumn_timeline) +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(theme = theme(plot.background = element_rect(fill = "#0B1320", colour = NA), plot.margin = margin(0, 0, 0, 0)))
  )

  combined <- (p / timeline_plot + plot_layout(heights = c(4.6, 1.65))) &
    theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "#0B1320", colour = NA))
  file_name <- paste0("map_02_", str_replace_all(str_to_lower(species$scientific_name), " ", "_"))
  ngulia_save(file.path(output_dir, paste0(file_name, ".png")), combined, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
}
