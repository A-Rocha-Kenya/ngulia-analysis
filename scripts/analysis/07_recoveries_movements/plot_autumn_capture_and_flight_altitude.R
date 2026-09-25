# Compare autumn capture timing and flight altitude near Ngulia -----------

library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(rnaturalearth)
library(sf)
library(stringr)

project_dir <- here::here()
source(file.path(project_dir, "scripts", "helpers", "plot_style.R"))
source(file.path(project_dir, "scripts", "helpers", "data_paths.R"))
paths <- get_data_paths(project_dir)
trajectory_path <- file.path(paths$geolocator_dir, "ngulia_most_likely_paths.csv")
capture_path <- file.path(paths$curated_dir, "daily_counts.csv")
source_dir <- Sys.getenv("NGULIA_GEOLOCATOR_DIR", unset = file.path(paths$dataset_dir, "data", "01_raw", "external", "geolocator"))
pressure_path <- file.path(source_dir, "pressurepaths.csv")
edges_path <- file.path(source_dir, "edges.csv")
output_dir <- ngulia_figure_dir(file.path(project_dir, "outputs", "analysis", "07_recoveries_movements", "figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ngulia_latitude <- -3.01361
ngulia_longitude <- 38.21098
dark <- c(ink = "#EAF2F7", grid = "#26384A", autumn = "#D1495B", capture = "#E9C46A")
species_labels <- c("Marsh Warbler", "Eurasian Nightjar", "European Roller", "Barred Warbler", "Red-backed Shrike", "Thrush Nightingale")
ngulia_altitude <- 920
ngulia_ridge_altitude <- 1821

trajectory_paths <- read_csv(trajectory_path, show_col_types = FALSE)
species_lookup <- trajectory_paths |> distinct(tag_id, scientific_name)

trajectory_points <- trajectory_paths |>
  mutate(point_time = as.POSIXct(start, tz = "UTC") + as.numeric(difftime(as.POSIXct(end, tz = "UTC"), as.POSIXct(start, tz = "UTC"), units = "secs")) / 2) |>
  group_by(datapackage_id, tag_id) |>
  arrange(stap_id, .by_group = TRUE) |>
  mutate(next_lat = lead(lat), next_lon = lead(lon), next_time = lead(point_time)) |>
  ungroup()

winter_anchors <- trajectory_points |>
  filter(stopover_days > 30) |>
  group_by(datapackage_id, tag_id) |>
  arrange(lat, desc(stopover_days), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  transmute(datapackage_id, tag_id, winter_anchor_time = point_time)

autumn_crossings <- trajectory_points |>
  left_join(winter_anchors, by = c("datapackage_id", "tag_id")) |>
  filter(!is.na(next_lat), !is.na(next_lon), !is.na(next_time), abs(next_lat - lat) > 0.05, (lat - ngulia_latitude) * (next_lat - ngulia_latitude) <= 0) |>
  mutate(fraction = (ngulia_latitude - lat) / (next_lat - lat), crossing_time = point_time + (next_time - point_time) * fraction, crossing_lon = lon + (next_lon - lon) * fraction, passage_season = if_else(crossing_time < winter_anchor_time, "Autumn passage", "Spring passage")) |>
  filter(fraction >= 0, fraction <= 1, passage_season == "Autumn passage", next_lat < lat) |>
  mutate(
    common_name = recode(scientific_name, "Acrocephalus palustris" = "Marsh Warbler", "Caprimulgus europaeus" = "Eurasian Nightjar", "Coracias garrulus" = "European Roller", "Curruca nisoria" = "Barred Warbler", "Lanius collurio" = "Red-backed Shrike", "Luscinia luscinia" = "Thrush Nightingale"),
    common_name = factor(common_name, levels = species_labels),
    seasonal_date = as.Date(paste0("2000-", format(crossing_time, "%m-%d")))
  )

captures <- read_csv(capture_path, show_col_types = FALSE) |>
  filter(common_name %in% species_labels) |>
  mutate(ringing_date = as.Date(ringing_date), common_name = factor(common_name, levels = species_labels), seasonal_date = as.Date(paste0("2000-", format(ringing_date, "%m-%d")))) |>
  filter(month(ringing_date) >= 8, month(ringing_date) <= 12)

capture_plot <- ggplot(captures, aes(seasonal_date, weight = n_records)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = dark[["capture"]], colour = NA, alpha = 0.7) +
  geom_vline(data = autumn_crossings, aes(xintercept = as.numeric(seasonal_date)), colour = dark[["autumn"]], alpha = 0.62, linewidth = 0.9) +
  facet_wrap(~common_name, ncol = 2, scales = "free_y", drop = FALSE) +
  scale_x_date(limits = as.Date(c("2000-10-15", "2001-01-01")), date_breaks = "1 month", date_labels = "%b", expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "Autumn capture timing and inferred Ngulia passage", subtitle = "Gold: Ngulia capture distribution; red lines: individual autumn passage dates", x = NULL, y = "Density") +
  ngulia_theme(dark = TRUE, base_size = 10) +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA), strip.background = element_rect(fill = "#0B1320", colour = dark[["grid"]]), strip.text = element_text(colour = dark[["ink"]], face = "bold"), plot.title = element_text(colour = dark[["ink"]], face = "bold"), plot.subtitle = element_text(colour = dark[["ink"]]), axis.text = element_text(colour = dark[["ink"]]), axis.title = element_text(colour = dark[["ink"]]), panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25))

pressure_columns <- c("tag_id", "datetime", "stap_id", "type", "lat", "lon", "j", "pressure_tag", "label", "altitude", "surface_pressure", "surface_pressure_norm", "sunset", "sunrise", "location_name", "life_stage", "ind", "nb_sample", "temperature_2m", "dewpoint_temperature_2m", "geopotential", "land_sea_mask", "u_component_of_wind_10m", "v_component_of_wind_10m")
rg_bin <- if (file.exists("/Applications/ChatGPT.app/Contents/Resources/rg")) "/Applications/ChatGPT.app/Contents/Resources/rg" else "rg"
tag_pattern <- paste0("^(", paste(unique(trajectory_paths$tag_id), collapse = "|"), "),")
pressure_cmd <- paste(shQuote(rg_bin), "-N", shQuote(tag_pattern), shQuote(pressure_path))
pressure_raw <- fread(cmd = pressure_cmd, header = FALSE, col.names = pressure_columns, showProgress = FALSE) |> as_tibble() |>
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"), sunset = as.POSIXct(sunset, tz = "UTC"), sunrise = as.POSIXct(sunrise, tz = "UTC")) |>
  filter(type == "most_likely", is.finite(altitude), !is.na(stap_id)) |>
  left_join(species_lookup, by = "tag_id") |>
  mutate(distance_km = 2 * 6371 * asin(sqrt(sin((lat - ngulia_latitude) * pi / 360)^2 + cos(lat * pi / 180) * cos(ngulia_latitude * pi / 180) * sin((lon - ngulia_longitude) * pi / 360)^2)))

pressurepaths <- pressure_raw |>
  filter(stap_id %% 1 != 0, distance_km <= 200)

flight_altitude <- pressurepaths |>
  inner_join(autumn_crossings |> select(tag_id, crossing_time), by = "tag_id", relationship = "many-to-many") |>
  distinct(tag_id, datetime, .keep_all = TRUE)

write_csv(flight_altitude |> select(tag_id, scientific_name, datetime, stap_id, lat, lon, altitude, distance_km, crossing_time), file.path(output_dir, "autumn_flight_altitude_near_ngulia.csv"))

altitude_panel_labels <- flight_altitude |>
  group_by(scientific_name) |>
  summarise(panel_label = paste0(first(scientific_name), "\n", dplyr::n(), " flight points; ", n_distinct(tag_id), " individuals"), .groups = "drop")

flight_altitude <- flight_altitude |>
  left_join(altitude_panel_labels, by = "scientific_name") |>
  mutate(panel_label = factor(panel_label, levels = altitude_panel_labels$panel_label))

altitude_max <- ceiling(max(flight_altitude$altitude, na.rm = TRUE) / 500) * 500

altitude_plot <- ggplot(flight_altitude, aes(altitude)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = dark[["autumn"]], colour = NA, alpha = 0.78) +
  geom_vline(xintercept = ngulia_altitude, colour = dark[["ink"]], linewidth = 0.8) +
  geom_vline(xintercept = ngulia_ridge_altitude, colour = dark[["capture"]], linewidth = 0.8) +
  facet_wrap(~panel_label, nrow = 1, scales = "free_y") +
  scale_x_continuous(limits = c(0, altitude_max), labels = scales::label_number(suffix = " m")) +
  coord_flip() +
  labs(title = "Flight altitude near Ngulia during autumn passage", subtitle = "Fractional STAPs only; pressure-path altitude above mean sea level; points within 200 km of Ngulia. White: Ngulia (920 m); gold: Ngulia ridge (1,821 m)", x = "Pressure-derived altitude above mean sea level", y = "Density") +
  ngulia_theme(dark = TRUE, base_size = 10) +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA), strip.background = element_rect(fill = "#0B1320", colour = dark[["grid"]]), strip.text = element_text(colour = dark[["ink"]], face = "bold"), plot.title = element_text(colour = dark[["ink"]], face = "bold"), plot.subtitle = element_text(colour = dark[["ink"]]), axis.text = element_text(colour = dark[["ink"]]), axis.title = element_text(colour = dark[["ink"]]), panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25))

# Flight timing relative to local sunrise and sunset ---------------------
edge_passages <- fread(edges_path, showProgress = FALSE) |> as_tibble() |>
  filter(type == "most_likely", tag_id %in% species_lookup$tag_id) |>
  left_join(species_lookup, by = "tag_id") |>
  mutate(
    start_time = as.POSIXct(start, tz = "UTC"),
    end_time = as.POSIXct(end, tz = "UTC"),
    longitude_scale = cos(ngulia_latitude * pi / 180),
    dx = (lon_t - lon_s) * longitude_scale,
    dy = lat_t - lat_s,
    closest_fraction = pmax(0, pmin(1, ((ngulia_longitude - lon_s) * longitude_scale * dx + (ngulia_latitude - lat_s) * dy) / (dx^2 + dy^2))),
    closest_latitude = lat_s + closest_fraction * (lat_t - lat_s),
    closest_longitude = lon_s + closest_fraction * (lon_t - lon_s),
    closest_distance_km = 2 * 6371 * asin(sqrt(sin((closest_latitude - ngulia_latitude) * pi / 360)^2 + cos(closest_latitude * pi / 180) * cos(ngulia_latitude * pi / 180) * sin((closest_longitude - ngulia_longitude) * pi / 360)^2)),
    passage_time = start_time + as.numeric(difftime(end_time, start_time, units = "secs")) * closest_fraction,
    common_name = recode(scientific_name, "Acrocephalus palustris" = "Marsh Warbler", "Caprimulgus europaeus" = "Eurasian Nightjar", "Coracias garrulus" = "European Roller", "Curruca nisoria" = "Barred Warbler", "Lanius collurio" = "Red-backed Shrike", "Luscinia luscinia" = "Thrush Nightingale"),
    common_name = factor(common_name, levels = species_labels)
  ) |>
  filter(lat_t < lat_s, closest_distance_km <= 100, month(start_time) %in% c(10, 11, 12, 1), !(month(start_time) == 10 & day(start_time) < 15))

passage_solar <- edge_passages |>
  select(tag_id, stap_s, stap_t, lat_s, lon_s, lat_t, lon_t, closest_latitude, closest_longitude, start_time, end_time, passage_time, common_name, closest_distance_km) |>
  left_join(pressure_raw |> select(tag_id, datetime, sunrise, sunset), by = "tag_id", relationship = "many-to-many") |>
  mutate(time_difference = abs(as.numeric(difftime(datetime, passage_time, units = "secs")))) |>
  group_by(tag_id, stap_s, stap_t, start_time) |>
  slice_min(time_difference, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    start_local = with_tz(start_time, "Africa/Nairobi"),
    end_local = with_tz(end_time, "Africa/Nairobi"),
    passage_local = with_tz(passage_time, "Africa/Nairobi"),
    sunrise_local = with_tz(sunrise, "Africa/Nairobi"),
    sunset_local = with_tz(sunset, "Africa/Nairobi"),
    start_hour = hour(start_local) + minute(start_local) / 60 + second(start_local) / 3600,
    duration_hours = as.numeric(difftime(end_time, start_time, units = "hours")),
    start_hour_noon = if_else(start_hour < 12, start_hour + 24, start_hour),
    end_hour_noon = start_hour_noon + duration_hours,
    passage_hour_noon = start_hour_noon + as.numeric(difftime(passage_time, start_time, units = "hours")),
    day_night = if_else(passage_local >= sunrise_local & passage_local < sunset_local, "Day", "Night")
  )

write_csv(passage_solar |> select(tag_id, common_name, stap_s, stap_t, start_time, end_time, passage_time, passage_local, closest_distance_km, sunrise_local, sunset_local, day_night), file.path(output_dir, "autumn_passage_timing_day_night.csv"))

flight_levels <- passage_solar |>
  arrange(tag_id, start_time) |>
  transmute(flight_label = paste(tag_id, format(start_local, "%d %b %Y"), sep = " · ")) |>
  pull(flight_label)

passage_solar <- passage_solar |>
  mutate(flight_label = factor(paste(tag_id, format(start_local, "%d %b %Y"), sep = " · "), levels = flight_levels))
sunrise_hour <- passage_solar |> summarise(value = mean(hour(sunrise_local) + minute(sunrise_local) / 60, na.rm = TRUE)) |> pull(value)
sunset_hour <- passage_solar |> summarise(value = mean(hour(sunset_local) + minute(sunset_local) / 60, na.rm = TRUE)) |> pull(value)
species_colours <- c("Marsh Warbler" = "#E76F51", "Eurasian Nightjar" = "#F4A261", "European Roller" = "#2A9D8F", "Barred Warbler" = "#8AB17D", "Red-backed Shrike" = "#577590", "Thrush Nightingale" = "#C77DFF")

timing_plot <- ggplot() +
  annotate("rect", xmin = 12, xmax = sunset_hour, ymin = -Inf, ymax = Inf, fill = "#1D2A35") +
  annotate("rect", xmin = sunrise_hour + 24, xmax = 36, ymin = -Inf, ymax = Inf, fill = "#1D2A35") +
  geom_vline(xintercept = c(sunset_hour, sunrise_hour + 24), colour = dark[["grid"]], linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = passage_solar, aes(x = start_hour_noon, xend = end_hour_noon, y = flight_label, yend = flight_label, colour = common_name), linewidth = 3.2, lineend = "round", alpha = 0.95) +
  geom_point(data = passage_solar, aes(passage_hour_noon, flight_label, colour = common_name), shape = 8, size = 4.4, stroke = 0.8) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c("12", "15", "18", "21", "00", "03", "06", "09", "12"), expand = c(0, 0)) +
  coord_cartesian(xlim = c(12, 36), clip = "on") +
  scale_colour_manual(values = species_colours, name = "Species", drop = TRUE) +
  labs(title = "Flight timing near Ngulia", subtitle = "Southbound autumn flight edges passing within 100 km. Stars mark the projected closest approach to Ngulia; darker central band is night (mean local sunset/sunrise).", x = "Local time of day (hour; noon–noon axis)", y = "Flight (tag ID · date)") +
  ngulia_theme(dark = TRUE, base_size = 10) +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA), plot.title = element_text(colour = dark[["ink"]], face = "bold"), plot.subtitle = element_text(colour = dark[["ink"]]), axis.text = element_text(colour = dark[["ink"]]), axis.title = element_text(colour = dark[["ink"]]), panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25), legend.background = element_rect(fill = "#0B1320", colour = NA), legend.key = element_rect(fill = "#0B1320", colour = NA), legend.text = element_text(colour = dark[["ink"]]), legend.title = element_text(colour = dark[["ink"]]))

# Duration of flight STAPs entering the 200-km zone -----------------------
flight_segments <- pressure_raw |>
  filter(stap_id %% 1 != 0) |>
  mutate(flight_id = floor(stap_id)) |>
  group_by(tag_id, scientific_name, flight_id) |>
  summarise(start_time = min(datetime, na.rm = TRUE), end_time = max(datetime, na.rm = TRUE), duration_hours = as.numeric(difftime(end_time, start_time, units = "hours")), within_200km = any(distance_km <= 200, na.rm = TRUE), .groups = "drop") |>
  filter(within_200km, is.finite(duration_hours), duration_hours > 0) |>
  mutate(common_name = recode(scientific_name, "Acrocephalus palustris" = "Marsh Warbler", "Caprimulgus europaeus" = "Eurasian Nightjar", "Coracias garrulus" = "European Roller", "Curruca nisoria" = "Barred Warbler", "Lanius collurio" = "Red-backed Shrike", "Luscinia luscinia" = "Thrush Nightingale"), common_name = factor(common_name, levels = species_labels))

write_csv(flight_segments, file.path(output_dir, "flight_duration_within_200km_ngulia.csv"))

duration_plot <- ggplot(flight_segments, aes(duration_hours, fill = common_name)) +
  geom_histogram(bins = 35, colour = NA, alpha = 0.9, position = "stack") +
  scale_fill_hue(h = c(15, 375), c = 100, l = 65, drop = TRUE, name = "Species") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Flight duration near Ngulia", subtitle = "Fractional STAPs entering the 200-km zone; duration from first to last pressure-path record", x = "Flight duration (hours)", y = "Number of flight STAPs") +
  ngulia_theme(dark = TRUE, base_size = 10) +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA), plot.title = element_text(colour = dark[["ink"]], face = "bold"), plot.subtitle = element_text(colour = dark[["ink"]]), axis.text = element_text(colour = dark[["ink"]]), axis.title = element_text(colour = dark[["ink"]]), panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25), legend.background = element_rect(fill = "#0B1320", colour = NA), legend.key = element_rect(fill = "#0B1320", colour = NA), legend.text = element_text(colour = dark[["ink"]]), legend.title = element_text(colour = dark[["ink"]]))

# Autumn flight edges and stationary periods around Ngulia ----------------
map_limits <- list(x = c(34.5, 42), y = c(-6.5, 0.5))
map_tags <- passage_solar |> distinct(tag_id) |> pull(tag_id) |> as.character()

nearby_staps <- trajectory_paths |>
  filter(tag_id %in% map_tags, !is.na(lat), !is.na(lon)) |>
  mutate(
    stopover_days_plot = pmin(stopover_days, 30),
    common_name = recode(scientific_name, "Acrocephalus palustris" = "Marsh Warbler", "Caprimulgus europaeus" = "Eurasian Nightjar", "Coracias garrulus" = "European Roller", "Curruca nisoria" = "Barred Warbler", "Lanius collurio" = "Red-backed Shrike", "Luscinia luscinia" = "Thrush Nightingale"),
    common_name = factor(common_name, levels = species_labels)
  ) |>
  filter(between(lon, map_limits$x[[1]], map_limits$x[[2]]), between(lat, map_limits$y[[1]], map_limits$y[[2]]))

world <- ne_countries(scale = "medium", returnclass = "sf")

flight_map <- ggplot() +
  geom_sf(data = world, fill = "#1D2A35", colour = dark[["grid"]], linewidth = 0.18) +
  geom_segment(data = passage_solar, aes(lon_s, lat_s, xend = lon_t, yend = lat_t, colour = common_name), linewidth = 0.9, alpha = 0.88) +
  geom_point(data = nearby_staps, aes(lon, lat, size = stopover_days_plot, fill = common_name), shape = 21, colour = dark[["ink"]], stroke = 0.25, alpha = 0.9) +
  geom_point(data = passage_solar, aes(closest_longitude, closest_latitude, colour = common_name), shape = 8, size = 2.7, stroke = 0.75) +
  geom_point(data = tibble(lon = ngulia_longitude, lat = ngulia_latitude), aes(lon, lat), shape = 8, size = 4.3, colour = dark[["ink"]]) +
  annotate("text", x = ngulia_longitude + 0.08, y = ngulia_latitude - 0.1, label = "Ngulia", hjust = 0, colour = dark[["ink"]], fontface = "bold", size = 3.6) +
  coord_sf(xlim = map_limits$x, ylim = map_limits$y, expand = FALSE) +
  scale_colour_manual(values = species_colours, name = "Species", drop = TRUE) +
  scale_fill_manual(values = species_colours, name = "Species", drop = TRUE) +
  scale_size_continuous(name = "Stationary period (days)", range = c(1.5, 7), trans = "sqrt", breaks = c(1, 7, 30), labels = c("1", "7", "30+")) +
  guides(colour = guide_legend(order = 1), fill = "none", size = guide_legend(order = 2)) +
  labs(title = "Autumn flight edges and stationary periods near Ngulia", subtitle = "Southbound flight edges passing within 100 km of Ngulia; stars mark their projected closest approach. Points are stationary periods, scaled by duration.", x = NULL, y = NULL) +
  ngulia_theme(dark = TRUE, base_size = 11) +
  theme(plot.background = element_rect(fill = "#0B1320", colour = NA), panel.background = element_rect(fill = "#0B1320", colour = NA), plot.title = element_text(colour = dark[["ink"]], face = "bold"), plot.subtitle = element_text(colour = dark[["ink"]]), axis.text = element_text(colour = dark[["ink"]]), panel.grid = element_line(colour = dark[["grid"]], linewidth = 0.25), legend.position = "right", legend.background = element_rect(fill = "#0B1320", colour = NA), legend.key = element_rect(fill = "#0B1320", colour = NA), legend.text = element_text(colour = dark[["ink"]]), legend.title = element_text(colour = dark[["ink"]], face = "bold"))

ngulia_save(file.path(output_dir, "autumn_capture_vs_ngulia_passage.png"), capture_plot, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
ngulia_save(file.path(output_dir, "autumn_flight_altitude_near_ngulia.png"), altitude_plot, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
ngulia_save(file.path(output_dir, "autumn_passage_timing_day_night.png"), timing_plot, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
ngulia_save(file.path(output_dir, "flight_duration_within_200km_ngulia.png"), duration_plot, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
ngulia_save(file.path(output_dir, "autumn_flight_edges_and_staps_near_ngulia.png"), flight_map, width = 11.69, height = 8.27, dpi = 320, bg = "#0B1320")
