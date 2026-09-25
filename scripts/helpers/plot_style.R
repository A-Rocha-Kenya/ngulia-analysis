# Shared figure defaults; keep plot-specific scales, labels and layout local.
#
# Set NGULIA_FIGURE_VARIANT=dark_ppt to create a presentation-ready dark
# counterpart of an output.  It is written below a dark_ppt/ subdirectory, so
# the paper figures are never overwritten.
ngulia_figure_variant <- function() {
  Sys.getenv("NGULIA_FIGURE_VARIANT", unset = "paper")
}

ngulia_is_dark <- function() {
  ngulia_figure_variant() == "dark_ppt"
}

ngulia_figure_dir <- function(path) {
  if (ngulia_is_dark()) file.path(path, "dark_ppt") else path
}
ngulia_colours <- c(
  teal = "#176D7A", blue = "#1F5F8B", gold = "#D5942A",
  green = "#287D6B", red = "#B24C3D", purple = "#8B6BAE",
  muted = "#7A8C93", pale_teal = "#8BBFC4"
)

ngulia_palette <- function(dark = ngulia_is_dark()) {
  if (dark) {
    c(ink = "#EAF2F7", muted = "#B9C8D4", grid = "#26384A", paper = "#0B1320",
      teal = "#58D6E7", blue = "#7ABBEA", gold = "#F3C45B", no_data = "#172538")
  } else {
    c(ink = "#24343A", muted = "#52666F", grid = "#DCE4E5", paper = "#FFFFFF",
      teal = ngulia_colours[["teal"]], blue = ngulia_colours[["blue"]],
      gold = ngulia_colours[["gold"]], no_data = "#EEF4F6")
  }
}

ngulia_theme <- function(dark = ngulia_is_dark(), base_size = 11, base_family = "sans") {
  colours <- ngulia_palette(dark)
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = colours[["ink"]]),
      plot.background = ggplot2::element_rect(fill = colours[["paper"]], colour = NA),
      panel.background = ggplot2::element_rect(fill = colours[["paper"]], colour = NA),
      panel.grid.major = ggplot2::element_line(colour = colours[["grid"]], linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(colour = colours[["muted"]]),
      plot.caption = ggplot2::element_text(colour = colours[["muted"]], hjust = 0),
      axis.text = ggplot2::element_text(colour = colours[["muted"]]),
      axis.title = ggplot2::element_text(colour = colours[["ink"]]),
      strip.text = ggplot2::element_text(face = "bold", colour = colours[["ink"]]),
      legend.text = ggplot2::element_text(colour = colours[["ink"]]),
      legend.title = ggplot2::element_text(colour = colours[["ink"]]),
      legend.background = ggplot2::element_rect(fill = colours[["paper"]], colour = NA)
    )
}

ngulia_save <- function(filename, plot, width, height, dpi = 220, bg = ngulia_palette()[["paper"]], ...) {
  ggplot2::ggsave(
    filename, plot, width = width, height = height, dpi = dpi, bg = bg, ...
  )
}
