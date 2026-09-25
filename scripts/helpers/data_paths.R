# Locate the published dataset checkout used by these analyses.
get_data_paths <- function(project_dir = here::here()) {
  dataset_dir <- Sys.getenv("NGULIA_DATASET_DIR", unset = file.path(dirname(project_dir), "ngulia-dataset"))
  list(
    dataset_dir = dataset_dir,
    curated_dir = file.path(dataset_dir, "data", "04_curated"),
    reference_dir = file.path(dataset_dir, "data", "02_reference"),
    geolocator_dir = file.path(dataset_dir, "data", "03_intermediate", "geolocator_paths"),
    analysis_output_dir = file.path(project_dir, "outputs", "analysis")
  )
}
