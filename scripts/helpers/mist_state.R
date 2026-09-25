# Apply the same observed-state constraints in curation and held-out models.
constrain_mist_probabilities <- function(probability, observation) {
  exact <- match(observation, c("none", "light_patchy", "good"))
  rows <- which(!is.na(exact))
  probability[rows, ] <- 0
  probability[cbind(rows, exact[rows])] <- 1
  present <- observation %in% "present_unspecified"
  probability[present, "none"] <- 0
  probability[present, ] <- probability[present, , drop = FALSE] /
    rowSums(probability[present, , drop = FALSE])
  probability
}

predict_mist_probabilities <- function(model, data) {
  data <- dplyr::mutate(data, cloud_base_height_km = cloud_base_height_mean_m / 1000)
  requireNamespace("nnet")
  probability <- as.matrix(stats::predict(model, newdata = data, type = "probs"))
  probability <- probability[, c("none", "light_patchy", "good"), drop = FALSE]
  constrain_mist_probabilities(probability, data$mist_observation)
}

draw_mist_state <- function(data, probability = NULL) {
  if (is.null(probability)) {
    probability <- as.matrix(data[, c("mist_probability_none", "mist_probability_light_patchy", "mist_probability_good")])
  }
  random_value <- runif(nrow(data))
  state <- ifelse(random_value <= probability[, 1], "none",
    ifelse(random_value <= probability[, 1] + probability[, 2], "light_patchy", "good"))
  dplyr::mutate(data, mist_state = factor(state, levels = c("none", "light_patchy", "good")))
}
