# First observed value whose cumulative count reaches the requested fraction.
# This is the empirical inverse-CDF convention used by the timing analyses.
weighted_quantile <- function(day, count, probability) {
  ordered <- order(day)
  day <- day[ordered]
  count <- count[ordered]
  day[which(cumsum(count) >= probability * sum(count))[1]]
}
