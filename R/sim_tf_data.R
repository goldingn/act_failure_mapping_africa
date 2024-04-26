#' Given a raster of treatment failure rates, simulate observations of treatment
#' failure.
#'
#' @title Simulate Treatment Failure Observations
#' @param tf_frequency a SpatRaster of treatment failure rates
#' @param n_tf_obs number of observations of treatment failure to simulate
#' @param tf_samples sample size of each treatment failure observation
#' @return a data.frame of the coordinates, number of attempted treatments at
#'   each coordinate, and number of treatment failures among the samples
#' @importFrom terra spatSample
#' @importFrom tidyverse `%*%` bind_cols pivot_longer mutate select
#' @author Nick Golding
sim_tf_data <- function(tf_frequency = tf_frequency,
                        n_tf_obs = n_tf_obs,
                        tf_samples = tf_samples) {

  # randomly sample locations in the raster and extract the expected treatment
  # failure rate
  coords <- terra::spatSample(tf_frequency,
                              size = n_tf_obs,
                              method = "random",
                              xy = TRUE,
                              na.rm = TRUE,
                              replace = FALSE)
  
  # add the number of samples and simulate the numbers of treatment failures
  coords %>%
    mutate(
      treated = tf_samples) %>%
    mutate(
      failed = rbinom(nrow(coords),
                          size = treated,
                          prob = failure_frequency)) %>%
    select(-failure_frequency)
}
