#' Given a multilayer SpatRaster of SNP frequencies, simulate observations of
#' SNP frequencies, possibly biased towards areas of high treatment failure
#'
#' @title Simulate SNP Observations
#' @param snp_frequency a multilayer SpatRaster of SNP frequencies
#' @param tf_frequency a SpatRaster of treatment failure rates
#' @param n_snp_obs number of observations of SNPs to simulate
#' @param snp_samples sample size of each SNP observation
#' @param snp_biased whether the locations of SNP samples should be biased
#' @return a data.frame of the coordinates, SNP, number of SNP samples at each
#'   coordinate, and number that had the SNP
#' @importFrom terra spatSample extract
#' @importFrom tidyverse `%*%` bind_cols pivot_longer mutate select
#' @author Nick Golding
sim_snp_data <- function(snp_frequency,
                         tf_frequency,
                         n_snp_obs,
                         snp_samples,
                         snp_biased) {

  # sample SNP points, possibly weighted towards higher TF frequency
  method <- ifelse(snp_biased, "weights", "random")
  
  # randomly sample locations in the raster
  coords <- terra::spatSample(tf_frequency,
                              size = n_snp_obs,
                              method = method,
                              xy = TRUE,
                              values = FALSE,
                              na.rm = TRUE,
                              replace = FALSE)
  
  # extract the expected SNP frequencies
  snp_freq_mat <- terra::extract(snp_frequency,
                                 coords,
                                 ID = FALSE)
  
  coords %>%
    bind_cols(snp_freq_mat) %>%
    pivot_longer(
      cols = starts_with("snp_"),
      names_to = "SNP",
      values_to = "frequency"
    ) %>%
    mutate(
      tested = snp_samples
    ) %>%
    mutate(
      positive = rbinom(n(), size = tested, prob = frequency)
    ) %>%
    select(
      -frequency
    )
    
}
