#' Simulate a synthetic spatial dataset of ACT treatment failure and genetic
#' markers of resistance
#'
#' Simulates data according to the described model, and provides both simulated
#' data dnt he true parameters and maps to compare with
#'
#' @title Simulate Synthetic Data
#' @param n_snp number of SNPs
#' @param n_snp_obs number of observations of SNPs to simulate
#' @param snp_samples sample size of each SNP observation
#' @param n_tf_obs number of observations of treatment failure to simulate
#' @param tf_samples sample size of each treatment failure observation
#' @param resolution vector of length two giving the dimensions of the simulated
#'   grid of values
#' @param n_latent the number of latent factors to simulate
#' @param snp_biased whether the locations of SNP samples should be biased
#'   towards areas with higher rates of treatment failure
#' @return
#' @author Nick Golding
#' @export
sim_dataset <- function(n_snp = 10,
                        n_snp_obs = 100,
                        snp_samples = 100,
                        n_tf_obs = 100,
                        tf_samples = 100,
                        dim = c(100, 100),
                        n_latent = 3,
                        snp_biased = TRUE) {
  
  # simulate a mask to represent land vs sea, and encode the resolution
  mask <- sim_mask(dim = dim)
  
  # simulate covariates as GPs
  covariates <- sim_covariates(mask = mask)
  
  # simulate latent spatial factors
  latent_factors <- sim_latent_factors(mask = mask,
                                       n_latent = n_latent)
  
  # simulate model parameters
  parameters <- sim_parameters(n_snp = n_snp,
                               n_latent = n_latent)
  
  # calculate true SNP frequencies
  snp_frequency <- calculate_snp_frequency(latent_factors = latent_factors,
                                           covariates = covariates,
                                           parameters = parameters)
  
  # calculate true treatment failure frequencies
  tf_frequency <- calculate_tf_frequency(snp_frequency = snp_frequency,
                                         parameters = parameters)

  # simulate treatment failure data
  tf_data <- sim_tf_data(tf_frequency = tf_frequency,
                         n_tf_obs = n_tf_obs,
                         tf_samples = tf_samples)

  # simulate SNP data (possibly biased to high TF areas)
  snp_data <- sim_snp_data(snp_frequency = snp_frequency,
                           tf_frequency = tf_frequency,
                           n_snp_obs = n_snp_obs,
                           snp_samples = snp_samples,
                           snp_biased = snp_biased)

  # return all data and true values
  list(
    data = list(
      tf_data = tf_data,
      snp_data = snp_data
    ),
    truth = list(
      parameters = parameters,
      covariates = covariates,
      latent_factors = latent_factors,
      snp_frequency = snp_frequency,
      tf_frequency = tf_frequency
    )
  )
  
}
