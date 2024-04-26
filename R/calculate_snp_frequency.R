#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title Create SNP Frequency Rasters
#' @param latent_factors A multilayer SpatRaster with the values of the latent
#'   factors
#' @param covariates A multilayer SpatRaster with the values of the covariates
#' @param parameters A named list of model parameters
#' @return
#' @author Nick Golding
#' @export
calculate_snp_frequency <- function(latent_factors,
                                    covariates,
                                    parameters) {

  # Extract the latent factors and covariates into matrices of non-missing
  # values
  cell_ids <- terra::cells(latent_factors)

  # epsilon matrix  
  latents <- terra::extract(latent_factors, cell_ids)
  
  # design matrix
  covs <- terra::extract(covariates, cell_ids)
  X <- cbind(1, as.matrix(covs))
  
  # combine these with parameters to get matrix of logit SNP frequencies
  logit_snp_freq_mat <- X %*% t(parameters$beta) +
    t(parameters$loadings %*% t(latents))
  
  # make a SpatRaster in which to store the results
  n_snp <- nrow(parameters$q)
  logit_snp_frequency <- terra::rast(replicate(n_snp, latent_factors[[1]])) * 0
  names(logit_snp_frequency) <- paste0("snp_", seq_len(n_snp))
  for (snp in seq_len(n_snp)) {
    logit_snp_frequency[[snp]][cell_ids] <- logit_snp_freq_mat[, snp]
  }

  # convert to frequency scale and return
  snp_frequency <- app(logit_snp_frequency, plogis)
  snp_frequency

}

