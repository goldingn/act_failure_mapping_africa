#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param snp_frequency
#' @param parameters
#' @return
#' @author Nick Golding
#' @export
calculate_tf_frequency <- function(snp_frequency,
                                   parameters) {
  
  # Extract the SNP frequencies into a matrix of non-missing values
  cell_ids <- terra::cells(snp_frequency)
  snp_freqs <- terra::extract(snp_frequency, cell_ids)
  
  # combine with the q parameters to get treatment failure probability
  individual_prob_failure <- sweep(snp_freqs, 2, parameters$q, FUN = "*") 
  tf_freq <- 1 - apply(1 - individual_prob_failure, 1, prod)
  
  # create a SpatRaster for the results
  tf_frequency <- snp_frequency[[1]] * 0
  names(tf_frequency) <- "failure_frequency"
  tf_frequency[cell_ids] <- tf_freq
  
  # return
  tf_frequency
  
}
