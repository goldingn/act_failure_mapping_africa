#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title Simulate Model Parameters
#' @param n_snp number of SNPs
#' @param n_latent the number of latent factors to simulate
#' @return
#' @importFrom greta calculate, cauchy
#' @author Nick Golding
sim_parameters <- function(n_snp = 10,
                           n_latent = n_latent) {

  # set up greta arrays for parameters
  params <- define_greta_parameters(n_snp = n_snp,
                                    n_latent = n_latent)
  
  # simulate parameters with greta
  calculate_args <- c(params, list(nsim = 1))
  parameters <- do.call(greta::calculate, calculate_args)

  # drop simulation dimension from greta params
  lapply(parameters, function(x) {
    dim(x) <- dim(x)[-1]
    x})
  
}
