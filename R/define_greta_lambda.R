#' Create a greta array for the loadings between SNPs and latent factors
#'
#' @title Define 'Lambda' Parameter
#' @param n_snp number of SNPs
#' @param n_latent number of latent factors
#' @return a greta array of dimension `n_snp` x `n_latent`
#' @importFrom greta zeros normal
#' @author Nick Golding
define_greta_lambda <- function(n_snp = n_snp, n_latent = n_latent) {

  # create lambda, subject to constraints to prevent rotational invariance when
  # fitting
  # upper triangular zero
  loadings <- greta::zeros(n_snp, n_latent)
  # diagonal positive
  diag(loadings) <- greta::normal(0, 1, truncation = c(0, Inf), dim = n_latent)
  # lower triangular unconstrained
  lower <- lower.tri(loadings)
  loadings[lower] <- greta::normal(0, 1, dim = sum(lower))

  loadings
}
