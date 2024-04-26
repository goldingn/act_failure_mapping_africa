#' Create a greta array for the hierarchical regression coefficients between
#' SNPs and covariates
#'
#' @title Define 'beta' Parameter
#' @param n_snp number of SNPs
#' @return a greta array of dimension `n_snp` x 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta <- function(n_snp,
                              mu,
                              sd) {
  
  # use hierarchical decentring specification, to remove prior correlation from
  # the posterior
  beta_raw <- greta::normal(0, 1, dim = c(n_snp, 3))
  beta_raw_sd <- greta::sweep(beta_raw, 2, sd, FUN = "*")
  beta <- greta::sweep(beta_raw_sd, 2, mu, FUN = "+")
  beta
  
}
