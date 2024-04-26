#' Create a greta array for the hierarchical mean of the regression coefficients
#' between SNPs and covariates
#'
#' @title Define 'mu' Parameter
#' @return a greta array of dimension 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta_mu <- function() {

  greta::normal(0, 3, dim = 3)

}
