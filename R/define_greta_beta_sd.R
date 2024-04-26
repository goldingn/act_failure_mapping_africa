#' Create a greta array for the hierarchical standard deviation of the
#' regression coefficients between SNPs and covariates
#'
#' @title Define 'sd' Parameter
#' @return a greta array of dimension 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta_sd <- function() {

  greta::normal(0, 1, dim = 3, truncation = c(0, Inf))

}
