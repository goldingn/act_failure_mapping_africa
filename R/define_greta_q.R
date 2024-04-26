#' Create a greta array for the individual probabilities of treatment failure
#' per SNP
#'
#' @title Define 'q' Parameter
#' @param n_snp number of SNPs
#' @return a greta array of dimension `n_snp`
#' @author Nick Golding
define_greta_q <- function(n_snp) {

  greta::cauchy(0, 0.1,
                dim = n_snp,
                truncation = c(0, 1))

}
