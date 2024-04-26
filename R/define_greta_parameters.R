#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
define_greta_parameters <- function(n_snp,
                                    n_latent) {
  

  beta_mu <- define_greta_beta_mu()
  beta_sd <- define_greta_beta_sd()
  
  list(
    # individual probabilities of failure per SNP
    q = define_greta_q(n_snp = n_snp),
    # SNP - latent factor loadings
    loadings = define_greta_lambda(n_snp = n_snp,
                                   n_latent = n_latent),
    # hierarchical regression coefficients
    beta_mu = beta_mu,
    beta_sd = beta_sd,
    beta = define_greta_beta(n_snp = n_snp,
                             mu = beta_mu,
                             sd = beta_sd)
  )

}
