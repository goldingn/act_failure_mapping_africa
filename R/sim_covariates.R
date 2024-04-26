#' Simulate synthetic covariate values as a SpatRaster object
#'
#' @title Simulate Synthetic Covariates
#' @param mask a SpatRaster mask layer within which to simulate the values
#' @return
#' @author Nick Golding
sim_covariates <- function(mask) {

  dim <- dim(mask)[1:2]
  
  # simulate on logit scale
  logit_PfPR <- -3 + 0.5 * sim_gp(dim = dim)
  logit_ACT_use <- sim_gp(dim = dim)

  # transform to unit interval  
  PfPR <- terra::app(logit_PfPR, plogis)
  ACT_use <- terra::app(logit_ACT_use, plogis)

  # combine and return
  covariates <- terra::rast(c(PfPR = PfPR,
                              ACT_use = ACT_use))
  
  covariates <- terra::mask(covariates, mask)
  covariates
}
