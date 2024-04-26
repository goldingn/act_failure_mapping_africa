#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title Simulate a grid from a Gaussian process
#' @param dim vector of length two giving the dimensions of the simulated
#'   grid of values
#' @param smooth whether to simulate a smooth (vs rough) GP
#' @return
#' @importFrom fields sim.rf matern.image.cov
#' @importFrom terra rast
#' @author Nick Golding
sim_gp <- function(dim = dim,
                   smooth = FALSE,
                   range = 0.07) {

  # define a grid
  grid <- list(x = seq(0, 1, length.out = dim[1]),
               y = seq(0, 1, length.out = dim[2]))
  
  # efficiently simulate a GP on this grid, using a block circulant embedding
  # using fields package
  
  # determine smoothness
  smoothness <- ifelse(smooth, 2, 0.5)
  
  # create the appropriate object
  obj <- fields::matern.image.cov(grid = grid,
                                  aRange = range,
                                  smoothness = smoothness,
                                  setup = TRUE)
  gp_mat <- fields::sim.rf(obj)

  # return as a SpatRaster object  
  gp_rast <- terra::rast(gp_mat)
  gp_rast
  
}

