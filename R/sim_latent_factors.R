#' Simulate random latent factors for use in data simulation.
#' 
#' @title Simulate Latent Factors
#' @param mask a SpatRaster mask layer within which to simulate the values
#' @param n_latent the number of latent factors to simulate
#' @return SpatRaster object with `n_latent` layers containing random latent
#'   factors
#' @importFrom terra rast
#' @author Nick Golding
sim_latent_factors <- function(mask,
                               n_latent = 3) {

  # simulate multiple latents and combine in a SpatRaster
  dim <- dim(mask)[1:2]
  latents <- replicate(n_latent,
                       sim_gp(dim = dim,
                              range = 0.15))
  latents <- terra::rast(latents)
  names(latents) <- paste0("latent_",
                           seq_len(n_latent))

  # mask out non-land bits
  latents <- terra::mask(latents, mask)
  latents
  
}
