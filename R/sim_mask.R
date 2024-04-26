#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param resolution
#' @return
#' @author Nick Golding
#' @export
sim_mask <- function(dim) {

  # simulate a random, smooth, GP surface
  mask <- sim_gp(dim = dim,
                 smooth = TRUE)
  
  # mask out the lowest quantile
  quant <- terra::global(mask, function(x) quantile(x, 0.85))[[1]]
  mask[mask > quant] <- NA
  
  # tidy up and return
  mask <- mask * 0
  names(mask) <- "mask"
  
  mask
  
}
