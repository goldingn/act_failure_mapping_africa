#' Build up the design matrix from covariate rasters and an optiona set of
#' coordinates
#'
#' @title Build Covariate Design Matrix
#' @param covariates A multilayer SpatRaster with the values of the covariates
#' @param coords An optional matrix of coordinates at which to extract the
#'   covariate values and construct the design matrix
#' @param scale Should the covariates be scaled (to mean zero, standard
#'   deviation 1) before extraction?
#' @return
#' @author Nick Golding
#' @export
build_design_matrix <- function(covariates, coords = NULL, scale = TRUE) {

  # optionally scale covariates
  if (scale) {
    covariates <- scale(covariates)
  }
  
  # if coords is provided use those cells, otherwise use all cells
  if (is.null(coords)) {
    cell_ids <- terra::cells(covariates)
  } else {
    cell_ids <- terra::cellFromXY(covariates, coords)
  }
  
  # extract, pad with intercept dummy, and return
  covs <- terra::extract(covariates, cell_ids)
  
  cbind(Intercept = 1, as.matrix(covs))

}
