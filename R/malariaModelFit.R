#------------------------------------------------
#' @title Fits the Griffin et. al. (2014) malaria transmission model
#'
#' @description Contains functions for re-fitting the Griffin et al. (2014)
#'   model to data using MCMC.
#'
#' @references Griffin et. al. (2014). Estimates of the changing age-burden of
#'   Plasmodium falciparum malaria disease in sub-Saharan Africa.
#'   doi:10.1038/ncomms4136
#'
#' @name malariaModelFit
#' @docType package
NULL

#------------------------------------------------
#' @useDynLib malariaModelFit, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("malariaModelFit", libpath)  # nocov
}