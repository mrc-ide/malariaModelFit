#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib malariaModelFit
#' @import compiler
#' @import statmod
#' @importFrom Rcpp evalCpp
NULL

#------------------------------------------------
#' Dummy function
#'
#' This is a dummy function
#'
#' @param x Some parameter
#'
#' @export
#' @examples
#' dummy1()

dummy1 <- function() {
    
    print("R function working!")
    
    dummy1_cpp()
    
}
