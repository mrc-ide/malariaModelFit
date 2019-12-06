#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
#' @title Load sytem file for this package
#'
#' @description Load and return file from within the inst folder of this
#'   package.
#'
#' @param name name of file
#'
#' @export

mmfit_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'malariaModelFit', mustWork = TRUE)
  ret <- readRDS(name_full)
}

# -----------------------------------
#' @title Break coverage
#' 
#' @description Feed in a vector of breaks x, and min and max values of a range
#'   that falls within x. Output the proportion of each slice that is covered by
#'   the range. Useful for calculating quantities like prevalence in a given age
#'   range.
#'
#' @param breaks vector of breaks
#' @param range_min minimum value of range of interest
#' @param range_max maximum value of range of interest
#' 
#' @importFrom stats punif
#' @export

break_coverage <- function(breaks, range_min, range_max) {
  
  # get lower and upper breaks
  breaks0 <- breaks[-length(breaks)]
  breaks1 <- breaks[-1]
  
  # get total proportion of each break covered
  p0 <- punif(range_min, breaks0, breaks1, lower.tail = FALSE)
  p1 <- punif(range_max, breaks0, breaks1, lower.tail = FALSE)
  ret <- p0 - p1
  
  return(ret)
}

# -----------------------------------
#' @title Gaussian quadrature of normal density
#'
#' @description Return node values and weights from Gaussian quadrature of
#'   normal distribution with n nodes.
#'
#' @param n number of nodes
#'
#' @export

gq_normal <- function(n) {
    statmod::gauss.quad.prob(n, dist = "normal")
}


