
#------------------------------------------------
# check that all required_names are in names
# (not exported)
#' @noRd
check_names <- function(names, required_names) {
  if (!all(required_names %in% names)) {
    missing_params <- setdiff(required_names, names)
    stop(sprintf("missing parameters: {%s}", paste(missing_params, collapse = ", ")))
  }
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
#' @export

break_coverage <- function(breaks, range_min, range_max) {
  
  # get lower and upper breaks
  breaks0 <- breaks[-length(breaks)]
  breaks1 <- breaks[-1]
  
  # get total proportion of each break covered
  ret <- punif(range_min, breaks0, breaks1, lower.tail=F) - punif(range_max, breaks0, breaks1, lower.tail=F)
  
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
    ret <- statmod::gauss.quad.prob(n, dist="normal")
    return(ret)
}

#------------------------------------------------
#' @title Load sytem file for this package
#'
#' @description Load and return file from within the inst folder of this
#'   package.
#'
#' @param name name of file
#'
#' @export

malariaModelFit_file <- function(name) {
    system.file(name, package="malariaModelFit", mustWork=TRUE)
}

