
# -----------------------------------
#' breakCoverage
#'
#' Feed in a vector of breaks x, and min and max values of a range that falls within x. Output the proportion of each slice that is covered by the range. Useful for calculating things like prevalence in a given age range.
#'
#' @param breaks vector of breaks
#' @param range_min minimum value of range of interest
#' @param range_max maximum value of range of interest
#'
#' @export

breakCoverageTest <- function(breaks, range_min, range_max) {
    # get lower and upper breaks
    breaks0 <- breaks[-length(breaks)]
    breaks1 <- breaks[-1]
    
    # get total proportion of each break covered
    ret <- punif(range_min, breaks0, breaks1, lower.tail=F) - punif(range_max, breaks0, breaks1, lower.tail=F)
    
    return(ret)
}
