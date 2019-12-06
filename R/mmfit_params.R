#------------------------------------------------
# overload print() function for mmfit_params
#' @noRd
print.mmfit_params <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for mmfit_params
#' @noRd
summary.mmfit_params <- function(x, ...) {
  
  # print as dataframe
  x_print <- as.data.frame(cbind(unclass(x)))
  names(x_print) <- "value"
  print(x_print)
  
}

