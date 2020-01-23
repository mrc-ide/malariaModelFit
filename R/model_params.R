#------------------------------------------------
# overload print() function for model_params
#' @method print model_params
#' @noRd
print.model_params <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for model_params
#' @method summary model_params
#' @noRd
summary.model_params <- function(x, ...) {
  
  # print as dataframe
  x_print <- as.data.frame(cbind(unclass(x)))
  names(x_print) <- "value"
  print(x_print)
  
}

