#' Extract function as text from .cpp file
#' 
#' Conversion of .cpp loglikelihood file to string for use with DrJacoby. The 
#' file provided should contain a single function which returns SEXP and has a name
#' starting "loglikelihood". The end of the function must be flagged with 
#' "// loglikelihood_end" on the final line"
#'
#' @param address File address
#'
#' @return function as string
#' @example An example .cpp file may look like:
#' # // Model likelihood file
#' # #include <Rcpp.h>
#' # [[Rcpp::export]]
#' # SEXP loglikelihood(std::vector<double> params, std::vector<double> x){
#' #    likelihood calculations
#' # }
#' # // loglikelihood_end
#' @export
function_to_string <- function(address){
  assert_string(address)
  if(!grepl(".cpp", address)){
    stop("The file at address should be a .cpp file")
  }
  
  # Read in file line by line
  function_text <- readr::read_lines(address)
  # Isolate where the function starts
  f_start <- which(grepl("SEXP loglikelihood", function_text))
  if(length(f_start) == 0){
    stop("The string: 'SEXP loglikelihood' could not be found in the specified file")
  }
  if(length(f_start) > 1){
    stop("The string: SEXP loglikelihood was found more than once in the specified file")
  }
  
  # Isolate where function ends
  f_end <- which(grepl("loglikelihood_end", function_text))
  if(length(f_end) == 0){
    stop("The string: '// loglikelihood_end' could not be found in the specified file. This
         must be included on the final line to flag the function end")
  }
  
  # Create string of function
  paste(t1[f_start:(f_end - 1)], collapse =  "\n")
}


function_to_string("src/loglikelihood_diarrhoea.cpp")
