
#pragma once

#include "misc_v9.h"

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List human_equilibrium_cpp(double EIR, double ft, const Rcpp::List& p, const Rcpp::NumericVector& age_,
                                 const Rcpp::List& h);
