
#ifndef __malariaModelFit__main__
#define __malariaModelFit__main__

#include <RcppArmadillo.h>

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model. Equivalent to the R-version function "human_equilibrium" by Jamie
Rcpp::List human_equilibrium_fast(double EIR, double ft, const Rcpp::List& p, const arma::vec& age0, const arma::vec& ghnodes, const arma::vec& ghweights);

//------------------------------------------------
// convert arma::vec object to NumericVector object
Rcpp::NumericVector tonv(arma::vec& v1);

#endif
