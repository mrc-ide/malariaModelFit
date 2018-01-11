
#pragma once

#include <RcppArmadillo.h>

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model. Equivalent to the R-version function "human_equilibrium" by Jamie
Rcpp::List human_equilibrium_cpp(double EIR, double ft, const Rcpp::List& p, const arma::vec& age0, const arma::vec& ghnodes, const arma::vec& ghweights);

//------------------------------------------------
// run the Metropolis MCMC fitting for nrep iterations
Rcpp::List runMCMC(int nrep, Rcpp::DataFrame& data_key, Rcpp::List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0, const Rcpp::NumericVector& ghnodes1, const Rcpp::NumericVector& ghweights1, const std::vector<Rcpp::IntegerVector>& update_blocks);


