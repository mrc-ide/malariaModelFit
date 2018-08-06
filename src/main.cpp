
#include "main.h"
#include "misc.h"
#include "MCMC.h"
#include "parameters.h"
#include "population.h"

using namespace std;

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model.
// Equivalent to the R-version "human_equilibrium()"
// [[Rcpp::export]]
Rcpp::List human_equilibrium_cpp(double EIR, double ft, const Rcpp::List& p, const Rcpp::NumericVector& age_,
           const Rcpp::List& h) {
  
  // process inputs
  vector<double> ghnodes = rcpp_to_vector_double(h["nodes"]);
  vector<double> ghweights = rcpp_to_vector_double(h["weights"]);
  vector<double> age = rcpp_to_vector_double(age_);
  int na = age.size();
  
  // initialise parameters and population
  parameters par(p);
  population pop(age, ghnodes, ghweights);
  
  // set to equilibrium
  pop.set_equilibrium(EIR, ft, par);
  
  // process output
  Rcpp::NumericMatrix states(na,13);
  for (int i=0; i<na; i++) {
    states(i,0) = age[i];
    states(i,1) = pop.S_sum[i];
    states(i,2) = pop.T_sum[i];
    states(i,3) = pop.D_sum[i];
    states(i,4) = pop.A_sum[i];
    states(i,5) = pop.U_sum[i];
    states(i,6) = pop.P_sum[i];
    states(i,7) = pop.inf_sum[i];
    states(i,8) = pop.prop[i];
    states(i,9) = pop.psi[i];
    states(i,10) = pop.pos_M_sum[i];
    states(i,11) = pop.pos_PCR_sum[i];
    states(i,12) = pop.inc_sum[i];
  }
  
  Rcpp::List result = Rcpp::List::create(Rcpp::_["states"]=states, Rcpp::_["FOIM"]=pop.FOIM);
  return result;
}
/*
//------------------------------------------------
// run the Metropolis MCMC fitting for nrep iterations
// 1. The parameters are updated in blocks as specified by update_blocks
// 2. update_blocks do not need to contain all the parameter indices, i.e. partial fitting is okay
// 3. Proposal distribution for each block is tuned adaptively to achieve reasonable acceptance ratio

// [[Rcpp::export]]
Rcpp::List runMCMC(int nrep, Rcpp::DataFrame& data_key, Rcpp::List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0, const Rcpp::NumericVector& ghnodes1, const Rcpp::NumericVector& ghweights1, const std::vector<Rcpp::IntegerVector>& update_blocks) {
    
    // read in the data
    all_data alldata(data_key, datasets, ghnodes1, ghweights1);
    
    return Rcpp::List::create(Rcpp::Named("foobar")=-9);
    
}
*/