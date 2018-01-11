
#include "main.h"
#include "misc.h"
#include "MCMC.h"
#include "parameters.h"
#include "population.h"

using namespace std;

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model. Equivalent to the R-version function "human_equilibrium" by Jamie

// [[Rcpp::export]]
Rcpp::List human_equilibrium_cpp(double EIR, double ft, const Rcpp::List& p, arma::vec& age0,
           const arma::vec& ghnodes, const arma::vec& ghweights) {
    
    // initialise parameters and population
    parameters par(p);
    population pop(age0, ghnodes, ghweights);
    
    // set to equilibrium
    pop.set_equilibrium(EIR, ft, par);
    
    // process output
    Rcpp::NumericMatrix E = Rcpp::cbind(tonv(age0), tonv(pop.S), tonv(pop.T), tonv(pop.D), tonv(pop.A), tonv(pop.U), tonv(pop.P), tonv(pop.inf), tonv(pop.prop), tonv(pop.psi), tonv(pop.pos_M), tonv(pop.pos_PCR), tonv(pop.inc));
    colnames(E) = Rcpp::CharacterVector::create("age","S","T","D","A","U","P","inf","prop","psi","pos_M","pos_PCR","inc");
    Rcpp::List result = Rcpp::List::create(Rcpp::_["states"]=E, Rcpp::_["FOIM"]=pop.FOIM);
    
    // return
    return(result);
}

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
