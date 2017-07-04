
#include "main.h"
#include "misc.h"
#include "parameters.h"
#include "population.h"

using namespace std;

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model. Equivalent to the R-version function "human_equilibrium" by Jamie
// [[Rcpp::export]]
Rcpp::List human_equilibrium_fast(double EIR, double ft, const Rcpp::List& p, const arma::vec& age0,
           const arma::vec& ghnodes, const arma::vec& ghweights) {
    
    // initialise parameters and population
    parameters par(p);
    population pop(age0, ghnodes, ghweights);
    
    // set to equilibrium
    pop.set_equilibrium(EIR, ft, par);
    
    // process output
    Rcpp::NumericMatrix E = Rcpp::cbind(tonv(pop.S), tonv(pop.T), tonv(pop.D), tonv(pop.A), tonv(pop.U), tonv(pop.P), tonv(pop.inf),
                            tonv(pop.prop), tonv(pop.psi), tonv(pop.pos_M), tonv(pop.pos_PCR), tonv(pop.inc));
    colnames(E) = Rcpp::CharacterVector::create("S","T","D","A","U","P","inf","prop","psi","pos_M","pos_PCR","inc");
    Rcpp::List result = Rcpp::List::create(Rcpp::_["states"]=E, Rcpp::_["FOIM"]=pop.FOIM);
    
    // return
    return(result);
}

//------------------------------------------------
// convert arma::vec object to NumericVector object
Rcpp::NumericVector tonv(arma::vec& v1) {
    Rcpp::NumericVector v(v1.begin(),v1.end());
    return(v);
}
