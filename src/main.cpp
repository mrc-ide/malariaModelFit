
#include "main.h"
#include "parameters.h"
#include "population.h"
#include "probability_v10.h"

using namespace std;

//------------------------------------------------
// Exportable function to find equilibrium solution of the transmission model.
// Equivalent to the R-version "human_equilibrium()"
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
