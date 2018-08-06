
#pragma once

#include <Rcpp.h>
#include "misc.h"

//------------------------------------------------
// class containing all model parameters
class parameters {
  
public:
  // PUBLIC OBJECTS
  
  // age, heterogeneity in exposure
  double eta;
  double rho;
  double a0;
  double s2;
  
  // rate of leaving infection states
  double rA;
  double rT;
  double rD;
  double rU;
  double rP;
  
  // human latent period and time lag from asexual parasites to infectiousness
  double dE;
  double tl;
  
  // infectiousness to mosquitoes
  double cD;
  double cT;
  double cU;
  double g_inf;
  
  // anti-parasite immunity
  double d1;
  double dd;
  double ID0;
  double kd;
  double ud;
  double ad0;
  double fd0;
  double gd;
  double aA;
  double aU;
  
  // anti-infection immunity
  double b0;
  double b1;
  double db;
  double IB0;
  double kb;
  double ub;
  
  // clinical immunity
  double phi0;
  double phi1;
  double dc;
  double IC0;
  double kc;
  double uc;
  double PM;
  double dm;
  
  // mosquito parameters
  double tau;
  double mu;
  double f;
  double Q0;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  parameters(); // default constructor
  parameters(const Rcpp::List& p); // import from an R list
  
};

