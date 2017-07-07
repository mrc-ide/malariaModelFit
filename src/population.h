
#ifndef __malariaModelFit__population__
#define __malariaModelFit__population__

#include <RcppArmadillo.h>
#include "parameters.h"


//------------------------------------------------
// population::
// class of population states
class population {
    
public:
    
    // PUBLIC OBJECTS
    
    double EIR; // annualised EIR
    double ft;
    
    arma::vec age; // age groups in years
    arma::vec dage; // age groups in days
    
    // different states
    arma::vec S;
    arma::vec T;
    arma::vec D;
    arma::vec A;
    arma::vec U;
    arma::vec P;
    arma::vec inf;
    arma::vec prop;
    arma::vec psi;
    arma::vec pos_M;
    arma::vec pos_PCR;
    arma::vec inc;
    
    double FOIM;
    
    // intermediate quantities needed to set equilibrium solution
    int nh;
    int na;
    arma::vec ghnodes;
    arma::vec ghweights;
    
    arma::vec zeta;
    
    arma::mat mICA;
    arma::mat mICM;
    
    arma::mat mS;
    arma::mat mA;
    arma::mat mT;
    arma::mat mD;
    arma::mat mU;
    arma::mat mP;
    
    arma::vec r;
    
    arma::mat mFOI;
    arma::mat mphi;
    arma::mat mq;
    arma::mat mcA;
    arma::mat mpos_M;
    arma::mat mpos_PCR;
    arma::mat minc;
    arma::mat minf;
    
    arma::vec IB;
    arma::vec IC;
    arma::vec ID;
    arma::vec b;
    
    double betaT, betaD, betaP;
    arma::vec betaS;
    arma::vec betaA;
    arma::vec betaU;
    
    arma::vec aT;
    arma::vec bT;
    arma::vec aD;
    arma::vec bD;
    arma::vec aP;
    arma::vec bP;
    
    arma::vec Y;
    
    
    // PUBLIC FUNCTIONS
    
    // constructors
    population(); // default constructor
    population(const arma::vec& age0, const arma::vec& ghn0, const arma::vec& ghw0);
    
    // find equilibrium solution given EIR (EIR0), proportion treated (ft) and model parameters. Translated from Jamie's R code.
    void set_equilibrium(double EIR0, double ft0, const parameters& p);
    // get the annual incidence rate between age0 and age1
    //double get_inc(double age0, double age1);
    // get the prevalence between age0 and age1 by Microscopy or PCR
    //double get_prev_M(double age0, double age1);
    //double get_prev_PCR(double age0, double age1);
    
};

#endif
