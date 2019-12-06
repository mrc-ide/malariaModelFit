
#pragma once

#include "parameters.h"
#include "misc_v9.h"

//------------------------------------------------
// class of population states
class population {
  
public:
  // PUBLIC OBJECTS
  
  // basic inputs
  std::vector<double> age;
  std::vector<double> ghnodes;
  std::vector<double> ghweights;
  int na;
  int nh;
  int age20;
  
  // ageing
  std::vector<double> age_days;
  std::vector<double> age_days_midpoint;
  std::vector<double> r;
  
  // EIR scaling
  std::vector<double> zeta;
  
  // force of infection on mosquitoes
  double FOIM;
  
  // model states
  std::vector<double> S;
  std::vector<double> T;
  std::vector<double> D;
  std::vector<double> A;
  std::vector<double> U;
  std::vector<double> P;
  
  std::vector<double> prop;
  std::vector<double> psi;
  std::vector<double> phi;
  std::vector<double> inf;
  std::vector<double> pos_M;
  std::vector<double> pos_PCR;
  std::vector<double> inc;
  
  std::vector<double> ICA;
  std::vector<double> FOI;
  std::vector<double> q;
  std::vector<double> cA;
  std::vector<double> ICM;
  
  // sums over nodes
  std::vector<double> S_sum;
  std::vector<double> T_sum;
  std::vector<double> D_sum;
  std::vector<double> A_sum;
  std::vector<double> U_sum;
  std::vector<double> P_sum;
  
  std::vector<double> inf_sum;
  std::vector<double> pos_M_sum;
  std::vector<double> pos_PCR_sum;
  std::vector<double> inc_sum;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  population() {};
  population(const std::vector<double>& age_, const std::vector<double>& ghnodes_, const std::vector<double>& ghweights_);
  
  // find equilibrium solution given EIR (EIR), proportion treated (ft) and
  // model parameters (p)
  void set_equilibrium(double EIR, double ft, const parameters& p);
  
  // get the incidence or prevalence between age0 and age1
  double get_inc(double age0, double age1);
  double get_prev_M(double age0, double age1);
  double get_prev_PCR(double age0, double age1);
    
};
