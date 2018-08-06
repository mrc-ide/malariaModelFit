
#include "population.h"

using namespace std;

//------------------------------------------------
// initialise a population given age groups (age) and Gauss-Hermit nodes (ghnodes) and
// the associated weights (ghweights)
population::population(const vector<double>& age_, const vector<double>& ghnodes_, const vector<double>& ghweights_) {
  
  // basic inputs
  age = age_;
  na = age.size();
  ghnodes = ghnodes_;
  ghweights = ghweights_;
  nh = ghnodes.size();
  
  // get age in days
  age_days = vector<double>(na);
  for (int i=0; i<na; ++i) {
    age_days[i] = age[i]*365;
  }
  
  // get midpoint of age range
  age_days_midpoint = vector<double>(na);
  for (int i=0; i<na; ++i) {
    age_days_midpoint[i] = (i==(na-1)) ? age_days[i] : (age_days[i]+age_days[i+1])*0.5;
  }
  
  // calculate age group of 20 year old
  age20 = 0;
  for (int i=0; i<na; ++i) {
    if (i<(na-1) && age_days_midpoint[i]<=20*365 && age_days_midpoint[i+1]>20*365) {
      age20 = i;
      break;
    }
  }
  
  // get rate of ageing in each group
  r = vector<double>(na);
  for (int i=0; i<na; ++i) {
    r[i] = (i==(na-1)) ? 0 : 1/double(age_days[i+1]-age_days[i]);
  }
  
  // EIR scaling
  zeta = vector<double>(nh);
  
  // force of infection on mosquitoes
  FOIM = 0;
  
  // model states
  S = vector<double>(na);
  T = vector<double>(na);
  D = vector<double>(na);
  A = vector<double>(na);
  U = vector<double>(na);
  P = vector<double>(na);
  
  prop = vector<double>(na);
  psi = vector<double>(na);
  phi = vector<double>(na);
  inf = vector<double>(na);
  pos_M = vector<double>(na);
  pos_PCR = vector<double>(na);
  inc = vector<double>(na);
  
  ICA = vector<double>(na);
  FOI = vector<double>(na);
  q = vector<double>(na);
  cA = vector<double>(na);
  ICM = vector<double>(na);
  
  // sums over nodes
  S_sum = vector<double>(na);
  T_sum = vector<double>(na);
  D_sum = vector<double>(na);
  A_sum = vector<double>(na);
  U_sum = vector<double>(na);
  P_sum = vector<double>(na);
  
  inf_sum = vector<double>(na);
  pos_M_sum = vector<double>(na);
  pos_PCR_sum = vector<double>(na);
  inc_sum = vector<double>(na);
  
}

//------------------------------------------------
// find equilibrium solution given EIR (EIR0), proportion treated (ft) and model
// parameters.
void population::set_equilibrium(double EIR, double ft, const parameters& p) {
  
  // calculate proportion and relative biting rate in each age group
  for (int i=0; i<na; ++i) {
    prop[i] = (i==0) ? p.eta/(r[i]+p.eta) : r[i-1]*prop[i-1]/(r[i]+p.eta);
    psi[i] = 1 - p.rho*exp(-age_days_midpoint[i]/p.a0);
  }
  
  // calculate EIR scaling factor over Gaussian quadrature nodes
  for (int i=0; i<nh; ++i) {
    zeta[i] = exp(-p.s2*0.5 + sqrt(p.s2)*ghnodes[i]);
  }
  
  // zero vectors for storing sums over Gausian quadrature nodes
  fill(S_sum.begin(), S_sum.end(), 0);
  fill(T_sum.begin(), T_sum.end(), 0);
  fill(D_sum.begin(), D_sum.end(), 0);
  fill(A_sum.begin(), A_sum.end(), 0);
  fill(U_sum.begin(), U_sum.end(), 0);
  fill(P_sum.begin(), P_sum.end(), 0);
  
  fill(inf_sum.begin(), inf_sum.end(), 0);
  fill(pos_M_sum.begin(), pos_M_sum.end(), 0);
  fill(pos_PCR_sum.begin(), pos_PCR_sum.end(), 0);
  fill(inc_sum.begin(), inc_sum.end(), 0);
  
  // loop through Gaussian quadrature nodes
  FOIM = 0;
  for (int i=0; i<nh; ++i) {
    
    // calculate immunity functions and onward infectiousness at equilibrium for 
    // all age groups
    double IB = 0, IC = 0, ID = 0;
    for (int j=0; j<na; ++j) {
      
      // rate of ageing plus death
      double re = r[j] + p.eta;
      
      // update pre-erythrocytic immunity IB
      double eps = zeta[i] * EIR/365 * psi[j];
      IB = (eps/(eps*p.ub+1) + re*IB)/(1/p.db+re);
      
      // calculate probability of infection from pre-erythrocytic immunity IB via
      // Hill function
      double b = p.b0*(p.b1 + (1-p.b1)/(1+pow(IB/p.IB0,p.kb)));
      
      // calculate force of infection (lambda)
      FOI[j] = b*eps;
      
      // update clinical immunity IC
      IC = (FOI[j]/(FOI[j]*p.uc+1) + re*IC)/(1/p.dc + re);
      ICA[j] = IC;
      
      // update detection immunity ID
      ID = (FOI[j]/(FOI[j]*p.ud+1) + re*ID)/(1/p.dd + re);
      
      // calculate probability that an asymptomatic infection (state A) will be
      // detected by microscopy
      double fd = 1 - (1-p.fd0)/(1+pow(age_days_midpoint[j]/p.ad0, p.gd));
      q[j] = p.d1 + (1-p.d1)/(1+pow(ID/p.ID0,p.kd)*fd);
      
      // calculate onward infectiousness to mosquitoes
      cA[j] = p.cU + (p.cD-p.cU)*pow(q[j],p.g_inf);
    }
    
    // calculate maternal clinical immunity, assumed to be at birth a proportion
    // of the acquired immunity of a 20 year old
    double IM0 = ICA[age20]*p.PM;
    for (int j=0; j<na; ++j) {
      
      // rate of ageing plus death
      double re = r[j] + p.eta;
      
      // maternal clinical immunity decays from birth
      double ICM_prev = (j==0) ? IM0 : ICM[j-1];
      ICM[j] = ICM_prev*re/(1/p.dm + re);
    }
    
    // calculate probability of acquiring clinical disease as a function of 
    // different immunity types
    for (int j=0; j<na; ++j) {
      phi[j] = p.phi0*(p.phi1 + (1-p.phi1)/(1+pow((ICA[j]+ICM[j])/p.IC0, p.kc)));
    }
    
    // calculate equilibrium solution of all model states. Again, see references
    // above for details
    for (int j=0; j<na; ++j) {
      
      // rate of ageing plus death
      double re = r[j] + p.eta;
      
      // calculate beta values
      double betaT = p.rT + re;
      double betaD = p.rD + re;
      double betaA = FOI[j]*phi[j] + p.rA + re;
      double betaU = FOI[j] + p.rU + re;
      double betaP = p.rP + re;
      
      // calculate a and b values
      double aT = ft*phi[j]*FOI[j]/betaT;
      double bT = (j==0) ? 0 : r[j-1]*T[j-1]/betaT;
      double aD = (1-ft)*phi[j]*FOI[j]/betaD;
      double bD = (j==0) ? 0 : r[j-1]*D[j-1]/betaD;
      double aP = p.rT*aT/betaP;
      double bP = (p.rT*bT + ((j==0) ? 0 : r[j-1]*P[j-1]))/betaP;
      
      // calculate Y
      double Y = (prop[j] - (bT+bD+bP))/(1+aT+aD+aP);
      
      // calculate final {T,D,P} solution
      T[j] = aT*Y+bT;
      D[j] = aD*Y+bD;
      P[j] = aP*Y+bP;
                
      // calculate final {A, U, S} solution
      double rA = 0, rU = 0;
      if (j>0) {
        rA = r[j-1]*A[j-1];
        rU = r[j-1]*U[j-1];
      }
      A[j] = (rA + (1-phi[j])*Y*FOI[j] + p.rD*D[j])/(betaA + (1-phi[j])*FOI[j]);
      U[j] = (rU + p.rA*A[j])/betaU;
      S[j] = Y - A[j] - U[j];
                
      // calculate proportion detectable by mocroscopy and PCR
      pos_M[j] = D[j] + T[j] + A[j]*q[j];
      pos_PCR[j] = D[j] + T[j] + A[j]*pow(q[j], p.aA) + U[j]*pow(q[j], p.aU);
      
      // calculate clinical incidence
      inc[j] = Y*FOI[j]*phi[j];
      
      // calculate incidence of infection
      inf[j] = p.cD*D[j] + p.cT*T[j] + cA[j]*A[j] + p.cU*U[j];
    }
    
    // add to sums over nodes
    double delta_FOIM = 0;
    for (int j=0; j<na; ++j) {
      S_sum[j] += ghweights[i]*S[j];
      T_sum[j] += ghweights[i]*T[j];
      D_sum[j] += ghweights[i]*D[j];
      A_sum[j] += ghweights[i]*A[j];
      U_sum[j] += ghweights[i]*U[j];
      P_sum[j] += ghweights[i]*P[j];
      
      inf_sum[j] += ghweights[i]*inf[j];
      pos_M_sum[j] += ghweights[i]*pos_M[j];
      pos_PCR_sum[j] += ghweights[i]*pos_PCR[j];
      inc_sum[j] += ghweights[i]*inc[j];
      
      delta_FOIM += inf[j]*psi[j];
    }
    FOIM += delta_FOIM*ghweights[i]*zeta[i];
    
  }  // end loop over nodes
  
  // complete overall force of infection on mosquitoes
  FOIM *= p.f*p.Q0/(1 - p.rho*p.eta/(p.eta+1/p.a0));
  
}
