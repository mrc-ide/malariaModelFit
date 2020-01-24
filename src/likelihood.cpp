#include <Rcpp.h>
// [[Rcpp::export]]
SEXP loglikelihood(std::vector<double> params, std::vector<double> x){
  
  // Unpack data ///////////////////////////////////////////////////////////////
    // Data input vector format:
      // Output type
      // Site number
      // Vector of site group number
        // List of age prop for each site
      // List of age r for each site
      // List of age age days midpoint for each site
      // List of age psi for each site
      // Vector of age 20 year old for each site
      // Number of heterogeneity classes
      // Vector of Gaussian quadrature nodes
      // Vector of Gaussian quadrature weights
      
      // List of site data denominators (Number of persons or person years)
      // List of site data numerators (Number +ve or cases)

    //Rcpp::Rcout << "Unpacking data" << std::endl;
    // Index of data vector
    int di = 0;
    
    // output_type
    int output_type = x[di];
    di++;
    // Number of sites
    int site_n = x[di];
    di++;
    // Site group number
    std::vector<int> ng(site_n);
    for(int i = 0; i < site_n; ++i){
      ng[i] = x[di];
      di++;
    }

    // Empty template site_n X group_n
    std::vector<std::vector<double> > template_dbl(site_n);
    for(int i = 0; i < site_n; ++i){
      template_dbl[i].resize(ng[i]);
    }

    // Age prop for each site
    std::vector<std::vector<double> > prop = template_dbl;
    for(int i = 0; i < site_n; ++i){
      for(int j = 0; j < ng[i]; ++j){
        prop[i][j] = x[di];
        di++;
      } 
    }
    // Age r for each site
    std::vector<std::vector<double> > r = template_dbl;
    for(int i = 0; i < site_n; ++i){
      for(int j = 0; j < ng[i]; ++j){
        r[i][j] = x[di];
        di++;
      } 
    }
    // Age day midpoint for each site
    std::vector<std::vector<double> > age_days_midpoint = template_dbl;
    for(int i = 0; i < site_n; ++i){
      for(int j = 0; j < ng[i]; ++j){
        age_days_midpoint[i][j] = x[di];
        di++;
      } 
    }
    // Age psi for each site
    std::vector<std::vector<double> > psi = template_dbl;
    for(int i = 0; i < site_n; ++i){
      for(int j = 0; j < ng[i]; ++j){
        psi[i][j] = x[di];
        di++;
      } 
    }
    // Index of 20 year old age for each site
    std::vector<int> age20(site_n);
    for(int i = 0; i < site_n; ++i){
      age20[i] = x[di];
      di++;
    }
    // Gaussian quadrature nodes and weights
    int nh = x[di];
    di++;
    std::vector<double> nodes(nh);
    std::vector<double> weights(nh);
    for(int i = 0; i < nh; ++i){
      nodes[i] = x[di];
      di++;
    }
    for(int i  = 0; i < nh; ++i){
      weights[i] = x[di];
      di++;
    }
    // Data numerators
    std::vector<std::vector<double> > numer = template_dbl;
    for(int i = 0; i < site_n; ++i){
      for(int j = 0; j < ng[i]; ++j){
        numer[i][j] = x[di];
        di++;
      } 
    }
    // Data denominators
    std::vector<std::vector<double> > denom = template_dbl;
    for(int i = 0; i < site_n; ++i){
     for(int j = 0; j < ng[i]; ++j){
       denom[i][j] = x[di];
       di++;
     } 
    }
    // Prevalence or incidence
    std::vector<int> type(site_n);
    for(int i = 0; i < site_n; ++i){
      type[i] = x[di];
      di++;
    }
    // Case detection
    std::vector<int> case_detection(site_n);
    for(int i = 0; i < site_n; ++i){
      case_detection[i] = x[di];
      di++;
    }
  //////////////////////////////////////////////////////////////////////////////
  
  // Unpack parameters /////////////////////////////////////////////////////////
     //Rcpp::Rcout << "Unpacking parameters" << std::endl;
    // Index of parameter vectors
    int pi = 0;
    // age, heterogeneity in exposure
    double eta = params[pi];
    pi++;
    double rho = params[pi];
    pi++;
    double a0 = params[pi];
    pi++;
    double s2 = params[pi];
    pi++;
    
    // rate of leaving infection states
    double prA = params[pi];
    pi++;
    double prT = params[pi];
    pi++;
    double prD = params[pi];
    pi++;
    double prU = params[pi];
    pi++;
    double prP = params[pi];
    pi++;
    
    // human latent period and time lag from asexual parasites to infectiousness
    double dE = params[pi];
    pi++;
    double tl = params[pi];
    pi++;
    
    // infectiousness to mosquitoes
    double cD = params[pi];
    pi++;
    double cT = params[pi];
    pi++;
    double cU = params[pi];
    pi++;
    double g_inf  = params[pi];
    pi++;
    
    // anti-parasite immunity
    double d1 = params[pi];
    pi++;
    double dd = params[pi];
    pi++;
    double ID0 = params[pi];
    pi++;
    double kd = params[pi];
    pi++;
    double ud = params[pi];
    pi++;
    double ad0 = params[pi];
    pi++;
    double fd0 = params[pi];
    pi++;
    double gd = params[pi];
    pi++;
    double aA = params[pi];
    pi++;
    double aU = params[pi];
    pi++;
    
    // anti-infection immunity
    double b0 = params[pi];
    pi++;
    double b1 = params[pi];
    pi++;
    double db = params[pi];
    pi++;
    double IB0 = params[pi];
    pi++;
    double kb = params[pi];
    pi++;
    double ub = params[pi];
    pi++;
    
    // clinical immunity
    double phi0 = params[pi];
    pi++;
    double phi1 = params[pi];
    pi++;
    double dc = params[pi];
    pi++;
    double IC0 = params[pi];
    pi++;
    double kc = params[pi];
    pi++;
    double uc = params[pi];
    pi++;
    double PM = params[pi];
    pi++;
    double dm = params[pi];
    pi++;
    
    // mosquito parameters
    double tau = params[pi];
    pi++;
    double mu = params[pi];
    pi++;
    double f = params[pi];
    pi++;
    double Q0 = params[pi];
    pi++;
    
    // EIR
    std::vector<double> EIR(site_n);
    for(int i = 0; i < site_n; ++i){
      EIR[i] = params[pi];
      pi++;
    }
    // Treatment coverage
    std::vector<double> ft(site_n);
    for(int i = 0; i < site_n; ++i){
      ft[i] = params[pi];
      pi++;
    }
    // Hyper parameters
    std::vector<double> alpha_hyper(site_n);
    for(int i = 0; i < site_n; ++i){
      alpha_hyper[i] = params[pi];
      pi++;
    }
    std::vector<double> beta_hyper(site_n);
    for(int i = 0; i < site_n; ++i){
      beta_hyper[i] = params[pi];
      pi++;
    }
  //////////////////////////////////////////////////////////////////////////////
  
  // Initialise output variables for all sites /////////////////////////////////
    std::vector<std::vector<double> > S_out = template_dbl;
    std::vector<std::vector<double> > T_out = template_dbl;
    std::vector<std::vector<double> > D_out = template_dbl;
    std::vector<std::vector<double> > A_out = template_dbl;
    std::vector<std::vector<double> > U_out = template_dbl;
    std::vector<std::vector<double> > P_out = template_dbl;
    std::vector<std::vector<double> > pos_M_out = template_dbl;
    std::vector<std::vector<double> > pos_PCR_out = template_dbl;
    std::vector<std::vector<double> > inf_out = template_dbl;
    std::vector<std::vector<double> > inc_out = template_dbl;
    std::vector<double> FOIM(site_n);
  //////////////////////////////////////////////////////////////////////////////
  
  for(int s = 0; s < site_n; ++s){
  // Initialise output variables for a single site /////////////////////////////
    std::vector<std::vector<double> > pos_M(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > pos_PCR(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > inc(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > inf(nh, std::vector<double>(ng[s], 0)); 
    std::vector<std::vector<double> > S(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > T(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > P(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > D(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > A(nh, std::vector<double>(ng[s], 0));
    std::vector<std::vector<double> > U(nh, std::vector<double>(ng[s], 0));
    double zeta;
    double rA;
    double rU;
  //////////////////////////////////////////////////////////////////////////////
  
  // Equilbrium ////////////////////////////////////////////////////////////////
     //Rcpp::Rcout << "Running EQ" << std::endl;
    // loop through all Gaussian quadrature nodes
    for(int n = 0; n < nh; ++n){
      zeta = exp(-s2 * 0.5 + std::sqrt(s2) * nodes[n]);
      double EIR_cur = EIR[s] / 365 * zeta;
      // Human EQ no-het:
      double IB = 0;
      double IC = 0;
      double ID = 0;
      std::vector<double> ICA(ng[s], 0);
      std::vector<double> FOI(ng[s], 0);
      std::vector<double> q(ng[s], 0);
      std::vector<double> cA(ng[s], 0);
      double re;
      double eps;
      double b;
      double fd;
      double IM0;
      std::vector<double> ICM(ng[s], 0);
      double ICM_prev;
      std::vector<double> phi(ng[s], 0);
      
       //Rcpp::Rcout << "Starting NA" << std::endl;
      for(int a = 0; a < ng[s]; ++a){
        re = r[s][a] + eta;
        eps = EIR_cur * psi[s][a];
        IB = (eps / (eps * ub + 1) + re * IB) / (1 / db + re);
        b = b0 * (b1 + (1 - b1) / (1 + std::pow((IB / IB0), kb)));
        FOI[a] = b * eps;
        IC = (FOI[a] / (FOI[a] * uc + 1) + re * IC) / (1 / dc + re);
        ICA[a] = IC;
        ID = (FOI[a] / (FOI[a] * ud + 1) + re * ID) / (1 / dd + re);
        fd = 1 - (1 - fd0) / (1 + std::pow((age_days_midpoint[s][a] / ad0), gd));
        q[a] = d1 + (1 - d1) / (1 + std::pow((ID / ID0), kd) * fd);
        cA[a] = cU + (cD - cU) * std::pow(q[a], g_inf);
      }
      
      // Rcpp::Rcout << "Starting prev1" << std::endl;
      IM0 = ICA[age20[s]] * PM;
      for(int a = 0; a < ng[s]; ++a){
        re = r[s][a] + eta;
        if (a == 0) {
          ICM_prev = IM0;
        } else {
          ICM_prev = ICM[a - 1];
        }
        ICM[a] = ICM_prev * re / (1 / dm + re);
      }
      
      // Rcpp::Rcout << "Starting phi" << std::endl;
      for(int a = 0; a < ng[s]; ++a){
        phi[a] = phi0 * (phi1 + (1 - phi1) / (1 + (std::pow((ICA[a] + ICM[a]) / IC0, kc))));
      }

      double betaT;
      double betaD;
      double betaA;
      double betaU;
      double betaP;
      double aT;
      double aP;
      double aD;
      double bT;
      double bP;
      double bD;
      double Y;
      
      // Rcpp::Rcout << "Starting states" << std::endl;
      for(int a = 0; a < ng[s]; ++a){
        re = r[s][a] + eta;
        betaT = prT + re;
        betaD = prD + re;
        betaA = FOI[a]*phi[a] + prA + re;
        betaU = FOI[a] + prU + re;
        betaP = prP + re;
        aT = ft[s] * phi[a] * FOI[a] / betaT;
        aP = prT * aT / betaP;
        aD = (1 - ft[s]) * phi[a] * FOI[a] / betaD;
        if (a == 0){
          bT = 0;
          bD = 0;
          bP = 0;
        } else {
          bT = r[s][a - 1] * T[n][a - 1] / betaT;
          bD = r[s][a - 1] * D[n][a - 1] / betaD;
          bP = prT * bT + r[s][a - 1] * P[n][a - 1] / betaP;
        }
        Y = (prop[s][a] - (bT + bD + bP)) / (1 + aT + aD + aP);

        T[n][a] = aT * Y + bT;
        D[n][a] = aD * Y + bD;
        P[n][a] = aP * Y + bP;
        if (a == 0){
          rA = 0;
          rU = 0;
        } else {
          rA = r[s][a - 1] * A[n][a - 1];
          rU = r[s][a - 1] * U[n][a - 1];
        }
        A[n][a] = (rA + (1 - phi[a]) * Y * FOI[a] + prD * D[n][a]) / (betaA + (1 - phi[a]) * FOI[a]);
        U[n][a] = (rU + prA * A[n][a]) / betaU;
        S[n][a] = Y - A[n][a] - U[n][a];
        pos_M[n][a] = D[n][a] + T[n][a] + A[n][a] * q[a];
        pos_PCR[n][a] = D[n][a] + T[n][a] + A[n][a] * std::pow(q[a], aA) + U[n][a] * std::pow(q[a], aU);
        inc[n][a] = Y * FOI[a] * phi[a];
      }
      
      for(int a = 0; a < ng[s]; ++a){
        inf[n][a] = cD * D[n][a] + cT * T[n][a] + cA[a] * A[n][a] + cU * U[n][a];
      }
    }
    
    // Average over nodes
    for(int n = 0; n < nh; ++n){
      zeta = exp(-s2 * 0.5 + std::sqrt(s2) * nodes[n]);
      for(int a = 0; a < ng[s]; ++a){
        S_out[s][a] += S[n][a] * weights[n];
        T_out[s][a] += T[n][a] * weights[n];
        D_out[s][a] += D[n][a] * weights[n];
        A_out[s][a] += A[n][a] * weights[n];
        U_out[s][a] += U[n][a] * weights[n];
        P_out[s][a] += P[n][a] * weights[n];
        pos_M_out[s][a] += pos_M[n][a] * weights[n];
        pos_PCR_out[s][a] += pos_PCR[n][a] * weights[n];
        inf_out[s][a] += inf[n][a] * weights[n];
        inc_out[s][a] += inc[n][a] * weights[n];
        FOIM[s] += inf[n][a] * psi[s][a] * weights[n] * zeta;
      }
    }
    
    double omega = 1 - rho * eta / (eta + 1 / a0);
    double alpha = f * Q0;
    FOIM[s] *= alpha/omega;
  }
  //Rcpp::Rcout << "Outer" << std::endl;
    if(output_type == 1){
       return Rcpp::List::create(
                         Rcpp::Named("S") = S_out,
                         Rcpp::Named("T") = T_out,
                         Rcpp::Named("D") = D_out,
                         Rcpp::Named("A") = A_out,
                         Rcpp::Named("U") = U_out,
                         Rcpp::Named("P") = P_out,
                         Rcpp::Named("inf") = inf_out,
                         Rcpp::Named("prop") = prop,
                         Rcpp::Named("psi") = psi,
                         Rcpp::Named("pos_M") = pos_M_out,
                         Rcpp::Named("pos_PCR") = pos_PCR_out,
                         Rcpp::Named("inc") = inc_out,
                         Rcpp::Named("FOIM") = FOIM);
    }
  //////////////////////////////////////////////////////////////////////////////
 
  // Likelihood ////////////////////////////////////////////////////////////////
  //Rcpp::Rcout << "Likelihood" << std::endl;
  double lL = 0;
  // For each site
  for(int s = 0; s < site_n; ++s){
    // For each age group
    for(int a = 0; a < ng[s]; ++a){
      switch(type[s]) {
      // Incidence calculation
      case 1:
        lL += 0;
        break;
      // Prevalence calculation
      case 2:
        lL += 0;
        break;
      default:
        break;
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  
  return Rcpp::wrap(lL);
}
// loglikelihood_end"