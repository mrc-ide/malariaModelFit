#include <Rcpp.h>
// [[Rcpp::export]]
SEXP loglikelihood(std::vector<double> params, std::vector<double> x){
  
  // Unpack data ///////////////////////////////////////////////////////////////
    // Rcpp::Rcout << "Unpacking data" << std::endl;
    // Index of data vector
    int di = 0;
    
    // output_type
    int output_type = x[di];
    di++;
    // Number of studies
    int study_n = x[di];
    di++;
    // EIR
    double EIR = x[di];
    di++;
    double ft = x[di];
    di++;
    // Age
    int na = x[di];
    di++;
    std::vector<double> age(na);
    for(int i = 0; i < na; ++i){
      age[i] = x[di];
      di++;
    }
    std::vector<double> prop(na);
    for(int i = 0; i < na; ++i){
      prop[i] = x[di];
      di++;
    }
    std::vector<double> r(na);
    for(int i = 0; i < na; ++i){
      r[i] = x[di];
      di++;
    }
    std::vector<double> age_days_midpoint(na);
    for(int i = 0; i < na; ++i){
      age_days_midpoint[i] = x[di];
      di++;
    }
    std::vector<double> psi(na);
    for(int i = 0; i < na; ++i){
      psi[i] = x[di];
      di++;
    }
    int age20 = x[di++];
    // Gaussian quadrature nodes
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
  
    // ETC
  //////////////////////////////////////////////////////////////////////////////
  
  // Unpack parameters /////////////////////////////////////////////////////////
    // Rcpp::Rcout << "Unpacking parameters" << std::endl;
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
  //////////////////////////////////////////////////////////////////////////////
  
  // Initialise output variables ///////////////////////////////////////////////
    std::vector<std::vector<double>> pos_M(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> pos_PCR(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> inc(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> inf(nh, std::vector<double>(na, 0)); 
    std::vector<std::vector<double>> S(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> T(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> P(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> D(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> A(nh, std::vector<double>(na, 0));
    std::vector<std::vector<double>> U(nh, std::vector<double>(na, 0));
    double zeta;
    double rA;
    double rU;
  //////////////////////////////////////////////////////////////////////////////
  
  // Equilbrium ////////////////////////////////////////////////////////////////
    // Rcpp::Rcout << "Running EQ" << std::endl;
    // loop through all Gaussian quadrature nodes
    for(int n = 0; n < nh; ++n){
      zeta = exp(-s2 * 0.5 + std::sqrt(s2) * nodes[n]);
      double EIR_cur = EIR / 365 * zeta;
      // Human EQ no-het:
      double IB = 0;
      double IC = 0;
      double ID = 0;
      std::vector<double> ICA(na, 0);
      std::vector<double> FOI(na, 0);
      std::vector<double> q(na, 0);
      std::vector<double> cA(na, 0);
      double re;
      double eps;
      double b;
      double fd;
      double IM0;
      std::vector<double> ICM(na, 0);
      double ICM_prev;
      std::vector<double> phi(na, 0);
      
      // Rcpp::Rcout << "Starting NA" << std::endl;
      for(int a = 0; a < na; ++a){
        re = r[a] + eta;
        eps = EIR_cur * psi[a];
        IB = (eps / (eps * ub + 1) + re * IB) / (1 / db + re);
        b = b0 * (b1 + (1 - b1) / (1 + std::pow((IB / IB0), kb)));
        FOI[a] = b * eps;
        IC = (FOI[a] / (FOI[a] * uc + 1) + re * IC) / (1 / dc + re);
        ICA[a] = IC;
        ID = (FOI[a] / (FOI[a] * ud + 1) + re * ID) / (1 / dd + re);
        fd = 1 - (1 - fd0) / (1 + std::pow((age_days_midpoint[a] / ad0), gd));
        q[a] = d1 + (1 - d1) / (1 + std::pow((ID / ID0), kd) * fd);
        cA[a] = cU + (cD - cU) * std::pow(q[a], g_inf);
      }
      
      // Rcpp::Rcout << "Starting prev1" << std::endl;
      IM0 = ICA[age20] * PM;
      for(int a = 0; a < na; ++a){
        re = r[a] + eta;
        if (a == 0) {
          ICM_prev = IM0;
        } else {
          ICM_prev = ICM[a - 1];
        }
        ICM[a] = ICM_prev * re / (1 / dm + re);
      }
      
      // Rcpp::Rcout << "Starting phi" << std::endl;
      for(int a = 0; a < na; ++a){
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
      for(int a = 0; a < na; ++a){
        re = r[a] + eta;
        betaT = prT + re;
        betaD = prD + re;
        betaA = FOI[a]*phi[a] + prA + re;
        if(a == 3) Rcpp::Rcout << "::: " << rA <<std::endl;
        betaU = FOI[a] + prU + re;
        betaP = prP + re;
        aT = ft * phi[a] * FOI[a] / betaT;
        aP = prT * aT / betaP;
        aD = (1 - ft) * phi[a] * FOI[a] / betaD;
        if (a == 0){
          bT = 0;
          bD = 0;
          bP = 0;
        } else {
          bT = r[a - 1] * T[n][a - 1] / betaT;
          bD = r[a - 1] * D[n][a - 1] / betaD;
          bP = prT * bT + r[a - 1] * P[n][a - 1] / betaP;
        }
        Y = (prop[a] - (bT + bD + bP)) / (1 + aT + aD + aP);

        T[n][a] = aT * Y + bT;
        D[n][a] = aD * Y + bD;
        P[n][a] = aP * Y + bP;
        if (a == 0){
          rA = 0;
          rU = 0;
        } else {
          rA = r[a - 1] * A[n][a - 1];
          rU = r[a - 1] * U[n][a - 1];
        }
        A[n][a] = (rA + (1 - phi[a]) * Y * FOI[a] + prD * D[n][a]) / (betaA + (1 - phi[a]) * FOI[a]);
        U[n][a] = (rU + prA * A[n][a]) / betaU;
        S[n][a] = Y - A[n][a] - U[n][a];
        pos_M[n][a] = D[n][a] + T[n][a] + A[n][a] * q[a];
        pos_PCR[n][a] = D[n][a] + T[n][a] + A[n][a] * std::pow(q[a], aA) + U[n][a] * std::pow(q[a], aU);
        inc[n][a] = Y * FOI[a] * phi[a];
      }
      
      for(int a = 0; a < na; ++a){
        inf[n][a] = cD * D[n][a] + cT * T[n][a] + cA[a] * A[n][a] + cU * U[n][a];
      }
    }
    
    // Average over nodes
    std::vector<double> S_out(na, 0);
    std::vector<double> T_out(na, 0);
    std::vector<double> D_out(na, 0);
    std::vector<double> A_out(na, 0);
    std::vector<double> U_out(na, 0);
    std::vector<double> P_out(na, 0);
    std::vector<double> pos_M_out(na, 0);
    std::vector<double> pos_PCR_out(na, 0);
    std::vector<double> inf_out(na, 0);
    std::vector<double> inc_out(na, 0);
    double FOIM = 0;
    
    for(int n = 0; n < nh; ++n){
      zeta = exp(-s2 * 0.5 + std::sqrt(s2) * nodes[n]);
      for(int a = 0; a < na; ++a){
        S_out[a] += S[n][a] * weights[n];
        T_out[a] += T[n][a] * weights[n];
        D_out[a] += D[n][a] * weights[n];
        A_out[a] += A[n][a] * weights[n];
        U_out[a] += U[n][a] * weights[n];
        P_out[a] += P[n][a] * weights[n];
        pos_M_out[a] += pos_M[n][a] * weights[n];
        pos_PCR_out[a] += pos_PCR[n][a] * weights[n];
        inf_out[a] += inf[n][a] * weights[n];
        inc_out[a] += inc[n][a] * weights[n];
        FOIM += inf[n][a] * psi[a] * weights[n] * zeta;
      }
    }
    
    double omega = 1 - rho * eta / (eta + 1 / a0);
    double alpha = f * Q0;
    FOIM *= alpha/omega;
    
    if(output_type == 1){
       return Rcpp::List::create(Rcpp::Named("age") = age,
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
  
  //////////////////////////////////////////////////////////////////////////////
  return 0;
}
// loglikelihood_end"