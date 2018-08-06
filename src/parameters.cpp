
#include "parameters.h"

using namespace std;

//------------------------------------------------
// default constructor
parameters::parameters() {
  
  // default constructor with values from Griffin et al. 2014 - SI
  eta = 0.0001305;
  rho = 0.85;
  a0 = 2920;
  s2 = 1.67;
  rA = 0.00512821;
  rT = 0.2;
  rD = 0.2;
  rU = 0.00906627;
  rP = 0.2;
  dE = 12;
  tl = 12.5;
  cD = 0.0676909;
  cT = 0.0034482;
  cU = 0.006203;
  g_inf =	1.82425;
  d1 = 0.160527;
  dd = 3650;
  ID0	= 1.577533;
  kd = 0.476614;
  ud = 9.44512;
  ad0 = 8001.99;
  fd0 = 0.007055;
  gd = 4.8183;
  aA = 0.757;
  aU = 0.186;
  b0 = 0.590076;
  b1 = 0.5;
  db = 3650;
  IB0 = 43.8787;
  kb = 2.15506;
  ub = 7.19919;
  phi0 = 0.791666;
  phi1 = 0.000737;
  dc = 10950;
  IC0	= 18.02366;
  kc = 2.36949;
  uc = 6.06349;
  PM = 0.774368;
  dm = 67.6952;
  tau	= 10;
  mu = 0.132;
  f = 0.33333333;
  Q0 = 0.92;
}

//------------------------------------------------
// informed constructor, taking values from an R list
parameters::parameters(const Rcpp::List& p) {
  
  // take parameter values from an R list
  eta = rcpp_to_double(p["eta"]);
  rho = rcpp_to_double(p["rho"]);
  a0 = rcpp_to_double(p["a0"]);
  s2 = rcpp_to_double(p["s2"]);
  rA = rcpp_to_double(p["rA"]);
  rT = rcpp_to_double(p["rT"]);
  rD = rcpp_to_double(p["rD"]);
  rU = rcpp_to_double(p["rU"]);
  rP = rcpp_to_double(p["rP"]);
  dE = rcpp_to_double(p["dE"]);
  tl = rcpp_to_double(p["tl"]);
  cD = rcpp_to_double(p["cD"]);
  cT = rcpp_to_double(p["cT"]);
  cU = rcpp_to_double(p["cU"]);
  g_inf =	rcpp_to_double(p["g_inf"]);
  d1 = rcpp_to_double(p["d1"]);
  dd = rcpp_to_double(p["dd"]);
  ID0	= rcpp_to_double(p["ID0"]);
  kd = rcpp_to_double(p["kd"]);
  ud = rcpp_to_double(p["ud"]);
  ad0 = rcpp_to_double(p["ad0"]);
  fd0 = rcpp_to_double(p["fd0"]);
  gd = rcpp_to_double(p["gd"]);
  aA = rcpp_to_double(p["aA"]);
  aU = rcpp_to_double(p["aU"]);
  b0 = rcpp_to_double(p["b0"]);
  b1 = rcpp_to_double(p["b1"]);
  db = rcpp_to_double(p["db"]);
  IB0 = rcpp_to_double(p["IB0"]);
  kb = rcpp_to_double(p["kb"]);
  ub = rcpp_to_double(p["ub"]);
  phi0 = rcpp_to_double(p["phi0"]);
  phi1 = rcpp_to_double(p["phi1"]);
  dc = rcpp_to_double(p["dc"]);
  IC0	= rcpp_to_double(p["IC0"]);
  kc = rcpp_to_double(p["kc"]);
  uc = rcpp_to_double(p["uc"]);
  PM = rcpp_to_double(p["PM"]);
  dm = rcpp_to_double(p["dm"]);
  tau	= rcpp_to_double(p["tau"]);
  mu = rcpp_to_double(p["mu"]);
  f = rcpp_to_double(p["f"]);
  Q0 = rcpp_to_double(p["Q0"]);
}
