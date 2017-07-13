
#include "parameters.h"

using namespace std;


//------------------------------------------------
// parameters::
// default constructor for parameters class
parameters::parameters()
{
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
// parameters::
// informed constructor for parameters class, taking values from an R list
parameters::parameters(const Rcpp::List& p)
{
    // take parameter values from an R list
    eta = Rcpp::as<double>(p["eta"]);
    rho = Rcpp::as<double>(p["rho"]);
    a0 = Rcpp::as<double>(p["a0"]);
    s2 = Rcpp::as<double>(p["s2"]);
    rA = Rcpp::as<double>(p["rA"]);
    rT = Rcpp::as<double>(p["rT"]);
    rD = Rcpp::as<double>(p["rD"]);
    rU = Rcpp::as<double>(p["rU"]);
    rP = Rcpp::as<double>(p["rP"]);
    dE = Rcpp::as<double>(p["dE"]);
    tl = Rcpp::as<double>(p["tl"]);
    cD = Rcpp::as<double>(p["cD"]);
    cT = Rcpp::as<double>(p["cT"]);
    cU = Rcpp::as<double>(p["cU"]);
    g_inf =	Rcpp::as<double>(p["g_inf"]);
    d1 = Rcpp::as<double>(p["d1"]);
    dd = Rcpp::as<double>(p["dd"]);
    ID0	= Rcpp::as<double>(p["ID0"]);
    kd = Rcpp::as<double>(p["kd"]);
    ud = Rcpp::as<double>(p["ud"]);
    ad0 = Rcpp::as<double>(p["ad0"]);
    fd0 = Rcpp::as<double>(p["fd0"]);
    gd = Rcpp::as<double>(p["gd"]);
    aA = Rcpp::as<double>(p["aA"]); // TODO find
    aU = Rcpp::as<double>(p["aU"]); // TODO find
    b0 = Rcpp::as<double>(p["b0"]);
    b1 = Rcpp::as<double>(p["b1"]);
    db = Rcpp::as<double>(p["db"]);
    IB0 = Rcpp::as<double>(p["IB0"]);
    kb = Rcpp::as<double>(p["kb"]);
    ub = Rcpp::as<double>(p["ub"]);
    phi0 = Rcpp::as<double>(p["phi0"]);
    phi1 = Rcpp::as<double>(p["phi1"]);
    dc = Rcpp::as<double>(p["dc"]);
    IC0	= Rcpp::as<double>(p["IC0"]);
    kc = Rcpp::as<double>(p["kc"]);
    uc = Rcpp::as<double>(p["uc"]);
    PM = Rcpp::as<double>(p["PM"]);
    dm = Rcpp::as<double>(p["dm"]);
    tau	= Rcpp::as<double>(p["tau"]);
    mu = Rcpp::as<double>(p["mu"]);
    f = Rcpp::as<double>(p["f"]);
    Q0 = Rcpp::as<double>(p["Q0"]);
}
