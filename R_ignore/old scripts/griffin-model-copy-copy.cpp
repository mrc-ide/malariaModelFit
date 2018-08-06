/*
Implementation of the Griffin et al. 2014 model
 - deterministic equilibrium solution for the transmission model
 - data fitting using MCMC
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rmath.h>
#include <cmath>
#include <fstream>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////////////
/////////////// Declarations of classes and functions
////////////////////////////////////////////////////////////////////////////////////////


////////////// Classes for the transmission model //////////////
////// NOTE: The code for the equilibrium solution is translated from Jamie's R code
////// But it is vectorised over the Gauss-Hermit numerical integration
//////

// parameters for the transmission model
class parameters
{
public:
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

    // constructors
    parameters(); // default parameters
    parameters(const List& p); // import from a R list

    // method
    //void get_par(const vp_eta& vp); // get model parameters from vp_eta
};

// class of population states
class population
{
public:
    double EIR; // annualised EIR
    //double dEIR; // daily EIR
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

    // Below are some intermediate quantities needed to set equilibrium solution
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

    // constructors
    //population();
    population(const arma::vec& age0, const arma::vec& ghn0, const arma::vec& ghw0);

    // find the equilibrium given EIR (EIR0), proportion treated (ft) and model parameters
    void set_equilibrium(double EIR0, double ft0, const parameters& p);
    // get the annual incidence rate between age0 and age1
    double get_inc(double age0, double age1);
    // get the prevalence between age0 and age1 by Microscopy or PCR
    double get_prev_M(double age0, double age1);
    double get_prev_PCR(double age0, double age1);

//private:

};


////////////// Classes for the MCMC fitting incl. data/parameters //////////////

// class of the Eta: the variable parameters in the MCMC fitting that could be sampled
// theta: original parameters for the transmission model, EIR and ft
// some theta's have restricted support - need to be transformed to eta with infinity support
class Eta
{
public:
    int ngp;         // number of global parameters for equilibrium solution
    int nhp;         // number of other hyper-parameters in the likelihood
    int nlogeir;     // number of site-specific logeir parameters
    int nft;         // number of site-specific ft parameters
    NumericVector eta;   // a vector containing the ngp + nlogeir + nft parameters

    NumericVector logeir_mean; // Normal prior mean for the log-EIR
    NumericVector logeir_sd;   // Normal prior sd for the log-EIR
    NumericVector ft_a;        // Beta prior for ft, Beta(a,b)
    NumericVector ft_b;

    // Constructor:
    // without initialising eta values
    Eta(const DataFrame& data_key);
    // initialise eta values
    // those with indices in init_eta_ids will be initialised randomly if r=true or
    // pre-specified default values if r=false
    Eta(const DataFrame& data_key, bool r, const IntegerVector& init_eta_ids);

    // method to return the log-likelihood (density) of the prior at current eta values
    double loglik();

    // method to return transmission model parameters corresponding to current eta values
    parameters get_parameters();

    // method to update the values of eta with indices in eta_ids using Metropolis random walk update
    // eta_new = N(eta_current, Sigma)
    void update(IntegerVector eta_ids, arma::mat& Sigma);
};

// class for incidence data at one site
class inc_data
{
public:
    // data
    IntegerVector clin_k; // number of clinical events
    NumericVector clin_p_years; // number of person years for clinical events
    NumericVector age0; // lower end of age group
    NumericVector age1; // upper end of age group
    IntegerVector clin_group; // group for clinical data random effects by age

    // for likelihood calculation
    NumericVector clin_mu; // model-predicted incidence of clinical malaria
    NumericVector clin_mu_copy; // a copy of clin_mu for temporary update

    // Constructors:
    inc_data();
    inc_data(IntegerVector& clink, NumericVector& clinpy, NumericVector& age0, NumericVector& age1,
             IntegerVector& clingroup, NumericVector& clinmu);

    // method to print out to visualise the data
    void print_out();
};

// class for prevalence data at one site
class prev_data
{
public:
    // data
    IntegerVector slide_k; // number positive by slide (Microscopy)
    IntegerVector slide_n; // number sampled for slide (and for PCR)
    IntegerVector pcr_k; // number positive by PCR
    NumericVector age0; // lower end of age group
    NumericVector age1; // upper end of age group

    // for likelihood calculation
    NumericVector slide_mu; // model-predicted prevalence by slide (Microscopy)
    NumericVector pcr_mu; // model-predicted prevalence by PCR
    NumericVector slide_mu_copy; // a copy of slide_mu
    NumericVector pcr_mu_copy; // a copy of pcr_mu

    // Constructors:
    prev_data();
    prev_data(IntegerVector& slidek, IntegerVector& sliden, IntegerVector& pcrk,
              NumericVector& age0, NumericVector& age1,
              NumericVector& slidemu, NumericVector& pcrmu);

    // method to print out to visualise the data
    void print_out();
};

// class for the data including incidence and/or prevalence data at one site
class site_data
{
public:
    int site_id; // index of this site
    int study_id; // index of the study to which this site belongs
    int logeir_id; // index of the logeir of this site in the logeir_mean/sd vectors
    int ft_id; // index of the ft of this site in the ft_a/b vectors
    int case_detection;
    bool for_inc; // if this site has data for inc
    bool for_prev; // if this site has data for prev

    // incidence data at this site if there is one
    inc_data inc;
    // prevalence data at this site if there is one
    prev_data prev;

    // Constructor:
    site_data(inc_data& incdata, prev_data& prevdata, int siteid, int studyid, int logeirid, int ftid,
              int casedetection, bool forinc, bool forprev);

    // method to evaluate the likelihood contribution from this site
    void loglik(bool update_data, Eta& eta, parameters& par, population& pop, const NumericVector& ghnodes, const NumericVector& ghweights,
                std::vector<NumericVector>& loglikincgh, std::vector<NumericVector>& loglikprevgh);

    // method to update the model-predicted incidence/prevalence if necessary
    // ONLY when the proposed eta values are accepted
    void update_data();

    // method to print out to visualise the data
    void print_out();
};

// class for all the data used in Griffin et al. 2014
class all_data
{
public:
    std::vector<site_data> sites; // in the order of site in data_key
    int nstudies; // number of studies
    int nsites;   // number of sites
    IntegerVector study_ids; // study ids
    std::vector<IntegerVector> study_site_ids; // site ids for each study

    // a nsites by ngh matrix for the contribution of incidence data to likelihood from each site
    std::vector<NumericVector> loglik_inc_gh;
    // a nsites by ngh matrix for the contribution of prevalence data to likelihood from each site
    std::vector<NumericVector> loglik_prev_gh;
    // a nstudies vector for the contribution to likelihood from each study
    NumericVector loglik_studies;

    // copies of the above contribution components
    // this is for the temporary storage since an update may be accepted or rejected
    // hence we need to store both the current and updated likelihood at the same time
    std::vector<NumericVector> loglik_inc_gh_copy;
    std::vector<NumericVector> loglik_prev_gh_copy;
    NumericVector loglik_studies_copy;

    // Constructors:
    // initialise the loglik_inc/prev_gh; note that ghnodes and ghweights are different from those in population
    all_data(DataFrame& data_key, List& datasets, const NumericVector& ghnodes, const NumericVector& ghweights);

    // Not initialise the Gauss-Hermit relevant matrices;
    // used for inference on u since the likelihood evaluation doesn't need G-H quadrature
    all_data(DataFrame& data_key, List& datasets);

    // method to evaluate the log-likelihood contribution from study with index of studyid
    void loglik_study(int studyid, const NumericVector& ghweights);

    // method to evaluate the log-likelihood from all data at the parameters eta
    double loglik(const IntegerVector& eta_ids, Eta& eta, const IntegerVector& eta_ids_site_ids,
                  const IntegerVector& eta_ids_study_ids, population& pop,
                  const NumericVector& ghnodes1, const NumericVector& ghweights1, bool initialise);

    // update the data depending on the acceptance or rejection
    void update_data(const IntegerVector& eta_ids_site_ids);

    // update the loglik_inc_gh, loglik_prev_gh, loglik_studies if necessary
    void update_loglik_gh();
};


// some simple functions
double t1(double eta)
{
    // eta to theta for theta > 0
    return(exp(eta));
}

double t2(double eta)
{
    // eta to theta for 0 < theta < 1
    return(1/(1+exp(-eta)));
}

double tb1(double theta)
{
    // theta to eta for theta > 0
    return(log(theta));
}

double tb2(double theta)
{
    // theta to eta for 0 < theta < 1
    return(log(theta/(1-theta)));
}

double dt1(double eta)
{
    // derivative of t1
    return(exp(eta));
}

double dt2(double eta)
{
    // derivative of t2
    return(exp(eta)/pow(1+exp(eta),2));
}

double ldt1(double eta)
{
    // log of the derivative of t1
    return(eta);
}

double ldt2(double eta)
{
    // log of the derivative of t2
    return(eta+2*log(1+exp(eta)));
}

double logdhnorm(double x, double mean, double sd, double scale)
{
    // log-density of half normal x/scale ~ N(mean, sd) with x>0
    double pab = 1 - R::pnorm(0,mean,sd,1,0);
    return(-log(scale) + R::dnorm(x/scale,mean,sd,1)-log(pab));
}

double logdht(double x, double df, double scale)
{
    // log-density of half t x/scale ~ t(df) with x>0
    double pab = 1 - R::pt(0,df,1,0);
    return(-log(scale) + R::dt(x/scale,df,1)-log(pab));
}

NumericVector logdbb(int k, int n, NumericVector a, NumericVector b)
{
    // log-density of Beta-binomial
    return(R::lchoose(n,k)+lbeta(a+k,n-k+b)-lbeta(a,b));
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma)
{
   // Draw n samples from multivariate Normal distribution N(mu, sigma)
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
}

arma::vec mvrnormArma2(arma::vec mu, arma::mat sigma)
{
    // Draw One samples from multivariate Normal distribution N(mu, sigma)
    int ncols = sigma.n_cols;
    arma::vec Y = arma::randn(ncols);
    return(mu + arma::chol(sigma).t() * Y);
}

NumericVector tonv(arma::vec& v1)
{
    // convert arma::vec object to NumericVector object
    NumericVector v(v1.begin(),v1.end());
    return(v);
}

IntegerVector which(LogicalVector& v)
{
    // Return the index of 'true' elements in v
    IntegerVector ind = seq(0, v.size()-1);
    return(ind[v]);
}


////////////////////////////////////////////////////////////////////////////////////////
/////////////// Definitions of classes and functions
////////////////////////////////////////////////////////////////////////////////////////

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

parameters::parameters(const List& p)
{
    // Return parameters object from a R list
    eta = as<double>(p["eta"]);
    rho = as<double>(p["rho"]);
    a0 = as<double>(p["a0"]);
    s2 = as<double>(p["s2"]);
    rA = as<double>(p["rA"]);
    rT = as<double>(p["rT"]);
    rD = as<double>(p["rD"]);
    rU = as<double>(p["rU"]);
    rP = as<double>(p["rP"]);
    dE = as<double>(p["dE"]);
    tl = as<double>(p["tl"]);
    cD = as<double>(p["cD"]);
    cT = as<double>(p["cT"]);
    cU = as<double>(p["cU"]);
    g_inf =	as<double>(p["g_inf"]);
    d1 = as<double>(p["d1"]);
    dd = as<double>(p["dd"]);
    ID0	= as<double>(p["ID0"]);
    kd = as<double>(p["kd"]);
    ud = as<double>(p["ud"]);
    ad0 = as<double>(p["ad0"]);
    fd0 = as<double>(p["fd0"]);
    gd = as<double>(p["gd"]);
    aA = as<double>(p["aA"]);
    aU = as<double>(p["aU"]);
    b0 = as<double>(p["b0"]);
    b1 = as<double>(p["b1"]);
    db = as<double>(p["db"]);
    IB0 = as<double>(p["IB0"]);
    kb = as<double>(p["kb"]);
    ub = as<double>(p["ub"]);
    phi0 = as<double>(p["phi0"]);
    phi1 = as<double>(p["phi1"]);
    dc = as<double>(p["dc"]);
    IC0	= as<double>(p["IC0"]);
    kc = as<double>(p["kc"]);
    uc = as<double>(p["uc"]);
    PM = as<double>(p["PM"]);
    dm = as<double>(p["dm"]);
    tau	= as<double>(p["tau"]);
    mu = as<double>(p["mu"]);
    f = as<double>(p["f"]);
    Q0 = as<double>(p["Q0"]);
}

population::population(const arma::vec& age0, const arma::vec& ghn0, const arma::vec& ghw0)
{
    // initialise a population given age groups (age0) and Gauss-Hermit nodes (ghn0) and
    // the associated weights (ghw0)

    na = age0.n_rows;
	age = age0;
    ghnodes = ghn0;
	ghweights = ghw0;

	S = arma::vec(na); S.fill(0.0);
	T = arma::vec(na); T.fill(0.0);
    D = arma::vec(na); D.fill(0.0);
    A = arma::vec(na); A.fill(0.0);
    U = arma::vec(na); U.fill(0.0);
    P = arma::vec(na); P.fill(0.0);
    inf = arma::vec(na); inf.fill(0.0);
    prop = arma::vec(na); prop.fill(0.0);
    psi = arma::vec(na); psi.fill(0.0);
    pos_M = arma::vec(na); pos_M.fill(0.0);
    pos_PCR = arma::vec(na); pos_PCR.fill(0.0);
    inc = arma::vec(na); inc.fill(0.0);

    FOIM = 0.0;

    nh = ghnodes.n_rows;
	zeta = arma::vec(nh); zeta.fill(0.0);

	mICA = arma::mat(nh, na); mICA.fill(0.0);
	mICM = arma::mat(nh, na); mICM.fill(0.0);

	mS = arma::mat(nh, na); mS.fill(0.0);
	mA = arma::mat(nh, na); mA.fill(0.0);
	mT = arma::mat(nh, na); mT.fill(0.0);
	mD = arma::mat(nh, na); mD.fill(0.0);
	mU = arma::mat(nh, na); mU.fill(0.0);
	mP = arma::mat(nh, na); mP.fill(0.0);

	r = arma::vec(na); r.fill(0.0);

	mFOI = arma::mat(nh, na); mFOI.fill(0.0);
	mphi = arma::mat(nh, na); mphi.fill(0.0);
	mq = arma::mat(nh, na); mq.fill(0.0);

	mcA = arma::mat(nh, na); mcA.fill(0.0);
	mpos_M = arma::mat(nh, na); mpos_M.fill(0.0);
	mpos_PCR = arma::mat(nh, na); mpos_PCR.fill(0.0);
	minc = arma::mat(nh, na); minc.fill(0.0);
	minf = arma::mat(nh, na); minf.fill(0.0);

	IB = arma::vec(nh); IB.fill(0.0);
    IC = arma::vec(nh); IC.fill(0.0);
    ID = arma::vec(nh); ID.fill(0.0);
    b = arma::vec(nh); b.fill(0.0);

    //double betaT, betaD, betaP;
    betaS = arma::vec(nh); betaS.fill(0.0);
    betaA = arma::vec(nh); betaA.fill(0.0);
    betaU = arma::vec(nh); betaU.fill(0.0);

    aT = arma::vec(nh); aT.fill(0.0);
    bT = arma::vec(nh); bT.fill(0.0);
    aD = arma::vec(nh); aD.fill(0.0);
    bD = arma::vec(nh); bD.fill(0.0);
    aP = arma::vec(nh); aP.fill(0.0);
    bP = arma::vec(nh); bP.fill(0.0);

    Y = arma::vec(nh); Y.fill(0.0);
}

void population::set_equilibrium(double EIR0, double ft0, const parameters& p)
{
    // find the equilibrium solution
    // translated from Jamie's R code
    EIR = EIR0;
    ft = ft0;

    FOIM = 0.0;

	zeta = exp(-p.s2*0.5 + sqrt(p.s2)*ghnodes);

	arma::vec dage = age*365; // age in days
	arma::vec dEIR = zeta*EIR/365; // daily EIR

	mICA.fill(0.0);
	mICM.fill(0.0);

	mS.fill(0.0);
	mA.fill(0.0);
	mT.fill(0.0);
	mD.fill(0.0);
    mU.fill(0.0);
	mP.fill(0.0);

	prop.fill(0.0);
	r.fill(0.0);

	mFOI.fill(0.0);
	mphi.fill(0.0);
	mq.fill(0.0);

	mcA.fill(0.0);
	mpos_M.fill(0.0);
	mpos_PCR.fill(0.0);
	minc.fill(0.0);
	minf.fill(0.0);


	for(int i = 0; i < na; ++i)
    {
        r[i] = (i==(na-1) ? 0 : 1/(dage[i+1]-dage[i]));
        prop[i] = (i==0 ? p.eta : r[i-1]*prop[i-1])/(r[i]+p.eta);
        if(i < (na-1))
        {
            dage[i] = (dage[i]+dage[i+1])*0.5;
        }
    }


    psi = 1-p.rho*exp(-dage/p.a0);
    IB.fill(0.0);
    IC.fill(0.0);
    ID.fill(0.0);
    b.fill(0.0);
    int age20 = 0;
    double re;

    for(int i = 0; i < na; ++i)
    {
        if(i < (na-1) && dage[i] <= 20*365 && dage[i+1] > 20*365)
        {
            age20 = i;
        }


        re = r[i] + p.eta;
        IB = (dEIR*psi[i]/(dEIR*psi[i]*p.ub+1) + re*IB)/(1/p.db+re);
        b = p.b0*(p.b1+(1-p.b1)/(1+pow(IB/p.IB0,p.kb)));
		mFOI.col(i) = dEIR%b*psi[i];


		IC = (mFOI.col(i)/(mFOI.col(i)*p.uc+1) + re*IC)/(1/p.dc+re);
		mICA.col(i) = IC;
		ID = (mFOI.col(i)/(mFOI.col(i)*p.ud+1) + re*ID)/(1/p.dd+re);
		mq.col(i)= p.d1+(1-p.d1)/(1+pow(ID/p.ID0,p.kd)*(1-(1-p.fd0)/(1+pow(dage[i]/p.ad0,p.gd))));
		mcA.col(i) = p.cU + (p.cD-p.cU)*pow(mq.col(i),p.g_inf);

    }

    arma::vec IM0 = mICA.col(age20)*p.PM;
    mICM.col(0) = IM0*(r[0] + p.eta)/(1/p.dm + r[0] + p.eta);
	for(int i = 1; i < na; ++i)
    {
        mICM.col(i) = (r[i] + p.eta)*mICM.col(i-1)/(1/p.dm + r[i] + p.eta);
    }

    mphi = p.phi0*(p.phi1+(1-p.phi1)/(1+pow((mICA+mICM)/p.IC0,p.kc)));

    betaS.fill(0.0);
    betaA.fill(0.0);
    betaU.fill(0.0);

    aT.fill(0.0);
    bT.fill(0.0);
    aD.fill(0.0);
    bD.fill(0.0);
    aP.fill(0.0);
    bP.fill(0.0);

    Y.fill(0.0);

	for(int i = 0; i < na; ++i)
    {
		re = r[i] + p.eta;
		betaS = mFOI.col(i) + re;
		betaT = p.rT + re;
		betaD = p.rD + re;
		betaA = mphi.col(i)%mFOI.col(i)+ p.rA + re;
		betaU = mFOI.col(i)+ p.rU + re;
		betaP = p.rP + re;

		aT = ft*mphi.col(i)%mFOI.col(i)/betaT;
		if(i!=0) bT = r[i-1]*mT.col(i-1)/betaT;
		aD = (1-ft)*mphi.col(i)%mFOI.col(i)/betaD;
		if(i!=0) bD = r[i-1]*mD.col(i-1)/betaD;
		aP = p.rT*aT/betaP;
		if(i==0)
        {
            bP = p.rT*bT/betaP;
        }
        else
        {
            bP = (p.rT*bT + r[i-1]*mP.col(i-1))/betaP;
        }

		Y = (prop[i]-(bT+bD+bP))/(1+aT+aD+aP);

		mT.col(i) = aT%Y + bT;
		mD.col(i) = aD%Y + bD;
		mP.col(i) = aP%Y + bP;

		if(i==0)
        {
            mA.col(i) = ((1-mphi.col(i))%mFOI.col(i)%Y + p.rD*mD.col(i))/(betaA+(1-mphi.col(i))%mFOI.col(i));
        }
        else
        {
            mA.col(i) = (r[i-1]*mA.col(i-1) + (1-mphi.col(i))%mFOI.col(i)%Y + p.rD*mD.col(i))/(betaA+(1-mphi.col(i))%mFOI.col(i));
        }

		if(i==0)
        {
            mU.col(i) = (p.rA*mA.col(i))/betaU;
        }
        else
        {
            mU.col(i) = (p.rA*mA.col(i)+r[i-1]*mU.col(i-1))/betaU;
        }
		mS.col(i) = Y-mA.col(i)-mU.col(i);

		mpos_M.col(i) = mD.col(i) + mT.col(i) + mA.col(i)%mq.col(i); // Microsopy
		mpos_PCR.col(i) = mD.col(i) + mT.col(i) + mA.col(i)%pow(mq.col(i),p.aA) + mU.col(i)%pow(mq.col(i),p.aU); // PCR
		minc.col(i) = mphi.col(i)%mFOI.col(i)%Y;
	}

	minf = p.cD*mD + p.cT*mT + mcA%mA + p.cU*mU;

	arma::mat mFOIM = (psi.t()*minf.t())*(ghweights%zeta);
	double omega = 1-p.rho*p.eta/(p.eta+1/p.a0);
	double alpha = p.f*p.Q0;

	S = mS.t()*ghweights;
	T = mT.t()*ghweights;
	D = mD.t()*ghweights;
	A = mA.t()*ghweights;
	U = mU.t()*ghweights;
	P = mP.t()*ghweights;
	inf = minf.t()*ghweights;
	pos_M = mpos_M.t()*ghweights;
	pos_PCR = mpos_PCR.t()*ghweights;
	inc = minc.t()*ghweights;
	FOIM = mFOIM(0,0)*alpha/omega;
}

double population::get_inc(double age0, double age1)
{
    // Find the incidence rate between age0 and age1 in a equilibrium population
    double sum_inc = 0.0;
    double sum_prop = 0.0;
    for(int i = 0; i < na; ++i)
    {
        sum_inc += (age[i]>=age0 && age[i]<age1) ? inc[i] : 0.0;
        sum_prop += (age[i]>=age0 && age[i]<age1) ? prop[i] : 0.0;
    }
    return(365*sum_inc/sum_prop);
}

double population::get_prev_M(double age0, double age1)
{
    // Find the prevalence by Microscopy between age0 and age1 in a equilibrium population
    double sum_pos = 0.0;
    double sum_prop = 0.0;
    for(int i = 0; i < na; ++i)
    {
        sum_pos += (age[i]>=age0 && age[i]<age1) ? pos_M[i] : 0.0;
        sum_prop += (age[i]>=age0 && age[i]<age1) ? prop[i] : 0.0;
    }
    return(sum_pos/sum_prop);
}

double population::get_prev_PCR(double age0, double age1)
{
    // Find the prevalence by PCR between age0 and age1 in a equilibrium population
    double sum_pos = 0.0;
    double sum_prop = 0.0;
    for(int i = 0; i < na; ++i)
    {
        sum_pos += (age[i]>=age0 && age[i]<age1) ? pos_PCR[i] : 0.0;
        sum_prop += (age[i]>=age0 && age[i]<age1) ? prop[i] : 0.0;
    }
    return(sum_pos/sum_prop);
}


////////////////////////////////////////////////////////////////////////////////////////////////

////// An exportable function to find the equilibrium solution of the transmission model
////// Equivalent to the R-version function "human_equilibrium" by Jamie

// [[Rcpp::export]]
List equi2(double EIR, double ft, const List& p, const arma::vec& age0,
          const arma::vec& ghnodes, const arma::vec& ghweights)
{
    parameters par(p);
    population pop(age0, ghnodes, ghweights);
    pop.set_equilibrium(EIR, ft, par);
    NumericMatrix E = cbind(tonv(pop.S), tonv(pop.T), tonv(pop.D), tonv(pop.A), tonv(pop.U), tonv(pop.P), tonv(pop.inf),
                                 tonv(pop.prop), tonv(pop.psi), tonv(pop.pos_M), tonv(pop.pos_PCR), tonv(pop.inc));
    colnames(E) = CharacterVector::create("S","T","D","A","U","P","inf","prop","psi","pos_M","pos_PCR","inc");
    List result = List::create(_["states"]=E, _["FOIM"]=pop.FOIM);

    return(result);
}


////////////////////////////////////////////////////////////////////////////////////////////////

// A blank default constructor for inc_data object
inc_data::inc_data(){};

// A customised constructor to create an inc_data object with user-input values
inc_data::inc_data(IntegerVector& clink, NumericVector& clinpy, NumericVector& a0, NumericVector& a1,
                   IntegerVector& clingroup, NumericVector& clinmu)
{
    clin_k = clink;
    clin_p_years = clinpy;
    age0 = a0;
    age1 = a1;
    clin_group = clingroup;
    clin_mu = clinmu;
    clin_mu_copy = clone(clin_mu);
}

// A blank default constructor for prev_data object
prev_data::prev_data(){};

// A customised constructor to create an prev_data object with user-input values
prev_data::prev_data(IntegerVector& slidek, IntegerVector& sliden, IntegerVector& pcrk,
                     NumericVector& a0, NumericVector& a1,
                     NumericVector& slidemu, NumericVector& pcrmu)
{
    slide_k = slidek;
    slide_n = sliden;
    pcr_k = pcrk;
    age0 = a0;
    age1 = a1;
    slide_mu = slidemu;
    pcr_mu = pcrmu;
    slide_mu_copy = clone(slide_mu);
    pcr_mu_copy = clone(pcr_mu);
}

// Print out an inc_data, prev_data or site_data object
void inc_data::print_out()
{
    Rcout << std::setw(8) << "clin_k" << std::setw(8) << "clin_py" << std::setw(8) << "age0" << std::setw(8) << "age1" << std::setw(8) << "clin_gr" << "\n";
    for(int i = 0; i < clin_k.size(); ++i)
    {
        Rcout << std::setw(8) << clin_k[i] << std::setw(8) << clin_p_years[i] << std::setw(8) << age0[i] << std::setw(8) << age1[i] << std::setw(8) << clin_group[i] << "\n";
    }
}

void prev_data::print_out()
{
    Rcout << std::setw(8) << "slide_k" << std::setw(8) << "slide_n" << std::setw(8) << "pcr_k" << std::setw(8) << "age0" << std::setw(8) << "age1" << "\n";
    for(int i = 0; i < slide_k.size(); ++i)
    {
        Rcout << std::setw(8) << slide_k[i] << std::setw(8) << slide_n[i] << std::setw(8) << pcr_k[i] << std::setw(8) << age0[i] << std::setw(8) << age1[i] << "\n";
    }
}

void site_data::print_out()
{
    Rcout << "Site - " << site_id << " in Study - " << study_id << " data \n";
    Rcout << "logeir_id: " << logeir_id << "ft_id: " << ft_id << "case detection: " << case_detection << "\n";
    if(for_inc)
    {
        Rcout << "    inc data:\n";
        inc.print_out();
    }
    if(for_prev)
    {
        Rcout << "    prev data:\n";
        prev.print_out();
    }
}

// A constructor to create a site_data object from user-input values
site_data::site_data(inc_data& incdata, prev_data& prevdata, int siteid, int studyid, int logeirid, int ftid,
              int casedetection, bool forinc, bool forprev)
{
    site_id = siteid;
    inc = incdata;
    prev = prevdata;
    logeir_id = logeirid;
    ft_id = ftid;
    for_inc = forinc;
    for_prev = forprev;
    study_id = studyid;
    case_detection = casedetection;
}


// A method to evaluate the log-likelihood contribution from one site
// Note: This is for the marginalised likelihood over random effects
//       Study level random effect u is integrated out numerically using adaptive GH quadrature
//       Site-age level random effect v_jk is integrated out analytically
void site_data::loglik(bool update_data, Eta& eta, parameters& par, population& pop,
                       const NumericVector& ghnodes, const NumericVector& ghweights,
                       std::vector<NumericVector>& loglikincgh, std::vector<NumericVector>& loglikprevgh)
{
    // Calculate the loglik contribution to inc and prev of this site
    // The log-likelihood is evaluated on each adaptive GH quadrature nodes
    double eir = exp(eta.eta[eta.ngp + eta.nhp + logeir_id]); // eir for the site
    double ft = t2(eta.eta[eta.ngp + eta.nhp + eta.nlogeir + ft_id]); // ft for the site

    // find the equilibrium for this site if required
    if(update_data)
    {
        pop.set_equilibrium(eir, ft, par);
    }

    // calculate the likelihood then...
    double a = 1/exp(eta.eta[24]); // a = 1/alphac; alphac/5 ~ half t-6
    double b = a;
    double theta = 1/exp(eta.eta[26]); // theta = 1/alphap; alphap/5 ~ half t-6
    double sigmac = exp(eta.eta[23]);
    double sigmap = exp(eta.eta[25]);

    int ngh = ghnodes.size(); // number of adaptive Gauss Hermit quadrature nodes

    // clinical incidence data
    NumericVector loglik_inc(ngh, 0.0);
    NumericVector loglik_prev(ngh, 0.0);
    if(for_inc)
    {
        // work on the temporary clin_mu_copy
        inc.clin_mu_copy = clone(inc.clin_mu);

        NumericVector u = sigmac*ghnodes - 0.5*sigmac*sigmac; // Adaptive to the mu/sigma - did Jamie do differently?
        LogicalVector in_groupk;
        IntegerVector k_ind;
        NumericVector sum1(ngh, 0.0);
        double sum2 = 0.0;
        NumericVector sum3(ngh, 0.0);
        double Tji, yji, muji, tempji, r;
        if(case_detection == 2)
        {
            r = 1/(1+exp(-eta.eta[21]));  // drW
        }
        else if(case_detection == 6)
        {
            r = 1/(1+exp(-eta.eta[22]));  // drP
        }
        else
        {
            r = 1.0;
        }
        for(int k = 0; k < 4; ++k)
        {
            // k for index of age group
            in_groupk = (inc.clin_group == k);
            k_ind = which(in_groupk);
            if(k_ind.size()!=0)
            {
                sum1.fill(0.0);
                sum2 = 0.0;
                sum3.fill(0.0);
                for(IntegerVector::iterator it = k_ind.begin(); it != k_ind.end(); ++it)
                {
                    if(update_data)
                    {
                        inc.clin_mu_copy[*it] = pop.get_inc(inc.age0[*it], inc.age1[*it]);
                    }
                    Tji = inc.clin_p_years[*it];
                    yji = inc.clin_k[*it];
                    muji = inc.clin_mu_copy[*it];
                    tempji = r*Tji*muji;
                    sum1 = sum1 + tempji*exp(u);
                    sum2 = sum2 + yji;
                    sum3 = sum3 + yji*(log(tempji)+u) - R::lgammafn(yji+1);
                }
                loglik_inc = loglik_inc + R::lgammafn(a+sum2) - (a+sum2)*log(b+sum1) + sum3 + a*log(b) - R::lgammafn(a);
            }
        }
    }

    // prevalence data
    if(for_prev)
    {
        // work on the temporary slide_mu_copy and pcr_mu_copy
        prev.slide_mu_copy = clone(prev.slide_mu);
        prev.pcr_mu_copy = clone(prev.pcr_mu);

        NumericVector w = sigmap*ghnodes;
        int nji, xji;
        double pji;
        NumericVector qji;
        for(int i = 0; i < prev.slide_n.size(); ++i)
        {
            nji = prev.slide_n[i];

            // Microscopy
            if(prev.slide_k[i]!=-1)
            {
                if(update_data)
                {
                    prev.slide_mu_copy[i] = pop.get_prev_M(prev.age0[i], prev.age1[i]);
                }
                xji = prev.slide_k[i];
                pji = prev.slide_mu_copy[i];
                qji = 1/(1+exp(-(log(pji/(1-pji))+w)));
                loglik_prev = loglik_prev + logdbb(xji, nji, qji*theta, (1-qji)*theta);
            }

            // PCR
            if(prev.pcr_k[i]!=-1)
            {
                if(update_data)
                {
                    prev.pcr_mu_copy[i] = pop.get_prev_PCR(prev.age0[i], prev.age1[i]);
                }
                xji = prev.pcr_k[i];
                pji = prev.pcr_mu_copy[i];
                qji = 1/(1+exp(-(log(pji/(1-pji))+w)));
                loglik_prev = loglik_prev + logdbb(xji, nji, qji*theta, (1-qji)*theta);
            }
        }
    }
    loglikincgh[site_id] = loglik_inc;
    loglikprevgh[site_id] = loglik_prev;
}

// method to update model-predicted incidence/prevalence to the accepted proposed values
void site_data::update_data()
{
    // incidence data
    if(for_inc)
    {
        inc.clin_mu = clone(inc.clin_mu_copy); // update to the accepted proposed values
    }

    // prevalence data
    if(for_prev)
    {
        prev.slide_mu = clone(prev.slide_mu_copy);
        prev.pcr_mu = clone(prev.pcr_mu_copy);
    }
}

// A constructor to create an all_data object from user-input values
all_data::all_data(DataFrame& data_key, List& datasets, const NumericVector& ghnodes, const NumericVector& ghweights)
{
    nsites = datasets.size();
    IntegerVector cd = as<IntegerVector>(data_key[7]); // case detection id
    IntegerVector si = as<IntegerVector>(data_key[8]); // study id
    study_ids = unique(si);
    nstudies = study_ids.size();
    IntegerVector pi = as<IntegerVector>(data_key[4]); // place id for ft
    study_site_ids = std::vector<IntegerVector>(nstudies);
    int ngh = ghnodes.size();
    for(int i = 0; i < nsites; ++i)
    {
        // site i
        int siteid_i = i;
        int logeirid_i = i;
        int ftid_i = pi[i];
        int casedetection_i = cd[i];
        int studyid_i = si[i];
        LogicalVector temp = study_ids==studyid_i;
        int temp_id = which(temp)[0];
        study_site_ids[temp_id].push_back(siteid_i);
        DataFrame data_i = datasets[i];
        IntegerVector sk = as<IntegerVector>(data_i[0]);
        IntegerVector sn = as<IntegerVector>(data_i[1]);
        IntegerVector pk = as<IntegerVector>(data_i[2]);
        IntegerVector ck = as<IntegerVector>(data_i[3]);
        NumericVector cpy = as<NumericVector>(data_i[4]);
        NumericVector a0 = as<NumericVector>(data_i[7]);
        NumericVector a1 = as<NumericVector>(data_i[8]);
        IntegerVector cg = as<IntegerVector>(data_i[9]);
        inc_data incdatai;
        prev_data prevdatai;
        bool forinci;
        bool forprevi;
        // clinical data
        if(is_true(all(ck == -1)))
        {
            forinci = false;
            //incdatai = inc_data();
        }
        else
        {
            forinci = true;
            LogicalVector temp = (ck!=-1);
            IntegerVector cki = ck[temp];
            NumericVector cpyi = cpy[temp];
            NumericVector a0i = a0[temp];
            NumericVector a1i = a1[temp];
            IntegerVector cgi = cg[temp];
            NumericVector ckmui(cki.size(), 0.0);
            incdatai = inc_data(cki, cpyi,a0i,a1i,cgi,ckmui);
        }
        // prevalence data
        if(is_true(all(sn == -1)))
        {
            forprevi = false;
            //prevdatai = prev_data();
        }
        else
        {
            forprevi = true;
            LogicalVector temp = (sn!=-1);
            IntegerVector ski = sk[temp];
            IntegerVector sni = sn[temp];
            IntegerVector pki = pk[temp];
            NumericVector a0i = a0[temp];
            NumericVector a1i = a1[temp];
            NumericVector skmui(ski.size(), 0.0);
            NumericVector pkmui(ski.size(), 0.0);
            prevdatai = prev_data(ski, sni, pki, a0i, a1i, skmui, pkmui);
        }
        site_data sitei = site_data(incdatai, prevdatai, siteid_i, studyid_i, logeirid_i, ftid_i, casedetection_i, forinci, forprevi);
        sites.push_back(sitei);
        NumericVector loglik_inc_gh_i(ngh, 0.0);
        NumericVector loglik_prev_gh_i(ngh, 0.0);
        loglik_inc_gh.push_back(loglik_inc_gh_i);
        loglik_prev_gh.push_back(loglik_prev_gh_i);
    }

    loglik_studies = NumericVector(nstudies, 0.0);

    loglik_inc_gh_copy = loglik_inc_gh;
    loglik_prev_gh_copy = loglik_prev_gh;
    loglik_studies_copy = NumericVector(nstudies, 0.0);
}

// A constructor to create an all_data object from user-input values
// without the Gauss-Hermit initialisation i.e. number of GH nodes = 0
all_data::all_data(DataFrame& data_key, List& datasets)
{
    nsites = datasets.size();
    IntegerVector cd = as<IntegerVector>(data_key[7]); // case detection id
    IntegerVector si = as<IntegerVector>(data_key[8]); // study id
    study_ids = unique(si);
    nstudies = study_ids.size();
    IntegerVector pi = as<IntegerVector>(data_key[4]); // place id for ft
    study_site_ids = std::vector<IntegerVector>(nstudies);
    int ngh = 0;
    for(int i = 0; i < nsites; ++i)
    {
        // site i
        int siteid_i = i;
        int logeirid_i = i;
        int ftid_i = pi[i];
        int casedetection_i = cd[i];
        int studyid_i = si[i];
        LogicalVector temp = study_ids==studyid_i;
        int temp_id = which(temp)[0];
        study_site_ids[temp_id].push_back(siteid_i);
        DataFrame data_i = datasets[i];
        IntegerVector sk = as<IntegerVector>(data_i[0]);
        IntegerVector sn = as<IntegerVector>(data_i[1]);
        IntegerVector pk = as<IntegerVector>(data_i[2]);
        IntegerVector ck = as<IntegerVector>(data_i[3]);
        NumericVector cpy = as<NumericVector>(data_i[4]);
        NumericVector a0 = as<NumericVector>(data_i[7]);
        NumericVector a1 = as<NumericVector>(data_i[8]);
        IntegerVector cg = as<IntegerVector>(data_i[9]);
        inc_data incdatai;
        prev_data prevdatai;
        bool forinci;
        bool forprevi;
        // clinical data
        if(is_true(all(ck == -1)))
        {
            forinci = false;
            //incdatai = inc_data();
        }
        else
        {
            forinci = true;
            LogicalVector temp = (ck!=-1);
            IntegerVector cki = ck[temp];
            NumericVector cpyi = cpy[temp];
            NumericVector a0i = a0[temp];
            NumericVector a1i = a1[temp];
            IntegerVector cgi = cg[temp];
            NumericVector ckmui(cki.size(), 0.0);
            incdatai = inc_data(cki, cpyi,a0i,a1i,cgi,ckmui);
        }
        // prevalence data
        if(is_true(all(sn == -1)))
        {
            forprevi = false;
            //prevdatai = prev_data();
        }
        else
        {
            forprevi = true;
            LogicalVector temp = (sn!=-1);
            IntegerVector ski = sk[temp];
            IntegerVector sni = sn[temp];
            IntegerVector pki = pk[temp];
            NumericVector a0i = a0[temp];
            NumericVector a1i = a1[temp];
            NumericVector skmui(ski.size(), 0.0);
            NumericVector pkmui(ski.size(), 0.0);
            prevdatai = prev_data(ski, sni, pki, a0i, a1i, skmui, pkmui);
        }
        site_data sitei = site_data(incdatai, prevdatai, siteid_i, studyid_i, logeirid_i, ftid_i, casedetection_i, forinci, forprevi);
        sites.push_back(sitei);
        NumericVector loglik_inc_gh_i(ngh, 0.0);
        NumericVector loglik_prev_gh_i(ngh, 0.0);
        loglik_inc_gh.push_back(loglik_inc_gh_i);
        loglik_prev_gh.push_back(loglik_prev_gh_i);
    }

    loglik_studies = NumericVector(nstudies, 0.0);

    loglik_inc_gh_copy = loglik_inc_gh;
    loglik_prev_gh_copy = loglik_prev_gh;
    loglik_studies_copy = NumericVector(nstudies, 0.0);
}

void all_data::loglik_study(int studyid, const NumericVector& ghweights)
{
    LogicalVector studyid_in_study_ids = (study_ids==studyid);
    int ind_study = which(studyid_in_study_ids)[0];
    IntegerVector studyid_site_ids = study_site_ids[ind_study];
    double loglik_studyid_inc = 0.0;
    double loglik_studyid_prev = 0.0;
    IntegerVector::iterator it;
    NumericVector colsum(ghweights.size(), 0.0);

    if(sites[studyid_site_ids[0]].for_inc)
    {
        colsum.fill(0.0);
        for(it = studyid_site_ids.begin(); it != studyid_site_ids.end(); ++it)
        {
            colsum = colsum + loglik_inc_gh_copy[*it];
        }
        loglik_studyid_inc = log(sum(ghweights*exp(colsum)));
    }

    if(sites[studyid_site_ids[0]].for_prev)
    {
        colsum.fill(0.0);
        for(it = studyid_site_ids.begin(); it != studyid_site_ids.end(); ++it)
        {
            colsum = colsum + loglik_prev_gh_copy[*it];
        }
        loglik_studyid_prev = log(sum(ghweights*exp(colsum)));
    }

    loglik_studies_copy[ind_study] = loglik_studyid_inc + loglik_studyid_prev;
}

double all_data::loglik(const IntegerVector& eta_ids, Eta& eta, const IntegerVector& eta_ids_site_ids,
                        const IntegerVector& eta_ids_study_ids, population& pop,
                        const NumericVector& ghnodes1, const NumericVector& ghweights1,
                        bool initialise)
{
    // copy the previous result to _copy for temporary use
    loglik_inc_gh_copy = loglik_inc_gh;
    loglik_prev_gh_copy = loglik_prev_gh;
    loglik_studies_copy = clone(loglik_studies);
    double loglik_total = eta.loglik();
    parameters par = eta.get_parameters();

    if(is_true(all(eta_ids >= eta.ngp)) && is_true(all(eta_ids < eta.ngp + eta.nhp)) && !initialise)
    {
        // if only hyperparameters are updated and it is not an initialisation step,
        // there is no need to update the model-predicted data
        for(IntegerVector::const_iterator it = eta_ids_site_ids.begin(); it != eta_ids_site_ids.end(); it++)
        {
            sites[*it].loglik(false, eta, par, pop, ghnodes1, ghweights1, loglik_inc_gh_copy, loglik_prev_gh_copy);
        }
    }
    else
    {
        for(IntegerVector::const_iterator it = eta_ids_site_ids.begin(); it != eta_ids_site_ids.end(); it++)
        {
            sites[*it].loglik(true, eta, par, pop, ghnodes1, ghweights1, loglik_inc_gh_copy, loglik_prev_gh_copy);
        }
    }

    for(IntegerVector::const_iterator it = eta_ids_study_ids.begin(); it != eta_ids_study_ids.end(); it++)
    {
        loglik_study(*it, ghweights1);
    }

    loglik_total += sum(loglik_studies_copy);
    return(loglik_total);
}

void all_data::update_data(const IntegerVector& eta_ids_site_ids)
{
    for(IntegerVector::const_iterator it = eta_ids_site_ids.begin(); it != eta_ids_site_ids.end(); it++)
    {
        sites[*it].update_data();
    }
}

void all_data::update_loglik_gh()
{
    loglik_inc_gh = loglik_inc_gh_copy;
    loglik_prev_gh = loglik_prev_gh_copy;
    loglik_studies = clone(loglik_studies_copy);
}

Eta::Eta(const DataFrame& data_key)
{
    ngp = 21;
    nhp = 6;
    logeir_mean = as<NumericVector>(data_key[2]);
    logeir_sd = as<NumericVector>(data_key[3]);
    nlogeir = logeir_mean.size();
    IntegerVector place_id = as<IntegerVector>(data_key[4]);
    //IntegerVector uni_place_id = unique(place_id);
    nft =  unique(place_id).size();
    NumericVector fa = as<NumericVector>(data_key[5]);
    NumericVector fb = as<NumericVector>(data_key[6]);
    //nplace =
    ft_a = NumericVector(17);
    ft_b = NumericVector(17);
    NumericVector fai, fbi;
    for(int i = 0; i < nft; ++i)
    {
        fai = fa[place_id == i];
        fbi = fb[place_id == i];
        ft_a[i] = fai[0];
        ft_b[i] = fbi[0];
    }
}

Eta::Eta(const DataFrame& data_key, bool r, const IntegerVector& init_eta_ids)
{
    ngp = 21;
    nhp = 6;
    logeir_mean = as<NumericVector>(data_key[2]);
    logeir_sd = as<NumericVector>(data_key[3]);
    nlogeir = logeir_mean.size();
    IntegerVector place_id = as<IntegerVector>(data_key[4]);
    //IntegerVector uni_place_id = unique(place_id);
    nft = unique(place_id).size();
    NumericVector fa = as<NumericVector>(data_key[5]);
    NumericVector fb = as<NumericVector>(data_key[6]);
    ft_a = NumericVector(17);
    ft_b = NumericVector(17);
    NumericVector fai, fbi;
    for(int i = 0; i < nft; ++i)
    {
        fai = fa[place_id == i];
        fbi = fb[place_id == i];
        ft_a[i] = fai[0];
        ft_b[i] = fbi[0];
    }
    eta = NumericVector(ngp+nhp+nlogeir+nft);
    // 21 gp: 0 - 20
    eta[0] = tb1(1/0.00906627); //rU = 0.00906627 // dU > 0 rU = 0.00906627
    eta[1] = tb2(0.160527); //d1 = 0.160527 // 0 < d1 < 1
    eta[2] = tb1(1.577533); //ID0 = 1.577533 // ID0 > 0
    eta[3] = tb1(0.476614); //kd = 0.476614 // kd > 0
    eta[4] = tb1(9.44512); //ud = 9.44512 // ud > 0
    eta[5] = tb1(8001.99/365); //ad0 = 8001.99 // ad0 > 0
    eta[6] = tb2(0.007055); //fd0 = 0.007055 // 0 < fd0 < 1
    eta[7] = tb1(4.8183); //gd = 4.8183 // gd > 0
    eta[8] = tb2(0.757); //aA = 0.757 // 0 < aA < 1
    eta[9] = tb2(0.186); //aU = 0.186 // 0 < aU < 1
    eta[10] = tb2(0.590076); //b0 = 0.590076 // 0 < b0 < 1
    eta[11] = tb1(43.8787); //IB0 = 43.8787 // IB0 > 0
    eta[12] = tb1(2.15506); //kb = 2.15506 // kb > 0
    eta[13] = tb1(7.19919); //ub = 7.19919 // ub > 0
    eta[14] = tb2(0.791666); //phi0 = 0.791666 // 0 < phi0 < 1
    eta[15] = tb2(0.000737); //phi1 = 0.000737 // 0 < phi1 < 1
    eta[16] = tb1(18.02366); //IC0	= 18.02366 // IC0 > 0
    eta[17] = tb1(2.36949); //kc = 2.36949 // kc > 0
    eta[18] = tb1(6.06349); //uc = 6.06349 // uc > 0
    eta[19] = tb2(0.774368); //PM = 0.774368 // 0 < PM < 1
    eta[20] = tb1(67.6952); //dm = 67.6952 // dm > 0
    // 6 hp: 21 - 26
    eta[21] = tb2(0.723); //drW = 0.723 // 0 < drW < 1
    eta[22] = tb2(0.342); //drP = 0.342 // 0 < drP < 1
    eta[23] = -0.7; //0.5 // sigmac > 0
    eta[24] = 0.7; //2 // alphac > 0
    eta[25] = -0.7; //0.5 // sigmap > 0
    eta[26] = 0.7; //2 // alphap > 0
    // eir
    for(int i = 0; i < nlogeir; ++i)
    {
        eta[ngp+nhp+i] = logeir_mean[i];
    }
    // ft
    for(int i = 0; i < nft; ++i)
    {
        eta[ngp+nhp+nlogeir+i] = tb2(ft_a[i]/(ft_a[i]+ft_b[i]));
    }
    // if random initial values for those with init_eta_ids are required with r=true
    if(r)
    {
        for(IntegerVector::const_iterator it = init_eta_ids.begin(); it != init_eta_ids.end(); it++)
        {
            eta[*it] = rnorm(1,0.0,2.0)[0];
        }
    }
}


double Eta::loglik()
{
    double loglik;
    loglik = R::dlnorm(exp(eta[0]), 4.00, 0.31, 1) + eta[0] +                                       // dU
             R::dbeta(1/(1+exp(-eta[1])), 10.1, 30.6, 1) + eta[1] - 2*log(1+exp(eta[1])) +          // d1
             R::dlnorm(exp(eta[2]), 2.38, 1.08, 1) + eta[2] +
             R::dlnorm(exp(eta[3]), 0.50, 0.32, 1) + eta[3] +
             R::dlnorm(exp(eta[4]), 1.40, 0.80, 1) + eta[4] +
             R::dlnorm(exp(eta[5]), 2.02, 0.43, 1) + eta[5] +
             R::dbeta(1/(1+exp(-eta[6])), 1, 1, 1) + eta[6] - 2*log(1+exp(eta[6])) +
             R::dlnorm(exp(eta[7]), 0.22, 0.64, 1) + eta[7] +
             R::dbeta(1/(1+exp(-eta[8])), 1, 1, 1) + eta[8] - 2*log(1+exp(eta[8])) +
             R::dbeta(1/(1+exp(-eta[9])), 4, 1, 1) + eta[9] - 2*log(1+exp(eta[9])) +
             R::dbeta(1/(1+exp(-eta[10])), 1.3, 1.3, 1) + eta[10] - 2*log(1+exp(eta[10])) +
             R::dlnorm(exp(eta[11]), 3.07, 0.93, 1) + eta[11] +
             R::dlnorm(exp(eta[12]), 0.50, 0.32, 1) + eta[12] +
             R::dlnorm(exp(eta[13]), 1.40, 0.80, 1) + eta[13] +
             R::dbeta(1/(1+exp(-eta[14])), 8.3, 2.1, 1) + eta[14] - 2*log(1+exp(eta[14])) +
             R::dbeta(1/(1+exp(-eta[15])), 1.1, 2.0, 1) + eta[15] - 2*log(1+exp(eta[15])) +
             R::dlnorm(exp(eta[16]), 3.76, 1.02, 1) + eta[16] +
             R::dlnorm(exp(eta[17]), 0.50, 0.32, 1) + eta[17] +
             R::dlnorm(exp(eta[18]), 1.40, 0.80, 1) + eta[18] +
             R::dbeta(1/(1+exp(-eta[19])), 1, 1, 1) + eta[19] - 2*log(1+exp(eta[19])) +
             R::dlnorm(exp(eta[20]), 5.26, 0.33, 1) + eta[20] +
             R::dbeta(1/(1+exp(-eta[21])), 8.9, 3, 1) + eta[21] - 2*log(1+exp(eta[21])) +
             R::dbeta(1/(1+exp(-eta[22])), 1, 1, 1) + eta[22] - 2*log(1+exp(eta[22])) +
             logdhnorm(exp(eta[23]), 0, 1, 5) + eta[23] +
             logdht(exp(eta[24]), 6, 5) + eta[24] +
             logdhnorm(exp(eta[25]), 0, 1, 5) + eta[25] +
             logdht(exp(eta[26]), 6, 5) + eta[26];
    for(int i = 0; i < nlogeir; ++i)
    {
        loglik += R::dnorm(eta[ngp+nhp+i],logeir_mean[i],logeir_sd[i],1);
    }
    for(int i = 0; i < nft; ++i)
    {
        loglik += R::dbeta(1/(1+exp(-eta[ngp+nhp+nlogeir+i])),ft_a[i],ft_b[i],1) + eta[ngp+nhp+nlogeir+i] -
                  2*log(1+exp(eta[ngp+nhp+nlogeir+i]));
    }
    return(loglik);
}

parameters Eta::get_parameters()
{
    parameters par;
    par.rU = 1/exp(eta[0]); // dU > 0; rU = 1/dU
    par.d1 = 1/(1+exp(-eta[1])); // 0 < d1 < 1
    par.ID0 = exp(eta[2]); // ID0 > 0
    par.kd = exp(eta[3]); // kd > 0
    par.ud = exp(eta[4]); // ud > 0
    par.ad0 = exp(eta[5])*365; // ad0 > 0; in days
    par.fd0 = 1/(1+exp(-eta[6])); // 0 < fd0 < 1
    par.gd = exp(eta[7]); // gd > 0
    par.aA = 1/(1+exp(-eta[8])); // 0 < aA < 1
    par.aU = 1/(1+exp(-eta[9])); // 0 < aU < 1
    par.b0 = 1/(1+exp(-eta[10])); // 0 < b0 < 1
    par.IB0 = exp(eta[11]); // IB0 > 0
    par.kb = exp(eta[12]); // kb > 0
    par.ub = exp(eta[13]); // ub > 0
    par.phi0 = 1/(1+exp(-eta[14])); // 0 < phi0 < 1
    par.phi1 = 1/(1+exp(-eta[15])); // 0 < phi1 < 1
    par.IC0 = exp(eta[16]); // IC0 > 0
    par.kc = exp(eta[17]); // kc > 0
    par.uc = exp(eta[18]); // uc > 0
    par.PM = 1/(1+exp(-eta[19])); // 0 < PM < 1
    par.dm = exp(eta[20]); // dm > 0
    //par.drW = 1/(1+exp(-eta[21])); // 0 < drW < 1
    //par.drP = 1/(1+exp(-eta[22])); // 0 < drP < 1
    return(par);
}

void Eta::update(IntegerVector eta_ids, arma::mat& Sigma)
{
    if(eta_ids.size()==1)
    {
        double eta0 = eta[eta_ids[0]];
        eta[eta_ids[0]] = rnorm(1,eta0,Sigma.col(0)[0])[0];
    }
    else
    {
        NumericVector eta0 = eta[eta_ids];
        arma::vec eta1 = mvrnormArma2(as<arma::vec>(eta0), Sigma);
        eta[eta_ids] = tonv(eta1);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
////// run the Metropolis MCMC fitting for nrep iterations
////// 1. The parameters are updated in blocks as specified by update_blocks
////// 2. update_blocks do not need to contain all the parameter indices, i.e. partial fitting is okay
////// 3. Proposal distribution for each block is tuned adaptively to achieve reasonable acceptance ratio

// [[Rcpp::export]]
List runMCMC(int nrep, DataFrame& data_key, List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0,
          const NumericVector& ghnodes1, const NumericVector& ghweights1, const std::vector<IntegerVector>& update_blocks)
{
    // read in the data
    all_data alldata(data_key, datasets, ghnodes1, ghweights1);
    // get a population object to work on
    population pop(age0, ghnodes0, ghweights0);

    // set up the updating blocks
    int nblocks = update_blocks.size();
    //std::vector<IntegerVector> update_blocks;
    std::vector<arma::mat> Sigma_blocks;
    std::vector<arma::vec> Mu_blocks;
    NumericVector lambda_blocks;

    IntegerVector update_eta_ids;
    for(int i = 0; i < nblocks; ++i)
    {
        for(int j = 0; j < update_blocks[i].size(); ++j)
        {
            update_eta_ids.push_back(update_blocks[i][j]);
        }
    }
    // initialise the eta parameters
    Eta eta(data_key, false, update_eta_ids);

    // the associated lambda, Mu, Sigma for proposal distribution
    for(int i = 0; i < nblocks; ++i)
    {
        int bi_len = update_blocks[i].size();
        Sigma_blocks.push_back(arma::eye(bi_len,bi_len));
        Mu_blocks.push_back(arma::ones(bi_len));
        lambda_blocks.push_back(1.0);
    }


    std::vector<arma::mat> proposal_Sigma_blocks = Sigma_blocks;
    std::vector<arma::vec> proposal_Mu_blocks = Mu_blocks;
    NumericVector proposal_lambda_blocks = clone(lambda_blocks);

    // update blocks the associated site or study ids
    std::vector<IntegerVector> update_blocks_site_ids(nblocks);
    std::vector<IntegerVector> update_blocks_study_ids(nblocks);


    for(int i = 0; i < nblocks; ++i)
    {
        IntegerVector eta_ids = update_blocks[i];
        if(is_true(any(eta_ids < (eta.ngp+eta.nhp))))
        {
            update_blocks_site_ids[i] = seq_len(alldata.nsites)-1;
            update_blocks_study_ids[i] = alldata.study_ids;
        }
        else
        {
            IntegerVector eta_ids_logeir = eta_ids[eta_ids<(eta.ngp+eta.nhp+eta.nlogeir)];
            eta_ids_logeir = eta_ids_logeir - eta.ngp - eta.nhp;
            IntegerVector eta_ids_ft = eta_ids[eta_ids>=(eta.ngp+eta.nhp+eta.nlogeir)];
            eta_ids_ft = eta_ids_ft - (eta.ngp+eta.nhp+eta.nlogeir);
            IntegerVector eta_ids_site_ids;
            IntegerVector eta_ids_study_ids;
            for(int j = 0; j < alldata.nsites; ++j)
            {
                if(is_true(any(eta_ids_logeir==alldata.sites[j].logeir_id))||is_true(any(eta_ids_ft==alldata.sites[j].ft_id)))
                {
                    eta_ids_site_ids.push_back(j);
                    eta_ids_study_ids.push_back(alldata.sites[j].study_id);
                }
            }
            update_blocks_site_ids[i] = unique(eta_ids_site_ids);
            update_blocks_study_ids[i] = unique(eta_ids_study_ids);
        }
    }

    // some useful variables to control the MCMC process
    double apopt_mul = 0.234;
    double apopt_uni = 0.44;
    int winsize = 100;
    int n0 = 50;

    // store the paths with arma::mat
    int neta = update_eta_ids.size();
    arma::mat eta_mcmc(neta, nrep);
    arma::mat accept_mcmc(nblocks, nrep);
    NumericVector ar_mcmc(nblocks, 0.0); // only record the current ar(winsize)
    int tune_iter = nrep + 1;
    int tune_done = 0;
    int update_proposal = 1;
    int num_good_proposal = 0;
    NumericVector accept_sum(nblocks);

    // initial proposal
    double ll = alldata.loglik(update_eta_ids, eta, seq_len(alldata.nsites)-1, alldata.study_ids, pop, ghnodes1, ghweights1, true);
    eta_mcmc.col(0) = as<arma::colvec>(eta.eta[update_eta_ids]);
    accept_mcmc.col(0) = arma::ones(nblocks);

    ar_mcmc.fill(1.0);

    accept_sum.fill(1.0);

    alldata.update_data(seq_len(alldata.nsites)-1);
    alldata.update_loglik_gh();


    // start the MCMC
    NumericVector eta_copy;
    double ll_test, ap;
    int nar0 = 0; // number of all-zero acceptance ratio
    IntegerVector eta_ids;// = seq_len(neta)-1;
    /*// write to file
    std::ofstream myfile("eta_mcmc.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < neta; ++i)
        {
            myfile << eta_mcmc.col(0)[i] << "\t";
        }
        myfile << "\n";
        myfile.close();
    }
    else
    {
        Rcout << "Unable to open file\n";
    }

    myfile.open("accept_mcmc.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < nblocks; ++i)
        {
            myfile << accept_mcmc.col(0)[i] << "\t";
        }
        myfile << "\n";
        myfile.close();
    }
    else
    {
        Rcout << "Unable to open file\n";
    }
    */

    for(int iter = 1; iter < nrep; ++iter)
    {
        // if tune not finished, update the proposal every winsize iterations
        if((update_proposal==1) && (iter % winsize ==0))
        {
             proposal_Sigma_blocks = Sigma_blocks;
             proposal_Mu_blocks = Mu_blocks;
             proposal_lambda_blocks = clone(lambda_blocks);
        }

        // work on each updating block
        for(int i = 0; i < nblocks; ++i)
        {
            eta_ids = update_blocks[i];
            eta_copy = clone(eta.eta); // make a copy in case rejected
            arma::mat Sig_i = proposal_lambda_blocks[i]*proposal_Sigma_blocks[i];
            eta.update(eta_ids, Sig_i);
            ll_test = alldata.loglik(eta_ids, eta, update_blocks_site_ids[i], update_blocks_study_ids[i], pop, ghnodes1, ghweights1, false);
            ap = (ll_test - ll > 0)? 1 : exp(ll_test - ll);
            if(runif(1)[0] < ap)
            {
                // accept the proposal
                ll = ll_test;
                accept_mcmc.col(iter)[i] = 1;
                alldata.update_data(update_blocks_site_ids[i]);
                alldata.update_loglik_gh();
            }
            else
            {
                // reject the proposal
                eta.eta = clone(eta_copy);
                accept_mcmc.col(iter)[i] = 0;
            }

            // smoothing acceptance ratio over the last winsize iterations
            accept_sum[i] = accept_sum[i] + accept_mcmc.col(iter)[i] - ((iter > winsize - 1)? accept_mcmc.col(iter - winsize)[i] : 0);
            ar_mcmc[i] = accept_sum[i]/((iter > winsize - 1)? winsize : (iter + 1));

            // update the Sigma and Mu for the proposal distributions if necessary
            // using the generated sample information to approximate the posterior adaptively
            if((tune_done == 0) && iter > n0)
            {
                double stepsize = 1.0/((iter > 200)? (200+(iter - 200)/20.0) : (iter));
                if(update_blocks[i].size()==1)
                {
                    Sigma_blocks[i] = exp(log(Sigma_blocks[i])+stepsize*((accept_mcmc.col(iter)[i]==1)? 1/apopt_uni : -1/(1-apopt_uni)));
                }
                else
                {
                    lambda_blocks[i] = exp(log(lambda_blocks[i])+stepsize*(ap - apopt_mul));
                    arma::vec temp = as<arma::vec>(eta.eta[eta_ids]) - Mu_blocks[i];
                    Mu_blocks[i] = Mu_blocks[i] + stepsize*temp;
                    Sigma_blocks[i] = Sigma_blocks[i] + stepsize*(temp*temp.t() - Sigma_blocks[i]);
                }
            }

            if((iter+1) % 100==0)
            {
                Rcout << (iter+1) << "-th iteration - tune done: " << tune_done << " - block - " << i <<
                " with current MA(" << winsize << ") acceptance ratio = " << ar_mcmc[i] << "\n";
            }
        }
        // store the eta path
        eta_mcmc.col(iter) = as<arma::vec>(eta.eta[update_eta_ids]);

        // Check whether the tuning process for proposal distribution has finished
        if((tune_done == 0) && iter > 299 && ((iter+1) % winsize == 0))
        {
            // see if the current proposal is good depending on how far the current acceptance ratio (ar_mcmc)
            // is away from the optimal ones
            int current_proposal_is_good = 1;
            for(int i = 0; i < nblocks; ++i)
            {
                if(std::abs(ar_mcmc[i] - ((update_blocks[i].size()==1)? apopt_uni : apopt_mul)) > ((update_blocks[i].size()==1)? 0.16 : 0.08))
                {
                    current_proposal_is_good = 0;
                    Rcout << "block - " << i << " failed the test... \n";
                    break;
                }
            }

            // if the current proposal is good - stop updating the proposal, increase the number of successful consecutive tests by 1,
            // and move onto the next 100 iterations for double check
            if(current_proposal_is_good == 1)
            {
                update_proposal = 0;
                num_good_proposal += 1;
                Rcout << "Current proposal passed " << num_good_proposal << " consecutive tests\n";
            }
            // if the current proposal is not good - continue updating the proposal and re-set the number of successful consecutive tests to 0
            else
            {
                update_proposal = 1;
                num_good_proposal = 0;
                Rcout << "Current proposal passed " << num_good_proposal << " consecutive tests\n";
            }
            // if the proposal has passed 5 successful consecutive tests, set the tuning procedure as "done" and stop updating the proposal any more
            if(num_good_proposal == 5)
            {
                tune_done = 1;
                update_proposal = 0;
                tune_iter = iter;
            }
        }

        /*// write to file
        myfile.open("eta_mcmc.txt", std::ios::app);
        if (myfile.is_open())
        {
            for(int i = 0; i < neta; ++i)
            {
                myfile << eta_mcmc.col(iter)[i] << "\t";
            }
            myfile << "\n";
            myfile.close();
        }
        else
        {
            Rcout << "Unable to open file\n";
        }

        myfile.open("accept_mcmc.txt", std::ios::app);
        if (myfile.is_open())
        {
            for(int i = 0; i < nblocks; ++i)
            {
                myfile << accept_mcmc.col(iter)[i] << "\t";
            }
            myfile << "\n";
            myfile.close();
        }
        else
        {
            Rcout << "Unable to open file\n";
        }
        */

        // If it happens that all the proposals are rejected for many iterations
        // There might be some problem with the code or getting stuck at some local optima
        // Return and check
        if(is_true(all(ar_mcmc == 0)))
        {
            nar0 += 1;
            if(nar0 == 100) // consecutive 100 ar all zeros... return
            {
                return(List::create(_["eta_mcmc"]=eta_mcmc, _["accept_mcmc"]=accept_mcmc, _["tune_iter"]=tune_iter,
                        _["ar_mcmc"]=ar_mcmc, _["Sigma"]=Sigma_blocks, _["lambda"]=lambda_blocks));
            }
        }
        else
        {
            nar0 = 0;
        }
    }


    return(List::create(_["eta_mcmc"]=eta_mcmc, _["accept_mcmc"]=accept_mcmc, _["tune_iter"]=tune_iter,
                        _["Sigma"]=Sigma_blocks, _["lambda"]=lambda_blocks));
}

///////////////////////////////////////////////////////////////////////////////////////////
////// Continue the MCMC fitting if necessary starting from the parameters given by eta0
////// Proposal distributions are given and do not need to be tuned as in runMCMC

// [[Rcpp::export]]
List contMCMC(int nrep, DataFrame& data_key, List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0,
          const NumericVector& ghnodes1, const NumericVector& ghweights1, const std::vector<IntegerVector>& update_blocks,
          const NumericVector& eta0, const std::vector<arma::mat>& Sigma_blocks, const NumericVector& lambda_blocks)
{
    // read in the data
    all_data alldata(data_key, datasets, ghnodes1, ghweights1);
    // get a population object to work on
    population pop(age0, ghnodes0, ghweights0);

    // set up the updating blocks
    int nblocks = update_blocks.size();
    //std::vector<IntegerVector> update_blocks;
    //std::vector<arma::mat> Sigma_blocks;
    //std::vector<arma::vec> Mu_blocks;
    //NumericVector lambda_blocks;

    IntegerVector update_eta_ids;
    for(int i = 0; i < nblocks; ++i)
    {
        for(int j = 0; j < update_blocks[i].size(); ++j)
        {
            update_eta_ids.push_back(update_blocks[i][j]); ///////// NEEDS TO BE CHANGED!!!
        }
    }
    // initialise the eta parameters
    Eta eta(data_key, false, update_eta_ids);

    std::vector<arma::mat> proposal_Sigma_blocks = Sigma_blocks;
    //std::vector<arma::vec> proposal_Mu_blocks = Mu_blocks;
    NumericVector proposal_lambda_blocks = clone(lambda_blocks);

    // update blocks the associated site or study ids
    std::vector<IntegerVector> update_blocks_site_ids(nblocks);
    std::vector<IntegerVector> update_blocks_study_ids(nblocks);


    for(int i = 0; i < nblocks; ++i)
    {
        IntegerVector eta_ids = update_blocks[i];
        IntegerVector eta_ids_logeir = eta_ids[eta_ids<(eta.ngp+eta.nhp+eta.nlogeir)];
        eta_ids_logeir = eta_ids_logeir - eta.ngp - eta.nhp;
        IntegerVector eta_ids_ft = eta_ids[eta_ids>=(eta.ngp+eta.nhp+eta.nlogeir)];
        eta_ids_ft = eta_ids_ft - (eta.ngp+eta.nhp+eta.nlogeir);
        IntegerVector eta_ids_site_ids;
        IntegerVector eta_ids_study_ids;
        for(int j = 0; j < alldata.nsites; ++j)
        {
            if(is_true(any(eta_ids_logeir==alldata.sites[j].logeir_id))||is_true(any(eta_ids_ft==alldata.sites[j].ft_id)))
            {
                eta_ids_site_ids.push_back(j);
                eta_ids_study_ids.push_back(alldata.sites[j].study_id);
            }
        }
        update_blocks_site_ids[i] = unique(eta_ids_site_ids);
        update_blocks_study_ids[i] = unique(eta_ids_study_ids);
    }

    // some useful variables to control the MCMC process
    //double apopt_mul = 0.234;
    //double apopt_uni = 0.44;
    int winsize = 100;

    // store the paths with arma::mat
    int neta = update_eta_ids.size();
    arma::mat eta_mcmc(neta, nrep);
    arma::mat accept_mcmc(nblocks, nrep);
    //arma::mat ar_mcmc(nblocks, nrep);
    //std::vector<NumericVector> ar_mcmc(nrep);
    NumericVector ar_mcmc(nblocks, 0.0); // only record the current ar(winsize)
    IntegerVector tune_iters(nblocks, 0);
    //IntegerVector tune_done(nblocks, 0);
    //int tune_done = 0;
    NumericVector accept_sum(nblocks);

    // initial proposal

    eta.eta[update_eta_ids] = clone(eta0);

    double ll = alldata.loglik(update_eta_ids, eta, seq_len(alldata.nsites)-1, alldata.study_ids, pop, ghnodes1, ghweights1, true);
    eta_mcmc.col(0) = as<arma::colvec>(eta.eta[update_eta_ids]);
    accept_mcmc.col(0) = arma::ones(nblocks);
    //ar_mcmc.col(0) = arma::ones(nblocks);
    ar_mcmc.fill(1.0);

    accept_sum.fill(1.0);

    alldata.update_data(seq_len(alldata.nsites)-1);
    alldata.update_loglik_gh();

    // start the MCMC
    NumericVector eta_copy;
    double ll_test, ap;
    int nar0 = 0; // number of all-zero acceptance ratio
    IntegerVector eta_ids; //= seq_len(neta)-1;

    proposal_Sigma_blocks = Sigma_blocks;
    proposal_lambda_blocks = clone(lambda_blocks);
    for(int iter = 1; iter < nrep; ++iter)
    {
        // work on each updating block
        for(int i = 0; i < nblocks; ++i)
        {
            eta_ids = update_blocks[i];
            eta_copy = clone(eta.eta); // make a copy in case rejected
            arma::mat Sig_i = proposal_lambda_blocks[i]*proposal_Sigma_blocks[i];
            eta.update(eta_ids, Sig_i);
            ll_test = alldata.loglik(eta_ids, eta, update_blocks_site_ids[i], update_blocks_study_ids[i], pop, ghnodes1, ghweights1, false);
            ap = (ll_test - ll > 0)? 1 : exp(ll_test - ll);
            if(runif(1)[0] < ap)
            {
                // accept the proposal
                ll = ll_test;
                accept_mcmc.col(iter)[i] = 1;
                alldata.update_data(update_blocks_site_ids[i]);
                alldata.update_loglik_gh();
            }
            else
            {
                // reject the proposal
                eta.eta = clone(eta_copy);
                accept_mcmc.col(iter)[i] = 0;
            }

            // smoothing acceptance ratio over the last winsize iterations
            accept_sum[i] = accept_sum[i] + accept_mcmc.col(iter)[i] - ((iter > winsize - 1)? accept_mcmc.col(iter - winsize)[i] : 0);
            ar_mcmc[i] = accept_sum[i]/((iter > winsize - 1)? winsize : (iter + 1));

            if((iter+1) % 100==0)
            {
                Rcout << (iter+1) << "-th iteration - block - " << i <<
                " with current MA(" << winsize << ") acceptance ratio = " << ar_mcmc[i] << "\n";
            }
        }
        // store the eta path
        eta_mcmc.col(iter) = as<arma::vec>(eta.eta[update_eta_ids]);

        if(is_true(all(ar_mcmc == 0)))
        {
            nar0 += 1;
            if(nar0 == 100) // consecutive 100 ar all zeros... return
            {
                return(List::create(_["eta_mcmc"]=eta_mcmc, _["accept_mcmc"]=accept_mcmc, _["tune_iters"]=tune_iters,
                        _["ar_mcmc"]=ar_mcmc, _["Sigma"]=Sigma_blocks, _["lambda"]=lambda_blocks));
            }
        }
        else
        {
            nar0 = 0;
        }
    }

    return(List::create(_["eta_mcmc"]=eta_mcmc, _["accept_mcmc"]=accept_mcmc, //_["tune_iters"]=tune_iters,
                        _["Sigma"]=Sigma_blocks, _["lambda"]=lambda_blocks));
}

///////////////////////////////////////////////////////////////////////////////////////////
////// Evaluate the log-likelihood contributed from the data (no parameters prior)
////// at the parameters given by eta0

// [[Rcpp::export]]
List runloglik(DataFrame& data_key, List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0,
          const NumericVector& ghnodes1, const NumericVector& ghweights1, const NumericVector& eta0, const IntegerVector& eta_ids)
{
    // read in the data
    all_data alldata(data_key, datasets, ghnodes1, ghweights1);
    // get a population object to work on
    population pop(age0, ghnodes0, ghweights0);
    // initialise the eta parameters
    Eta eta(data_key, false, eta_ids);
    eta.eta[eta_ids] = clone(eta0);

    //int neta = eta.ngp + eta.nlogeir + eta.nft;

    //IntegerVector eta_ids = seq_len(neta)-1;
    double ll = alldata.loglik(eta_ids, eta, seq_len(alldata.nsites)-1, alldata.study_ids, pop, ghnodes1, ghweights1, true);

    return List::create(_["ll"]=ll, _["ll_c_gh"]=alldata.loglik_inc_gh_copy, _["ll_p_gh"]=alldata.loglik_prev_gh_copy, _["ll_eta"]=eta.loglik(), _["ll_s_gh"]=alldata.loglik_studies_copy);
}

///////////////////////////////////////////////////////////////////////////////////////////
////// Print all the data information that are input to an all_data object
////// to see if the creation is correct and complete

// [[Rcpp::export]]
void check_alldata(DataFrame& data_key, List& datasets,
          const NumericVector& ghnodes1, const NumericVector& ghweights1)
{
    // read in the data
    all_data alldata(data_key, datasets, ghnodes1, ghweights1);
    for(int i = 0; i < alldata.nsites; ++i)
    {
        site_data sitei = alldata.sites[i];
        sitei.print_out();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
////// MCMC fitting for the study level random effect u at study - studyid
////// conditional upon the other varying parameters given by eta0
////// NOTE: the site-age random effect v_jk is integrated out analytically from the likelihood
//////       whose posterior can be analytically obtained conditional upon
//////       all the other parameters and study level random effect u

// [[Rcpp::export]]
List runMCMC_for_u(int nrep, int studyid, DataFrame& data_key, List& datasets, const arma::vec& age0, const arma::vec& ghnodes0, const arma::vec& ghweights0,
          const IntegerVector& update_eta_ids, const NumericVector& eta0)
{
    // run MCMC procedure for the study level random effect u for the study with studyid conditional upon fixed parameters eta0
    // read in the data
    all_data alldata(data_key, datasets);
    // get a population object to work on
    population pop(age0, ghnodes0, ghweights0);

    // initialise the eta and parameters
    Eta eta(data_key, false, update_eta_ids);
    eta.eta[update_eta_ids] = clone(eta0);
    parameters par = eta.get_parameters();

    double sigmac = exp(eta.eta[23]);

    // some useful variables to control the MCMC process
    double apopt_uni = 0.44;
    int winsize = 100;
    int n0 = 50;

    // store the paths with vectors
    NumericVector u_mcmc(nrep, 0.0);
    IntegerVector accept_mcmc(nrep, 0);
    NumericVector ll_mcmc(nrep, 0.0);
    NumericVector sigma_mcmc;
    double ar_mcmc = 0.0; // only record the current ar(winsize)
    int tune_done = 0;
    int update_proposal = 1;
    int num_good_proposal = 0;
    double accept_sum = 0;
    int tune_iter = 0;

    // prepare some quantities that are required repeatedly in likelihood calculation
    LogicalVector studyid_in_study_ids = (alldata.study_ids==studyid);
    int ind_study = which(studyid_in_study_ids)[0];
    IntegerVector studyid_site_ids = alldata.study_site_ids[ind_study];
    NumericVector components1;
    NumericVector components2;

    double a = 1/exp(eta.eta[24]); // eta[24]: alphac
    double b = a;

    for(IntegerVector::iterator it = studyid_site_ids.begin(); it != studyid_site_ids.end(); ++it)
    {
        double eir = exp(eta.eta[eta.ngp + eta.nhp + alldata.sites[*it].logeir_id]); // eir for the site_id
        double ft = t2(eta.eta[eta.ngp + eta.nhp + eta.nlogeir + alldata.sites[*it].ft_id]);

        // find the equilibrium status for this site
        pop.set_equilibrium(eir, ft, par);

        // store the likelihood components
        LogicalVector in_groupk;
        IntegerVector k_ind;
        double comp1 = 0.0;
        double comp2 = 0.0;
        double Tji, yji, muji, r;
        if(alldata.sites[*it].case_detection == 2)
        {
            r = 1/(1+exp(-eta.eta[21]));  // drW
        }
        else if(alldata.sites[*it].case_detection == 6)
        {
            r = 1/(1+exp(-eta.eta[22]));  // drP
        }
        else
        {
            r = 1.0;
        }
        for(int k = 0; k < 4; ++k)
        {
            // k for age group
            in_groupk = (alldata.sites[*it].inc.clin_group == k);
            k_ind = which(in_groupk);
            if(k_ind.size()!=0)
            {
                comp1 = 0.0;
                comp2 = 0.0;
                for(IntegerVector::iterator it1 = k_ind.begin(); it1 != k_ind.end(); ++it1)
                {
                    Tji = alldata.sites[*it].inc.clin_p_years[*it1];
                    yji = alldata.sites[*it].inc.clin_k[*it1];
                    muji = pop.get_inc(alldata.sites[*it].inc.age0[*it1], alldata.sites[*it].inc.age1[*it1]);
                    comp1 += yji;
                    comp2 += r*Tji*muji;

                }
                components1.push_back(comp1);
                components2.push_back(comp2);
            }
        }
    }

    // MCMC sampling starts from here...
    double update_sigma = 1.0;//sigma0; // initial sigma for normal proposal
    double proposal_sigma = update_sigma;
    double ll, ll_test, ap;
    // initial sample
    double u, u_test;
    u = R::rnorm(0, proposal_sigma);
    ll = -pow(u+0.5*sigmac*sigmac, 2.0)/(2*sigmac*sigmac);
    for(int l = 0; l < components1.size(); ++l)
    {
        ll += (-(a+components1[l])*log(b + exp(u)*components2[l])+u*components1[l]);
    }
    u_mcmc[0] = u;
    accept_mcmc[0] = 1;
    accept_sum += 1.0;
    ll_mcmc[0] = ll;

    for(int iter = 1; iter < nrep; ++iter)
    {
        if((update_proposal==1) && (iter % winsize ==0))
        {
             proposal_sigma = update_sigma;
        }
        // make a proposal
        u_test = R::rnorm(u, proposal_sigma);
        ll_test = -pow(u_test+0.5*sigmac*sigmac, 2.0)/(2*sigmac*sigmac);

        for(int l = 0; l < components1.size(); ++l)
        {
            ll_test += (-(a+components1[l])*log(b + exp(u_test)*components2[l])+u_test*components1[l]);
        }
        ll_mcmc[iter] = ll_test;
        ap = (ll_test - ll > 0)? 1 : exp(ll_test - ll);
        if(runif(1)[0] < ap)
        {
            // accept the proposal
            u = u_test;
            ll = ll_test;
            accept_mcmc[iter] = 1;
        }
        else
        {
            // reject the proposal
            accept_mcmc[iter] = 0;
        }
        accept_sum = accept_sum + accept_mcmc[iter] - ((iter > winsize - 1)? accept_mcmc[iter - winsize] : 0);
        ar_mcmc = accept_sum/((iter > winsize - 1)? winsize : (iter + 1));

        // store the path
        u_mcmc[iter] = u;

        if((iter+1) % 100==0)
        {
            Rcout << (iter+1) << "-th iteration - tune done: " << tune_done << " and current MA(" << winsize << ") acceptance ratio = " << ar_mcmc << "\n";
        }

        // update the proposal if needed to achieve reasonable acceptance ratio
        if((tune_done == 0) && iter > n0)
        {
            double stepsize = 1.0/((iter > 200)? (200+(iter - 200)/20.0) : (iter));
            update_sigma = exp(log(update_sigma)+stepsize*((accept_mcmc[iter]==1)? 1/apopt_uni : -1/(1-apopt_uni)));
            sigma_mcmc.push_back(update_sigma);
        }
        // see if the tuning has done
        if((tune_done == 0) && iter > 299 && ((iter+1) % winsize == 0))
        {
            // see if the current proposal is good
            int current_proposal_is_good = 1;
            if(std::abs(ar_mcmc - apopt_uni) > 0.16)
            {
                current_proposal_is_good = 0;
                Rcout << "current proposal failed the test... \n";
            }

            // if the current proposal is good - stop updating the proposal and increase the number of successful consecutive tests by 1
            if(current_proposal_is_good == 1)
            {
                update_proposal = 0;
                num_good_proposal += 1;
                Rcout << "Current proposal passed " << num_good_proposal << " consecutive tests\n";
            }
            // if the current proposal is not good - continue updating the proposal and re-set the number of successful consecutive tests to 0
            else
            {
                update_proposal = 1;
                num_good_proposal = 0;
                Rcout << "Current proposal passed " << num_good_proposal << " consecutive tests\n";
            }
            // if there are 5 successful consecutive tests, set the tuning procedure to be done and stop updating the proposal any more
            if(num_good_proposal == 5)
            {
                tune_done = 1;
                update_proposal = 0;
                tune_iter = iter;
            }
        }

    }
    return(List::create(_["u_mcmc"]=u_mcmc, _["accept_mcmc"]=accept_mcmc, _["tune_iter"]=tune_iter, _["sigma"]=proposal_sigma, _["ll"]=ll_mcmc, _["sigma_path"]=sigma_mcmc));
}
