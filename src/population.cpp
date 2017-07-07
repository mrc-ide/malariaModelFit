
#include "population.h"

using namespace std;


//------------------------------------------------
// population::
// default constructor for population class
population::population() {}

//------------------------------------------------
// population::
// informed constructor for population class
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

//------------------------------------------------
// population::
// find equilibrium solution given EIR (EIR0), proportion treated (ft) and model parameters. Translated from Jamie's R code.
void population::set_equilibrium(double EIR0, double ft0, const parameters& p)
{
    // define starting parameters
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
