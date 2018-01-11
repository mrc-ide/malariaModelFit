
#ifndef __malariaModelFit__MCMC__
#define __malariaModelFit__MCMC__

#include <RcppArmadillo.h>

//------------------------------------------------
// Eta::
// class containing the variable parameters in the MCMC fitting that could be sampled. theta is defined as the original parameters for the transmission model, EIR and ft, but some theta's have restricted support and need to be transformed to eta, which has infinite support
class Eta {
    
public:
    
    int ngp;         // number of global parameters for equilibrium solution
    int nhp;         // number of other hyper-parameters in the likelihood
    int nlogeir;     // number of site-specific logeir parameters
    int nft;         // number of site-specific ft parameters
    Rcpp::NumericVector eta;   // a vector containing the ngp + nlogeir + nft parameters
    
    Rcpp::NumericVector logeir_mean; // Normal prior mean for the log-EIR
    Rcpp::NumericVector logeir_sd;   // Normal prior sd for the log-EIR
    Rcpp::NumericVector ft_a;        // Beta prior for ft, Beta(a,b)
    Rcpp::NumericVector ft_b;
    
    // constructors
    Eta(const Rcpp::DataFrame& data_key); // without initialising eta values
    Eta(const Rcpp::DataFrame& data_key, bool r, const Rcpp::IntegerVector& init_eta_ids); // initialise eta values. those with indices in init_eta_ids will be initialised randomly if r=true or pre-specified default values if r=false
    
    // method to return the log-likelihood (density) of the prior at current eta values
    //double loglik();
    
    // method to return transmission model parameters corresponding to current eta values
    //parameters get_parameters();
    
    // method to update the values of eta with indices in eta_ids using Metropolis random walk update
    // eta_new = N(eta_current, Sigma)
    //void update(IntegerVector eta_ids, arma::mat& Sigma);
};

//------------------------------------------------
// inc_data::
// class for incidence data at one site
class inc_data {
    
public:
    
    // data
    Rcpp::IntegerVector clin_k; // number of clinical events
    Rcpp::NumericVector clin_p_years; // number of person years for clinical events
    Rcpp::NumericVector age0; // lower end of age group
    Rcpp::NumericVector age1; // upper end of age group
    Rcpp::IntegerVector clin_group; // group for clinical data random effects by age
    
    // for likelihood calculation
    Rcpp::NumericVector clin_mu; // model-predicted incidence of clinical malaria
    Rcpp::NumericVector clin_mu_copy; // a copy of clin_mu for temporary update
    
    // Constructors:
    inc_data();
    inc_data(Rcpp::IntegerVector& clink, Rcpp::NumericVector& clinpy, Rcpp::NumericVector& age0, Rcpp::NumericVector& age1, Rcpp::IntegerVector& clingroup, Rcpp::NumericVector& clinmu);
    
    // method to print out to visualise the data
    void print_out();
};

//------------------------------------------------
// prev_data::
// class for prevalence data at one site
class prev_data {

public:
    
    // data
    Rcpp::IntegerVector slide_k; // number positive by slide (Microscopy)
    Rcpp::IntegerVector slide_n; // number sampled for slide (and for PCR)
    Rcpp::IntegerVector pcr_k; // number positive by PCR
    Rcpp::NumericVector age0; // lower end of age group
    Rcpp::NumericVector age1; // upper end of age group
    
    // for likelihood calculation
    Rcpp::NumericVector slide_mu; // model-predicted prevalence by slide (Microscopy)
    Rcpp::NumericVector pcr_mu; // model-predicted prevalence by PCR
    Rcpp::NumericVector slide_mu_copy; // a copy of slide_mu
    Rcpp::NumericVector pcr_mu_copy; // a copy of pcr_mu
    
    // Constructors:
    prev_data();
    prev_data(Rcpp::IntegerVector& slidek, Rcpp::IntegerVector& sliden, Rcpp::IntegerVector& pcrk,
              Rcpp::NumericVector& age0, Rcpp::NumericVector& age1,
              Rcpp::NumericVector& slidemu, Rcpp::NumericVector& pcrmu);
    
    // method to print out to visualise the data
    void print_out();
};

//------------------------------------------------
// site_data::
// class for the data including incidence and/or prevalence data at one site
class site_data {
    
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
    //void loglik(bool update_data, Eta& eta, parameters& par, population& pop, const NumericVector& ghnodes, const NumericVector& ghweights, std::vector<NumericVector>& loglikincgh, std::vector<NumericVector>& loglikprevgh);
    
    // method to update the model-predicted incidence/prevalence if necessary
    // ONLY when the proposed eta values are accepted
    //void update_data();
    
    // method to print out to visualise the data
    void print_out();
};

//------------------------------------------------
// all_data::
// class for all the data used in Griffin et al. 2014
class all_data {
    
public:
    
    std::vector<site_data> sites; // in the order of site in data_key
    int nstudies; // number of studies
    int nsites;   // number of sites
    Rcpp::IntegerVector study_ids; // study ids
    std::vector<Rcpp::IntegerVector> study_site_ids; // site ids for each study
    
    // a nsites by ngh matrix for the contribution of incidence data to likelihood from each site
    std::vector<Rcpp::NumericVector> loglik_inc_gh;
    // a nsites by ngh matrix for the contribution of prevalence data to likelihood from each site
    std::vector<Rcpp::NumericVector> loglik_prev_gh;
    // a nstudies vector for the contribution to likelihood from each study
    Rcpp::NumericVector loglik_studies;
    
    // copies of the above contribution components
    // this is for the temporary storage since an update may be accepted or rejected
    // hence we need to store both the current and updated likelihood at the same time
    std::vector<Rcpp::NumericVector> loglik_inc_gh_copy;
    std::vector<Rcpp::NumericVector> loglik_prev_gh_copy;
    Rcpp::NumericVector loglik_studies_copy;
    
    // Constructors:
    // initialise the loglik_inc/prev_gh; note that ghnodes and ghweights are different from those in population
    all_data(Rcpp::DataFrame& data_key, Rcpp::List& datasets, const Rcpp::NumericVector& ghnodes, const Rcpp::NumericVector& ghweights);
    
    // Not initialise the Gauss-Hermit relevant matrices;
    // used for inference on u since the likelihood evaluation doesn't need G-H quadrature
    all_data(Rcpp::DataFrame& data_key, Rcpp::List& datasets);
    
    // method to evaluate the log-likelihood contribution from study with index of studyid
    //void loglik_study(int studyid, const NumericVector& ghweights);
    
    // method to evaluate the log-likelihood from all data at the parameters eta
    //double loglik(const IntegerVector& eta_ids, Eta& eta, const IntegerVector& eta_ids_site_ids, const IntegerVector& eta_ids_study_ids, population& pop, const NumericVector& ghnodes1, const NumericVector& ghweights1, bool initialise);
    
    // update the data depending on the acceptance or rejection
    //void update_data(const IntegerVector& eta_ids_site_ids);
    
    // update the loglik_inc_gh, loglik_prev_gh, loglik_studies if necessary
    //void update_loglik_gh();
};

#endif
