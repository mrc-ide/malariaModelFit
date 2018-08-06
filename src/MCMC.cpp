/*
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// some simple transformation functions
double t1(double eta) {
    // eta to theta for theta > 0
    return(exp(eta));
}

double t2(double eta) {
    // eta to theta for 0 < theta < 1
    return(1/(1+exp(-eta)));
}

double tb1(double theta) {
    // theta to eta for theta > 0
    return(log(theta));
}

double tb2(double theta) {
    // theta to eta for 0 < theta < 1
    return(log(theta/(1-theta)));
}

double dt1(double eta) {
    // derivative of t1
    return(exp(eta));
}

double dt2(double eta) {
    // derivative of t2
    return(exp(eta)/pow(1+exp(eta),2));
}

double ldt1(double eta) {
    // log of the derivative of t1
    return(eta);
}

double ldt2(double eta) {
    // log of the derivative of t2
    return(eta+2*log(1+exp(eta)));
}

//------------------------------------------------
// Eta::
// constructur without initialising eta values
Eta::Eta(const Rcpp::DataFrame& data_key) {
    ngp = 21;
    nhp = 6;
    logeir_mean = Rcpp::as<Rcpp::NumericVector>(data_key[2]);
    logeir_sd = Rcpp::as<Rcpp::NumericVector>(data_key[3]);
    nlogeir = logeir_mean.size();
    Rcpp::IntegerVector place_id = Rcpp::as<Rcpp::IntegerVector>(data_key[4]);
    //IntegerVector uni_place_id = unique(place_id); // TODO - Bob delete if not needed?
    nft =  unique(place_id).size();
    Rcpp::NumericVector fa = Rcpp::as<Rcpp::NumericVector>(data_key[5]);
    Rcpp::NumericVector fb = Rcpp::as<Rcpp::NumericVector>(data_key[6]);
    //nplace = // TODO - Bob delete if not needed?
    ft_a = Rcpp::NumericVector(17);
    ft_b = Rcpp::NumericVector(17);
    Rcpp::NumericVector fai, fbi;
    for (int i = 0; i < nft; ++i) {
        fai = fa[place_id == i];
        fbi = fb[place_id == i];
        ft_a[i] = fai[0];
        ft_b[i] = fbi[0];
    }
}

//------------------------------------------------
// Eta::
// constructor to initialise eta values. those with indices in init_eta_ids will be initialised randomly if r=true or pre-specified default values if r=false
Eta::Eta(const Rcpp::DataFrame& data_key, bool r, const Rcpp::IntegerVector& init_eta_ids) {
    ngp = 21;
    nhp = 6;
    logeir_mean = Rcpp::as<Rcpp::NumericVector>(data_key[2]);
    logeir_sd = Rcpp::as<Rcpp::NumericVector>(data_key[3]);
    nlogeir = logeir_mean.size();
    Rcpp::IntegerVector place_id = Rcpp::as<Rcpp::IntegerVector>(data_key[4]);
    //IntegerVector uni_place_id = unique(place_id); // TODO - Bob delete if not needed?
    nft = unique(place_id).size();
    Rcpp::NumericVector fa = Rcpp::as<Rcpp::NumericVector>(data_key[5]);
    Rcpp::NumericVector fb = Rcpp::as<Rcpp::NumericVector>(data_key[6]);
    ft_a = Rcpp::NumericVector(17);
    ft_b = Rcpp::NumericVector(17);
    Rcpp::NumericVector fai, fbi;
    
    for(int i = 0; i < nft; ++i) {
        fai = fa[place_id == i];
        fbi = fb[place_id == i];
        ft_a[i] = fai[0];
        ft_b[i] = fbi[0];
    }
    eta = Rcpp::NumericVector(ngp+nhp+nlogeir+nft);
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
    for (int i = 0; i < nlogeir; ++i) {
        eta[ngp+nhp+i] = logeir_mean[i];
    }
    
    // ft
    for (int i = 0; i < nft; ++i) {
        eta[ngp+nhp+nlogeir+i] = tb2(ft_a[i]/(ft_a[i]+ft_b[i]));
    }
    
    // if random initial values for those with init_eta_ids are required with r=true
    if (r) {
        for (Rcpp::IntegerVector::const_iterator it = init_eta_ids.begin(); it != init_eta_ids.end(); it++) {
            eta[*it] = Rcpp::rnorm(1,0.0,2.0)[0];
        }
    }
}

//------------------------------------------------
// inc_data::
// default constructor
inc_data::inc_data(){};

// informed constructor
inc_data::inc_data(Rcpp::IntegerVector& clink, Rcpp::NumericVector& clinpy, Rcpp::NumericVector& a0, Rcpp::NumericVector& a1, Rcpp::IntegerVector& clingroup, Rcpp::NumericVector& clinmu) {
    
    clin_k = clink;
    clin_p_years = clinpy;
    age0 = a0;
    age1 = a1;
    clin_group = clingroup;
    clin_mu = clinmu;
    clin_mu_copy = clone(clin_mu);
}

//------------------------------------------------
// inc_data::
// print out an inc_data object
void inc_data::print_out() {
    
    Rcpp::Rcout << setw(8) << "clin_k" << setw(8) << "clin_py" << setw(8) << "age0" << setw(8) << "age1" << setw(8) << "clin_gr" << "\n";
    for (int i = 0; i < clin_k.size(); ++i) {
        Rcpp::Rcout << setw(8) << clin_k[i] << setw(8) << clin_p_years[i] << setw(8) << age0[i] << setw(8) << age1[i] << setw(8) << clin_group[i] << "\n";
    }
}

//------------------------------------------------
// prev_data::
// default constructor
prev_data::prev_data(){};

// informed constructor
prev_data::prev_data(Rcpp::IntegerVector& slidek, Rcpp::IntegerVector& sliden, Rcpp::IntegerVector& pcrk, Rcpp::NumericVector& a0, Rcpp::NumericVector& a1, Rcpp::NumericVector& slidemu, Rcpp::NumericVector& pcrmu) {
    
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

//------------------------------------------------
// prev_data::
// print out a prev_data object
void prev_data::print_out() {
    
    Rcpp::Rcout << setw(8) << "slide_k" << setw(8) << "slide_n" << setw(8) << "pcr_k" << setw(8) << "age0" << setw(8) << "age1" << "\n";
    for (int i = 0; i < slide_k.size(); ++i) {
        Rcpp::Rcout << setw(8) << slide_k[i] << setw(8) << slide_n[i] << setw(8) << pcr_k[i] << setw(8) << age0[i] << setw(8) << age1[i] << "\n";
    }
}

//------------------------------------------------
// site_data::
// constructor to create a site_data object from user-input values
site_data::site_data(inc_data& incdata, prev_data& prevdata, int siteid, int studyid, int logeirid, int ftid, int casedetection, bool forinc, bool forprev) {
    
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

//------------------------------------------------
// site_data::
// print out a site_data object
void site_data::print_out()
{
    Rcpp::Rcout << "Site - " << site_id << " in Study - " << study_id << " data \n";
    Rcpp::Rcout << "logeir_id: " << logeir_id << "ft_id: " << ft_id << "case detection: " << case_detection << "\n";
    if (for_inc) {
        Rcpp::Rcout << "    inc data:\n";
        inc.print_out();
    }
    if (for_prev) {
        Rcpp::Rcout << "    prev data:\n";
        prev.print_out();
    }
}

//------------------------------------------------
// all_data::
// constructor to create an all_data object from user-input values
all_data::all_data(Rcpp::DataFrame& data_key, Rcpp::List& datasets, const Rcpp::NumericVector& ghnodes, const Rcpp::NumericVector& ghweights) {
    // TODO - Bob, errors that need looking into
    
    nsites = datasets.size();
    Rcpp::IntegerVector cd = Rcpp::as<Rcpp::IntegerVector>(data_key[7]); // case detection id
    
    Rcpp::IntegerVector si = Rcpp::as<Rcpp::IntegerVector>(data_key[8]); // study id
    study_ids = unique(si);
    nstudies = study_ids.size();
    Rcpp::IntegerVector pi = Rcpp::as<Rcpp::IntegerVector>(data_key[4]); // place id for ft
    study_site_ids = std::vector<Rcpp::IntegerVector>(nstudies);
    int ngh = ghnodes.size();
    
    for (int i = 0; i < nsites; ++i) {
        // site i
        int siteid_i = i;
        int logeirid_i = i;
        int ftid_i = pi[i];
        int casedetection_i = cd[i];
        int studyid_i = si[i];
        Rcpp::LogicalVector temp = study_ids==studyid_i;
        int temp_id = which(temp)[0];
        study_site_ids[temp_id].push_back(siteid_i);
        Rcpp::DataFrame data_i = datasets[i];
        Rcpp::IntegerVector sk = Rcpp::as<Rcpp::IntegerVector>(data_i[0]);
        Rcpp::IntegerVector sn = Rcpp::as<Rcpp::IntegerVector>(data_i[1]);
        Rcpp::IntegerVector pk = Rcpp::as<Rcpp::IntegerVector>(data_i[2]);
        Rcpp::IntegerVector ck = Rcpp::as<Rcpp::IntegerVector>(data_i[3]);
        Rcpp::NumericVector cpy = Rcpp::as<Rcpp::NumericVector>(data_i[4]);
        Rcpp::NumericVector a0 = Rcpp::as<Rcpp::NumericVector>(data_i[7]);
        Rcpp::NumericVector a1 = Rcpp::as<Rcpp::NumericVector>(data_i[8]);
        Rcpp::IntegerVector cg = Rcpp::as<Rcpp::IntegerVector>(data_i[9]);
        inc_data incdatai;
        prev_data prevdatai;
        bool forinci;
        bool forprevi;
        
        // clinical data
        if (is_true(all(ck == -1))) {
            forinci = false;
            //incdatai = inc_data();
        }
        else {
            forinci = true;
            Rcpp::LogicalVector temp = (ck!=-1);
            Rcpp::IntegerVector cki = ck[temp];
            Rcpp::NumericVector cpyi = cpy[temp];
            Rcpp::NumericVector a0i = a0[temp];
            Rcpp::NumericVector a1i = a1[temp];
            Rcpp::IntegerVector cgi = cg[temp];
            Rcpp::NumericVector ckmui(cki.size(), 0.0);
            incdatai = inc_data(cki, cpyi,a0i,a1i,cgi,ckmui);
        }
        // prevalence data
        if (is_true(all(sn == -1))) {
            forprevi = false;
            //prevdatai = prev_data();
        }
        else {
            forprevi = true;
            Rcpp::LogicalVector temp = (sn!=-1);
            Rcpp::IntegerVector ski = sk[temp];
            Rcpp::IntegerVector sni = sn[temp];
            Rcpp::IntegerVector pki = pk[temp];
            Rcpp::NumericVector a0i = a0[temp];
            Rcpp::NumericVector a1i = a1[temp];
            Rcpp::NumericVector skmui(ski.size(), 0.0);
            Rcpp::NumericVector pkmui(ski.size(), 0.0);
            prevdatai = prev_data(ski, sni, pki, a0i, a1i, skmui, pkmui);
        }
        site_data sitei = site_data(incdatai, prevdatai, siteid_i, studyid_i, logeirid_i, ftid_i, casedetection_i, forinci, forprevi);
        sites.push_back(sitei);
        Rcpp::NumericVector loglik_inc_gh_i(ngh, 0.0);
        Rcpp::NumericVector loglik_prev_gh_i(ngh, 0.0);
        loglik_inc_gh.push_back(loglik_inc_gh_i);
        loglik_prev_gh.push_back(loglik_prev_gh_i);
    }
    
    loglik_studies = Rcpp::NumericVector(nstudies, 0.0);
    
    loglik_inc_gh_copy = loglik_inc_gh;
    loglik_prev_gh_copy = loglik_prev_gh;
    loglik_studies_copy = Rcpp::NumericVector(nstudies, 0.0);
    
}
*/