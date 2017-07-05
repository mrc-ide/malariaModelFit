// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// human_equilibrium_fast
Rcpp::List human_equilibrium_fast(double EIR, double ft, const Rcpp::List& p, arma::vec& age0, const arma::vec& ghnodes, const arma::vec& ghweights);
RcppExport SEXP malariaModelFit_human_equilibrium_fast(SEXP EIRSEXP, SEXP ftSEXP, SEXP pSEXP, SEXP age0SEXP, SEXP ghnodesSEXP, SEXP ghweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type EIR(EIRSEXP);
    Rcpp::traits::input_parameter< double >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type age0(age0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ghnodes(ghnodesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ghweights(ghweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(human_equilibrium_fast(EIR, ft, p, age0, ghnodes, ghweights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"malariaModelFit_human_equilibrium_fast", (DL_FUNC) &malariaModelFit_human_equilibrium_fast, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_malariaModelFit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
