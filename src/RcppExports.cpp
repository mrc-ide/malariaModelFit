// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// loglikelihood
SEXP loglikelihood(std::vector<double> params, std::vector<double> x);
RcppExport SEXP _malariaModelFit_loglikelihood(SEXP paramsSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood(params, x));
    return rcpp_result_gen;
END_RCPP
}
// human_equilibrium_cpp
Rcpp::List human_equilibrium_cpp(double EIR, double ft, const Rcpp::List& p, const Rcpp::NumericVector& age_, const Rcpp::List& h);
RcppExport SEXP _malariaModelFit_human_equilibrium_cpp(SEXP EIRSEXP, SEXP ftSEXP, SEXP pSEXP, SEXP age_SEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type EIR(EIRSEXP);
    Rcpp::traits::input_parameter< double >::type ft(ftSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type age_(age_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(human_equilibrium_cpp(EIR, ft, p, age_, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_malariaModelFit_loglikelihood", (DL_FUNC) &_malariaModelFit_loglikelihood, 2},
    {"_malariaModelFit_human_equilibrium_cpp", (DL_FUNC) &_malariaModelFit_human_equilibrium_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_malariaModelFit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
