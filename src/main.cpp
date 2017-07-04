
#include "main.h"
#include "misc.h"

using namespace std;

// [[Rcpp::export]]
Rcpp::List dummy1_cpp() {

    print("Rcpp link working!");

    return Rcpp::List::create(Rcpp::Named("foo")=-9);
}
