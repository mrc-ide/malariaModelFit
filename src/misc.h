
#ifndef __MRCfellowship__misc__
#define __MRCfellowship__misc__

#include <RcppArmadillo.h>

//------------------------------------------------
// define very small number for catching underflow problems
#define UNDERFLO   1e-100

//------------------------------------------------
// dummy function to check Rcpp working correctly
void dummy1_cpp();

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
    TYPE output = 0;
    for (int i=0; i<int(x.size()); i++)
        output += x[i];
    return(output);
}

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(std::vector<TYPE> &x) {
    return(sum(x)/double(x.size()));
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB);

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    Rcpp::Rcout << x << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void printVector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << x[i] << " ";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void printMatrix(std::vector< std::vector<TYPE> > &M) {
    for (int i=0; i<M.size(); i++) {
        for (int j=0; j<M[i].size(); j++) {
            Rcpp::Rcout << M[i][j] << " ";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void printArray(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                Rcpp::Rcout << x[i][j][k] << " ";
            }
            Rcpp::Rcout << "\n";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void printStars(std::string title, int n);

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by=1);

//------------------------------------------------
// searches for integer item in vector. Returns true if found, otherwise false
bool item_in_vector(std::vector<int> &vec, int item);

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int Rcpp_to_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double Rcpp_to_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> Rcpp_to_vector_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> Rcpp_to_vector_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
std::vector<std::string> Rcpp_to_vector_string(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector< std::vector<int> > Rcpp_to_mat_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector< std::vector<double> > Rcpp_to_mat_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector< std::vector< std::vector<double> > > Rcpp_to_array_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector< std::vector< std::vector<int> > > Rcpp_to_array_int(Rcpp::List x);

//------------------------------------------------
// convert arma::vec object to NumericVector object
Rcpp::NumericVector tonv(arma::vec& v1);

#endif
