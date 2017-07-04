
#include "misc.h"

using namespace std;

//------------------------------------------------
// define very small number for catching underflow problems
// DEFINED IN HEADER

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// mean of vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB) {
    if (logA-logB > 100) {
        return(logA);
    } else if (logB-logA > 100) {
        return(logB);
    }
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return(output);
}

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void printStars(string title, int n) {
    Rcpp::Rcout << title << " ";
    for (int i=0; i<n; i++) {
        Rcpp::Rcout << "*";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
    int n = floor((to-from)/double(by)) + 1;
    vector<int> ret(n,from);
    for (int i=1; i<n; i++) {
        from += by;
        ret[i] = from;
    }
    return(ret);
}

//------------------------------------------------
// searches for integer item in vector. Returns true if found, otherwise false
bool item_in_vector(vector<int> &vec, int item) {
    bool output = find(vec.begin(), vec.end(), item)!=vec.end();
    return(output);
}

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int Rcpp_to_int(SEXP x) {
    return(Rcpp::as<int>(x));
}

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double Rcpp_to_double(SEXP x) {
    return(Rcpp::as<double>(x));
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
vector<int> Rcpp_to_vector_int(SEXP x) {
    return(Rcpp::as<vector<int> >(x));
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
vector<double> Rcpp_to_vector_double(SEXP x) {
    return(Rcpp::as<vector<double> >(x));
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
vector<string> Rcpp_to_vector_string(SEXP x) {
    return(Rcpp::as<vector<string> >(x));
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
vector< vector<int> > Rcpp_to_mat_int(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<int> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<int> >(x[i]);
    }
    return(x_mat);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
vector< vector<double> > Rcpp_to_mat_double(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<double> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<double> >(x[i]);
    }
    return(x_mat);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector< vector< vector<int> > > Rcpp_to_array_int(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<int> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<int> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
        }
    }
    return(ret);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector< vector< vector<double> > > Rcpp_to_array_double(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<double> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<double> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
        }
    }
    return(ret);
}
