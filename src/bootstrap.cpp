#include <Rcpp.h>
using namespace Rcpp;


//' execute the bootstrap procedure
//' 
//' @param x A single integer.
// [[Rcpp::export]]
int C_bootstrap(int x) {
   return x * 2;
}
