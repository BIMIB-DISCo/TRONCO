#include <Rcpp.h>
using namespace Rcpp;


//' execute the bootstrap procedure
//' 
//' @param x A single integer.
// [[Rcpp::export]]
int C_bootstrap(NumericMatrix A, int x) {
    //Rcout << A.size();
    //Rcout << '\n';
    Rcout << A[0,0];
    Rcout << '\n';
    Rcout << x;
    Rcout << '\n';
    return x * 2;
}

// [[Rcpp::export]]
Rcpp::List C_get_dag_scores(Rcpp::NumericMatrix dataset,
                    Rcpp::NumericMatrix adj_matrix, 
                    double epos, 
                    double eneg ) {
    Rcout << "dataset\n";
    Rcout << dataset;
    Rcout << "\n adj.matrix \n";
    Rcout << adj_matrix;
    Rcout << "\n epos \n";
    Rcout << epos;
    Rcout << "\n eneg \n";
    Rcout << eneg;
    Rcout << "\n";
    Rcpp::NumericMatrix marginal_probs(3, 3);
    Rcpp::NumericMatrix joint_probs(3, 3);
    Rcpp::NumericMatrix prima_facie_model(3, 3);
    Rcpp::NumericMatrix prima_facie_null(3, 3);
    Rcpp::List scores = Rcpp::List::create(Rcpp::Named("marginal.probs") = marginal_probs,
        Rcpp::Named("joint.probs") = joint_probs,
        Rcpp::Named("prima.facie.model") = prima_facie_model,
        Rcpp::Named("prima.facie.null") = prima_facie_null);
    return scores;
}
