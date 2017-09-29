# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Function to update posterior for dirichlet concentration parameters.
//'
//' @param phi Current value of posterior dirichlet weights.
//' @param i Current iteration of the SUGS algorithm.
//' @param clustnew The cluster to which the \code{i}th observation has allocated.
//' @param n A vector indicating the number of observations allocated to each cluster.
//' @param betaHat The grid prior for the dirichlet concentration parameter.
//'
//' @return An update value for phi.


// [[Rcpp::export]]
NumericMatrix phiUpdatecpp( NumericMatrix phi, int i, int clustnew, NumericVector n, NumericVector betaHat){

  int L = betaHat.size();

  for(int jj=0; jj<L; jj++){
    phi(i-1, jj) = phi(i-2, jj) * n(clustnew - 1)/(i + betaHat(jj) - 1) ;
  }
  phi.row(i-1) = phi.row(i-1)/sum(phi.row(i-1)) ;

  return phi ;
}

