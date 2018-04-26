# include <RcppArmadillo.h>
# include "dtscaled_h.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' A C++ accelerated version of sugscompPml
//'
//' @param X The data matrix with rows as observations
//' @inheritParams sugsclustMargcpp
//' @inheritParams sugs
//' @param N The total number of people to be clustered
//' 
//' @return The log PML.
//' @export
// [[Rcpp::export]]
arma::vec sugscompPmlcpp(NumericMatrix X,
                         int K,
                         int N,
                         int D,
                         NumericVector n,
                         NumericMatrix phi,
                         NumericVector betaHat,
                         NumericMatrix m,
                         NumericVector nu,
                         NumericVector lambda,
                         NumericMatrix S,
                         NumericVector mu_0,
                         double nu_0,
                         double lambda_0,
                         NumericVector S_0){


  int L = betaHat.size() ;
  arma::mat predclust(N, K+1) ;
  arma::mat scaledPi(K, L) ;
  arma::vec CRPz(K) ;
  arma::vec pred(N) ;
  arma::vec PML(1) ;
  arma::vec CRPzz(1) ;
  NumericVector scalenew(D);

  for(int j = 0; j<K; j++)
  {
    for(int l=0; l<L; l++)
    {
      scaledPi(j,l) = phi(N-1, l) *  n(j)/(N + betaHat(l) - 1) ; //prior clust via chinese rest process non empty clusters
    }
  }

  CRPz = arma::sum(scaledPi, 1) ;
  CRPzz = sum(phi.row(N-1) * (betaHat)/(N + betaHat - 1)) ;  //prior clust via CRP new cluster

  if(K > 1){
    NumericMatrix scale(K, D) ;
    for (int i = 0; i<N; i++)
    {
      for (int j = 0; j<K; j++){
        scale.row(j) = sqrt((1 + lambda(j)) * S.row(j)/lambda(j));
        predclust(i,j) = exp(log(CRPz(j)) + arma::sum(dtscaled(X.row(i), rep(nu(j),D), m.row(j), scale.row(j), 1)));

      }
    }
  } else{
    NumericVector scale(D) ;
    scale = sqrt(((1 + lambda(0)) * S.row(0)/lambda(0))) ;
    for (int i = 0; i<N; i++)
    {
      predclust(i,0) = exp(log(CRPz(0)) + arma::sum(dtscaled(X.row(i), rep(nu(0),D), m, scale, 1))) ;

    }
  }
  //empty cluster
  scalenew = sqrt(((1 + lambda_0) * S_0/lambda_0)) ;

  for (int i = 0; i<N; i++)
  {
    predclust(i, K) = exp(log(CRPzz(0)) + arma::sum(dtscaled(X.row(i), rep(nu_0,D), mu_0, scalenew, 1))) ; //predictive of new cluster
  }

  pred = arma::sum(predclust, 1) ; //sum along rows, that is sum across clusters
  PML = arma::sum(log(pred), 0) ;

  return PML ;
}
