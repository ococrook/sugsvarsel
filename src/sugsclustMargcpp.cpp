#include <RcppArmadillo.h>
#include "dtscaled_h.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' A C++ accelerated version of sugsclusterMarg
//'
//' @param x A numeric vector containing data of the observation
//' @param K An integer specifiying the number of clusters
//' @param i An integer specifying the current iteration of SUGS algorithm
//' @param D An integer specifiying the number of variables
//' @param n A numeric vector contain the number of people already allocated to each cluster
//' @param betaHat A numeric vector containing the grid of hyperpriors for the dirichlet
//' concentration parameter
//' @param phi A numeric matrix containing the weights for the dirichlet hyperpiors
//' @param m A numeric matrix containing the means for each component
//' @param nu A numeric vector containg the degrees of freedom for each component
//' @param S A numeric matrix containing the scale prior for each component
//' @param lambda A numeric vector containing the mean variance hyperparamter for each component
//'
//' @return An arma::vec containg the probability that observation \code{x} belongs to each component
//' @export
// [[Rcpp::export]]
arma::vec sugsclustMargcpp(NumericVector x,
                                         int K,
                                         int i,
                                         int D,
                                         NumericVector n,
                                         NumericVector betaHat,
                                         NumericMatrix phi,
                                         NumericMatrix m,
                                         NumericVector nu,
                                         NumericMatrix S,
                                         NumericVector lambda,
                                         arma::rowvec intfeature)
{

  int L = betaHat.size();

  arma::vec probz(K);
  arma::vec probmember(K);
  arma::mat probxFeat(K, D);
  arma::mat probx(K, D);
  arma::rowvec dtscaledRes;



  if(K>1){
  NumericMatrix scale(K,D);
  arma::mat scaledPi(K,L);

    for(int j = 0; j < K; j++)
    {
      for(int l = 0; l < L; l++)
      {
        scaledPi(j, l) = (phi(i-2, l) * (n(j))/(i + betaHat(l) - 1)); //discrete gamma prior for dirichlet concentration
      }
    }

  //Note that the 1 in the below ensures we get row sums
  probz  = arma::sum(scaledPi, 1);   //marginal for z
    for(int j = 0; j < K; j++)
    {
      //scale[j,]<-((1+lambda[j])*S[j,]/lambda[j])^(1/2)
      scale.row(j)  = sqrt((1 + lambda(j)) * S.row(j)/lambda(j));
      probx.row(j)  = dtscaled(x, rep(nu(j),D), m.row(j), scale.row(j), 1);
      probxFeat.row(j) = intfeature % probx.row(j);
      probmember(j) = exp(log(probz(j)) + arma::sum(probxFeat.row(j))); //multiply, normalise later , exp(sum(log)) for stability
    }
  }
  else{
    NumericVector scale(D);
    arma::vec scaledPi(L);

      for(int l = 0; l < L; l++)
      {
        scaledPi(l) = (phi(i-2, l) * (n(0))/(i + betaHat(l) - 1)); //discrete gamma prior for dirichlet concentration
      }
    probz = arma::sum(scaledPi);

  //compute scale matrix
    scale      = sqrt((1 + lambda(0)) * S/lambda(0)) ;
    probx      = dtscaled(x, rep(nu,D), m, scale, 1) ;
    probxFeat  = intfeature % probx ;
    probmember = exp(log(probz) + arma::sum(probx, 1)) ; //multiply, normalise later , exp(sum(log)) for stability

  }
  return(probmember);
}

//' A C++ accelerated version of sugsnewclustMarg
//'
//' @inheritParams sugsclustMargcpp
//' @inheritParams sugsclustMarg
//' @inheritParams sugs
//' @return The unnormalised probability of belong to any of the currently unoccupied clusters.
// [[Rcpp::export]]
arma::vec sugsnewclustMargcpp(NumericVector x,
                              int i,
                              int D,
                              NumericMatrix phi,
                              NumericVector betaHat,
                              NumericVector mu_0,
                              double nu_0,
                              NumericVector S_0,
                              double lambda_0,
                              arma::rowvec intfeature)
{
  int L = betaHat.size();

  NumericVector scale(D);
  arma::vec newProb(1);
  arma::vec newClustprob(1);
  arma::rowvec newPredictFeat(D);
  arma::rowvec newPredict(D);
  arma::vec scaledPi(L);

  for(int l = 0; l < L; l++)
  {
    scaledPi(l) = (phi(i-2, l) * betaHat(l)/(i + betaHat(l) - 1)); //discrete gamma prior for dirichlet concentration
  }
  newProb = arma::sum(scaledPi); //marginal for z

 //compute scale matrix
    scale = sqrt(((1 + lambda_0) * S_0/lambda_0)) ;
    newPredict = dtscaled(x, rep(nu_0,D), mu_0, scale, 1) ;
    newPredictFeat = intfeature % newPredict ;
    newClustprob = exp(log(newProb) + sum(newPredictFeat)) ; // product

  return newClustprob ;
}
