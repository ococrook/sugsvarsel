#include <RcppArmadillo.h>
using namespace Rcpp;

arma::rowvec dtscaled(NumericVector x, NumericVector df, NumericVector mean, NumericVector sd, int lg);
