#' Computes the log marginal likelihood of clustering for vanilla SUGS
#'
#' @inheritParams addStatsDiag
#' @param D The number of variable in the data
#'
#' @return The log Marginal likelihood of the clustering

sugscompMl <- function(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S){

  Marg  <- matrix(0, D)
  scale <- matrix(0, K, D)
  logMarg <- matrix(0, K, D)

  priorMarg <- log(lambda_0)/2 + nu_0 * log((nu_0 * S_0))/2 - lgamma(nu_0/2) #compute normalising constant of prior

    for(j in 1:K){
      #compute normalising constant of the jth cluster
      logMarg[j,] <- lgamma(nu[j]/2) - log(pi^(n[j]/2))/2 - log(lambda[j])/2 - nu[j]*log((nu[j]*S[j,]))/2
    }
  MLSUGS <- sum(colSums(logMarg) + priorMarg)



  return(ML=MLSUGS)
}
