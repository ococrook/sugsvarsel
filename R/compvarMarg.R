#' A function to compute the marginal probability that of a variable of being relevant.
#'
#' @inheritParams addStatsDiag
#' @inheritParams sugsComp
#' @param D The number of variables in the data matrix
#'
#' @return The log marginal probability of a variable being relevant.

compvarMarg <- function(K, D, n, lambda_0, nu_0, S_0, m, nu, lambda, S){

  lgammaMarg <- matrix(0, D)
  scale <- matrix(0, K, D)
  logclusterMarg <- matrix(0, K, D)
  priorMarg <- matrix(0, D)

  priorMarg <- log(lambda_0)/2 + nu_0 * log((nu_0 * S_0))/2 - lgamma(nu_0/2) #compute normalising constant of prior

  for (j in 1:K) {

      #compute normalising constant of the kth cluster
      logclusterMarg[j, ] <- lgamma(nu[j]/2) - (n[j]/2) * log(pi)/2 - log(lambda[j])/2 - (nu[j] * log((nu[j] * S[j,]))/2)

  }
  lgammaMarg <- colSums(logclusterMarg) + priorMarg

  return(lgammaMarg=lgammaMarg)

}

