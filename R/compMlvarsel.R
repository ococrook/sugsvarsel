#'Computes the log marginal likelihood of clustering and variable selection
#'
#' @inheritParams sugsComp
#' @inheritParams addStatsDiag
#' @inheritParams sugsclustMarg
#' @inheritParams sugsvarAlloc
#'
#' @return The log marginal likelhood

compMlvarsel<-function(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S, lognullMarg, intfeature){

  Marg  <- matrix(0, D)
  scale <- matrix(0, K, D)
  logMarg <- matrix(0, K, D)

  priorMarg <- log(lambda_0)/2 + nu_0 * log((nu_0 * S_0))/2 - lgamma(nu_0/2) #compute normalising constant of prior

  for(j in 1:K){
    #compute normalising constant of the jth cluster
    logMarg[j, ] <- lgamma(nu[j]/2) - (n[j]/2) * log(pi)/2 - log(lambda[j])/2 - nu[j] * log((nu[j] * S[j, ]))/2

  }
  ML <- sum(colSums(apply(logMarg, 1, function(x) intfeature * x))) + sum(colSums(apply(priorMarg, 1, function(x) intfeature * x))) + sum(apply(lognullMarg, 1, function(x) (1 - intfeature) * x))



  return(ML=ML)
}
