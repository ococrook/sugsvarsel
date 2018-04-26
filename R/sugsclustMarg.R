#' Computes the marginal probability of an observation belonging to a currently occupied cluster.
#'
#' @param x The observation currently underconsideration.
#' @param K The number of clusters.
#' @param i The current iteration state of the algorithm.
#' @param D The number of variable in the data matrix.
#' @param n The vector indicating the number of observations in each cluster
#' @param betaHat The grid of concentration parameters.
#' @param phi The probability weights of \code{beteHat}.
#' @param m The current posterior mean.
#' @param nu The current posterior degrees of freedom.
#' @param S The current posterior scale vector.
#' @param lambda The current posterior mean variance.
#' @param intfeature A binary vector indicating whether features are irrelevant (0) or relevant (1).
#'
#' @return The unnormalised probability of belong to any of the currently occupied clusters.

sugsclustMarg<-function(x, K, i, D, n, betaHat, phi, m, nu, S, lambda, intfeature){

  L <- length(betaHat)

  probz <- matrix(0, K)
  probx <- matrix(0, K, D)
  probmember <- matrix(0, K)
  probxFeat <- matrix(0, K, D)
  scaledPi <- matrix(0, K, L)
  scale <- matrix(0, K, D)

  if (K>1) {
    for (j in 1:K) {
      for (l in 1:L) {
        scaledPi[, l] <- phi[i-1, l] * (n[j])/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
      }
      probz[j] <- sum(scaledPi[j, ])  # marginal for z

      #compute scale matrix
      scale[j,] <- ((1 + lambda[j]) * S[j, ]/lambda[j])^(1/2)

      probx[j,] <- dt.scaled(x, mean = m[j, ], df = nu[j], sd = scale[j, ], log = TRUE)
      probxFeat[j,] <- intfeature * probx[j, ] #remove features that are switched off by multiplcation by 0
      probmember[j] <- exp(log(probz[j]) + sum(probxFeat[j, ])) #multiply, normalise later , exp(sum(log)) for stability
    }
  } else{
    for (l in 1:L) {
      scaledPi[, l] <- phi[i-1, l]*(n)/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
    }

    probz <- sum(scaledPi)

    #compute scale matrix
    scale <- ((1 + lambda) * S/lambda)^(1/2)

    probx <- dt.scaled(x, mean = m, df = nu, sd = scale, log = TRUE)
    probxFeat <- intfeature * probx
    probmember <- exp(log(probz) + sum(probxFeat)) #multiply, normalise later , exp(sum(log)) for stability
  }

  #unit test
  if (sum(is.na(probmember)) > 0) {
    print("marginals error")
    print(S)
    print(probx)
    stop(sum(is.na(probmember)) > 1)
  }

  return(probmember)
}

#' Computes the marginal probability of an observation belonging to a currently unoccupied cluster.
#'
#' @param x The observation currently underconsideration.
#' @param i The current iteration state of the algorithm.
#' @param D The number of variable in the data matrix.
#' @param betaHat The grid of concentration parameters.
#' @param phi The probability weights of \code{beteHat}.
#' @param mu_0 The prior mean.
#' @param nu_0 The prior degrees of freedom.
#' @param S_0 The prior scale vector.
#' @param lambda_0 The prior mean variance.
#' @param intfeatures A binary vector indicating whether features are irrelevant (0) or relevant (1).
#'
#' @return The unnormalised probability of belong to any of the currently unoccupied clusters.

sugsnewclustMarg <- function(x, i, D, phi, betaHat, mu_0, lambda_0, nu_0, S_0, intfeatures){

  L <- length(betaHat)
  scaledPi <- matrix(0,L)

  for(l in 1:L){
    scaledPi[l] <- phi[i-1, l] * (betaHat[l])/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
  }
  newProb <- sum(scaledPi) #marginal for z

  #compute scale matrix
  scale <- ((1 + lambda_0) * S_0/lambda_0)^(1/2)

  newPredict <- dt.scaled(x, mean = mu_0, df = nu_0, sd = scale, log = TRUE)
  newPredictFeat <- intfeatures * newPredict
  newClustprob <- exp(log(newProb) + sum(newPredictFeat)) # product

  return(newClustProb=newClustprob)
}

