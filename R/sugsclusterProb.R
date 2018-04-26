#' A function to compute the marginal probability of belong to an occupied cluster for the SUGS algorithm
#'
#' @param x The current observation under consideration.
#' @param K The number of occupied clusters
#' @param i The current iteration of the SUGS algorithm
#' @param D The number of variables
#' @param n A numeric vector containing the number of observations allocated to each cluster
#' @param betaHat The grid of values for the dirichlet concentration parameter
#' @param phi The weights associated with the dirichlet grid prior
#' @param m The current posterior mean
#' @param nu The current posterior degrees of freedom
#' @param S The current posterior scale vector
#' @param lambda The current posterior for the mean variance
#'
#' @return The numeric vector with the unnormalised probability of beloging to each cluster.


sugsclusterProb<-function(x, K, i, D, n, betaHat, phi, m, nu, S, lambda){

  L <- length(betaHat)

  probz <- matrix(0, K)
  probx <- matrix(0, K, D)
  probmember <- matrix(0, K)
  scaledPi <- matrix(0, K, L)

  scale <- matrix(0, K, D)

  if (K > 1) {
    for (j in 1:K) {
      for (l in 1:L) {
        scaledPi[, l] <- phi[i-1, l] * (n[j])/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
      }
      probz[j] <- sum(scaledPi[j, ])  # marginal for z
      #compute scale matrix
      scale[j, ] <- ((1 + lambda[j])*S[j, ]/lambda[j])^(1/2)
      probx[j,]  <- dt.scaled(x, mean = m[j,], df = nu[j], sd = scale[j,], log = TRUE)
      probmember[j] <- exp(log(probz[j]) + sum(probx[j,])) #multiply, normalise later , exp(sum(log)) for stability
    }

  } else{
    for (l in 1:L) {
      scaledPi[, l] <- phi[i-1, l] * (n)/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
    }

    probz <- sum(scaledPi)
    #compute scale matrix
    scale <- ((1 + lambda) * S/lambda)^(1/2)
    probx <- dt.scaled(x, mean = m, df = nu, sd = scale, log = TRUE)
    probmember <- exp(log(probz) + sum(probx)) #multiply, normalise later , exp(sum(log)) for stability
  }

  #Test
  if(sum(is.na(probmember)) > 0){
    print("marginals error")
    print(S)
    print(probx)
    stop(sum(is.na(probmember)) > 1)
  }

  return(probmember)
}

#' A function to compute the probability of belonging to an unoccupied cluster for the SUGS algorithm
#'
#' @inheritParams sugsclusterProb
#' @inheritParams sugsComp
#'
#' @return The unnormlised probability of the observation belonging to a new cluster

sugsnewclusterProb <- function(x, i, D, phi, betaHat, mu_0, lambda_0, nu_0, S_0){

  L <- length(betaHat)
  scaledPi <- matrix(0, L)

  for(l in 1:L){
    scaledPi[l]<-phi[i-1, l] * (betaHat[l])/(i + betaHat[l] - 1) #discrete gamma prior for dirichlet concentration
  }
  newProb <- sum(scaledPi) #marginal for z
  #compute scale matrix
  scale <- ((1 + lambda_0) * S_0/lambda_0)^(1/2)
  newPredict <- dt.scaled(x, mean = mu_0, df = nu_0, sd = scale, log = TRUE)
  newClustprob <- exp(log(newProb) + sum(newPredict)) # product

  return(newClustProb=newClustprob)
}





