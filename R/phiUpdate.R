#' Function to update posterior for dirichlet concentration parameters.
#'
#' @param phi Current value of posterior dirichlet weights.
#' @param i Current iteration of the SUGS algorithm.
#' @param clustnew The cluster to which the \code{i}th observation has allocated.
#' @param n A vector indicating the number of observations allocated to each cluster.
#' @param betaHat The grid prior for the dirichlet concentration parameter.
#'
#' @return An update value for phi.

phiUpdate<-function(phi, i, clustnew, n, betaHat){

  L <- length(betaHat)

  for (jj in 1:L) {
    phi[i,jj] <- phi[i-1, jj] * (n[clustnew]/(i + betaHat[jj] - 1))
  }
  phi[i, ] <- phi[i, ]/sum(phi[i, ])

  return(phi=phi)
}
