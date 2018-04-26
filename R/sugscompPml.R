#' Function to compute the PML for vanilla SUGS
#'
#' @inheritParams addStatsDiag
#' @inheritParams sugsclusterProb
#' @inheritParams sugsnewclusterProb
#' @param X The data matrix with rows as observations
#' @param N The total number of people to be clustered
#'
#' @return The log PML.

sugscompPml <- function(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0){

  L <- length(betaHat)
  predclust <- matrix(0, N, K+1)
  scale <- matrix(0, K, D)
  scaledPi <- matrix(0, K, L)
  CRPz<-matrix(0,K)

  for (j in 1:K) {
    for (l in 1:L) {
      scaledPi[,l] <- phi[N, l] * (n[j])/(N + betaHat[l] - 1)  #prior clust via chinese rest process non empty clusters
    }
    CRPz[j] <- sum(scaledPi[j, ])
  }
  CRPzz <- sum(phi[N, ] * (betaHat)/(N + betaHat - 1)) #prior clust via CRP new cluster

  if (K > 1) {
    for (i in 1:N) {
      for (j in 1:K) { #non-empty clusters
        scale[j, ] <- ((1 + lambda[j]) * S[j, ]/lambda[j])^(1/2)
        predclust[i, j] <- exp(log(CRPz[j]) + sum(dt.scaled(X[i, ], df=nu[j], mean = m[j,], scale[j, ], log=TRUE)))
      }
    }
  } else{
     scale <- ((1 + lambda[1]) * S[1,]/lambda[1])^(1/2)
     for(i in 1:N){
       predclust[i, 1] <- exp(log(CRPz[1]) + sum(dt.scaled(X[i, ], df = nu[1], mean = m, scale, log=TRUE)))
     }
  }

  #empty cluster
  scalenew <- ((1 + lambda_0) * S_0/lambda_0)^(1/2)
  for(i in 1:N){
    predclust[i, K+1] <- exp(log(CRPzz) + sum(dt.scaled(X[i, ], mean = mu_0, df = nu_0, sd = scalenew[1,], log = TRUE))) #predictive of new cluster
  }

  pred <- rowSums(predclust) #sum along rows, that is sum across clusters
  PML <- sum(log(pred))

  return(PML=PML)
}
