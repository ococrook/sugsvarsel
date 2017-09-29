#' Computes the log PML for SUGS VarSel in the case of Gaussian mixtures
#'
#' @param X The data matrix with observations as rows
#' @param N The number of observations
#' @param D The number of variables
#' @inheritParams sugsvarAlloc
#' @inheritParams sugsComp
#' @inheritParams addStatsDiag
#' @inheritParams phiUpdate
#'
#' @return The log PML

compPmlvarsel <- function(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0, intfeature, w){

  L <- length(betaHat)
  predclust <- array(0, c(N, D, K+1))
  scale <- matrix(0, K, D)
  scaledPi <- matrix(0, K, L)
  CRPz <- matrix(0, K)

  #compute values for NULL view (no clustering structure)
  #functions as if no clustering structure, vectorised over D
  #statistics
  x_barNULL <- colSums(X)/N
  SCNULL <- colSums(X^2)

  #posterior updates
  nuNULL <- nu_0 + N
  lambdaNULL <- lambda_0 + N
  mNULL <- (lambda_0 * mu_0 + N * x_barNULL)/lambdaNULL
  SNULL <- SCNULL/nuNULL + nu_0*S_0/nuNULL - (lambdaNULL * mNULL^2/nuNULL) + (lambda_0 * mu_0^2/nuNULL)
  scaleNULL <- ((1 + lambdaNULL) * SNULL/lambdaNULL)^(1/2)

  for (j in 1:K) {
    for (l in 1:L) {
      scaledPi[, l] <- phi[N, l] * (n[j])/(N + betaHat[l] - 1)  #prior clust via chinese rest process non empty clusters
    }
    CRPz[j] <- sum(scaledPi[j, ])
  }
  CRPzz <- sum(phi[N, ] * (betaHat)/(N + betaHat - 1)) #prior clust via CRP new cluster

  if (K>1) {
      for (j in 1:K){ #non-empty clusters
        scale[j, ] <- ((1 + lambda[j]) * S[j, ]/lambda[j])^(1/2)
        predclust[, , j] <- CRPz[j] * dt.scaled(X, df = nu[j], mean = m[j, ], scale[j, ], log = FALSE)

      }
  } else{
    scale <- ((1 + lambda[1]) * S[1, ]/lambda[1])^(1/2)
    predclust[, , 1] <- CRPz[1] * dt.scaled(X, df = nu[1], mean = m, scale, log = FALSE)
  }

  #empty cluster
  scalenew <- ((1 + lambda_0) * S_0/lambda_0)^(1/2)
  predclust[, , K+1] <- CRPzz * dt.scaled(X, df = nu_0, mean = mu_0, scalenew[1, ], log=FALSE)

  pred <- rowSums(apply(log(apply(predclust, c(1, 2), sum)),1,function(x) intfeature * x))
    + log(w[1]) + log(w[2])
    + rowSums(apply(dt.scaled(X, df = nuNULL, mean = mNULL, scaleNULL[1, ], log=TRUE), 1, function(x) (1 - intfeature) * x))   #sum along rows, that is sum across clusters

  PML <- sum((pred))

  return(PML=PML)
}

#' Computes the log PML for SUGS VarSel in one dimension in the case of Gaussian mixtures
#'
#' @param X The data matrix with observations as rows
#' @param N The number of observations
#' @param D The number of variables
#' @inheritParams sugsvarAlloc
#' @inheritParams sugsComp
#' @inheritParams addStatsDiag
#' @inheritParams phiUpdate
#'
#' @return The log PML

compPmlvarsel1D <- function(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0){

  L <- length(betaHat)
  predclust <- matrix(0, N, K+1)
  scale <- matrix(0, K, D)
  scaledPi <- matrix(0, K, L)
  CRPz <- matrix(0, K)

  for (j in 1:K) {
    for (l in 1:L) {
      scaledPi[, l] <- phi[N, l] * (n[j])/(N + betaHat[l] - 1)  #prior clust via chinese rest process non empty clusters
    }
    CRPz[j] <- sum(scaledPi[j, ])
  }
  CRPzz <- sum(phi[N, ] * (betaHat)/(N + betaHat - 1)) #prior clust via CRP new cluster

  if (K > 1) {
    for (i in 1:N){
      for (j in 1:K){ #non-empty clusters
        scale[j] <- ((1 + lambda[j]) * S[j]/lambda[j])^(1/2)
        predclust[i, j] <- exp(log(CRPz[j]) + sum(dt.scaled(X[i], df = nu[j], mean = m[j], scale[j], log = TRUE)))
      }
    }
  } else {
    scale <- ((1 + lambda[1]) * S[1]/lambda[1])^(1/2)
     for (i in 1:N){
       predclust[i, 1]<-exp(log(CRPz[1]) + sum(dt.scaled(X[i], df = nu[1], mean = m, scale, log = TRUE)))
     }
  }

  #empty cluster
  scalenew <- ((1 + lambda_0) * S_0/lambda_0)^(1/2)
  for(i in 1:N){
    predclust[i, K+1]<-exp(log(CRPzz) + sum(dt.scaled(X[i], mean = mu_0, df = nu_0, sd = scalenew[1], log = TRUE))) #predictive of new cluster
  }
  pred <- rowSums(predclust) #sum along rows, that is sum across clusters
  PML  <- sum(log(pred))

  return(PML=PML)
}
