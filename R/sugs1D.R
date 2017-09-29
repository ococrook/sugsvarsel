#' The SUGS clustering algorithm by Wang and Dunson (2011) for 1 dimensional data
#'
#' @param mydata Data matrix with observations as rows
#' @param Model Character string indicating whether to use PML, ML or
#' both for model selection. "Both" defaults to PML.
#'
#' @return A matrix of cluster allocation, K the number of clusters, a matrix
#' indicating the number of observations allocated to each cluster. The value of
#' the model selection criteria either log PML, log ML or both.

sugs1D <- function(mydata, Model){

  #get data matrix
  X <- mydata

  #data dimensions
  N <- length(X)
  D <- 1
  K <- 1 #initially 1 cluster

  #INITIALISE HYPERPARAMETERs
  lambda_0 <- 0.01
  S_0 <- matrix(0, 1)
  S_0 <- var(X)/3
  S_0 <- rbind(c(S_0))
  nu_0 <- 1
  mu_0 <- mean(X)

  priors <- list(mu_0=mu_0, lambda_0=lambda_0, nu_0 =nu_0, S_0=S_0)

  betaHat<-c(0.1, 1, 5, 10) #grid of betas
  etaz<-dgamma(betaHat, 4, 2) #unnormalised discretisation of gamma dist
  eta<-etaz/sum(etaz) #normalise

  #Storage
  n <- matrix(0, K)
  x_bar <- matrix(0, K)
  m  <- matrix(0, K)
  nu <- matrix(0, K)
  lambda < -matrix(0, K)
  S   <- matrix(0, c(K))
  SCL <- matrix(0,K)
  cluster <- matrix(0, N)
  phi <- matrix(0, N, length(betaHat))

  #assigned first person to cluster 1
  cluster[1] <- 1

  #compute statistics for first cluster
  n[1]  <- 1
  x_bar <- X[1]
  SCL <- X[1]^2

  #compute first phi
  phi[1, ] <- (eta/betaHat)/sum(eta/betaHat)

  #posterior parameter updates after first cluster
  nu[1] <- nu_0 + n[1]
  lambda <-lambda_0 + n[1]
  m <- (lambda_0 * mu_0 + n[1] * x_bar)/lambda[1]
  S[1] <- SCL/nu[1] + nu_0 * S_0/nu[1] - (lambda[1] * m^2/nu[1]) + (lambda_0 * mu_0^2/nu[1])

  for (i in 2:N){ #loop over people

    #information for this iteration
    x <- X[i]

    clusterprob <- sugsclusterProb(x, K, i, D, n, betaHat, phi, m, nu, S, lambda)
    newclusterprob <- sugsnewclusterProb(x, i, D, phi, betaHat, mu_0, lambda_0, nu_0, S_0 )

    #greedly get new label
    prob <- c(clusterprob, newclusterprob) #additional cluster allowed no need to normalise
    clustnew <- which.max(prob) #sample new label
    cluster[i] <- clustnew #set new cluster

    #increase storage size
    if(clustnew > K){
      K <- K + 1
      n <- c(n, 0)
      x_bar <- rbind(x_bar, rep(0, D))
      m <- rbind(m, rep(0, D))
      nu <- c(nu,0)
      lambda <- c(lambda, 0)
      SCL <- rbind(SCL, rep(0, D))

      if(ncol(S)!=D){ #test
        print(S)
        print(dim(S))
      }
      S <- rbind(S, rep(0, D))
    }

    #update cluster with new data
    addStats <- addStatsDiag(x, K, clustnew, n, x_bar, SCL, m, lambda, S, nu, beta, S_0, lambda_0, mu_0, nu_0)

    #update parameters
    x_bar <- addStats$x_bar
    n <- addStats$n
    SCL <- addStats$SCL
    m <- addStats$m
    lambda <- addStats$lambda
    nu <- addStats$nu
    S <- addStats$S

    #function to update phi
    phi <- phiUpdate(phi, i, clustnew, n, betaHat)

  }  #i loop ends here

  posterior <- list(m, lambda, nu, S)

  if (Model=="PML") {
    LPML <- compPmlvarsel1D(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0) #calculate estimate of logPML
  } else if(Model=="Both"){
    LPML <- compPmlvarsel1D(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0) #calculate estimate of logPML
    ML <- sugscompMl(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S)
  } else{
    ML <- sugscompMl(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S)
  }
  member<-cluster

  if(Model=="PML"){
    return(list(member=member, K=K, n=n, LPML=LPML, priors=priors, posterior=posterior))
  }else if(Model=="Both"){
    return(list(member=member, K=K, n=n, LPML=LPML, ML=ML, priors=priors, posterior=posterior))
  } else{
    return(list(member=member, K=K, n=n, ML=ML, priors=priors, posterior=posterior))
  }

}
