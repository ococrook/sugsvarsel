#' The SUGS clustering algorithm by Wang and Dunson (2011)
#'
#'
#' @param mydata Data matrix with observations as rows.
#' @param Model Character string indicating whether to use PML, ML or
#' both for model selection. "Both" defaults to PML.
#' @param mu_0 The mean hyperparameter, default is the column means of the data matrix.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{2 * (D + 2)}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the default value is a tenth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(1, 5, 15, 30, 50, 100)}.
#' @param a The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.
#' @param b The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.
#'
#'
#' @return A matrix of cluster allocation, K the number of clusters, a matrix
#' indicating the number of observations allocated to each cluster. The value of
#' the model selection criteria either log PML, log ML or both.

sugs<-function(mydata,
               Model,
               mu_0 = NULL,
               lambda_0 = 0.01,
               nu_0 = NULL,
               S_0 = NULL,
               betaHat = c(1, 5, 15, 30, 50, 100),
               a = 10,
               b = 1
               ) {

  #get data matrix
  X <- mydata

  #data dimensions
  N <- nrow(X)
  D <- ncol(X)
  K <- 1 #initially 1 cluster


  #Hyperparameters
  if( is.null(mu_0)) {
    mu_0 <- colMeans(X)
  } else{
    stopifnot( length(mu_0)==D )
    mu_0 <- mu_0
  }
  if ( is.null(nu_0) ) {
    nu_0 <- 2 * (D + 2)
  } else{
    stopifnot(is.numeric(nu_0))
    nu_0 <- nu_0
  }
  if (is.null(S_0)) {
    S_0 <- ((colMeans(X * X) - colMeans(X)^2) * (N / (N - 1)))/10
  } else {
    stopifnot( length(S_0)==D )
    S_0 <- S_0
  }
  S_0  <- rbind(c(S_0))

  #dirichlet hyperparameter
  betaHat <- betaHat #grid of betas
  etaz <- dgamma(betaHat, shape = a, b) #unnormalised discretisation of gamma dist
  eta <- etaz/sum(etaz)

  #Sampler
  n  <- matrix(0, K)
  x_bar <- matrix(0, K, D)
  m  <- matrix(0, K, D)
  nu <- array(0, K)
  lambda <- matrix(0, K)
  S <- array(0, c(K, D))
  SCL <- matrix(0, K, D)
  cluster <- matrix(0, N)
  phi <- matrix(0,N ,length(betaHat))

  #assigned first person to cluster 1
  cluster[1] <- 1

  #compute statistics for first cluster
  n[1]  <- 1
  x_bar <- X[1, ]
  SCL <- X[1, ]^2

  #compute first phi
  phi[1, ] <- (eta/betaHat)/sum(eta/betaHat)

  #posteriorparameter updates after first cluster
    nu[1]  <- nu_0 + n[1]
    lambda <- lambda_0 + n[1]
    m <- (lambda_0*mu_0 + n[1] * x_bar)/lambda[1]
    S[1,] <- SCL/nu[1] + nu_0 * S_0/nu[1] - (lambda[1] * m^2/nu[1]) + (lambda_0 * mu_0^2/nu[1])

    for (i in 2:N) { #loop over people

      #information for this iteration
      x <- X[i,]

      #compute predicitve amongst current clusters
      #probcluster <- sugsclusterProb(x, K, i, D, n, betaHat, phi, m, nu, S, lambda) R version
      probcluster <- sugsclusterProbcpp(x, K, i, D, n, betaHat, phi, as.matrix(m), nu, as.matrix(S), lambda)

      #consider new cluster
      probclusternew <- sugsnewclusterProb(x, i, D, phi, betaHat, mu_0, lambda_0, nu_0, S_0)

      #greedly get new label
      prob <- c(probcluster, probclusternew) #additional cluster allowed no need to normalise
      clustnew <- which.max(prob) #sample (greedy) new label
      cluster[i]<-clustnew #set new cluster

      #increase storage size
      if (clustnew > K) {
        K <- K+1
        n <- c(n, 0)
        x_bar <- rbind(x_bar, rep(0, D))
        m <- rbind(m, rep(0, D))
        nu <- c(nu, 0)
        lambda <- c(lambda, 0)
        SCL <- rbind(SCL, rep(0, D))

        if (ncol(S) != D){ #test
          print(S)
          print(dim(S))
          print("Error with Scale matrix")
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

      #update weights for concentration parameter
      #phi <- phiUpdate(phi, i, clustnew, n, betaHat) R version
      phi <- phiUpdatecpp(phi, i, clustnew, n, betaHat)

    }  #i loop ends here

    if (Model=="PML"){
      #LPML <- sugscompPml(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0)  #R version
      LPML <- sugscompPmlcpp(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0) #calculate estimate of logPML
    } else if(Model=="Both"){
      LPML <- sugscompPmlcpp(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0) #calculate estimate of logPML
      ML <- sugscompMl(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S) #calculate ML
    } else{
      ML <- sugscompMl(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S) #calculate ML
    }
    member<-cluster #set member as cluster

    if(Model=="PML"){
      return(list(member=member, K=K, n=n, LPML=LPML))
    }else if(Model=="Both"){
      return(list(member=member, K=K, n=n, LPML=LPML, ML=ML))
    } else{
      return(list(member=member, K=K, n=n, ML=ML))
    }
}
