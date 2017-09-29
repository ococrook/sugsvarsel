#' The computational function for the SUGS with variable selection algorithm.
#'
#' @param mydata Data matrix with observations as rows.
#' @param intfeature A binary vector of feature which are parition as irrelevant (0)
#' or relevant (1).
#' @param Model A character string sating whether PML, ML or both are used for feature selection
#' and returned.
#' @param mu_0 The mean hyperparameter, default is the column means of the data matrix.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{D}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the deault value is a fifth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100)}.
#' @param a The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.
#' @param b The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.
#' @param w The prior probability of a variable belong to the irrelevant or relevant partition. The vector must
#' contain two entries the first entry being the probabiliy of being irreleavnt and the second being the probability of being relevant
#' The default value is \code{c(0.5,0.5)}.
#'
#' @return An allocation vector, called member, of observation to clusters, K the number of clusters,
#' n the number of observation in each cluster, the log PML or log ML or both as appropriate, a binary
#' vector of features that have been classified as irrelevant (0) or relevant (1).

sugsComp<-function(mydata,
                   intfeature,
                   Model,
                   mu_0 = NULL,
                   lambda_0 = 0.01,
                   nu_0 = NULL,
                   S_0 = NULL,
                   betaHat = c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100),
                   a = 10,
                   b = 1,
                   w = c(0.5, 0.5)
                  ){

  #get data matrix
  X <- mydata

  #data dimensions
  N <- nrow(X)
  D <- ncol(X)
  K <- 1 #initially 1 cluster

  #Hyperparameters
  if ( is.null(mu_0)) {
    mu_0 <- colMeans(X)
  } else{
    stopifnot( length(mu_0)==D )
    mu_0 <- mu_0
  }
  if ( is.null(nu_0) ) {
    nu_0 <- D
  } else{
    stopifnot(is.numeric(nu_0))
    nu_0 <- nu_0
  }
  if (is.null(S_0)) {
    S_0 <- ((colMeans(X * X) - colMeans(X)^2) * (N / (N - 1)))/5
  } else {
    stopifnot( length(S_0)==D )
    S_0 <- S_0
  }
  S_0  <- rbind(c(S_0))

  #dirichlet hyperparameter
  betaHat <- betaHat #grid of betas
  etaz <- dgamma(betaHat, shape = a, b) #unnormalised discretisation of gamma dist
  eta <- etaz/sum(etaz)

  lognullMarg <- nullView(X, D, N, lambda_0, nu_0, mu_0, S_0) #get null view for features

  #Storage
  n <- matrix(0, K)
  x_bar <- matrix(0, K, D)
  m <- matrix(0, K, D)
  nu <- array(0, K)
  lambda <- matrix(0, K)
  S <- array(0, c(K, D))
  SCL <- matrix(0, K, D)
  labelprob <- matrix(0, K, N)
  cluster <- matrix(0, N)
  phi <- matrix(0, N, length(betaHat))

  #assigned first person to cluster 1
  cluster[1] <- 1

  #compute statistics for first cluster
  n[1] <- 1
  x_bar <- X[1, ]
  SCL <- X[1, ]^2

  #compute first phi
  phi[1, ] <- (eta/betaHat)/sum(eta/betaHat)

  #posteriorparameter updates after first cluster
  nu[1] <- nu_0 + n[1]
  lambda <- lambda_0 + n[1]
  m <- (lambda_0 * mu_0 + n[1] * x_bar)/lambda[1]
  S[1,] <- SCL/nu[1] + nu_0 * S_0/nu[1] - (lambda[1] * m^2/nu[1]) +(lambda_0 * mu_0^2/nu[1])

  for(i in 2:N){ #loop over people
    #information for this iteration
    x <- X[i, ]

    #compute predicitve amongst current clusters
    #clustprob <- sugsclustMarg(x, K, i, D, n, betaHat, phi, m, nu, S, lambda, intfeature) #R function
    clustprob <- sugsclustMargcpp(x, K, i, D, n, betaHat, phi, as.matrix(m), nu, as.matrix(S), lambda, intfeature)

    #consider new label
    #newclustprob <- sugsnewclustMarg(x, i, D, phi, betaHat, mu_0, lambda_0, nu_0, S_0, intfeature) #R function
    newclustprob <- sugsnewclustMargcpp(x, i, D, phi, betaHat, mu_0, nu_0, S_0, lambda_0, intfeature)

    #greedly get new label
    prob <- c(clustprob, newclustprob) #additional cluster allowed no need to normalise
    clustnew <- which.max(prob) #sample new label
    cluster[i] <- clustnew #set new cluster

    #increase storage size
    if (clustnew > K){
      K <- K + 1
      n <- c(n, 0)
      x_bar <- rbind(x_bar, rep(0, D))
      m <- rbind(m, rep(0,D))
      nu <- c(nu, 0)
      lambda <- c(lambda, 0)
      SCL <- rbind(SCL, rep(0, D))

      if (ncol(S) != D){ #test
        print(S)
        print(dim(S))
        print("Error scale matrix incorrect size")
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

  #compute feature marginals
  varMarg <- compvarMarg(K, D, n, lambda_0, nu_0, S_0, m, nu, lambda, S)
  intfeature <- sugsvarAlloc(D, varMarg, lognullMarg, intfeature, w)

  if(Model == "PML"){
    LPML <- compPmlvarsel(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0, intfeature, w) #calculate estimate of logPML
  } else if(Model=="Both"){
    LPML <- compPmlvarsel(X, K, N, D, n, phi, betaHat, m, nu, lambda, S, mu_0, nu_0, lambda_0, S_0, intfeature, w) #calculate estimate of logPML
    ML <- compMlvarsel(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S, lognullMarg, intfeature)
  } else{
    ML <- compMlvarsel(K, D, n, nu_0, S_0, lambda_0, nu, lambda, S, lognullMarg, intfeature)
  }

  member <- cluster #assign members

  if (Model=="PML"){
    return(list(member=member, K=K, n=n, LPML=LPML, features=intfeature))
  } else if(Model=="Both"){
    return(list(member=member, K=K, n=n, LPML=LPML, ML=ML, features=intfeature))
  } else{
    return(list(member=member, K=K, n=n, ML=ML, features=intfeature))
  }

}
