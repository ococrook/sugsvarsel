#' Implements the p-variate selection strategy
#'
#' @param mydata A data matrix with observations as rows
#' @param iter The number of random orderingsdev for the SUGS algorithm
#' @param p The number of variables to select
#' @param numSelect The number of times to randomly select variables
#' @param Model Character String indicating whether PML, ML or (both) should be computed.
#' Model selection is performed using PML if PML or both is selected, else model
#' selection is performed using ML.
#' @param mu_0 The mean hyperparameter, default is the column means of the data matrix.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{D}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the deault value is a tenth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100)}.
#' @param a The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.
#' @param b The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.
#' @param w The prior probability of a variable belong to the irrelevant or relevant partition. The vector must
#' contain two entries the first entry being the probabiliy of being irreleavnt and the second being the probability of being relevant
#' The default value is \code{c(0.5,0.5)}.
#'
#' @param BPPARAM Support for parallel processing using the
#' \code{BiocParallel} infrastructure. When missing (default),
#'     the default registered \code{BiocParallelParam} parameters are
#'     used. Alternatively, one can pass a valid
#'     \code{BiocParallelParam} parameter instance: \code{SnowParam},
#'     \code{MulticoreParam}, \code{DoparParam}, \ldots see the
#'     \code{BiocParallel} package for details. To revert to the
#'     origianl serial implementation, use \code{NULL}.
#'
#' @return A matrix where each row is a binary vector indicating whether that variables belongs
#' to the relevant (1) or irrelevant (0) partition.

pSelect<-function(mydata,
                  iter,
                  p,
                  numSelect,
                  Model,
                  mu_0 = NULL,
                  lambda_0 = 0.01,
                  nu_0 = NULL,
                  S_0 = NULL,
                  betaHat = c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100),
                  a = 10,
                  b = 1,
                  w = c(0.5, 0.5),
                  BPPARAM = bpparam()
                  ){

  #defensive tests
  stopifnot(is.numeric(mydata))
  stopifnot(is.numeric(iter))
  if (!(Model %in% c("PML", "ML", "Both"))) {
    stop("Please state a valid model selection criteria")
  }
  stopifnot(is.numeric(numSelect))
  stopifnot(numSelect > 0)
  stopifnot(iter > 0)


  X <- mydata
  D <- ncol(mydata)
  N <- nrow(mydata)

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


  #compute null view
  lognullMarg <- nullView(X, D, N, lambda_0, nu_0, mu_0, S_0)

  ind <- seq(1:D) #variable index

  intfeature <- array(0, c(numSelect, D))

  for(r in 1:numSelect){

    index <- sample(ind, p)
    tempdata <- mydata[, index]

    SUGSres <- runsugsvarsel(tempdata, 1, iter, rep(1,p), Model,
                             mu_0 = NULL, lambda_0 = lambda_0, nu_0 = NULL,
                             S_0 = NULL, betaHat = betaHat, a = a, b = b,
                             w = w, BPPARAM = BPPARAM )

    #compute feature marginals for best clustering
    X <- mydata[SUGSres$ordering[which.max(SUGSres$ML),],]

    #get best cluster
    numclust <- max(unique(SUGSres$member[which.max(SUGSres$ML),]))
    n <- as.matrix(table(SUGSres$member[which.max(SUGSres$ML),]))
    clustlabels <- SUGSres$member[which.max(SUGSres$ML),]

    #compute posterior for best clusering
    posterior <- compsugsStats(n, nu_0, lambda_0, mu_0, S_0, clustlabels, X, numclust)

    lambda <- posterior$lambda
    m  <- posterior$m
    nu <- posterior$nu
    S  <- posterior$S

    varMarg <- compvarMarg(numclust, D, n, lambda_0, nu_0, S_0, m, nu, lambda, S)
    intfeature[r,] <- sugsvarAlloc(D, varMarg, lognullMarg, rep(1,D), w)

  }

  return(intfeature)
}
