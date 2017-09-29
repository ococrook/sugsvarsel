#' Perform clustering and variable selection using the SUGS algorithm.
#'
#' @param X A data matrix with rows as observations.
#' @param featiter The number of iterations of variable selection
#' @param clustiter The number of random ordering, for which to apply SUGS.
#' @param intfeatures A binary matrix of the intial variable set, probably chosen using function \code{pSelect}.
#' See documentation for \code{pSelect} for more details.
#' @param numSelect The total number of feature sets for the algorithm to intialise.
#' @param Model The method used for Model select, either PML, ML or Both. If you select both
#' the PML will be used to perform model selection.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{D}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the deault value is a tenth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(0.01, 0.1, 1, 5, 10, 15, 30, 50, 100)}.
#' @param a The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.
#' @param b The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.
#' @param w The prior probability of a variable belong to the irrelevant or relevant partition. The vector must
#' contain two entries the first entry being the probabiliy of being irreleavnt and the second being the probability of being relevant
#' The default value is \code{c(0.5,0.5)}.
#' @param BPPARAM Support for parallel processing using the
#' \code{BiocParallel} infrastructure. When missing (default),
#'     the default registered \code{BiocParallelParam} parameters are
#'     used. Alternatively, one can pass a valid
#'     \code{BiocParallelParam} parameter instance: \code{SnowParam},
#'     \code{MulticoreParam}, \code{DoparParam}, \ldots see the
#'     \code{BiocParallel} package for details. To revert to the
#'     origianl serial implementation, use \code{NULL}.
#'
#'
#' @return A vector of log marginal likelihoods, a matrix of memberships, the reorderings of the data and the associated feature sets.

sugsvarsel <- function(X,
                       featiter,
                       clustiter,
                       intfeatures,
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
                       BPPARAM = bpparam(),
                       Verbose = T
                       ){

  #defensive tests
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(featiter))
  stopifnot(is.numeric(clustiter))
  if (!(Model %in% c("PML", "ML", "Both"))) {
    stop("Please state a valid model selection criteria")
  }
  stopifnot(is.numeric(numSelect))
  stopifnot(numSelect > 0)
  stopifnot(intfeatures %in% c(0,1))
  stopifnot(!(is.null(ncol(X))))
  stopifnot(ncol(X)==ncol(intfeatures))

  SUGSfeatRes <- vector("list", 0)

  if (numSelect==1) {
    SUGSfeatRes <- runsugsvarsel(X, featiter, clustiter, intfeatures, Model,
                                 mu_0 = mu_0, lambda_0 = lambda_0, nu_0 = nu_0,
                                 S_0 = S_0, betaHat = betaHat, a = a, b = b,
                                 w = w, BPPARAM=BPPARAM)
    Result <- SUGSfeatRes
  } else {
     for (d in 1:numSelect) {
       if(isTRUE(Verbose)){
         cat("\n", "Starting initial iteration", d)
       }

      SUGSfeatRes[[d]] <- runsugsvarsel(X, featiter, clustiter, intfeatures[d,], Model,
                                        mu_0 = mu_0, lambda_0 = lambda_0, nu_0 = nu_0,
                                        S_0 = S_0, betaHat = betaHat, a = a, b = b,
                                        w = w, BPPARAM=BPPARAM)
     }
    Result <- extract(X, SUGSfeatRes)
  }

  return((Result=Result))
}

