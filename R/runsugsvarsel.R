#' Internal function for SUGS with variable selection algorithm. This function
#' manages the random orderings and loops.
#' @param mydata A data matrix with rows as observations.
#' @param featiter The number of iterations of variable selection
#' @param clustiter The number of random ordering, for which to apply SUGS.
#' @param intfeatures A binary matrix of the intial variable set, probably chosen using function \code{pSelect}.
#' See documentation for \code{pSelect} for more details.
#' @param Model The method used for Model select, either PML, ML or Both. If you select both
#' the PML will be used to perform model selection.
#' @param mu_0 The mean hyperparameter, default is the column means of the data matrix.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{2 * (D + 2)}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the deault value is a tenth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(1, 5, 15, 30, 50, 100)}.
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
#' @return A matrix of cluster allocations for the number of iterations, a vector of either log PML, log ML or both,
#' the number of clusters at each interations, the random ordering used and the last function output for these.

runsugsvarsel<-function(mydata,
                        featiter,
                        clustiter,
                        intfeatures,
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
                        ) {

  N <- nrow(mydata)
  D <- ncol(mydata)
  T <- featiter

  #produce random orderings of the data
  rand <- lapply(seq(1:(clustiter-1)), function(x) sample(nrow(mydata)))
  rand <- c(list(seq(1:N)), rand)

  res <- bplapply(rand, function(x, mydata, featInt, Model, T, mu_0, lambda_0, nu_0, S_0,
                                 betaHat, a, b, w) {

                          suppressMessages(library(sugsvarsel))
                          for(t in 1:T){
                            res <- sugsComp(mydata[x,], featInt, Model, mu_0 = mu_0,
                                          lambda_0 = lambda_0, nu_0 = nu_0, S_0 = S_0,
                                          betaHat = betaHat, a = a, b = b, w = w) #call sugs
                            featInt <- res$features
                          }
                          return(res)

  }, mydata = mydata, featInt = intfeatures, Model = Model, T = featiter,
     mu_0 = mu_0, lambda_0 = lambda_0, nu_0 = nu_0, S_0 = S_0,
     betaHat = betaHat, a = a, b = b, w = w,  BPPARAM = BPPARAM)

  member <- matrix(unlist(lapply(res, function(x) x$member)), ncol = N, byrow = TRUE)
  clusters <- unlist(lapply(res, function(x) x$K))
  ordering <- matrix(unlist(rand), ncol = N, byrow = TRUE)
  features <- matrix(unlist(lapply(res, function(x) x$features)), ncol = ncol(mydata), byrow = TRUE)

  ML <- unlist(lapply(res, function(x) x$ML))

  if (Model=="PML") {
       LPML <- unlist(lapply(res, function(x) x$LPML))
    } else if(Model=="Both"){
       LPML <- unlist(lapply(res, function(x) x$LPML))
       ML <- unlist(lapply(res, function(x) x$ML))
    } else{
       ML <- unlist(lapply(res, function(x) x$ML))
    }

  if(Model=="PML"){
    return(list(member=member, clusters=clusters, LPML=LPML, ordering = ordering, features=features))
  } else if(Model=="Both"){
    return(list(member=member, clusters=clusters, LPML=LPML, ML=ML, ordering = ordering, features=features))
  } else{
    return(list(member=member, clusters=clusters, ML=ML, ordering = ordering, features=features))

  }

}
