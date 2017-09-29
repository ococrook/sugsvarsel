#' The SUGS clustering algorithm by Wang and Dunson (2011)
#'
#' @import BiocParallel
#'
#' @param mydata Data matrix with observations as rows
#' @param Model Character string indicating whether to use PML, ML or
#' both for model selection. "Both" defaults to PML.
#' @param iter The number of random orderings of the observations for which to run
#' sugs.
#' @param mu_0 The mean hyperparameter, default is the column means of the data matrix.
#' @param lambda_0 The variance of the Guassian mean prior, the dafault value is \code{0.01}.
#' @param nu_0 The degrees of freedom hyperparameter, the default value is \code{2 * (D + 2)}, where \code{D} is the number of variables.
#' @param S_0 The scale hyperparamter, the deault value is a tenth of the column variance of the data matrix.
#' @param betaHat A grid of hyperparameters for the dirichlet concentration parameter, the default is \code{c(1, 5, 15, 30, 50, 100)}.
#' @param a The scale of the gamma prior for the dirichlet concentration parameter, the dafault value is \code{10}.
#' @param b The rate of the gamma prior for the dirichlet concentration parameter, the default value is \code{1}.
#' @param BPPARAM Support for parallel processing using the
#' \code{BiocParallel} infrastructure. When missing (default),
#'     the default registered \code{BiocParallelParam} parameters are
#'     used. Alternatively, one can pass a valid
#'     \code{BiocParallelParam} parameter instance: \code{SnowParam},
#'     \code{MulticoreParam}, \code{DoparParam}, \ldots see the
#'     \code{BiocParallel} package for details. To revert to the
#'     origianl serial implementation, use \code{serialParam}.
#'
#'
#' @return A matrix of cluster allocation, K the number of clusters, a matrix
#' indicating the number of observations allocated to each cluster. The value of
#' the model selection criteria either log PML, log ML or both and the random orderings used


runSugs <- function(iter,
                    mydata,
                    Model,
                    mu_0 = NULL,
                    lambda_0 = 0.01,
                    nu_0 = NULL,
                    S_0 = NULL,
                    betaHat = c(1, 5, 15, 30, 50, 100),
                    a = 10,
                    b = 1,
                    BPPARAM = bpparam()
                    ) {

  #defensive tests
  stopifnot(is.numeric(mydata))
  stopifnot(is.numeric(iter))
  if (iter <= 1) {
    stop("iter must be an integer of at least 2 ")
  }
  if (!(Model %in% c("PML", "ML", "Both"))) {
    stop("Please state a valid model selection criteria")
  }

  N <- nrow(mydata)
  T <- iter

  member <- matrix(0, T, N)
  clusters <- matrix(0, T)
  ordering <- matrix(0, T, N)

  if (Model=="PML") {
    LPML <- matrix(0, T)
  } else if(Model=="Both") {
    LPML <- matrix(0, T)
    ML <- matrix(0, T)
  } else{
    ML <- matrix(0, T)
  }

  rand <- lapply(seq(1:(T-1)), function(x) sample(nrow(mydata)))
  rand <- c(list(seq(1:N)), rand)

  res <- bplapply(rand, function(x, mydata, Model, mu_0, lambda_0, nu_0, S_0,
                                 betaHat, a, b){
                        suppressMessages(library(sugsvaRsel))
                        sugs(mydata[x,], Model, mu_0 = mu_0, lambda_0 = lambda_0,
                              nu_0 = nu_0, S_0 = S_0, betaHat = betaHat, a = a, b = b)
                          },
                  mydata = mydata, Model = Model, mu_0 = mu_0, lambda_0 = lambda_0,
                  nu_0 = nu_0, S_0 = S_0, betaHat = betaHat, a = a, b = b, BPPARAM = BPPARAM )

  member <- matrix(unlist(lapply(res, function(x) x$member)), ncol = N, byrow = TRUE)
  clusters <- unlist(lapply(res, function(x) x$K))
  ordering <- matrix(unlist(rand), ncol = N, byrow = TRUE)

  if (Model=="PML") {
    LPML <- unlist(lapply(res, function(x) x$LPML))
  } else if (Model=="Both"){
    LPML <- unlist(lapply(res, function(x) x$LPML))
    ML <- unlist(lapply(res, function(x) x$ML))
  } else{
    ML <- unlist(lapply(res, function(x) x$ML))
  }

  if(Model=="PML"){
    return(list(member=member, clusters=clusters, LPML=LPML, ordering = ordering))
  }else if(Model=="Both"){
    return(list(member=member, clusters=clusters, LPML=LPML, ML=ML, ordering = ordering))
  } else{
    return(list(member=member, clusters=clusters, ML=ML, ordering = ordering))
  }
}
