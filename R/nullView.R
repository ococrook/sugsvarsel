#'A function to compute the marginal probability of belonging to the irrelevant partition.
#'
#' @param X The data matrix with observations as rows.
#' @param D The number of variables of the data matrix.
#' @param N The number of observations.
#' @param lambda_0 The prior value of the mean variance.
#' @param nu_0 The prior number of degrees of freedom.
#' @param mu_0 The prior mean.
#' @param S_0 The prior scale for the inverse-chisquared distribution.
#'
#' @return Returns, in log form, the marginal probability of belong to the irrelevant
#' partition.


nullView<-function(X, D, N, lambda_0, nu_0, mu_0, S_0) {

    lognullMarg <- matrix(0, D)

    #functions as if no clustering structure, vectorised over D
    #statistics
    x_barNULL <- colSums(X)/N
    SCNULL <- colSums(X^2)

    #posterior updates
    nuNULL <- nu_0 + N
    lambdaNULL <- lambda_0 + N
    mNULL <- (lambda_0 * mu_0 + N * x_barNULL)/lambdaNULL
    SNULL <- SCNULL/nuNULL + nu_0 * S_0/nuNULL - (lambdaNULL * mNULL^2/nuNULL) + (lambda_0 * mu_0^2/nuNULL)

    #compute log marginal likeihood under NULL view (no clustering structure)
    priorMarg <- log(lambda_0)/2 + nu_0 * log(nu_0 * S_0)/2 - lgamma(nu_0/2) #compute normalising constant of prior

    logfullMarg <- lgamma(nuNULL/2) - (N/2) * log(pi)/2 - log(lambdaNULL)/2 - nuNULL * log((nuNULL*SNULL))/2 #compute normalising constant of the NULL view

    lognullMarg <- priorMarg + logfullMarg # marginal of NULL view


  return(lognullMarg=lognullMarg)
}

#' The one-dimensional version of the nullView function
#'
#' @param X The data matrix with observations as rows.
#' @param D The number of variables of the data matrix.
#' @param N The number of observations.
#' @param lambda_0 The prior value of the mean variance.
#' @param nu_0 The prior number of degrees of freedom.
#' @param mu_0 The prior mean.
#' @param S_0 The prior scale for the inverse-chisquared distribution.
#'
#' @return Returns, in log form, the marginal probability of belong to the irrelevant
#' partition.
#'
nullView1D <- function(X, D, N, lambda_0, nu_0, mu_0, S_0){

  lognullMarg <- matrix(0, D)

  #functions as if no clustering structure, vectorised over D
  #statistics
  x_barNULL <- sum(X)/N
  SCNULL <- sum(X^2)

  #posterior updates
  nuNULL <- nu_0 + N
  lambdaNULL <- lambda_0 + N
  mNULL <- (lambda_0 * mu_0 + N * x_barNULL)/lambdaNULL
  SNULL <- SCNULL/nuNULL + nu_0 * S_0/nuNULL - (lambdaNULL * mNULL^2/nuNULL) + (lambda_0 * mu_0^2/nuNULL)

  #compute log marginal likeihood under NULL view (no clustering structure)
  priorMarg <- log(lambda_0 * (nu_0 * S_0)^(nu_0))/2 - lgamma(nu_0/2) #compute normalising constant of prior

  logfullMarg <- lgamma(nuNULL/2) - log(pi^(N/2))/2 - log(lambdaNULL)/2 - nuNULL * log((nuNULL*SNULL))/2 #compute normalising constant of the NULL view

  lognullMarg <- priorMarg + logfullMarg # marginal of NULL view


  return(lognullMarg=lognullMarg)
}
