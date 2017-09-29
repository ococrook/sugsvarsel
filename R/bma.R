#' Bayesian model averaging for SUGS VarSel
#'
#' @import Matrix
#'
#' @param X Data matrix with observations as rows
#' @param MLnorm A vector of normalised marginal likelihood, most likely computed from \code{occamsWindow}
#' @param Result The output from the function \code{sugsvarsel}
#'
#' @return Returns a numerical bayesian model averaged coclustering matrix
#' and numerical vector of bayesian model averaged selected variables

bma<-function(X, MLnorm, Result){

  Asim <- matrix(0, nrow(X), nrow(X))

  for(j in 1:(length(MLnorm))){
    A <- matrix(0, nrow(X), nrow(X))
    a <- Result$member[, order(Result$ML,decreasing=TRUE)[j]][invPerm(Result$ordering[, order(Result$ML,decreasing=TRUE)][, j])]
    for(i in 1:nrow(X)){
      A[i, ] <- 1 * (a==a[i])
    }
    Asim<-MLnorm[j]*A+Asim
  }
  diag(Asim) <- rep(1, ncol(Asim))

  Asim[Asim>1] <- 1

  Fsim<-matrix(0, ncol(X))
  for(j in 1:(length(MLnorm))){
    Fsim <- Fsim + MLnorm[j]*Result$features[, order(Result$ML, decreasing=TRUE)[j]]
  }
    return(list(Asim=Asim, Fsim=Fsim))
}

#' Helper function to compute occams Window for SUGS VarSel
#'
#' @param ML The vector of log MLs returned by sugsvarsel
#' @param window A numerical value indicating the size of occams window
#'
#' @return A numerical vector containg the normalised MLs.

occamsWindow<-function(ML, window){

  MLorder  <- ML[order(ML, decreasing=TRUE)]
  MLwindow <- MLorder[(MLorder[1] - MLorder) <= window]
  MLnorm <- MLwindow-min(MLwindow)
  MLnorm <- exp(MLnorm)/sum(exp(MLnorm))

  return(MLnorm=MLnorm)
}
