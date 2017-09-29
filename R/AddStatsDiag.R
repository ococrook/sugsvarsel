#' Function to sequentially update statistics and posterior in SUGS (Gaussian mixtures)
#' @param x The observation currently under considetion
#' @param K The number of currently occupied clusters
#' @param n The vector indicating the number of observations in each cluster
#' @param x_bar The sample mean statistics
#' @param SCL The sample variable statistics
#' @param m The current posterior mean
#' @param lambda The current posterior mean variance
#' @param S The current posterior scale vector
#' @param nu The current posterior degrees of freedom
#' @param beta The current posterior concentration parameter (not currently used in this implementation)
#' @param S_0 The prior scale vector
#' @param lamda_0 The prior mean variance
#' @param mu_0 The prior mean
#' @param nu_0 The prior degrees of freedom
#'
#' @return Returns the updated statistics and posteriors.

addStatsDiag <- function(x, K, clustnew, n, x_bar, SCL, m, lambda, S, nu, beta, S_0, lambda_0, mu_0, nu_0){

  if (K > 1) {

    if (n[clustnew]==0){ #case new cluster
      x_bar[clustnew, ] <- x
      n[clustnew] <- n[clustnew] + 1
      SCL[clustnew, ] <- x^2
      m[clustnew, ] <- (lambda_0 * mu_0 + n[clustnew] * x_bar[clustnew, ])/(lambda_0 + 1)
      lambda[clustnew] <- lambda_0 + 1
      nu[clustnew] <- nu_0 + 1
      S[clustnew, ] <- (SCL[clustnew, ]/nu[clustnew] + nu_0 * S_0/nu[clustnew] - (lambda[clustnew]/nu[clustnew]) * m[clustnew,]^2 + (lambda_0/nu[clustnew]) * mu_0^2)

    } else{ #general case
      x_bar[clustnew, ] <- (n[clustnew] * x_bar[clustnew, ] + x)/(n[clustnew] + 1)
      n[clustnew] <- n[clustnew] + 1
      SCL[clustnew,] <- SCL[clustnew,]+x^2
      m[clustnew,] <- (lambda[clustnew] * m[clustnew, ] + x)/(lambda[clustnew] + 1)
      lambda[clustnew] <- lambda[clustnew] + 1
      nu[clustnew] <- nu[clustnew] + 1
      S[clustnew,] <- (SCL[clustnew,]/nu[clustnew] + nu_0 * S_0/nu[clustnew] - (lambda[clustnew]/nu[clustnew]) * m[clustnew, ]^2 +(lambda_0/nu[clustnew]) * mu_0^2)

    }
  } else { #case only 1 cluster

    x_bar <- (n[clustnew] * x_bar + x)/(n[clustnew] + 1)
    n <- n + 1
    SCL <- SCL + x^2
    m <- (lambda * m + x)/(lambda + 1)
    lambda <- lambda + 1
    nu <- nu + 1
    S <- SCL/nu[clustnew] + (nu_0 * S_0/nu[clustnew]) - (lambda/nu[clustnew]) * m^2 + (lambda_0/nu[clustnew]) * mu_0^2
  }




  return(list(x_bar=x_bar, n=n, SCL=SCL, m=m, lambda=lambda, S=S, nu=nu))

}
