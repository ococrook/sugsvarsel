#' A function to compute statistics and posterior given clustering and priors
#'
#' @inheritParams addStatsDiag
#' @inheritParams sugsComp
#' @param member A vector of cluster membership allocaations
#' @param numclust The total number of occupied clusters
#'
#' @return The posterior hyperparameters for the normal inverse chisquared distribution


compsugsStats<-function(n, nu_0, lambda_0, mu_0, S_0, member, X, numclust){

  #compute statistics
  X <- X
  K <- numclust
  D <- ncol(X)

  x_bar <- matrix(0, K, D)
  SCL <- matrix(0, K, D)
  nu  <- matrix(0, K)
  lambda <- matrix(0, K)
  m <- matrix(0 ,K, D)
  S <- matrix(0, K, D)

  for (j in 1:K){
    if (n[j]==1) {
      x_bar[j, ] <- X[member==j, ]    #colsums doesn't work if height 1
      SCL[j,] <- (X[member==j, ])^2  #case 1 member in cluster
    }else if(n[j]==0){
      x_bar[j, ] <- 0           #case where there are no points in cluster
      SCL[j, ] <- 0
    }else{
      x_bar[j, ] <- colSums(X[member==j, ])/n[j]
      for(i in 1:n[j]){
        SCL[j, ] <- colSums((X[member==j,]^2)) #colSums doesn't work in dim=1
      }
    }
  }
  #parameter updates
  for(j in 1:K){
    nu[j] <- nu_0 + n[j]
    lambda[j] <- lambda_0 + n[j]
    m[j,] <- ((lambda_0 * mu_0) + (n[j] * x_bar[j, ]))/lambda[j]
    S[j,] <- SCL[j, ]/nu[j] + nu_0 * S_0/nu[j] - (lambda[j] * m[j, ]^2/nu[j]) + (lambda_0 * mu_0^2/nu[j])
  }
  return(list(m=m, lambda=lambda, nu=nu, S=S))
}


#' A function to compute statistics and posterior given clustering and priors
#'
#' @inheritParams addStatsDiag
#' @inheritParams sugsComp
#' @param member A vector of cluster membership allocaations
#' @param numclust The total number of occupied clusters
#'
#' @return The posterior hyperparameters for the normal inverse chisquared distribution

compsugsStats1D<-function(n, nu_0, lambda_0, mu_0, S_0, member, mydata, numclust){

  #compute statistics
  X <- mydata
  K <- numclust

  x_bar <- matrix(0, K)
  SCL <- matrix(0, K)
  nu <- matrix(0, K)
  lambda <- matrix(0, K)
  m <- matrix(0, K)
  S <- matrix(0, K)

  for (j in 1:K) {
    if (n[j]==1) {
      x_bar[j] <- X[member==j]    #colsums doesn't work if height 1
      SCL[j] <- (X[member==j])^2  #case 1 member in cluster
    } else if (n[j]==0){
       x_bar[j] <- 0           #case where there are no points in cluster
       SCL[j] <- 0
    }else{
       x_bar[j] <- sum(X[member==j])/n[j]

        for(i in 1:n[j]){
         SCL[j] <- sum((X[member==j]^2)) #colSums doesn't work in dim=1
        }
    }
  }

  #parameter updates
  for (j in 1:K) {
    nu[j] <- nu_0 + n[j]
    lambda[j] <- lambda_0 + n[j]
    m[j] <- (lambda_0 * mu_0 + n[j] * x_bar[j])/lambda[j]
    S[j] <- SCL[j]/nu[j] + nu_0 * S_0/nu[j] - (lambda[j] * m[j]^2/nu[j]) + (lambda_0 * mu_0^2/nu[j])
  }

  return(list(m=m, lambda=lambda, nu=nu, S=S))
}
