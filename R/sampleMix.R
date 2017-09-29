#' Generate data from a mixture of normals
#'
#' @param N The number of points to be generated
#' @param weights A vector of weights for each mixture component
#' @param mean A list containing the means of each component
#' @param cov A list of covariance matricies for each component
#'
#' @return A list returning the weights, the sampled values and the component allocations.
#'
#'
sampleMix<- function(N, weights, mean, cov){

  components <- sample(1:length(weights), weights, size = N, replace=TRUE)


  data <- matrix(unlist(lapply(cov[components],  function(x)
    rmnorm(n = 1, mean = rep(0, length(mean[[1]])), varcov = x))),
                        ncol = length(mean[[1]]), byrow = TRUE) +
    matrix(unlist(mean[components]), ncol = length(mean[[1]]), byrow = TRUE)


  return(list(weights = weights, data = data, components = components))
}

#' Generate data from a mixture of normals with diagonal structure
#'
#' @param N The number of points to be generated
#' @param weights A vector of weights for each mixture component
#' @param mean A list containing the means of each component
#' @param cov A list of covariance as a vector for each component
#'
#' @return A list returning the weights, the sampled values and the component allocations.
#'
#'


samplemixDiag <- function(N, weights, mean, cov){

  components <- sample(1:length(weights), weights, size = N, replace=TRUE)


  data <- matrix(unlist(lapply(cov[components],  function(x)
    rnorm(n = length(mean[[1]]), mean = rep(0, length(mean[[1]])), sd = x))),
    ncol = length(mean[[1]]), byrow = TRUE) +
    matrix(unlist(mean[components]), ncol = length(mean[[1]]), byrow = TRUE)

  return(list(weights = weights, data = data, components = components))
}

#' Generate data from a mixture of normals with additionally noisy variables
#' @param N The number of observations to be generated
#' @param weights A vector of weights for each mixture component
#' @param mean A list containing the means for each component
#' @param cov A list containing covariance matrix or variance vectors. See paramater \code{struc}
#' @param p The number of noisy variables
#' @param meanNoise A vector of the means of each of the noisy variables
#' @param covNoise A vector of standard deviations
#' @param struc A character string indicating the structure of the covariance matrix. Diagonal is selected by default
#' and \code{cov} should be provided as a list of vectors. If "Full" is specified, then \code{cov} should be provided
#' as list of covariance matricies.
#'
#' @return A list of three of length three. The first the mixture weights,
#' the second a data matrix containing the observations, and the third a vector of the componenet
#' allocations.


samplemixNoise <- function(N, weights, mean, cov, p, meanNoise, covNoise , struc = "Diag"){

  if (struc == "Diag"){
    mix <- samplemixDiag(N, weights, mean, cov)

  } else if (struc == "Full"){
    mix <- sampleMix(N, weights, mean, cov)
  }

  noise <- matrix(rnorm(p * N, mean = meanNoise, sd = covNoise), ncol = p, byrow = TRUE)
  mix$data <-cbind(mix$data, noise)

  return(mix = mix)

}

