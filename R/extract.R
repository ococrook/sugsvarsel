#'This is a helper function for formatting data outputs from sugsvarsel
#'
#' @param X A data matrix with rows as observations.
#' @param results The output provided from runsugsvarsel
#'
#' @return A list including a vector of log marginal likelihoods, a matrix of memberships allocations, the random orderings
#' use in the algorithm a matrix of feature allocations.

extract<-function(X, results){

  ML <- matrix(0, 1)
  member <- matrix(0, nrow(X))
  ordering <- matrix(0, nrow(X))
  features <- matrix(0, ncol(X))
  D <- length(results)

  for(d in 1:D){
    ML <- c(ML, results[[d]]$ML)
    member <- cbind(member, t(results[[d]]$member))
    ordering <- cbind(ordering, t(results[[d]]$ordering))
    features <- cbind(features, t(results[[d]]$feat))
  }
  ML <- ML[-1]
  member <- member[, -1]
  ordering <- ordering[, -1]
  features <- features[, -1]


  return(list(ML=ML, member=member, ordering=ordering, features=features))
}
