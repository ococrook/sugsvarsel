#' A function to partition the variable into relevant and irrelevant
#'
#' @inheritParams compvarMarg
#' @inheritParams sugsclustMarg
#' @param varMarg The log marginal probability of a cluster being relevant
#' @param lognullMarg The log marginal probability of a variable belonging to the null model
#' @param w A numerical vector of length 2 giving the prior probability of a variable being
#' irrelevant or relevant. The first slot is irrelevant the second relevant.
#'
#' @return A binary vector indicating whether a variable is relevant (1) or irrelevant (0)

sugsvarAlloc <- function(D, varMarg ,lognullMarg, intfeature, w){

  if ( length(w) != 2) {
    print("Incorrect numbers supplied for feature priors")
  }

  featon  <- matrix(0, D)
  featoff <- matrix(0, D)

  for (d in 1:D){

    featon[d]  <- exp(log(w[2]) + varMarg[d])  #probability of feature being switch on
    featoff[d] <- exp(log(w[1]) + lognullMarg[d])    #probability of feature being switch off
  }

  #normlisation
  ff <- featon + featoff  #denominator

  featon <- featon/ff
  featoff <- featoff/ff

  #greedy assignment of indicators
  intfeature <- 1*(varMarg > lognullMarg)

  if(sum(is.na(intfeature) > 0)){
    print(intfeature)
    stop()
  }


  return(intfeature=intfeature)

}

#' A function to partition the variable into relevant and irrelevant in one dimension
#'
#' @inheritParams compvarMarg
#' @inheritParams sugsclustMarg
#' @param varMarg The log marginal probability of a cluster being relevant
#' @param lognullMarg The log marginal probability of a variable belonging to the null model
#' @param w A numerical vector of length 2 giving the prior probability of a variable being
#' irrelevant or relevant. The first slot is irrelevant the second relevant.
#'
#' @return A binary vector indicating whether a variable is relevant (1) or irrelevant (0)



sugsvarAlloc1D<-function(varMarg ,lognullMarg, w){

  if(length(w)!=2){
    print("Incorrect numbers supplied for feature priors")
  }

  featon  <- exp(log(w[2]) + varMarg)  #probability of feature being switch on
  featoff <- exp(log(w[1]) + lognullMarg)    #probability of feature being switch off

  #normlisation
  ff <- featon + featoff  #denominator


  #greedy assignment of indicators
  if (featon >= featoff) {
    intfeature <- 1
  } else {
    intfeature <- 0
  }

  if (sum(is.na(intfeature) > 0)) {
    print(intfeature)
    stop()
  }


  return(intfeature=intfeature)

}
