#' The SUGS clustering algorithm by Wang and Dunson (2011) for 1 dimensional data
#'
#' @param mydata Data matrix with observations as rows
#' @param Model Character string indicating whether to use PML, ML or
#' both for model selection. "Both" defaults to PML.
#' @param iter The number of random orderings of the observations for which to run
#' sugs
#'
#' @return A matrix of cluster allocation, K the number of clusters, a matrix
#' indicating the number of observations allocated to each cluster. The value of
#' the model selection criteria either log PML, log ML or both and the random orderings used

runsugs1D<-function(iter, mydata, Model){

  N <- length(mydata)
  T <- iter

  w <- c(0.5, 0.5) #prior for featuring being swithed off/on

  member <- matrix(0, T, N)
  clusters <- matrix(0, T)
  ordering <- matrix(0, T, N)

    if (Model=="PML"){
      LPML <- matrix(0, T)
    } else if(Model=="Both"){
      LPML <- matrix(0, T)
      ML <- matrix(0, T)
    } else{
      ML <- matrix(0, T)
    }
  SUGSres <- SUGS1D(mydata, Model=Model) #initial sugs run

  member[1, ] <- SUGSres$member
  clusters[1, ] <- SUGSres$K
  ordering[1, ] <- c(seq(1:N))

    if (Model=="PML") {
      LPML[1] <- SUGSres$LPML
    } else if (Model=="Both"){
      LPML[1] <- SUGSres$LPML
      ML[1] <- SUGSres$ML
    } else{
      ML[1] <- SUGSres$ML
    }

  for(t in 2:T){

    rand <- sample(length(mydata)) #reorder data

    SUGSres <- sugs1D(mydata[rand], Model) #call sugs

    #store
    member[t, ]   <- SUGSres$member
    clusters[t, ] <- SUGSres$K
    ordering[t, ] <- rand

      if (Model=="PML") {
        LPML[t] <- SUGSres$LPML
      } else if (Model=="Both"){
        LPML[t] <- SUGSres$LPML
        ML[t] <- SUGSres$ML
      } else{
        ML[t] <- SUGSres$ML
      }
  }

  #get priors which where used (all fixed)
  priors <- SUGSres$priors
  lambda_0 <- priors$lambda_0
  mu_0 <- priors$mu_0
  nu_0 <- priors$nu_0
  S_0 <- priors$S_0

  if (Model=="PML") {

    #compute feature marginals for best clustering
    X <- mydata[ordering[which.max(LPML), ]]
    D <- 1
    N <- length(mydata)

    #get best cluster
    numclust <- max(unique(member[which.max(LPML), ]))
    n <- as.matrix(table(member[which.max(LPML), ]))
    clustlabels <- member[which.max(LPML), ]
  } else if (Model=="Both"){

    #compute feature marginals for best clustering
    X <- mydata[ordering[which.max(LPML), ]]
    D <- 1
    N <- length(mydata)

    #get best cluster
    numclust <- max(unique(member[which.max(LPML), ]))
    n <- as.matrix(table(member[which.max(LPML), ]))
    clustlabels <- member[which.max(LPML),]
  } else { #marginal likelihod case

    #compute feature marginals for best clustering
    X <- mydata[ordering[which.max(ML), ]]
    D <- 1
    N <- length(mydata)

    #get best cluster
    numclust <- max(unique(member[which.max(ML), ]))
    n <- as.matrix(table(member[which.max(ML), ]))
    clustlabels <- member[which.max(ML), ]

  }

  #compute posterior for best clusering
  posterior <- compsugsStats1D(n, nu_0, lambda_0, mu_0, S_0, clustlabels, X, numclust)

  #update parameters
  lambda <- posterior$lambda
  m  <- posterior$m
  nu <- posterior$nu
  S  <- posterior$S


  #generate marginals from NULL view (no clustering structure)
  lognullMarg <- nullView1D(X, D, N, lambda_0, nu_0, mu_0, S_0)


  #feature selection step
  varMarg <- compvarMarg(X, numclust, D, n, lambda_0, nu_0, S_0, cluster,m,nu,lambda,S)
  intfeature <- sugsvarAlloc1D(varMarg ,lognullMarg, w)


  if(Model=="PML"){
    return(list(member=member, clusters=clusters, LPML=LPML, ordering = ordering, intfeature=intfeature))
  }else if(Model=="Both"){
    return(list(member=member, clusters=clusters, LPML=LPML, ML=ML, ordering = ordering, intfeature=intfeature))
  } else{
    return(list(member=member, clusters=clusters, ML=ML, ordering = ordering, intfeature=intfeature))
  }
}
