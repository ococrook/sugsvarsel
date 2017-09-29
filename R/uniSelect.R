#univariate selection strategy

uniSelect<-function(mydata, iter){
  
  X<-mydata
  D<-ncol(mydata)
  N<-nrow(mydata)
  
  #priors
  lambda_0<-0.01
  S_0<-matrix(0,D)
  for(d in 1:D){
    S_0[d]<-var(X[,d])/10
  }  
  S_0<-rbind(c(S_0))
  nu_0 <- D
  mu_0 <- colMeans(X) 
  alpha<-0.01
  
  
  #compute null view
  lognullMarg<-nullView(mydata,D,N, lambda_0, nu_0, mu_0, S_0)
  
 
  
  gammafeat<-array(0,c(D,D))
  
  for(d in 1:D){
    
    
    SUGSres<-runSUGS1D(iter,X[,d])
    
    
    #compute feature marginals for best clustering
    X<-mydata[SUGSres$ordering[which.max(SUGSres$LPML),],]
    D<-ncol(mydata)
    N<-nrow(mydata)
    
    #get best cluster
    numclust<-max(unique(SUGSres$member[which.max(SUGSres$LPML),]))
    n<-as.matrix(table(SUGSres$member[which.max(SUGSres$LPML),]))
    clustlabels<-SUGSres$member[which.max(SUGSres$LPML),]
    
    #compute posterior for best clusering
    
    posterior<-SUGSstatsfull(n,nu_0, lambda_0, mu_0,S_0,clustlabels,X, numclust)
    
    lambda<-posterior$lambda
    m<-posterior$m
    nu<-posterior$nu
    S<-posterior$S
    
    featMarg<-gammaMarg(X,numclust,D,n,lambda_0, nu_0, S_0, cluster,m,nu,lambda,S)
    gammafeat[d,]<-featureMarginalsSUGS(D,featMarg ,lognullMarg, rep(1,D), w)
    
  }
  
  return(gammafeat) 
}