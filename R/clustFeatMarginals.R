#cluster marginals with feature selection

ClustFeatMarginals<-function(x,K,N,D,n,beta, m, nu,lambda,S, gammaFeat){
  
  probz<-matrix(0,K)
  probx<-matrix(0,K,D)
  probxFeat<-matrix(0,K,D)
  probmember<-matrix(0,K)
  
  scale<-matrix(0,K,D)
  if(K>1){
    for(j in 1:K){
      probz[j]<-(n[j])/(N+beta-1)  # marginal for z
      
      #compute scale matrix
      scale[j,]<-((1+lambda[j])*S[j,]/lambda[j])^(1/2)
      
      probx[j,]<-dt.scaled(x, mean = m[j,], df = nu[j], sd = scale[j,], log = TRUE)
      
      probxFeat[j,]<-gammaFeat*probx[j,] #remove features that are switched off by multiplcation by 0
      probmember[j]<-exp(log(probz[j])+sum(probxFeat[j,])) #multiply, normalise later , exp(sum(log)) for stability
    }
  } else{
    probz<-(n)/(N+beta-1)
    
    #compute scale matrix
    scale<-((1+lambda)*S/lambda)^(1/2)
    
    probx<-dt.scaled(x, mean = m, df = nu, sd = scale, log = TRUE)
    probxFeat<-gammaFeat*probx #remove features that are switched off by multiplcation by 0
    probmember<-exp(log(probz)+sum(probxFeat)) #multiply, normalise later , exp(sum(log)) for stability
    
  }  
  
  #unit test
  if(sum(is.na(probmember))>0){ 
    print("marginals error")
    print(S)
    print(probx)
    stop(sum(is.na(probmember))>1)
  }
  
  return(list(probmember=probmember, logclstMarg=probx))
}