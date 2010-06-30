VECM.rankEstim<-function(data, lag.max=6, include = c( "const", "trend","none", "both")){
  include<-match.arg(include)
  models<-vector("list", lag.max)
  k<-ncol(data)
  T<-nrow(data)

  for(i in 1:lag.max) models[[i]]<-VECM(data, lag=i, include=include, estim="ML")
  allIC<-function(x,r) {
    A<-matrix(NA, ncol=1+k, nrow=4)
    for(j in 0:k)A[,j+1] <-c(AIC(x, r=j), AIC(x, r=j, k=log(T)), AIC(x, r=j, k=2*log(log(T))), AIC(x, r=j, k=(log(T)+2*log(log(T)))/2))
    dimnames(A)<-list(c("AIC", "BIC", "HQ", "LLIC"), 0:k)
    A
  }
  l<-lapply(models, allIC)

  which.min.index<-function(x) which(x==min(x), arr.ind = TRUE)
  AIC_best<-which.min.index(sapply(l, function(x) x["AIC",]))
  BIC_best<-which.min.index(sapply(l, function(x) x["BIC",]))
  HQ_best<-which.min.index(sapply(l, function(x) x["HQ",]))
  LLIC_best<-which.min.index(sapply(l, function(x) x["LLIC",]))
  
  ALL<-rbind(AIC_best, BIC_best, HQ_best, LLIC_best)
  dimnames(ALL)<-list(c("AIC", "BIC", "HQ", "LLIC"), c("Rank", "Lag"))

##Return results
  res<-list()
  res$lag.max<-lag.max
  res$bests<-ALL
  class(res)<-c("rankEstim")
  
  return(res)
}

print.rankEstim<-function(x,...){
cat("Rank and lag selection by IC\n\n")
print(t(x$bests))
}


if(FALSE){
library(tsDyn)
data(zeroyld)

est<-VECM.rankEstim(zeroyld, lag.max=6)
class(est)
print(est)

##Test with var
library(vars)
data(Canada)
VECM.rankEstim(Canada, lag.max=12)

## To figure out...
b<-VECM(data, lag=1, include=include, estim="ML")
b2<-VECM(data, lag=1, include=include, estim="2OLS")
AIC(b,r=0)
AIC(b2,r=0)
logLik(b)
logLik(b2)

##Check with vars
var_vars<-VAR(Canada, p=2)
var_my<-lineVar(Canada, lag=2)
logLik(var_vars)
logLik(var_my)

AIC(var_vars)
AIC(var_my)
logLik(lineVar(Canada, 4, include="none"))

logLik(lineVar(Canada, 3, I="diff", include="none"))

sum(lineVar(Canada, 4)$residuals[-c(1,2)]^ 2)
deviance(lineVar(Canada, 5, I="diff", include="none"))



summary(ca.jo(Canada))
}


