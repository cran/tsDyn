timeVECM<-function(data, lag,m=2, r=1, include=c("const","none") )
{

#check args, nane var
  y <- as.matrix(data)
  T <- nrow(y) 		#Size of original sample
  k <- ncol(y) 		#Number of variables
  p<-lag
  t<-T-p -1		#Size of end sample
  include<-match.arg(include)
#   if(m<1 | m>T) stop("M Should be greater than 0 and less than T\n")  
  if(is.null(colnames(data)))
    colnames(data)<-paste("Var", c(1:k), sep="")

mycheb<-function(y,p,chebdim){
  T<-nrow(y)
  if(chebdim==0){
    ret<-y[(p+1):(T-1),]
  }else {
    mat<-y[(p+1):(T-1),]
    n<-nrow(y)-p-1
    s<-seq(p+2, length.out=n)
    for(i in 1:chebdim)    mat<-cbind(mat, sqrt(2)*cos(i*pi*(s-0.5)/n)*y[(p+1):(T-1),])
    ret<-mat
  }
  return(ret)
}

#construct matrices
  Y <- y[(p+1):T,] #
  X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix
  DeltaY<-diff(y)[(p+1):(T-1),]
  DeltaX<-embed(diff(y),p+1)[,-(1:k)]
  if(include=="const") 
    DeltaX<-cbind(1,DeltaX)
  Xminus1<-mycheb(y,p,m)
  
##Estimation of cointegrating vector
#Auxiliary regression 1
  reg_res1<-lm.fit(DeltaX,DeltaY)
  u<-residuals(reg_res1)
#Auxiliary regression 2
  reg_res2<-lm.fit(DeltaX,Xminus1)
  v<-residuals(reg_res2)
#Moment matrices
  S00<-crossprod(u)
  S11<-crossprod(v)
  SSSS<-solve(S11)%*%(t(v)%*%u)%*%solve(S00)%*%(t(u)%*%v)
  eig<-eigen(SSSS)
  ve<-eig$vectors
  va<-eig$values

#normalize eigenvectors
  ve_no<-apply(ve,2, function(x) x/sqrt(t(x)%*%S11%*%x))
  ve_2<-t(t(ve_no)/diag(ve_no)) 
  ve_3<-ve_2[,1:r, drop=FALSE]
  C2 <- matrix(0, nrow = nrow(ve_2) - r, ncol = r)
  C <- rbind(diag(r), C2)
  ve_4 <- ve_3 %*% solve(t(C) %*% ve_3)

#compute A (speed adjustment)
  z0<-t(u)%*%v%*%ve_no[,1:r]%*%t(ve_no[,1:r])

  ###Slope parameters
  ECTminus1<-Xminus1%*%ve_4
  Z_final<-cbind(ECTminus1,DeltaX)
  B<-t(DeltaY)%*%Z_final%*%solve(t(Z_final)%*%Z_final)		#B: OLS parameters, dim 2 x npar

###naming of variables and parameters
  if(m>0)
    rownames(ve_4)<-c(colnames(data), paste(colnames(data), "Cheb", rep(1:m,each=k)))
  else
    rownames(ve_4)<-colnames(data)
  npar<-ncol(B)*nrow(B)
  rownames(B)<-paste("Equation",colnames(data))
  LagNames<-c(paste(rep(colnames(data),length(1:p)), -rep(1:p, each=k)))
  ECT<-paste("ECT", 1:r, sep="")

if(include=="const")
	Bnames<-c(ECT,"Intercept", LagNames)
else
	Bnames<-c(ECT,LagNames)
  colnames(B)<-Bnames

###fitted, residuals, regressors matrix
  fitted<-Z_final%*%t(B)
  res<-DeltaY-fitted
  naX<-rbind(matrix(NA, ncol=ncol(Z_final), nrow=T-t), Z_final)
  rownames(naX)<-rownames(data)
  colnames(naX)<-Bnames
  YnaX<-cbind(data, naX)
  

###Return outputs
  model.specific<-list()
  model.specific$nthresh<-0
  model.specific$r<-r
  model.specific$S00<-S00
  model.specific$lambda<-eig$values
  

  z<-list(residuals=res,  coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", fitted.values=fitted, model.x=Z_final,lag=lag, model=YnaX, df.residual=t-npar/k, model.specific=model.specific, coint=ve_4)
  class(z)<-c("VaryVECM","VECM", "nlVar")
  
  attr(z, "varsLevel")<-"diff"
  return(z)
}


if(FALSE) {
y<-log(barry)
timeVECM(y,1,0)
TiVaryVECM(y,1)
}
