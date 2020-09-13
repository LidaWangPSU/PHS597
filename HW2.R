library(glmnet)

####import data from the internet
X<-data.frame(X1 = c(7, 1, 11, 11, 7, 11, 3, 1, 2, 21, 1, 11, 10), 
              X2 = c(26,29, 56, 31, 52, 55, 71, 31, 54, 47, 40, 66, 68), 
              X3 = c(6, 15, 8, 8, 6,9, 17, 22, 18, 4, 23, 9, 8), 
              X4 = c(60, 52, 20, 47, 33, 22, 6, 44, 22, 26,34, 12, 12))
Y<-as.matrix((c(78.5, 74.3, 104.3, 87.6, 95.9, 109.2, 102.7, 72.5, 93.1, 115.9, 83.8, 113.3, 109.4)))
X<-as.matrix((X))


###use cyclic coordinate descent algorithm

###define soft threshold function
soft_threshold<-function(r,lambda){
  if(r>lambda){
    return(r-lambda)
  }
  if(r<(-lambda)){
    return(r+lambda)
  }
  if(r<=lambda&r>=(-lambda)){
    return(0)
  }
}

###define coordinate descent lasso
cd_lasso<-function(x,y,lambda,iteration,intercept){
  m<-nrow(x)
  n<-ncol(x)
  theta<-as.matrix(rep(1, times=n))
  if(intercept){
  intercept<-sum(y-x%*%theta)/m
  for(e in 1:iteration){
    for(i in 1:n){
      r<-matrix(nrow=m,ncol=1)
      z<-matrix(nrow=m,ncol=1)
      rk<-sum((y-x%*%theta+x[,i]*theta[i]-intercept)*x[,i])
      zk<-sum(x[,i]^2)
      theta[i]<-soft_threshold(rk,lambda)/zk
      intercept<-sum(y-x%*%theta)/m
    }
  }
  return(rbind(intercept,theta))}
  else{
    for(e in 1:iteration){
      for(i in 1:n){
        r<-matrix(nrow=m,ncol=1)
        z<-matrix(nrow=m,ncol=1)
        rk<-sum((y-x%*%theta+x[,i]*theta[i])*x[,i])
        zk<-sum(x[,i]^2)
        theta[i]<-soft_threshold(rk,lambda)/zk
      }
    }
    return(theta)}

}


###Compare with glmnet
fit1<-glmnet(X,Y,family='gaussian',standardize = FALSE,lambda=1,intercept=T)
coef(fit1)
fit2<-cd_lasso(X,Y,lambda=1,iteration=10000,intercept=T)
fit2





