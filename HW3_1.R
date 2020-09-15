library(pls)

#Construct a positive definite matrix as the covariance matrix
a<-matrix(c(1,-0.2,0.1,
            0,1,0.2,
            0,0,1),3,3)
Sigma<-a%*%t(a)

#Simulate x and y
set.seed(1)
mu<-c(1,1,1)
X<-mvrnorm(10,mu,Sigma)
beta<-c(1,2,3)
Y<-rnorm(10,mean=beta%*%mu,sd=1)
Y<-as.matrix(Y)

X<-scale(X)
Y<-scale(Y)

partial_least_square<-function(X,Y,l){
  n<-nrow(X)
  p<-ncol(X)
  Z<-matrix(nrow=n,ncol=l)
  Fi<-matrix(nrow=p,ncol=1)
  theta<-matrix(nrow=l,ncol=1)
  X0<-X
  #m<-1
  for(m in 1:l){
    Z[,m]<-rep(0, times=n)
    for(j in 1:p){
      Fi[j]<-t(X0[,j])%*%Y
      Z[,m]<-Z[,m]+Fi[j]%*%X0[,j]
      }
    for(j in 1:p){
      X0[,j]<-X0[,j]-(t(Z[,m])%*%X0[,j]/(t(Z[,m])%*%Z[,m]))*Z[,m]
    }
  }
  
  return(Z)
}

t<-partial_least_square(X,Y,3)
t


###use package pls
fit<-plsr(Y~X,ncomp=3,scale=F)
fit$scores
t2<-fit$scores


###compare the cosine between the two results, we can find that they are in the same direction. 
sum(t2[,1]*t[,1])/(sqrt(sum(t2[,1]*t2[,1]))*sqrt(sum(t[,1]*t[,1])))
sum(t2[,2]*t[,2])/(sqrt(sum(t2[,2]*t2[,2]))*sqrt(sum(t[,2]*t[,2])))
sum(t2[,3]*t[,3])/(sqrt(sum(t2[,3]*t2[,3]))*sqrt(sum(t[,3]*t[,3])))
