#####regression on indicator matrix
library(MASS)
library(mvtnorm)
#####Simulate three clusters
a<-matrix(c(0.5,0.1,
            0,0.5),2,2)
Sigma<-a%*%t(a)
set.seed(1)
mu1<-c(-2,-2)
mu2<-c(0,0)
mu3<-c(2,2)

a1<-mvrnorm(20,mu1,Sigma)
a2<-mvrnorm(20,mu2,Sigma)
a3<-mvrnorm(20,mu3,Sigma)

plot(a1,xlim=c(-3,3),ylim=c(-3,3),col="red")
points(a2,col="blue")
points(a3,col="green")

cluster<-rbind(a1,a2,a3)
x<-cbind(rep(1,60),cluster)

y<-matrix(0,nrow=60,ncol=3)
y[1:20,1]=1
y[21:40,2]=2
y[41:60,3]=3

b<-ginv(t(x)%*%x)%*%(t(x)%*%y)
yhat<-apply(x%*%b,1,which.max)

###Only get cluster1 and 3
abline
###
b0<-b[,1]-b[,3]
abline(a=-b0[1]/b0[2],b=-b0[1]/b0[3],col="grey",lty=2)

