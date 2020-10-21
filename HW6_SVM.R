library(MASS)
library(mvtnorm)
library(e1071)
###Simulate cluster
a<-matrix(c(0.5,0.1,
            0,0.5),2,2)
Sigma<-a%*%t(a)

set.seed(6)
mu1<-c(2,2)
mu2<-c(0,0)
a1<-mvrnorm(100,mu1,Sigma)
a2<-mvrnorm(100,mu2,Sigma)
data1<-rbind(cbind(rep(1,100),a1),cbind(rep(-1,100),a2))
colnames(data1)<-c("group","x1","x2")
plot(data1[1:100,2:3],xlim=c(-2,4),ylim=c(-2,4),col="red",xlab = "x1",ylab="x2")
points(data1[101:200,2:3],col="blue")


svm_0<-function(data){
  a0<-rep(1,nrow(data))
  y<-data[,1]
  x<-data[,2:ncol(data)]
  n<-nrow(x)
  x1<-data[1:100,2:ncol(data)]
  x2<-data[101:200,2:ncol(data)]
  
  beta<-matrix(nrow=2,ncol=1)
  a<-a0
  beta[1]<-sum(a*y*x[,1])/n
  beta[2]<-sum(a*y*x[,2])/n
  min(x1%*%beta)
  max(x2%*%beta)
  s<-2/(min(x1%*%beta)-max(x2%*%beta))
  beta[1]<-beta[1]*s
  beta[2]<-beta[2]*s
  beta0<-1-min(x1%*%beta)
  #abline(-beta0/beta[1],-beta[2]/beta[1])
  support<-matrix(nrow=n,ncol=3)
  for(j in 1:n){
    if(y[j]*(x[j,]%*%beta+beta0)<=1+10^(-5)){
      a[j]=a[j]
      support[j,]<-data[j,]
    }else{
      a[j]=0
      support[j,]<-NA
    }
  }
  support<-na.omit(support)
  
  beta<-(support[1,2:3]-support[2,2:3])
  s<-2/abs(support[1,2:3]%*%beta-support[2,2:3]%*%beta)
  beta[1]<-beta[1]*s
  beta[2]<-beta[2]*s
  beta0<-1-max(support[1,2:3]%*%beta,support[2,2:3]%*%beta)
  
if(min(x1%*%beta)+beta0>=1){  
  cat("beta0=",beta0,"\n")
  cat("beta=",beta,"\n")
  abline(-beta0/beta[2],-beta[1]/beta[2],col="green")
  abline((1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
  abline((-1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
}
  
if(min(x1%*%beta)+beta0<1){
  a0<-rep(1,nrow(data))
  y<-data[,1]
  x<-data[,2:ncol(data)]
  n<-nrow(x)
  x1<-data[1:100,2:ncol(data)]
  x2<-data[101:200,2:ncol(data)]
  
  beta<-matrix(nrow=2,ncol=1)
  a<-a0
  beta[1]<-sum(a*y*x[,1])/n
  beta[2]<-sum(a*y*x[,2])/n
  min(x1%*%beta)
  max(x2%*%beta)
  s<-2/(min(x1%*%beta)-max(x2%*%beta))
  beta[1]<-beta[1]*s
  beta[2]<-beta[2]*s
  beta0<-1-min(x1%*%beta)
  #abline(-beta0/beta[1],-beta[2]/beta[1])
  
  
  support<-matrix(nrow=n,ncol=4)
  for(j in 1:n){
    support[j,1]<-y[j]*(x[j,]%*%beta+beta0)
    support[j,2:4]<-data[j,]
  }
  support[order(support[,1],decreasing=FALSE)[1:3],]
  
  support1<-support[order(support[,1],decreasing=FALSE)[1:3],2:4]
  ss<-subset(support1,support1[,1]==sum(support1[,1]))
  ss0<-subset(support1,support1[,1]==-sum(support1[,1]))
  
  beta[1]<--(ss[1,2:3]-ss[2,2:3])[2]
  beta[2]<-(ss[1,2:3]-ss[2,2:3])[1]
  s<-2/abs(ss[2,2:3]%*%beta-ss0[1,2:3]%*%beta)
  beta[1]<-beta[1]*s
  beta[2]<-beta[2]*s
  beta0<-1-max(ss[1,2:3]%*%beta,ss0[1,2:3]%*%beta)
  
  if(min(abs(max(x1%*%beta)+beta0),abs(min(x1%*%beta)+beta0))>=1){  
    cat("beta0=",beta0,"\n")
    cat("beta=",beta,"\n")
    abline(-beta0/beta[2],-beta[1]/beta[2],col="green")
    abline((1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
    abline((-1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
  }  
}
  
}

svm_0(data1)
###compare with svm
m1<-svm(group~x1+x2,data=data1,kernel="linear",scale=FALSE,type="C-classification")
coef(m1)
cf <- coef(m1)
abline(-cf[1]/cf[3], -cf[2]/cf[3], col = "red")
#plot margin 
abline(-(cf[1] + 1)/cf[3], -cf[2]/cf[3], col = "blue",lty=2)
abline(-(cf[1] - 1)/cf[3], -cf[2]/cf[3], col = "blue",lty=2)
legend("bottomleft", legend=c("boundary by SVM package", "margin by SVM package","boundary by myself","margin by myself"),
       col=c("red", "blue","green","grey"), lty=c(1,2,1,2), cex=0.8)
title(main="SVM")