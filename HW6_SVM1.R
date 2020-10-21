library(MASS)
library(mvtnorm)
library(e1071)
library(quadprog)
###Simulate cluster
a<-matrix(c(0.2,0.1,
            0,0.2),2,2)
Sigma<-a%*%t(a)

set.seed(6)
mu1<-c(2,2)
mu2<-c(0,0)
a1<-mvrnorm(100,mu1,Sigma)
a2<-mvrnorm(100,mu2,Sigma)
data1<-rbind(cbind(rep(1,100),a1),cbind(rep(-1,100),a2))
colnames(data1)<-c("group","x1","x2")


svm_1<-function(data){
  #data<-data1
  err<-1e-10
  X<-data[,2:3]
  y<-data[,1]
  n<-nrow(data)
  Dmat<-diag(y)%*%X%*%t(X)%*%diag(y)+ err*diag(1,n)
  dvec<-as.matrix(rep(1,n))
  Amat<-cbind(as.matrix(y),diag(1,nrow=n))
  sol<-solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,meq=1,factorized = F)
  alpha<-sol$solution
  beta<-t(X)%*%matrix(alpha*y,nrow =n)
  beta0<--0.5*(min(X[y==1,]%*%beta)+max(X[y==-1,]%*%beta))
  
  abline(-beta0/beta[2],-beta[1]/beta[2],col="green")
  abline((1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
  abline((-1-beta0)/beta[2],-beta[1]/beta[2],col="grey",lty=2)
  cat("beta0 =",beta0,"\n")
  cat("beta = ",beta,"\n")
}

plot(data1[1:100,2:3],xlim=c(-2,4),ylim=c(-2,4),col="red",xlab = "x1",ylab="x2")
points(data1[101:200,2:3],col="blue")
svm_1(data1)

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



