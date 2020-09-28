###LDA
library(MASS)
library(mvtnorm)
###Simulate cluster
a<-matrix(c(0.5,0.1,
            0,0.5),2,2)
Sigma<-a%*%t(a)
set.seed(1)
mu1<-c(-2,-2)
mu2<-c(0,0)
mu3<-c(2,2)

a1<-mvrnorm(40,mu1,Sigma)
a2<-mvrnorm(40,mu2,Sigma)
a3<-mvrnorm(40,mu3,Sigma)

plot(a1,xlim=c(-4,4),ylim=c(-4,4),col="red")
points(a2,col="blue")
points(a3,col="green")


mu11=colMeans(a1)
mu22=colMeans(a2)
mu33=colMeans(a3)
mu<-rbind(mu11,mu22,mu33)
pi1=nrow(a1)/(nrow(a1)+nrow(a2)+nrow(a3))
pi2=nrow(a2)/(nrow(a1)+nrow(a2)+nrow(a3))
pi3=nrow(a3)/(nrow(a1)+nrow(a2)+nrow(a3))
pi<-c(pi1,pi2,pi3)
Sigma1<-1/(n-3)*((t(a1)-mu11)%*%(a1-mu11)+(t(a2)-mu22)%*%(a2-mu22)+(t(a3)-mu33)%*%(a3-mu33))


###standard methods
###discriminant function
#dis_func<-function(x,k,sigma,mu,pi){
#dk<-t(x)%*%ginv(sigma)%*%mu[k,]-0.5*t(mu[k,])%*%ginv(sigma)%*%mu[k,]+log(pi[k])
#return(dk)
#}
#x<-c(1,1)
#dis_func(x,1,Sigma1,mu,pi)
#dis_func(x,2,Sigma1,mu,pi)
#dis_func(x,3,Sigma1,mu,pi)


###standard decision boundary
#b12=ginv(sigma)%*%(mu[2,]-mu[1,])
#a12=-1/2*(mu[1,]+mu[2,])%*%ginv(sigma)%*%(mu[2,]-mu[1,])
#abline(-a12/b12[2],-b12[1]/b12[2],col="grey",lty=2)

#b13=ginv(sigma)%*%(mu[3,]-mu[1,])
#a13=-1/2*(mu[1,]+mu[3,])%*%ginv(sigma)%*%(mu[3,]-mu[1,])
#abline(-a13/b13[2],-b13[1]/b13[2],col="grey",lty=2)

#b23=ginv(sigma)%*%(mu[3,]-mu[2,])
#a23=-1/2*(mu[3,]+mu[2,])%*%ginv(sigma)%*%(mu[3,]-mu[2,])
#abline(-a23/b23[2],-b23[1]/b23[2],col="grey",lty=2)


###equal distance 
library(rmutil)
###sphering the data 
T<-svd(Sigma1)
D<-diag((T$d)^(-1/2))
s<-D%*%t(T$u)
aa1<-t(s%*%t(a1))
aa2<-t(s%*%t(a2))
aa3<-t(s%*%t(a3))
mu1s=colMeans(aa1)
mu2s=colMeans(aa2)
mu3s=colMeans(aa3)


###find the boundary and transfer it back
mu_12<-(mu1s+mu2s)/2
k12<-(mu2s[2]-mu1s[2])/(mu2s[1]-mu1s[1])

x_12<-seq(-12,12,length=300)
y_12<-mu_12[2]+mu_12[1]/k12+x_12*(-1/k12)

bd12<-cbind(x_12,y_12)
bd12_new<-t(ginv(s)%*%t(bd12))
lines(bd12_new,pch=16,cex=0.2,col="grey")


mu_13<-(mu1s+mu3s)/2
k13<-(mu3s[2]-mu1s[2])/(mu3s[1]-mu1s[1])

x_13<-seq(-12,12,length=300)
y_13<-mu_13[2]+mu_13[1]/k13+x_d*(-1/k13)

bd13<-cbind(x_13,y_13)
bd13_new<-t(ginv(s)%*%t(bd13))
lines(bd13_new,pch=16,cex=0.2,col="gold")



mu_23<-(mu2s+mu3s)/2
k23<-(mu3s[2]-mu2s[2])/(mu3s[1]-mu2s[1])

x_23<-seq(-12,12,length=300)
y_23<-mu_23[2]+mu_23[1]/k23+x_23*(-1/k23)

bd23<-cbind(x_23,y_23)
bd23_new<-t(ginv(s)%*%t(bd23))
lines(bd23_new,pch=16,cex=0.2,col="black")


legend("bottomright", legend=c("boundary between cluster red and blue", "boundary between cluster red and green","boundary between cluster blue and green"),
       col=c("grey", "gold","black"), lty=1, cex=0.8)
title(main="LDA Plots")
