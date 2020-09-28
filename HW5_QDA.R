###QDA
library(MASS)
library(mvtnorm)
###Simulate cluster
a1<-matrix(c(0.85,0.1,
            0,0.5),2,2)
Sigma1<-a1%*%t(a1)

a2<-matrix(c(1.5,0.2,
             0,1.5),2,2)
Sigma2<-a1%*%t(a2)

a3<-matrix(c(0.75,0.1,
             0,0.75),2,2)
Sigma3<-a1%*%t(a3)
set.seed(1)
mu1<-c(-3,-3)
mu2<-c(0,0)
mu3<-c(3,3)


a1q<-mvrnorm(50,mu1,Sigma1)
a2q<-mvrnorm(100,mu2,Sigma2)
a3q<-mvrnorm(70,mu3,Sigma3)

mu1q=colMeans(a1q)
mu2q=colMeans(a2q)
mu3q=colMeans(a3q)
mu_q<-rbind(mu1q,mu2q,mu3q)
pi1q=nrow(a1q)/(nrow(a1q)+nrow(a2q)+nrow(a3q))
pi2q=nrow(a2q)/(nrow(a1q)+nrow(a2q)+nrow(a3q))
pi3q=nrow(a3q)/(nrow(a1q)+nrow(a2q)+nrow(a3q))
piq<-c(pi1q,pi2q,pi3q)
Sigma1q<-1/(n-1)*((t(a1q)-mu_q[1,])%*%(a1q-mu_q[1,]))
Sigma2q<-1/(n-1)*((t(a2q)-mu_q[2,])%*%(a2q-mu_q[2,]))
Sigma3q<-1/(n-1)*((t(a3q)-mu_q[3,])%*%(a3q-mu_q[3,]))

###plot clusters
plot(a1q,xlim=c(-5,5),ylim=c(-5,5),col="red")
points(a2q,col="blue")
points(a3q,col="green")


###Sphering the data 
library(rmutil)
T1<-svd(Sigma1q)
T2<-svd(Sigma2q)
T3<-svd(Sigma3q)

D1<-diag((T1$d)^(-1/2))
s1<-D1%*%t(T1$u)
a1new<-t(s1%*%t(a1q))

D2<-diag((T2$d)^(-1/2))
s2<-D2%*%t(T2$u)
a2new<-t(s2%*%t(a2q))

D3<-diag((T3$d)^(-1/2))
s3<-D3%*%t(T3$u)
a3new<-t(s3%*%t(a3q))


mu1new=colMeans(a1new)
mu2new=colMeans(a2new)
mu3new=colMeans(a3new)
mu_new<-rbind(mu1new,mu2new,mu3new)


###Standard methods

###Plots in the sphering data
plot(a1new,xlim=c(-20,20),ylim=c(-20,20),col="red")
points(a2new,col="blue")
points(a3new,col="green")
title(main="Plots for the sphering data")

#b12=mu_new[2,]-mu_new[1,]
#a12=-1/2*(mu_new[1,]+mu_new[2,])%*%(mu_new[2,]-mu_new[1,])+log(pi1q)-log(pi2q)
#abline(-a12/b12[2],-b12[1]/b12[2],col="grey",lty=2)

#b13=mu_new[3,]-mu_new[1,]
#a13=-1/2*(mu_new[1,]+mu_new[3,])%*%(mu_new[3,]-mu_new[1,])+log(pi1q)-log(pi3q)
#abline(-a13/b13[2],-b13[1]/b13[2],col="grey",lty=2)

#b23=mu_new[2,]-mu_new[3,]
#a23=-1/2*(mu_new[2,]+mu_new[3,])%*%(mu_new[2,]-mu_new[3,])+log(pi3q)-log(pi2q)
#abline(-a23/b23[2],-b23[1]/b23[2],col="grey",lty=2)



###equal distance   
mu_12<-(mu1new+mu2new)/2
k12<-(mu2new[2]-mu1new[2])/(mu2new[1]-mu1new[1])

x_12<-seq(-20,20,length=300)
y_12<-mu_12[2]+mu_12[1]/k12+x_12*(-1/k12)

bd12<-cbind(x_12,y_12)
#bd12_new<-t(ginv(s)%*%t(bd12))
abline(mu_12[2]+mu_12[1]/k12,-1/k12,col="grey",lty=2)


mu_13<-(mu1new+mu3new)/2
k13<-(mu3new[2]-mu1new[2])/(mu3new[1]-mu1new[1])

x_13<-seq(-20,20,length=300)
y_13<-mu_13[2]+mu_13[1]/k13+x_12*(-1/k13)

bd13<-cbind(x_13,y_13)
#bd13_new<-t(ginv(s)%*%t(bd13))
abline(mu_13[2]+mu_13[1]/k13,-1/k13,col="gold",lty=2)


mu_23<-(mu2new+mu3new)/2
k23<-(mu3new[2]-mu2new[2])/(mu3new[1]-mu2new[1])

x_23<-seq(-20,20,length=300)
y_23<-mu_23[2]+mu_23[1]/k23+x_23*(-1/k23)

bd23<-cbind(x_23,y_23)
#bd12_new<-t(ginv(s)%*%t(bd12))
abline(mu_23[2]+mu_23[1]/k23,-1/k23,col="black",lty=2)

legend("bottomright", legend=c("boundary between cluster red and blue", "boundary between cluster red and green","boundary between cluster blue and green"),
       col=c("grey", "gold","black"), lty=1, cex=0.8)


###Boundary on the original data
plot(a1q,xlim=c(-5,5),ylim=c(-5,5),col="red")
points(a2q,col="blue")
points(a3q,col="green")
title(main="QDA Plots for the original data")

#b=(x,y)
#(b-mu_1)^t\sigma_1^{-1}(b-mu_1)-(b-mu_2)^t\sigma_2^{-1}(b-mu_2)=C
#By calcualtion we can find a relation ship between x and y, Ay^2+By+C=0, y=(-B+-\sqrt{B^2-4AC}/2A)

transfer<-function(mu1,mu2,sigma1,sigma2,X){
  y<-matrix(nrow = length(X),ncol=2)
  for(i in 1:length(x)){
  x<-X[i]
    
  mu00<-mu1[1]
  mu01<-mu1[2]
  mu10<-mu2[1]
  mu11<-mu2[2]
  a<-sigma1[1,1]
  b<-sigma1[1,2]
  c<-sigma1[2,1]
  d<-sigma1[2,2]
  p<-sigma2[1,1]
  q<-sigma2[1,2]
  r<-sigma2[2,1]
  s<-sigma2[2,2]
  A<-d-s
  B<--2*d*mu11+2*s*mu01+b*x-b*mu10+c*x-c*mu10-q*x+q*mu00-r*x+r*mu00
  C<-a*(x-mu10)^2-p*(x-mu00)^2-b*mu11*x-c*mu11*x+q*mu01*x+r*mu01*x+d*mu11^2-s*mu01^2+b*mu10*mu11+c*mu10*mu11-q*mu01*mu00-r*mu01*mu00
  y[i,1]<-(-B-(B^2-4*A*C)^(1/2))/(2*A)
  y[i,2]<-(-B+(B^2-4*A*C)^(1/2))/(2*A)
  }
  return(y)
}

x<-seq(-5,5,length=3000)
y12<-transfer(mu1q,mu2q,ginv(Sigma1q),ginv(Sigma2q),x)

lines(cbind(x,y12[,1]),col="grey",pch=16,cex=0.3)
lines(cbind(x,y12[,2]),col="grey",pch=16,cex=0.3)


y13<-transfer(mu1q,mu3q,ginv(Sigma1q),ginv(Sigma3q),x)

lines(cbind(x,y13[,1]),col="gold",pch=16,cex=0.3)
lines(cbind(x,y13[,2]),col="gold",pch=16,cex=0.3)


y23<-transfer(mu2q,mu3q,ginv(Sigma2q),ginv(Sigma3q),x)
lines(cbind(x,y23[,1]),col="black",pch=16,cex=0.3)
lines(cbind(x,y23[,2]),col="black",pch=16,cex=0.3)

legend("bottomright", legend=c("boundary between cluster red and blue", "boundary between cluster red and green","boundary between cluster blue and green"),
       col=c("grey", "gold","black"), lty=1, cex=0.8)
