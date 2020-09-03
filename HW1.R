library(MASS)

#####1. comapre one stage multivariate regression and 2-stage regression

#Construct a positive definite matrix as the covariance matrix
a<-matrix(c(1,-0.2,0.1,
            0,1,0.2,
            0,0,1),3,3)
Sigma<-a%*%t(a)

#Simulate x and y
set.seed(1)
mu<-c(1,1,1)
x<-mvrnorm(10,mu,Sigma)
beta<-c(1,2,3)
y<-rnorm(10,mean=beta%*%mu,sd=1)

#multivariate regression
fm1<-lm(y~x[,1]+x[,2]+x[,3])
summary(fm1)

#2-stage regression
fm2_1<-lm(y~x[,1]+x[,2])
res1<-residuals(fm2_1)

fm2_2<-lm(x[,3]~x[,1]+x[,2])
res2<-residuals(fm2_2)

fm2<-lm(res1~res2)
summary(fm2)

#The estimated parameter for x[,3] in fm1 and the estimated parameter for res2 in fm2 are equal which show thw equivalence of one stage multivariate regression and 2-stage regression.


#####2. multivariate regression to the effects obtained from successive orthogonalization
a<-rep(1, times=10)
Z<-
Z<-cbind(a,0,0,0)
Z[,2]<-x[,1]-Z[,1]*(x[,1]%*%Z[,1])/(Z[,1]%*%Z[,1])
Z[,3]<-x[,2]-Z[,1]*(x[,2]%*%Z[,1])/(Z[,1]%*%Z[,1])-Z[,2]*(x[,2]%*%Z[,2])/(Z[,2]%*%Z[,2])
Z[,4]<-x[,3]-Z[,1]*(x[,3]%*%Z[,1])/(Z[,1]%*%Z[,1])-Z[,2]*(x[,3]%*%Z[,2])/(Z[,2]%*%Z[,2])-Z[,3]*(x[,3]%*%Z[,3])/(Z[,3]%*%Z[,3])

fm3<-lm(y~Z[,4])
summary(fm3)
#The estimated parameter for x[,3] in fm1 and the estimated parameter for Z[,4] in fm3 are equal which show thw equivalence of multivariate regression to the effects obtained from successive orthogonalization.

