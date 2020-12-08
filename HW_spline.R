library(splines2)
#Simulate x and y
set.seed(1)
X<-as.matrix(rnorm(100,mean=0,sd=1))
Y<-exp(X)/(1+exp(X))+rnorm(100,mean=0,sd=1)
plot(X,Y)
X_cubic<-matrix(nrow=100,ncol=7)
X_cubic[,1]<-1
X_cubic[,2]<-X
X_cubic[,3]<-X^2
X_cubic[,4]<-X^3
X_cubic[,5]<-apply(X,1,function(x)ifelse((x>-1),1,0))*(X+1)^3
X_cubic[,6]<-apply(X,1,function(x)ifelse((x>0),1,0))*(X)^3
X_cubic[,7]<-apply(X,1,function(x)ifelse((x>1),1,0))*(X-1)^3

###Cubic splines
fm1<-lm(Y~X_cubic+0)
fm1$fitted.values
points(X,fm1$fitted.values,col="red")

###natural cubic spline
X_ncs<-matrix(nrow=100,ncol=5)
X_ncs[,1]<-1
X_ncs[,2]<-X

X_ncs[,3]<-(X_cubic[,5]-X_cubic[,7])/(2)-(X_cubic[,6]-X_cubic[,7])/(1)
X_ncs[,4]<-(X_cubic[,6]-X_cubic[,7])/(1)-(X_cubic[,6]-X_cubic[,7])/(1)
X_ncs[,5]<--(X_cubic[,6]-X_cubic[,7])/(1)

fm2<-lm(Y~X_ncs+0)
points(X,fm2$fitted.values,col="green")


###smoothing
fm3<-smooth.spline(X,Y)
fm3$lambda
points(fm3$x,fm3$y,col="yellow")

###Bspline
knots<-c(-1,0,1)
X_bsp<-bSpline(X,knots=knots)
fm4<-lm(Y~X_bsp+0)
points(X,fm4$fitted.values,col="blue")

legend("topleft", legend=c("cubic spline", "natural cubic spline","smoothing spline","B-spline"),
       col=c("red","green","yellow","blue"), pch=c(1,1,1,1), cex=0.8)
title(main="Spline")


