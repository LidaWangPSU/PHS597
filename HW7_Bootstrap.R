library(MASS)

compare_bootstraps<-function(mu,n0,B){
#Simulate X and Y
a1<-matrix(c(1,0.1,-0.1,
               0,1,0.1,
               0,0,1),3,3)
Sigma1<-a1%*%t(a1)
set.seed(5)
X<-mvrnorm(n0,mu,Sigma1)
Y<-as.matrix(rnorm(n0,10,sd=4))
fm<-lm(Y~X+0)
lm<-cbind(fm$coefficients,confint(fm))
colnames(lm)<-c("lse","2.5%","97.5%")

data<-cbind(Y,X)
colnames(data)<-c("y","x1","x2","x3")

###Nonparametric Bootstrap
index<-matrix(ncol=n0,nrow=B)
beta<-matrix(ncol=3,nrow=B)
for(i in 1:B){
index[i,]<-sample.int(n0, n0, replace =TRUE, prob = NULL)
data_b<-data[index[i,],]
beta[i,]<-t(ginv(t(data_b[,2:4])%*%data_b[,2:4])%*%t(data_b[,2:4])%*%data_b[,1])
}

ci<-cbind(quantile(beta[,1],c(0.025,0.975)),quantile(beta[,2],c(0.025,0.975)),quantile(beta[,3],c(0.025,0.975)))
mb<-t(rbind(c(mean(beta[,1]),mean(beta[,2]),mean(beta[,3])),ci))
colnames(mb)<-c("NP-bootstrap","2.5%","97.5%")


###Parametric bootstrap
beta_0<-t(ginv(t(data[,2:4])%*%data[,2:4])%*%t(data[,2:4])%*%data[,1])
beta_1<-matrix(ncol=3,nrow=B)
for(i in 1:B){
  data_b<-data[index[i,],]
  y_s<-data_b[,2:4]%*%t(beta_0)+rnorm(n0,0,1)
  beta_1[i,]<-t(ginv(t(data_b[,2:4])%*%data_b[,2:4])%*%t(data_b[,2:4])%*%y_s)
}
ci1<-cbind(quantile(beta_1[,1],c(0.025,0.975)),quantile(beta_1[,2],c(0.025,0.975)),quantile(beta_1[,3],c(0.025,0.975)))
mb1<-t(rbind(c(mean(beta_1[,1]),mean(beta_1[,2]),mean(beta_1[,3])),ci1))
colnames(mb1)<-c("P-bootstrap","2.5%","97.5%")
list(lm,mb,mb1)

}
mu<-c(3,4,5)
B<-100
compare_bootstraps(mu,5,B)
compare_bootstraps(mu,50,B)
compare_bootstraps(mu,500,B)
compare_bootstraps(mu,5000,B)

#For nonparametric bootstrap and parametric bootstrap, I chose to use mean of beta in B sample sets as the estimate, 
#and the 95% confidence interval were obtained using the 2.5% and 97.5% percentile of the estiamtes.
#Meanwhile, I compared the results with LSE.
#I found that when the original sample size was small, like less than 10. The estimate of non-parametic bootstrap was inaccurate. 
#In my opinion, when the sample size is small, there are only a few values and nonparametric bootstrap may decrease the amount of variation in the simulated data.



