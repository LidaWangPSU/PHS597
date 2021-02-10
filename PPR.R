library(MASS)
library(splines)

set.seed(123)
n<-2
N<-30

sigma<-diag(rep(1,n))
mu<-rep(0,n)

X<-mvrnorm(N,mu,sigma)
Y<-rnorm(N,mean=0,sd=10)


ppr_1<-function(X,Y){
d<-1
n<-ncol(X)
N<-nrow(X)
w<-rep(1/n,n)
w<-w/sqrt(sum(w^2))

while(d>=10^(-5)){
xw<-X%*%w
xw_s<-sort(xw)
X_cubic<-matrix(nrow=N,ncol=N+4)
X_cubic[,1]<-1
X_cubic[,2]<-xw
X_cubic[,3]<-xw^2
X_cubic[,4]<-xw^3
for(i in 1:N){
X_cubic[,4+i]<-apply(xw,1,function(x)ifelse((x>xw_s[i]),1,0))*(xw-xw_s[i])^3
}

X_ncs<-matrix(nrow=N,ncol=N)
X_ncs[,1]<-1
X_ncs[,2]<-xw
for(i in 1:(N-2)){
X_ncs[,i+2]<-(X_cubic[,i+4]-X_cubic[,N+4])/(xw_s[i]-xw_s[N])-(X_cubic[,N+3]-X_cubic[,N+4])/(xw_s[N-1]-xw_s[N])
}


#X_ncs
omega<- function(x) {
  N<-length(x)
  Omega <- matrix(rep(0, N*N), ncol = N)
  knot<-sort(x)
  Omega.sub<-matrix(rep(0, (N-2)*(N-2)), ncol = N-2)
 
  for (k in 1:(N-2)) {
    for (j in k:(N-2)) {
      K<-length(knot)
      Km1<-K-1
      C_k1 <- 1/(knot[k]-knot[K])
      C_k2 <- 1/(knot[Km1]-knot[K])
      C_j1 <- 1/(knot[j]-knot[K])
      C_j2 <- 1/(knot[Km1]-knot[K])
      integrand_p1<-function(x) {36*C_k1*(x-knot[k])*C_j1*(x-knot[j])}
      v1 <- integrate(integrand_p1, lower = knot[j], upper=knot[Km1])$value
      integrand_p2<- function(x) {(6*C_k1*(x-knot[k])- 6*C_k2*(x-knot[Km1]))*(6*C_j1*(x-knot[j])- 6*C_j2*(x-knot[Km1]))}
      v2<-integrate(integrand_p2, lower = knot[Km1], upper=knot[K])$value
      Omega.sub[k, j] <- v1+v2
      Omega.sub[j, k] <- Omega.sub[k, j]
    }
  }
  Omega[3:N,3:N]<-Omega.sub
  return(Omega)
}

lambda<-10
Omega<-omega(xw)
beta<-solve(t(X_ncs)%*%X_ncs+lambda*Omega)%*%t(X_ncs) %*%Y

y_hat<-X_ncs%*%beta


###calculate g'(wx)
de_mat<-matrix(nrow=N,ncol=N-2)
#rank(xw)
for(i in 1:N){
  for(j in 1:(N-2)){
    if(j<=rank(xw)[i]){
    de_mat[i,j]<-3*(xw[i]-xw_s[j])^2/(xw_s[j]-xw_s[N])
    }else{
    de_mat[i,j]<-0
    }
  }
}



de<-matrix(nrow=N,ncol=1)
de<-matrix(beta[2],nrow=N,ncol=1)+de_mat%*%beta[3:N,1]

s<-xw+(Y-y_hat)/de

W<-diag(c(de))

w_new<-ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%s

#w_new<-c()
#w_new[1]<-(B*C1-B*C2)/(A1*A2)
#w_new[1]<-(B*C1-B*C2)/(A1*A2)
#w_new[2]<-(C2-B*w_new[1])/A2
d<-sqrt(sum((w-w_new)^2))
w<-w_new
w<-w/sqrt(sum(w^2))
}
print(w)
}
ppr_1(X,Y)

f1<-ppr(X,Y,nterms=2)
f1$alpha[,1]


