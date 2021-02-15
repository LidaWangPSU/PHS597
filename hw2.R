###Hw2

##Foward
#z=wx+b
#y^hat=sigma(z)=exp(z)/(1+exp(z))
#L=(y-y^hat)^2


##Backward
#x(n,1) w(1,n) b(1,1) z(1,1) y^hat(1,1)

#dL/dw=(dL/dy^hat)*(dy^hat/dw)
#dL/dy^hat=2*(y^hat-y)
#dy^hat/dw=exp(wx+b)/(1+exp(wx+b))^2%*%t(x)
#dL/dw=2*(y^hat-y)*exp(wx+b)/(1+exp(wx+b))^2%*%t(x)

#dL/db=(dL/dy^hat)*(dy^hat/db)
#dL/dy^hat=2*(y^hat-y)
#dy^hat/db=exp(wx+b)/(1+exp(wx+b))^2
#dL/db=2*(y^hat-y)*exp(wx+b)/(1+exp(wx+b))^2
