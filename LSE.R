rm(list=ls())
lsefun<-function(z)
{
    ans<-rep(0,3)
    n<-length(z)
    z1<-sort(z)
    mu_lse<-z1[1]
    z1<-z1-mu_lse
    Fx=X=Y<-rep(0,n-1)
    z1<-z1[2:n]
    for (i in 1:n-1)    {
        Fx[i]=(i-0.3)/(n+0.4)
        Y[i]=log(-log(1-Fx[i]))
        X[i]=log(z1[i])
    }
    linearMod<-lm(Y~X)
    parameters<-linearMod$coefficients
    #print(parameters)
    lamb<-exp(parameters[1])
    alp<-parameters[2]
    ans<-c(alp, lamb, mu_lse)
    return (ans)
}

k=100
values<-matrix(0,k,3)
for (i in 1:k)  {
    source("generateweibull.R")
    ans<-lsefun(z)
    values[i,1]=ans[1]
    values[i,2]=ans[2]
    values[i,3]=ans[3]
}
for (i in i:k-1)  {
    if (values[i,2]<0.5)   {
        values[i,]=values[i+1,]
    }
}
if(values[k,2]<0.5) {
    values[k,]=values[1,]
}
mean_values<-rep(0,3)
for (i in 1:3)  {
    mean_values[i]=mean(values[,i])
}
true_alp=rep(alpha,k)
true_lam=rep(lambda,k)
true_mu=rep(mu,k)

MSE<-rep(0,3)
MSE[1]=(1/k)*sum((true_alp-mean_values[1])^2)
MSE[2]=(1/k)*sum((true_lam-mean_values[2])^2)
MSE[3]=(1/k)*sum((true_mu-mean_values[3])^2)
print(mean_values)
print(MSE)