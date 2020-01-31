options(warn=1)
library(cubature)
rm(list=ls())

posterior_upper<-function(x)
{
  f=1
  z=sort(z)
  for(i in 1:N-1)
  {
    f=f*x[1]*x[2]*((z[i]-x[3])^(x[1]-1))*exp((-x[2])*((z[i]-x[3])^x[1]))
  }
  num=(f)*(x[1]^(c-1))*exp((-d)*x[1])*(x[2]^(a-1))*exp((-b)*x[2])
  # deno=gam2ma(a)*gamma(c)*z1[1]
  # p=num/(deno)
  return (num)
}



hfun<-function(alphaf, lambdaf, muf, param, z1)
{
  lower<-c(0,0,0)
  upper<-c(3,3,z[1])
  g<-adaptIntegrate(posterior_upper, lower, upper, tol=1e-3)
  alpha_d <- d-sum(log(z1-muf))
  lambda_b <- b+sum((z1-muf)^alphaf)
  h1<-log(param)+log(gamma(a+N))+log(gamma(c+N))+log(z[1])
  P=as.numeric(g[1])
  #print(h1)
  h2<-log(P)+(c+N)*log(alpha_d)+(a+N)*log(lambda_b)+sum(log(z1-muf))
  #print(h2)
   
  return (exp(h1-h2))
}
#test<-hfun(a,b,z[1]/2)

bayes<-function(z1)
{
  #z1<-sort(z)
  m=100
  lambda_f=mu_f=alpha_f=rep(0,m)
  h_value=matrix(0,m,3)
  for(i in 1:m)
  {
    mu_f[i] = runif(1,0,z[1])
    z1<-z[2:N]
    #N<-N-1
    alpha_a=c+N
    alpha_b=d-sum(log(z1-mu_f[i]))
    alpha_f[i] = rgamma(1,alpha_a, rate=alpha_b)
    lambda_a=a+N
    lambda_b=b+sum((z1-mu_f[i])^alpha_f[i])
    lambda_f[i] = rgamma(1,lambda_a, rate=lambda_b)
    #print(hfun(alpha_f[i], lambda_f[i], mu_f[i]))
    h_value[i,1]=hfun(alpha_f[i], lambda_f[i], mu_f[i], alpha_f[i], z1)
    h_value[i,2]=hfun(alpha_f[i], lambda_f[i], mu_f[i], lambda_f[i], z1)
    h_value[i,3]=hfun(alpha_f[i], lambda_f[i], mu_f[i], mu_f[i], z1)
  }
  ans<-c(mean(h_value[,1]), mean(h_value[,2]), mean(h_value[,3]))
  #print(h_value)
  # means<-rep(0,3)
  # for (i in 1:3)  {
  #   means[i]=mean(h_value[,i])
  # }
  #print(means)
  return (ans)
  
}
ptm<-proc.time()
# source("generateweibull.R")
# z<-sort(z)
# a<-10.
# b<-10.
# c<-10.
# d<-10.

# print(bayes(z,N))
# print(proc.time()-ptm)
k<-40

values<-matrix(0,k,3)
for (i in 1:k)  {
  source("generateweibull.R")
  n<-length(z)
  z<-sort(z)
  a<-1.
  b<-1.
  c<-1.
  d<-1.
  ans<-bayes(z)
  values[i,1]=ans[1]
  values[i,2]=ans[2]
  values[i,3]=ans[3]
}
mean_values<-rep(0,3)
for (i in 1:3)  {
    mean_values[i]=mean(values[,i])
}
truev_alpha<-rep(alpha,k)
truev_lambda<-rep(lambda,k)
truev_mu<-rep(mu,k)
mse<-rep(0,3)
mse[1]<-(1/k)*sum((truev_alpha-mean_values[1])^2)
mse[2]<-(1/k)*sum((truev_lambda-mean_values[2])^2)
mse[3]<-(1/k)*sum((truev_mu-mean_values[3])^2)
print(mse)
print(mean_values)
print(proc.time()-ptm)