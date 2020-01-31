options(warn=1)
library(cubature)
rm(list=ls())

source('generateweibull.R')
z<-sort(z)
a<-100.
b<-100.
c<-100.
d<-100.
I1<-function(alpha) {
    num=(alpha^(c+N))*(exp(-d*alpha))*(prod(z^(alpha-1)))
    den=((b+sum(z^alpha))^(a+N))*(1e-300)
    return (num/den)
}

I2<-function(alpha) {
    num=(alpha^(c+N+1))*(exp(-d*alpha))*(prod(z^(alpha-1)))
    den=((b+sum(z^alpha))^(a+N))*(1e-300)
    return (num/den)
}

K<-function(alpha) {
    num=(alpha^(c+N-1))*(exp(-d*alpha))*(prod(z^(alpha-1)))
    den=((b+sum(z^alpha))^(a+N))*(1e-300)
    return (num/den)
}

lower<- 0
upper<- 5
gk<-adaptIntegrate(K, lower, upper, tol=1e-3)
i1<-adaptIntegrate(I1, lower, upper, tol=1e-3)
i2<-adaptIntegrate(I2, lower, upper, tol=1e-3)

print(c(gk, i1, i2))

p<- -1
q<- -1
i<- 0
while ((p<0) || (q<0))  {
    gk<-adaptIntegrate(K, lower, upper, tol=1e-3)
    i1<-adaptIntegrate(I1, lower, upper, tol=1e-3)
    i2<-adaptIntegrate(I2, lower, upper, tol=1e-3)
    gk<-as.numeric(gk[1])
    i1<-as.numeric(i1[1])
    i2<-as.numeric(i2[1])
    print(c(gk, i1, i2))
    print('************')
    p<-(i2-i1)/i1
    q<-(p*gk)/i1
    i=i+1
    source('generateweibull.R')
    z<-sort(z)
    if(i>1000)    {
        break
    }
} 

print(c(i, p, q))
print(p/q)
val1<-rgamma(1, p, rate=q)
val2<-rgamma(1, a+N, rate=b+sum(z^val1))
print(c(val1, val2))


# ptm<-proc.time()
# k<-100
# al<-rep(0,k)
# lam<-rep(0,k)
# for (i in 1:k)  {
#     p<- -1
#     q<- -1
#     i<- 0
#     while ((p<0) || (q<0))  {
#         gk<-adaptIntegrate(K, lower, upper, tol=1e-3)
#         i1<-adaptIntegrate(I1, lower, upper, tol=1e-3)
#         i2<-adaptIntegrate(I2, lower, upper, tol=1e-3)
#         gk<-as.numeric(gk[1])
#         i1<-as.numeric(i1[1])
#         i2<-as.numeric(i2[1])
        
#         p<-(i2-i1)/i1
#         q<-(p*gk)/i1
#         if(is.nan(p))   {
#             p<- -0.1600051
#         }
#         if(is.nan(q))   {
#             q<- -0.141291
#         }
#         i=i+1
#         source('generateweibull.R')
#         z<-sort(z)
#         if(i>1000)    {
#             break
#         }
#     } 
#     al[i]<-rgamma(1,p, rate=q)
#     if (is.nan(al[i]))  {
#         al[i]=.07428718
#     }
#     lam[i]<-rgamma(1, a+N, rate=b+sum(z^al[i]))
# }
# ALPHA<-mean(al)
# LAMBDA<-mean(lam)
# print(ALPHA)
# print(LAMBDA)
# print(proc.time()-ptm)