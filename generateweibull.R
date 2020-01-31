N<-25
R=runif(N,0,1)
lambda=1.
alpha=0.5
mu=0.001
z<-(((-1/lambda)*log(1-R))^(1/alpha))+mu
#z<-(((-lambda)*log(1-R))^(1/alpha))+mu
#print(z)
#write.table(z, file='weibdata.csv', na="", sep=",")
