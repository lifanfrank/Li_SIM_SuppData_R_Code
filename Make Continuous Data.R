#############################################################
# Simulate correlated Gaussian outcomes for CRXO trials

# by Fan Li
# June 2018

# INPUT
# n: Number of clusters
# m: Cluster-period size
# t: Number of periods
# delta: Effect size
# s2: Dispersion parameter or total variance (assume to be 1
#     and suppressed in the input)
# beta: Vector of period effects
# alpha: Vector of correlations 
#        alpha_0=within period correlation
#        alpha_1=inter-period correlation
#############################################################

contGEN<-function(n,m,t,delta,beta,alpha){
  
  require(mvtnorm)
  
  ########################################################
  # Create nested exchangeable correlation matrix.
  ########################################################
  
  nxch<-function(alpha){
    rho<-alpha[1]
    nu<-alpha[2]
    bm1<-(1-rho)*diag(t*m)
    bm2<-(rho-nu)*kronecker(diag(t),matrix(1,m,m))
    bm3<-nu*matrix(1,t*m,t*m)
    return(bm1+bm2+bm3)
  }
  
  ########################################################
  # returns variance matrix of Gaussian variables with 
  # dispersion parameter s2=1 and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,s2){
    return(s2*r)
  }
  
  # Create treatment sequences
  trtSeq<-matrix(c(1,0,
                   0,1),2,2,byrow=TRUE)
  g<-n/2 # number of clusters per group
  
  # Simulate correlated Gaussian outcomes
  s2<-1
  y<-NULL
  r<-nxch(alpha)
  v<-cor2var(r,s2)   # v <- cov matrix
  for(i in 1:2){
    u_c<-c(beta+delta*trtSeq[i,])
    u<-rep(u_c,each=m)
    y<-cbind(y,t(rmvnorm(g,u,v))) # simulate data matrix
  }
  
  # Return simulated data matrix
  return(y)
}

set.seed(062018)
n<-20
m<-65
t<-2
delta<- -0.25
beta<-c(0,-0.2)
alpha<-c(0.1,0.05)   # correlation parameters
# Generate outcome
y<-contGEN(n,m,t,delta,beta,alpha)
y<-c(y)

# marginal mean design matrix including period and treatment indicators
X<-NULL
trtSeq<-matrix(c(1,0,
                 0,1),2,2,byrow=TRUE)
g<-n/2 # number of clusters per group
for(i in 1:2){
  for(j in 1:g){
    X<-rbind(X,kronecker(cbind(diag(t),trtSeq[i,]),rep(1,m)))}
}

# column names of design matrix
colnames(X)<-c("period1","period2","treatment")
cluster<-rep(1:n,each=t*m)         # create cluster id
ind<-rep(rep(1:m,t),n)             # create individual id
period<-rep(rep(1:t,each=m),n)     # create period label

simdata_cont<-data.frame(cbind(y,ind,cluster,period,X))
setwd("D:/Research/CRT Methodology/CRXOSampleSize/Latex/Submission/R Code")
write.csv(simdata_cont, file = "simdata_cont.csv", row.names = FALSE)

