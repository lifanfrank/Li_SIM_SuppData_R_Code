#############################################################
# Simulate correlated binary outcomes for CRXO trials

# by Fan Li
# June 2018

# Ref: Qaqish, B. F. (2003). A family of multivariate binary distributions 
# for simulating correlated binary variables. Biometrika 90, 455-463.

# INPUT
# n: Number of clusters
# m: Cluster-period size
# t: Number of periods
# delta: Effect size in log odds ratio
# beta: Vector of period effects
# alpha: Vector of correlations 
#        alpha_0=within period correlation
#        alpha_1=inter-period correlation
#############################################################

binGEN<-function(n,m,t,delta,beta,alpha){
  
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
  # a[1:n, 1:n] is the input covariance matrix of Y[1:n].
  # Returns  b[1:n,1:n] such that b[1:t-1, t] are the 
  # slopes for regression of y[t] on y[1:t-1], for t=2:n.
  # Diagonals and lower half of b[,] are copied from a[,].
  # a[,] is assumed +ve definite symmetric, not checked.
  ########################################################
  
  allreg<-function(a){
    n<-nrow(a)
    b<-a
    for(t in 2:n){
      t1<-t-1
      gt<-a[1:t1,1:t1]
      st<-a[1:t1,t]
      bt<-solve(gt,st)
      b[1:t1,t]<-bt
    }
    return(b)
  }
  
  ########################################################
  # returns variance matrix of binary variables with mean
  # vector u[] and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,u){
    s<-diag(sqrt(u*(1-u)))
    return(s%*%r%*%s)
  }
  
  ########################################################
  # r[1:n, 1:n] is the corr mtx
  # u[1:n] is the mean of a binary vector
  # checks that pairwise corrs are in-range for the given u[]
  # only upper half of r[,] is checked.
  # return 0 if ok
  # return 1 if out of range
  ########################################################
  
  chkbinc<-function(r,u){
    n<-length(u)
    s<-sqrt(u*(1-u))
    for(i in 1:(n-1)){
      for(j in (i+1):n){
       uij<-u[i]*u[j]+r[i,j]*s[i]*s[j]
       ok<-((uij <= min(u[i], u[j])) & (uij >= max(0, u[i]+u[j]-1)))
       if(!ok) {return(1)}
      }
    }
    return(0)
  }
  
  ########################################################
  # Multivariate Binary Simulation by Linear Regression.
  # Simulate a single vector.
  # Returns a simulated binary random vector y[1:n] with mean 
  # u[1:n] and regression coefs matrix b[1:n,1:n] (obtained 
  # by calling allreg() above).
  # y[] and u[] are column vectors.
  # Returns -1 if the cond. linear family not reproducible
  ########################################################
  
  mbslr1<-function(b,u){
    n<-nrow(b)
    y<-rep(-1,n)
    y[1]<-rbinom(1,1,u[1])
    for(i in 2:n){
      i1<-i-1
      r<-y[1:i1]-u[1:i1]              # residuals
      ci<-u[i]+sum(r*b[1:i1,i])       # cond.mean
      if(ci < 0 | ci > 1){
        stop(paste("mbslr1: ERROR:",ci))
        return(-1)
      }
      y[i]<-rbinom(1,1,ci)
    }
    return(y)
  }
  
  ########################################################
  # Multivariate Binary Simulation by Linear Regression.
  # Simulate m independent columns (stored in columns).
  # Returns a simulated binary random matrix y[1:n, 1:m].
  # Each column will be a binary random vector y[1:n] with mean 
  # u[1:n] and regression coefs matrix b[1:n,1:n] (obtained 
  # by calling allreg() above).
  # u[] is a column vector.
  ########################################################
  
  mbslrm<-function(b,u,m){
    n<-nrow(b)
    y<-matrix(-1,n,m)
    for(i in 1:m){
      y[,i]<-mbslr1(b,u)
    }
    return(y)
  }
  
  # Create treatment sequences
  trtSeq<-matrix(c(1,0,
                   0,1),2,2,byrow=TRUE)
  g<-n/2 # number of clusters per group
  
  # Simulate correlated binary outcomes
  y<-NULL
  r<-nxch(alpha)
  for(i in 1:2){
    u_c<-c(plogis(beta+delta*trtSeq[i,]))
    u<-rep(u_c,each=m)
    v<-cor2var(r,u)   # v <- cov matrix
    oor<-chkbinc(r,u) # corr out of range?
    if(oor){
      stop("ERROR: Corr out of range for given mean")
    }
    b<-allreg(v)      # prepare coeffs
    y<-cbind(y,mbslrm(b,u,g)) # simulate data matrix
  }
  
  # Return simulated data matrix
  return(y)
}

set.seed(062018)
n<-20
m<-65
t<-2
delta<-log(0.6)
beta<-c(log(0.5/(1-0.5)),log(0.8))
alpha<-c(0.1,0.05)   # correlation parameters
# Generate outcome
y<-binGEN(n,m,t,delta,beta,alpha)
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

simdata_bin<-data.frame(cbind(y,ind,cluster,period,X))
setwd("D:/Research/CRT Methodology/CRXOSampleSize/Latex/Submission/R Code")
write.csv(simdata_bin, file = "simdata_bin.csv", row.names = FALSE)


