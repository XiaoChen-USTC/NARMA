library(MASS)
library("foreach")
library("doParallel")
library(plyr)

p=3
ran_num=100
A = read.csv("邻接矩阵A.csv",header=F)
SZ = read.csv("市值（一年100天）.csv",header=F)
Z = read.csv("Z.csv",header=T)
A=A[-1,-1]
A=as.matrix(A)
SZ=SZ[,2:262]
SZ=as.matrix(SZ)
SZ1=matrix(0,nrow=nrow(SZ),ncol=ncol(SZ)-1)
for(i in 1:nrow(SZ))
  for(j in 1:(ncol(SZ)-1))
  SZ1[i,j]=na.omit(diff(log(SZ[i,])))[j]
Z=Z[,3:5]
Z=as.matrix(Z)
N1=matrix(0,ncol=5,nrow=s_num)

i=1
while(i<=nrow(A)){
  q=0
  for(j in 1:ncol(A))
    if(A[i,j]==0) q=q+1
  if(q==nrow(A)){
    A=A[-i,]
    A=A[,-i]
    SZ=SZ[-i,]
    SZ1=SZ1[-i,]
    Z=Z[-i,]
  }
  else i=i+1
}
T=250

SZ=SZ[-85,]
SZ1=SZ1[-85,]
A=A[-85,]
A=A[,-85]
Z=Z[-85,]
N=nrow(A)
W=A/rowSums(A)
Y=SZ1[,1:T]


Q_rec <- function(para1){
  e_rec=matrix(0, nrow=N, ncol=T)
  e_rec[,1]=0
  for(i in 2:T){
    M=cbind(rep(1,N),as.vector(Y[,i-1]),as.vector(W%*%Y[,i-1]),as.vector(e_rec[,i-1]),as.vector(W%*%e_rec[,i-1]))
    a=M%*%para1
    e_rec[,i]= Y[,i]- a
  }
  return(sum(e_rec^2)/(N*(T-1)))
}
M=matrix(0,nrow=ran_num,ncol=6)
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)
N2 <- foreach(i=1:ran_num, .combine="rbind") %dopar% {
    fit=optim(par=runif(5),Q_rec)
    M[i,]=c(fit$par,fit$value)
}
a=N2[N2[,6]==min(N2[,6]),][1:5]
N2
a
stopCluster(cl)

Y1=matrix(0,nrow=N,ncol=T)
for(i in 1:T) Y1[,i]=SZ1[,i]
for(j in 1:t) Y1=cbind(Y1,0)
e_rec1=matrix(0,nrow=N,ncol=261)
for(i in 2:251){
  M1=cbind(rep(1,N),as.vector(Y1[,i-1]),as.vector(W%*%Y1[,i-1]),as.vector(e_rec1[,i-1]),as.vector(W%*%e_rec1[,i-1]))
  e_rec1[,i]=Y1[,i]-M1%*%a
}
for(i in 1:t){
  M=cbind(rep(1,N),as.vector(Y1[,T+i-1]),as.vector(W%*%Y1[,T+i-1]),as.vector(e_rec1[,i+T-1]),as.vector(W%*%e_rec1[,i+T-1]))
  Y1[,T+i]=M%*%a
}
D=matrix(0,nrow=N,ncol=10)
for(i in 1:t){
  for(j in 1:N){
    if(Y1[j,i+T]*SZ1[j,i+T]>0)
      D[j,i]=1
  }
}
acc=colSums(D)/N
plot(acc,type="b")

sum((Y1[,2]-SZ[,252])^2)



NAR_est <- function(Ymat,W,Z)
{
  Ymat1 = W%*%Ymat                                                                                             ### obtain WY
  Time = ncol(Ymat)-1
  X = cbind(rep(1, nrow(Ymat)*Time),                                                                           ### the intercept
            as.vector(Ymat1[,-ncol(Ymat)]),                                                                    ### WY_{t-1}
            as.vector(Ymat[,-ncol(Ymat)]))                                                                     ### Y_{t-1}
            #do.call("rbind", rep(list(Z), Time)))
  invXX = solve(crossprod(X))                                                                                  ### {t(X)X}^{-1}
  Yvec = as.vector(Ymat[,-1])                                                                                  ### the response vector
  thetaEst = invXX%*%colSums(X*Yvec) 
  return(thetaEst)
}

predict.NAR1<-function(SZ1,t,Z,W,N)
{
  L=250
  Z=as.matrix(Z)
  data1=matrix(0,nrow=N,ncol=L)
  for(i in 1:L) data1[,i]=SZ1[,i]
  for(j in 1:t) data1=cbind(data1,0)
  D=matrix(0,nrow=N,ncol=t) 
  for(k in 1:t){
    theta = NAR_est(data1[,k:(k+L-1)],W,Z)
    M=cbind(rep(1,N),as.vector(W%*%data1[,k+L-1]),as.vector(data1[,k+L-1]))
    data1[,k+L]=M%*%theta
    for(j in 1:N){
      if(data1[j,k+L]*SZ1[j,k+L]>0) 
        D[j,k]=1
      }
  }
  return(D)
}


D=predict.NAR1(SZ1,t=10,Z,W,N)
acc1=colSums(D)/N
plot(acc1,type="b")


