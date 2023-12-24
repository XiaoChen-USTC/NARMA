library(MASS)
library("foreach")
library("doParallel")
N=100
T=100
s_num=300
p=3
ran_num=100
para=c(0.3,-0.5,0.3,0.8,0.5,0.3,0.5,0.2)
N1=matrix(0,ncol=p+5,nrow=s_num)

get_Z <- function(N,p){     #generate Z
 set.seed(N*p)
 Sigma=matrix(0, nrow=p, ncol=p)
 for(i in 1:p)
  for(j in 1:p)
   Sigma[i,j] = 0.5^(abs(i-j))
 return(mvrnorm(N,rep(0,p),Sigma))
}
Z=get_Z(N,p)

get_A <- function(N){        #Dyad Independence Model
 A=matrix(0, nrow=N, ncol=N)
 m=matrix(0, nrow=N, ncol=N)
 for(i in 2:N)
  for(j in 1:(i-1)){
   m[i,j]=runif(1,0,1)
   if(m[i,j]<(20/N))
    A[i,j]=A[j,i]=1 
   else if(m[i,j]<(20/N+0.5*N^(-0.8)))
    A[i,j]=1
    else if(m[i,j]<(20/N+N^(-0.8)))
     A[j,i]=1
  }
 for(i in 1:N)
  if(sum(A[i,]==0)) A[i,(i+1)%%N]=1
 return(A)
}
A=get_A(N)
TNOE=0
for(i in 1:N)
  for(j in 1:N)
    if(A[i,j]==1) TNOE=TNOE+1
TNOE
W=A/rowSums(A)


epsilon=matrix(0, nrow=N, ncol=T)
Y=matrix(0, nrow=N, ncol=T)
Y[,1]=rnorm(N)
epsilon[,1]=rnorm(N)
for(i in 2:T){      
 set.seed(1000+200*i)
 epsilon[,i]=rnorm(N)
 for(j in 1:N)
  Y[j,i]= para[1] + Z[j,]%*% as.matrix(para[2:(p+1)])+ para[p+2]*Y[j,i-1]+ para[p+3]*(W[j,]%*%Y[,i-1]) + para[p+4]*epsilon[j,i-1] + para[p+5]*(W[j,]%*%epsilon[,i-1])+ epsilon[j,i]
}

no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)
N2 <- foreach(jj=1:s_num, .combine="rbind") %dopar% {
  set.seed(1000+200*jj)
  Q_rec <- function(para1,Y){
    e_rec=matrix(0, nrow=N, ncol=T)
    e_rec[,1]=rnorm(N)
    for(i in 2:T){
      M=cbind(rep(1,N),Z,as.vector(Y[,i-1]),as.vector(W%*%Y[,i-1]),as.vector(e_rec[,i-1]),as.vector(W%*%e_rec[,i-1]))
      a=M%*%para1
      e_rec[,i]= Y[,i]- a
    }
    return(sum(e_rec^2)/(N*(T-1)))
  }
  M=matrix(0,nrow=ran_num,ncol=(p+6))
  for(i in 1:ran_num){     #optimization
    fit=optim(par=runif(p+5),f=function(para) Q_rec(para,Y))
    M[i,]=c(fit$par,fit$value)
  }
  N1[jj,]=M[M[,p+6]==min(M[,p+6]),][1:(p+5)]
}
a=colMeans(N2)
a
sd1=NULL
for(i in 1: (p+5))
  sd1[i]=sd(N2[,i])
sd1
stopCluster(cl)