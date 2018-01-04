#l2
library(genlasso)
#library(comlasso)
library(Rglpk)
library("MethylCapSig")
setwd("D:/Jeon/rcode/ComLasso/test")
source('comlasso_funcs_l2.R')
source('comlasso_l2.R')

ns <- c(50)
ps <- c(50)
parset <- matrix(0,length(ns)*length(ps),2)
runtime.list <- vector(mode="list",length=nrow(parset))
cnt <- 0
for(ii in 1:length(ns)){
  for(jj in 1:length(ps)){
    cnt <- cnt+1
    parset[cnt,] <- c(ns[ii],ps[jj])
  }
}


for(ll in 1:nrow(parset)){
  n <- parset[ll,1] ; p <- parset[ll,2]
  sigma <- 0.2^(abs(outer(1:p,1:p,"-")))
  diag(sigma) <- 1
  mm<-rep(0,p);mm[1:5]<-p/2;mm[6:p]<-1
  
  Rnum <- 1
  runtime<-matrix(0,Rnum,2)
  for(r in 1:Rnum){
    set.seed(1)
    w <- mvlognormal(n=n, Mu=mm,Sigma=rep(1,p),R=sigma);
    z <- diag(1/rowSums(w)) %*% w
    
    b <- rep(0,p)
    b[1:16] <- c(1,-0.8,0.4,0,0,-0.6,0,0,0,0,-1.5,0,1.2,0,0,0.3)
    sum(b)
    y <- drop(z %*% b)+rnorm(n,0,0.5^2)
    
    Cm <- matrix(0,8,p)
    Cm[1,1:10]<-1;Cm[2,11:16]<-1;
    Cm[3,17:20]<-1;Cm[4,21:23]<-1;
    Cm[5,24:30]<-1;Cm[6,31:32]<-1;
    Cm[7,33:40]<-1;Cm[8,41:p]<-1;  
    drop(Cm%*%b)
    
    runtime[r,1]<-system.time(gfun <- genlasso(y=y,X=z,D=Cm,approx=FALSE,maxsteps=2000,minlam=0,
      rtol=1e-07,btol=1e-07,eps=1e-4,verbose=FALSE,svd=FALSE)
    )[3]
    
    runtime[r,2] <- system.time(cfun <- comlasso_l2(n=n,p=8,K=c(10,6,4,3,7,2,8,p-40),X.raw=z,y=y, 
      weights=NULL,max.steps=8*min(n,40)+1,lam.min=0,tol=1e-08,trace=FALSE)
    )[3]
  }
  runtime.list[[ll]] <- runtime
}
#source('./library_comlasso/comlasso_l2.R')
#source('./library_comlasso/comlasso_funcs_l2.R')
#cfun <- comlasso_l2(n=n,p=8,K=c(10,6,4,3,7,2,8,p-40),X.raw=z,y=y,
#  weights=NULL,max.steps=length(K)*min(n,sum(K))+1,lam.min=0,tol=1e-08, trace=FALSE)
#cfun <- comlasso_l2(n=n,p=1,K=p,X.raw=z,y=y,
#  weights=NULL,max.steps=length(K)*min(n,sum(K))+1,lam.min=0,tol=1e-08, trace=FALSE)

tot <- c()
for(ll in 1:nrow(parset)){
  tot <- rbind(tot,cbind(parset[ll,1],parset[ll,2],runtime.list[[ll]]))
}

aa <- matrix(unlist(lapply(runtime.list,colMeans)),nrow(parset),2,byrow=T)
library(xtable)
xtable(cbind(parset,aa))



pk = rep(c(3,3,4),1*200)
Cm <- matrix(0,K,p)
for ( i in 1:K)
{
  sidx <- idx_gs[i]:idx_ge[i]
  Cm[i,sidx]<-1
}
system.time(gfun <- genlasso(y=y,X=X,D=Cm,approx=FALSE,maxsteps=2000,minlam=0,
                             rtol=1e-07,btol=1e-07,eps=1e-4,verbose=FALSE,svd=FALSE)
)[3]




