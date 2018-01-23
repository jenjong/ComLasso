rm(list=ls())
library(Rcpp)
library(inline)
library(RcppArmadillo)
library(genlasso) #library(comlasso)
library(Rglpk); library("MethylCapSig")
sourceCpp('D:/Jeon/rcode/ComLasso/library/inner.cpp')
source("D:/Jeon/rcode/ComLasso/library/ComLassoR.R")
source("D:/Jeon/rcode/ComLasso/library/ComLassoC.R")

ns <- c(50,100,200)
ps <- c(6,11,21,51,101,201)
parset <- matrix(0,length(ns)*length(ps),2)
parset2 <- parset
runtime.list <- vector(mode="list",length=nrow(parset))
cnt <- 0
for(ii in 1:length(ns)){
  for(jj in 1:length(ps)){
    cnt <- cnt+1
    parset[cnt,] <- c(ns[ii],ps[jj])
  }
}
Rnum <- 20
ll = 9
for(ll in 12:nrow(parset)){
  cat(ll,"th learning!!!\n")
  n <- parset[ll,1] ; nk <- parset[ll,2]
  
  pk <- c(5,5,rep(c(3,3,4),nk-2))
  idx_gs <- cumsum(pk)-pk+1
  idx_ge <- cumsum(pk)   
  
  p <- sum(pk); parset2[ll,1] <-n; parset2[ll,2] <-p
  
  sigma <- 0.2^(abs(outer(1:p,1:p,"-")))
  diag(sigma) <- 1
  mm<-rep(0,p);mm[1:5]<-p/2;mm[6:p]<-1
  
  runtime<-matrix(0,Rnum,3); colnames(runtime) <-c("genlasso","comlasso_old","comlasso_new") 
  r = 1
  for(r in 1:Rnum){
    set.seed(r)
    w <- mvlognormal(n=n,Mu=mm,Sigma=rep(1,p),R=sigma);
    z <- diag(1/rowSums(w)) %*% w
    b <- rep(0,p)
    b[1:5] <- c(1,-0.8,0.4,0,-0.6)
    b[6:10] <- c(-1.5,0,1.2,0,0.3)
    #sum(b)
    y <- drop(z %*% b)+rnorm(n,0,0.5^2)
    
    # reparametrization for genlasso
    rX = X = log(z)
    Cm <- diag(1,p)
    for (i in 1:length(pk))
    {
      sidx = idx_gs[i]:idx_ge[i]
      for (j in sidx) rX[,j] = rX[,j] - rX[,idx_ge[i]]
      Cm[idx_ge[i],sidx] <- 1
    }
    rX <- rX[,-(idx_ge)]
    Cm <- Cm[,-(idx_ge)]

    for(qq in 1:(nk)){
      Cm[qq,(idx_gs[qq]:idx_ge[qq])] <- 1
    }

    runtime[r,1]<-system.time(gfun<-genlasso(y=y,X=rX,D=Cm,approx=FALSE,maxsteps=2000,minlam=0,
                                             rtol=1e-07,btol=1e-07,eps=1e-4,verbose=FALSE,svd=FALSE)
    )[3]

    runtime[r,2] <- system.time(cfun1<-comlasso_l2(n=n,p=length(pk),K=pk,X.raw=X,y=y,max.steps=1e5,lam.min=0,tol=1e-08)
    )[3]

    runtime[r,3] <- system.time(cfun2 <- comLassoC(X,y,pk=pk,lam_min=0,tol=1e-08,KKT_check=FALSE) # Prof. Jeon
    )[3]
    cat(runtime[r,],"\n")
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

aa <- matrix(unlist(lapply(runtime.list,colMeans)),nrow(parset2),3,byrow=T)
library(xtable)
tmp <- cbind(parset2,aa); colnames(tmp) <- c("n","p","genlasso","com_old","com_new")
xtable(tmp)





