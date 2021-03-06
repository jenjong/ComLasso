rm(list=ls())
setwd("~/GitHub/ComLasso/devp")
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(devtools)
library(Rglpk)
library(MethylCapSig)
install_github("glmgen/genlasso")
library(genlasso)
sourceCpp('inner.cpp')
source("ComLassoC.R")
# Set parameters in table 1
# para_vec[[1]] denotes that n = 50, p = 200
para_vec = list()
para_vec[[1]] <- c(50,rep(10,20))
para_vec[[2]] <- c(50,500)
para_vec[[3]] <- c(50,1000)
para_vec[[4]] <- c(100,200)
para_vec[[5]] <- c(100,500)
para_vec[[6]] <- c(100,1000)
runtime.list <- vector(mode="list",length=length(para_vec))
# number of repetitions
Rnum <- 20
ll = 1
# for(ll in 1:length(para_vec))
# {
  n = para_vec[[ll]][1]
  pk <- para_vec[[ll]][-1]
  idx_gs <- cumsum(pk)-pk+1
  idx_ge <- cumsum(pk)   
  idx_gs_r  = idx_gs - 0:(length(pk)-1)
  idx_ge_r  = idx_ge - 1:length(pk)
  
  p = sum(pk)
  B_list = list()
  
  for (j in 1:length(pk))
  {
    sigma <- 0.5^(abs(outer(1:(pk[j]-1),1:(pk[j]-1),"-")))  
    svdFit <- svd(sigma)
    B = diag(svdFit$d)^0.5 %*% t(svdFit$v)
    if ( j == 1) 
    {
      mm<-rep(0,pk[j]-1); mm[1:5]<- log(0.5*pk[j])
    }
    B_list[[j]] = B
  }
  
  runtime<-matrix(0,Rnum,2)
  colnames(runtime) <-c("comlasso", "genlasso") 
  r = 1
  
  #for(r in 1:Rnum)
  #{
    set.seed(r)
    w = NULL
    z = NULL
    for (j in 1:length(pk))
    {
      B <- B_list[[j]]
      tmp = matrix(rnorm(n*(pk[j]-1)), n, pk[j]-1) %*% B
      if (j == 1) tmp = tmp + mm
      tmpU = cbind(exp(tmp), 1)
      z = cbind(z, sweep(tmpU, 1, rowSums(tmpU), FUN = "/") )
      w <- cbind(w, tmp)
      
    }
    b <- rep(0,p-length(pk))    
    b[1:5] <- c(1,-0.8,0.6,0,0)
    b[6:8] <- c(-1.5, -0.5, 1.2)
    y <- drop(w %*% b)+rnorm(n,0,0.5^2)
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
    Cm <- cbind(0,Cm)
    rX <- cbind(1, rX)
    lam_min=0;tol=1e-08;KKT_check=FALSE    
    
    #if (p>n) Cm <- Matrix(Cm, sparse = TRUE) 
#     
#     runtime[r,1] <- system.time(cfun2 <- comLassoC(X,y,pk=pk,lam_min=0,
#                                                    tol=1e-08,KKT_check=FALSE) # Prof. Jeon
#     )[3]
#     
#     runtime[r,2]<-system.time(gfun<-genlasso(y=y,X=rX,D=Cm,approx=FALSE,
#                                              maxsteps=length(cfun2$lambda_vec)+1,
#                                              minlam=0,
#                                              rtol=1e-07,btol=1e-07,eps=1e-4,
#                                              verbose=FALSE,svd=FALSE))[3]
#     cat(runtime[r,],"\n")
#   }
#   runtime.list[[ll]] <- runtime
# }
# 
# lapply(runtime.list, colMeans)
# 
# 


