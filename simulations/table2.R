rm(list=ls())
library(Matrix)
library(Rcpp)
library(inline)
library(RcppArmadillo)
library(genlasso)
library(Rglpk); library("MethylCapSig")
sourceCpp('D:/Jeon/rcode/ComLasso/library/inner.cpp')
source("D:/Jeon/rcode/ComLasso/library/ComLassoC.R")
# table 2
para_vec = list()
para_vec[[1]] <- c(50,rep(10,20))
para_vec[[2]] <- c(50,rep(10,50))
para_vec[[3]] <- c(50,rep(10,100))
para_vec[[4]] <- c(100,rep(10,20))
para_vec[[5]] <- c(100,rep(10,50))
para_vec[[6]] <- c(100,rep(10,100))
runtime.list <- vector(mode="list",length=length(para_vec))
Rnum <- 20
ll = 1
for(ll in 1:length(para_vec))
{
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
  
  for(r in 1:Rnum)
  {
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
    #if (p>n) Cm <- Matrix(Cm, sparse = TRUE) 
    
    runtime[r,1] <- system.time(cfun2 <- comLassoC(X,y,pk=pk,lam_min=0,
                                                   tol=1e-08,KKT_check=FALSE) # Prof. Jeon
    )[3]
    
    runtime[r,2]<-system.time(gfun<-genlasso(y=y,X=rX,D=Cm,approx=FALSE,
                                             maxsteps=length(cfun2$lambda_vec)+1,
                                             minlam=0,
                                             rtol=1e-07,btol=1e-07,eps=1e-4,
                                             verbose=FALSE,svd=FALSE))[3]
    cat(runtime[r,],"\n")
  }
  runtime.list[[ll]] <- runtime
}




tot <- c()
for(ll in 1:nrow(parset)){
  tot <- rbind(tot,cbind(parset[ll,1],parset[ll,2],runtime.list[[ll]]))
}

aa <- matrix(unlist(lapply(runtime.list,colMeans)),nrow(parset2),3,byrow=T)
library(xtable)
tmp <- cbind(parset2,aa); colnames(tmp) <- c("n","p","genlasso","com_old","com_new")
xtable(tmp)




# gfun recover
gfun$beta
# figure
i = 18
setwd("C:/Users/Jeon/Documents/fig")
for (i in 1:200)
{
  jpeg(filename = paste0("fit-",i,'.jpg'))
  plot(cfun2$coefMat[,i+1], type = 'b')
  lines(gfun$beta[i+1,][-(1:2)])
  dev.off()  
}
gfun$lambda[-1]
cfun2$lambda_vec

betaMat_genlasso <- matrix(0, ncol(gfun$beta), p + 1)
betaMat_genlasso[,1] <- gfun$beta[1,]
j = 1
for (j in 1:length(pk))
{
  sdix <-idx_gs[j]:idx_ge[j] + 1
  sdix_r <-idx_gs_r[j]:idx_ge_r[j] + 1
  tmp = t(gfun$beta[sdix_r,])
  tmp = cbind(tmp,-rowSums(tmp))
  betaMat_genlasso[,sdix]  = tmp
}


# sensitivity to the threshold
df_mat = matrix(0,5, nrow(betaMat_genlasso))
df_mat[1,]<- c(1,rowSums((abs(cfun2$coefMat)!=0)))
thres_vec = c(1e-3,1e-4,1e-5,1e-6)

for (i in 1:length(thres_vec))
{
  thres = thres_vec[i]  
  for (j in 1:ncol(gfun$beta))
  {
    df_mat[i+1,j] <- sum(abs(gfun$beta[,j])>thres)  
  }
}
df_mat





