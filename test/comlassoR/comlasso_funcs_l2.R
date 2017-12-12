#library(Rglpk)
# For an intermediate lambda, 
# find new mu and maximum of step size corresponding to new category.
get.new.mu.step <- function(corr.B, d_corr.B, K.B, w.list.B, lam, d_lam)
{    
  tol <- 1e-8 
  feasibility <- TRUE
  num.B.categ <- length(corr.B)
  step.vec <- rep(Inf, num.B.categ)
  #categ <- matrix(0, num.B.categ, 2)
  categ1 <- vector(mode="list", num.B.categ)
  categ2 <- categ1
  
  cand.mu <- rep(0, num.B.categ)
  obj <- c(0,0,1)
  for(j in 1:num.B.categ){
    dir <- rep("<=", 2*K.B[j])
    rhs <- c(-corr.B[[j]]+lam*w.list.B[[j]], corr.B[[j]]+lam*w.list.B[[j]])
    Cmat <- matrix(0, 2*K.B[j], 3)
    Cmat[,1] <- rep(c(-1,1), each=K.B[j])
    Cmat[,2] <- -Cmat[,1]
    Cmat[,3] <- c(-d_corr.B[[j]]-d_lam*w.list.B[[j]],d_corr.B[[j]]-d_lam*w.list.B[[j]])
    
    insol <- Rglpk::Rglpk_solve_LP(obj,Cmat,dir,rhs,max=TRUE)
    step.vec[j] <- insol$solution[3]
    cand.mu[j] <- insol$solution[1]-insol$solution[2]

    tmp <- -corr.B[[j]]+cand.mu[j]

    j1 <- which(abs(tmp+lam*w.list.B[[j]]+ step.vec[j]*(d_corr.B[[j]]+d_lam*w.list.B[[j]]))<tol)
    j2 <- which(abs(tmp-lam*w.list.B[[j]]+ step.vec[j]*(d_corr.B[[j]]-d_lam*w.list.B[[j]]))<tol)
#cat("j1",j1,"\n")
#cat("j2",j2,"\n")
    if(length(j1)==length(j2)){
      if(sum(j1 == j2)==length(j1)){
        #j1 <- j1[1]; j2 <- j2[length(j1)]
        return(list(feasibility=FALSE))
      }
    }  
    if(length(j1)==0)
      j1 <- 0
    #categ[j,] <- c(j1, j2)
    categ1[[j]] <- j1
    categ2[[j]] <- j2
#cat("step.vec[j]", step.vec[j], cand.mu[j], "\n")
    if(step.vec[j] <= 1e-10 & j1==0)
      step.vec[j] <- Inf
  }
  return(list(step.vec=step.vec, categ1=categ1, categ2=categ2, cand.mu=cand.mu, feasibility=feasibility))
}
# For sufficiently large lambda, find lambda and mu.
get.init.solution_l2 <- function(corr, w.list, K, tol=1e-9)
{
  diff <- c(); cand.mu <- c()
  p <- length(corr)
  #A.cand <- matrix(0, p, 2)
  A.cand1 <- vector(mode="list",p)
  A.cand2 <- vector(mode="list",p)
  obj <- c(0,0,1)
  
  for(j in 1:p){
    dir <- rep("<=", 2*K[j])
    rhs <- c(corr[[j]], -corr[[j]])
    Cmat <- matrix(0, 2*K[j], 3)
    Cmat[1:K[j],1] <- 1; Cmat[(K[j]+1):(2*K[j]),1] <- -1;
    Cmat[,2] <- -Cmat[,1]
    Cmat[,3] <- c(-w.list[[j]],-w.list[[j]]) # rep(-w.list[[j]],times=2)
  
    insol <- Rglpk::Rglpk_solve_LP(obj,Cmat,dir,rhs,max=FALSE)
    diff[j] <- insol$solution[3]
    #cat("max lambda=", diff[j], "\n")
    cand.mu[j] <- insol$solution[1]-insol$solution[2]
    kkt1 <- -corr[[j]] + cand.mu[j] + diff[j]*w.list[[j]]
    kkt2 <- -corr[[j]] + cand.mu[j] - diff[j]*w.list[[j]]
#    print(min(abs(kkt1)))
#    print(min(abs(kkt2)))
    
    j1 <- which(abs(-corr[[j]] + cand.mu[j] + diff[j]*w.list[[j]])<tol)
    j2 <- which(abs(-corr[[j]] + cand.mu[j] - diff[j]*w.list[[j]])<tol)
#print(which(abs(-corr[[j]] + cand.mu[j] + diff[j]*w.list[[j]])<tol*1))
    A.cand1[[j]] <- j1
    A.cand2[[j]] <- j2
  }
  # Set lambda and mu accordingly
  act.categ <- which.max(diff)
  return(list(act.categ=act.categ, lambda=diff[act.categ], 
    mu=cand.mu[act.categ], A.cand1=A.cand1, A.cand2=A.cand2))
}
# sugar functions
Res.reg_l2 <- function(y,fx){
  drop(y-fx)
}
d_getgcorr_l2 <- function(tX, d_fx){
  return(drop(tX%*%d_fx))
}
getgcorr.reg_l2 <- function(y, tX, res){
  return(drop(tX%*%res))
}
# objetive value
getobjective.res_l2 <- function(res){
  1/2*sum(res^2)
}


### Function ends
obcoeff <- function(b1, k, p){
  if(k == 1) 
    b1 = c(-sum(b1),b1)
  if(k == p) 
    b1 = c(b1,-sum(b1))
  if(k>1 & k<p) 
    b1 = c(b1[1:(k-1)],-sum(b1), b1[(k):(p-1),])
  b1 
}
obmcoeff <- function(bm, k){
  b <- matrix(0, nrow(bm)+1, ncol(bm))
  for(j in 1:ncol(bm)){
    b[,j] <- obcoeff(bm[,j],k)  
  }
  return(b)
}  


