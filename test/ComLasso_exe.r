rm(list = ls()) ; gc()
# Rtools is required and Path for rtools should be set.
#library(Rcpp)
#library(inline)
#library(RcppArmadillo)
#setwd("D:/Jeon/rcode/ComLasso/test")
# sourceCpp('inner.cpp')

############### comlasso
### Input
# n = number of observation
# p = number of category
# K = for each category, vector of the numbers of dummies.
#     for example, 3 category, c(4, 5, 10) for each category.
# s0 = when s0 == sum(abs(beta)), stop

### Output
# s = solution beta's 1-norm
# beta0 = intercept
# weights = vector: adaptive weights
### Function start
#    comlasso_l2(n=n,p=8,K=c(10,6,4,3,7,2,8,p-40),X.raw=z,y=y, 
#                    weights=NULL,max.steps=8*min(n,40)+1,lam.min=0,tol=1e-08,trace=FALSE)
#### input
rm( list = ls())
n = 30
K = 3*1
pk = rep(c(3,3,4),1*1)
p = sum(pk)
set.seed(1)
X = matrix(rnorm(n*p), n, p)
y = rnorm(n)
weights=NULL
max.steps = 30
lam_min=0
tol=1e-08
trace=FALSE

# Setting index function #######################################################
# search
  # start 
  idx_gs <- cumsum(pk)-pk+1
  # end
  idx_ge <- cumsum(pk)
  # index of x_{kj} ---> single index
  k = 2; j = 1
  idx_gs[k] + (j-1)
  
# inverse search by matrix
  # group
  dict_idx_k <- rep(1:K, pk)
  # covariates  
  dict_idx_j <- unlist(lapply(pk ,FUN = seq, from = 1, by = 1))
  # single index ---> group index
  dict_idx_k[4]
  # single index ---> within group index
  dict_idx_j[4]
  
# active group  
  i_g_A <-  0L
  set_g_A <- rep(0,K)
  num_g_A <-  rep(0,K)
  

################################################################################    
# sugar function arg
  Res_l2 <- function(y,fx){
    drop(y-fx)
  }
  getobjective_l2 <- function(res){
    1/2*sum(res^2)
  }
  getgcorr_l2 <- function(tX, res){
    return(drop(tX%*%res))
  }
################################################################################
# comlasso  
  s0 <- 1e8
  y <- as.matrix(y,n,1)
  X_sum <- colSums(X)
  OUT = FALSE # TRUE if one of the KKT conditions is violated
  
  ### Variable definition
  mu <- rep(NA, K)
  lambda = 0
  beta_vec <- rep(0,p)
  beta_sign_vec<- rep(0,p)
  beta_vec_A<- rep(F,p)
  

  beta0 <- mean(y)
  res <- y-beta0
  
  # find the initial lambda and mu
  grad_vec <- -drop(t(X)%*%res)
  
  v <- -Inf
  for(istar in 1:length(idx_gs))
  {
    sidx <- idx_gs[istar]:idx_ge[istar]
    for (i in sidx)
    {
      for (j in sidx)
      {
        if (i == j ) next
        v0 <- -grad_vec[i]+grad_vec[j]
        if (v0 > v) 
        {
          v <- v0
          j1 <- i
          j2 <- j
        }
      }
    }
  }

  # j1 : + sign, # j2 : - sign
  beta_sign_vec[j1] <- 1
  beta_sign_vec[j2] <- -1
  beta_vec_A[j1] <- TRUE
  beta_vec_A[j2] <- TRUE
  # activated group index
  
  k_A <- dict_idx_k[j1]
  
  # Management of group index
    # intially set i_g_A is 1
    # if a new activated group is added then i_g_A <- i_g_A + 1;
    #                                        set_g_A[i_g_A] <- k_A;      
    # if group is deleted, then  set_g_A[i_g_A] <- 0;
    #                                         i_g_A <- i_g_A - 1;
    # active group can be returned by    set_g_A[1:i_g_A]
  i_g_A <-  i_g_A + 1L
  set_g_A[i_g_A] <- k_A
  
  # lambda
  lambda <- v/2
  # mu : mu[k_A]<- lambda  - corr_vec[j2] or the following
  mu[k_A] <- -lambda  - grad_vec[j1]
  # activated group index vector (bool) : 
      # redundancy of beta_gvec_A and num_g_A.
      # beta_gvec_A[k_A] <- 1
  # update num_g_A: 
  num_g_A[k_A] <- 2
  act_group <- sort(set_g_A[1:i_g_A])
  
# START LOOP !!  
  # Finding direction
  a <- sum(beta_vec_A)
  q <- i_g_A

  pA_vec <- c(1, a, 1, q)
  pA <- sum(pA_vec)
  XA <- X[,beta_vec_A,drop = F]
  XA1 <- X_sum[beta_vec_A]
  tXXA <- t(XA)%*%XA
  Amat <- matrix(0, a ,q)
  
  num_g_A0 <- num_g_A[num_g_A>0]
  k = 0
  for (i in 1:q)
  {
    Amat[(1:num_g_A0[i]) + k, i] <- 1
    k <- k + num_g_A0[i]
  }
  Bmat <- t(Amat)
  sign_vec <- beta_sign_vec[beta_vec_A]
  
  # 
  Dmat <- matrix(0, pA, pA)
  bvec <- rep(0, pA)
  bvec[pA] <- 1
  
  # bar
  Dmat[1,1] <- n
  Dmat[1,(1:a)+1] <- XA1 
  Dmat[(1:a)+1,1] <- XA1 
  Dmat[(1:a)+1, (1:a)+1] <- tXXA
  Dmat[(1:a + 1), (1:q + 2 + a)] <- Amat
  Dmat[(1:q + 1 + a), (1:a)+ 1] <- Bmat
  Dmat[1:a+1, 1+1+a] <- sign_vec
  Dmat[1+1+a+q , 1:a+1] <- sign_vec
  
  rderiv <- solve(Dmat,bvec)
  # delta
  delta <- rep(Inf,4)
  
  ### 2-1. Derive when active dummy variables become inactive
  step_vec <- -beta_vec[beta_vec_A]/rderiv[(1:a)+1]
  step_vec[step_vec < tol] = Inf
  delta[1] <- min(step_vec)

  ### 2-2-a. Compute the distance needed to activate a new categorical variable
  dy <- rderiv[1] + XA%*%rderiv[(1:a)+1]
  corr_vec <- drop(t(X)%*%dy)
  act_type <- "none"
  #(i,j) = (j,j')
  v <- Inf
  for(istar in 1:K)
  { 
    if ( istar %in% act_group) next
    sidx <- idx_gs[istar]:idx_ge[istar]
    for (i in sidx)
    {
      for (j in sidx)
      {
        if (i == j ) next
        cjj <- -corr_vec[i]+ corr_vec[j] - 2*rderiv[1+a+1]
        if (cjj <= 0) next
        djj <- grad_vec[i] - grad_vec[j]
        v0 <- (2*lambda+djj)/cjj
        if (v0<=0) next
        
        if (v0 < v) 
        {
          v <- v0
          j1 <- i
          j2 <- j
        }
      }
    }
  }
  delta[2]<- v
  # determine the type of violation : sign of active set
  j_vec <- c(j1, j2)
  j1_v <- grad_vec[j1] + corr_vec[j1]*v
  j2_v <-grad_vec[j2] + corr_vec[j2]*v
  j_order <-order(j_vec, decreasing = T)
  # jstar1: positive activation ; jstar2: negative activation
  jstar1 <- j_vec[j_order[1]]
  jstar2 <- j_vec[j_order[2]]
  mu_tmp <- lambda + rderiv[1+a+1]*v -grad_vec[jstar1]-corr_vec[jstar1]*v 
  # mu_tmp should be equal to the following value
  # mu_tmp <- -lambda -rderiv[1+a+1]*v -grad_vec[jstar2] -corr_vec[jstar2]*v 
  
  ### 2-2-b. activated individual variable in the active categorical variable
  v <- Inf
  q_istar = 0
  act_sign <- 0
  tmp_sign <- 0
  j1 <- NA
  for(istar in act_group)
  { 
    q_istar = q_istar + 1
    sidx <- idx_gs[istar]:idx_ge[istar]
    for (j in sidx)
    {
      if (beta_sign_vec[j]!=0) next
      c_j1 <- -rderiv[1+a+1]- corr_vec[j] - rderiv[1+a+1+q_istar]
      d_j1 <- lambda + grad_vec[j] + mu[istar]

      c_j2 <- -rderiv[1+a+1]+ corr_vec[j] + rderiv[1+a+1+q_istar]
      d_j2 <- lambda - grad_vec[j] - mu[istar]
      
      if (c_j1 <=0 & c_j2<=0) next
      if (c_j1 <=0 & c_j2>0) 
      {
        v0 <- d_j2/c_j2
        tmp_sign <- -1
      }
      if (c_j1 >0 & c_j2<=0) 
      {
        v0 <- d_j1/c_j1
        tmp_sign <- 1
      }
      if (c_j1 >0 & c_j2>0) 
      {
        v0 <- min(d_j1/c_j1,d_j2/c_j2)
        tmp_sign <- c(1,-1)[which.min(c(d_j1/c_j1,d_j2/c_j2))]
      }
      
      if(v0<v)
      {
        v <- v0
        act_sign <- tmp_sign
        j1 <- j
      }
    }  
  }
  delta[3]<- v
  jstar3 <- j1
  
  ### 2-3. Compute the distance needed for lambda to become 0
  if(rderiv[1+a+1]>0){
    delta[4] <- 0
    break
  }
  lam_final <- (lam_min-lambda)/rderiv[1+a+1]
  delta[4] = ifelse(lam_final < tol, Inf, lam_final)
  
  
################################################################################
  delta_f <- min(delta)
  # update beta_vec
  beta_vec[beta_vec_A] <- beta_vec[beta_vec_A]+delta_f*rderiv[(1:a)+1] 
  beta0 <- beta0 + delta_f*rderiv[1]
  lambda <- lambda + delta_f*rderiv[1+a+1]
  mu[act_group]   <- mu[act_group] + delta_f*rderiv[(1+a+1+1):pA] 
  
  # check vanishing type 
  if (which.min(delta) == 1) 
  {
    v_idx<- which(beta_vec_A==T)[which.min(step_vec)]
    k_A <- dict_idx_k[v_idx]
    if (num_g_A[k_A] == 2) 
      act_type <- "g_van" else act_type <- "i_van" 
  }
  # check activation type
  if (which.min(delta) == 2) act_type <- "g_act"  
  if (which.min(delta) == 3) act_type <- "i_act"  
  if (which.min(delta) == 4) act_type <- "t_act"

  cat("Type of update:", act_type, "\n")

  if (act_type == "g_van")
  {
    sdix <- idx_gs[k_A]:idx_ge[k_A]
    beta_sign_vec[sidx] <- 0
    beta_vec_A[sidx] <- FALSE
    mu[k_A] <- NA
    set_g_A[i_g_A] <- 0;
    i_g_A <- i_g_A - 1;
    num_g_A[k_A] <- 0
    act_group <- sort(set_g_A[1:i_g_A])
  }
  
  if (act_type == "i_van")
  {
    beta_sign_vec[v_idx] <- 0
    beta_vec_A[v_idx] <- FALSE
    mu[act_group]   <- mu[act_group] + delta_f*rderiv[(1+a+1+1):pA]     
    num_g_A[k_A] <- num_g_A[k_A] -1
  }
  
  if (act_type == "g_act")
  {
    beta_sign_vec[jstar1] <- 1
    beta_sign_vec[jstar2] <- -1
    beta_vec_A[jstar1] <- TRUE
    beta_vec_A[jstar2] <- TRUE
    k_A <- dict_idx_k[jstar1]
    mu[k_A] <- mu_tmp
    i_g_A <- i_g_A + 1
    set_g_A[i_g_A] <- k_A
    num_g_A[k_A] <- 2
    act_group <- sort(set_g_A[1:i_g_A])
  }
  
  if (act_type == "i_act")
  {
    beta_sign_vec[jstar3] <- act_sign
    beta_vec_A[jstar3] <- TRUE
    k_A <- dict_idx_k[jstar3]
    num_g_A[k_A] <- num_g_A[k_A] + 1
  }
  
  
  # intially set i_g_A is 1
  # if a new activated group is added then i_g_A <- i_g_A + 1;
  #                                        set_g_A[i_g_A] <- k_A;      
  # if group is deleted, then  set_g_A[i_g_A] <- 0;
  #                                         i_g_A <- i_g_A - 1;
  # active group can be returned by    set_g_A[1:i_g_A]
  
  