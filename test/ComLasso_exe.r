rm(list = ls()) ; gc()
library(Matrix)
# Rtools is required and Path for rtools should be set.
 library(Rcpp)
 library(inline)
 library(RcppArmadillo)
library(Rglpk); library("MethylCapSig")
 setwd("D:/Jeon/rcode/ComLasso/test")
# sourceCpp('inner.cpp')

############### comlasso
### Input
# n = number of observation
# p = number of category
# K = for each category, vector of the numbers of dummies.
#     for example, 3 category, c(4, 5, 10) for each category.
# s0 = when s0 == sum(abs(beta)), stop

#rm( list = ls())
 # n = 100
 # pk = rep(c(3,3,4),5000)
 # p = sum(pk)
 # set.seed(3)
 # X = matrix(rnorm(n*p), n, p)
 # y = rnorm(n)
 lam_min=0
 tol=1e-08
 max_iter = 1000
 KKT_check = TRUE

KKT_fun<- function(X,y, beta0, beta_vec, mu)
{
  res <- drop(y - (beta0 + X%*%beta_vec))
  mu_nonNA <- mu
  mu_nonNA[is.na(mu)] = 0
  (-t(X)%*%res + rep(mu_nonNA, pk))
}
################################################################################

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
ll <- 1
#for(ll in 1:nrow(parset)){
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
  
#  for(r in 1:Rnum){
    r = 2
    cat("r is", r, '\n')
    set.seed(r)
    w <- mvlognormal(n=n,Mu=mm,Sigma=rep(1,p),R=sigma);
    z <- diag(1/rowSums(w)) %*% w
    
    b <- rep(0,p)
    b[1:5] <- c(1,-0.8,0.4,0,-0.6)
    b[6:10] <- c(-1.5,0,1.2,0,0.3)
    #sum(b)
    y <- drop(z %*% b)+rnorm(n,0,0.5^2)
    X <- log(z)
    K = length(pk) 
    
    
  # Setting index function #######################################################
  # starting and ending index in each group
  idx_gs <- cumsum(pk)-pk+1
  idx_ge <- cumsum(pk)
  # inverse search by matrix
  # single index ---> group index
  dict_idx_k <- rep(1:K, pk)
  # single index ---> within group index
  dict_idx_j <- unlist(lapply(pk ,FUN = seq, from = 1, by = 1))
  
  # active group  
  i_g_A <-  0L
  set_g_A <- rep(0,K)
  num_g_A <-  rep(0,K)
  
  # comlasso  
  y <- as.matrix(y,n,1)
  X_sum <- colSums(X)
  
  ### Variable definition
  mu <- rep(NA, K)
  lambda = 0
  beta_vec <- rep(0,p)
  beta_sign_vec<- rep(0,p)
  beta_vec_A<- rep(F,p)
  beta_mat <- Matrix(data=0,nrow= 10*n, ncol= 1 + p + 1 + K, sparse=TRUE)
  rec_idx <- 1
  # Initialization start ---------------------------------------------------------  
  # find the initial lambda and mu
  beta0 <- mean(y)
  res <- y-beta0
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

  # C code 
  # init_fit <- init_fun_C(idx_gs, idx_ge, grad_vec)
  # v <- init_fit$v
  # j1 <- init_fit$jvec[1]
  # j2 <- init_fit$jvec[2]
  
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
  
  beta_mat[rec_idx, 1] <- beta0
  beta_mat[rec_idx, 1+p+1] <-lambda
  beta_mat[rec_idx, 1+p+1+act_group] <- mu[act_group]
  rec_idx <- rec_idx + 1
  iter  <- 0
# Initialization end -----------------------------------------------------------  
  while(T)
  {
    # LOOP procudure start ---------------------------------------------------------
    iter <- iter + 1
    #cat("iter:::",iter,'\n')


    # Finding direction start --->
    a <- sum(beta_vec_A)
    q <- i_g_A
    pA_vec <- c(1, a, 1, q)
    pA <- sum(pA_vec)
    XA <- X[, beta_vec_A,drop = F]
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
    #solve_C(Dmat,bvec)
    # Finding direction end ---------------------------------------------------->
    
    # Finding delta start ------------------------------------------------------>
    delta <- rep(Inf,4)
    ### 2-1. Derive when active dummy variables become inactive
    step_vec <- -beta_vec[beta_vec_A]/rderiv[(1:a)+1]
    step_vec[step_vec < tol] = Inf
    delta[1] <- min(step_vec)
    #load("test.rdata")
    ### 2-2-a. Compute the distance needed to activate a new categorical variable
    dy <- rderiv[1] + XA%*%rderiv[(1:a)+1]
    corr_vec <- drop(t(X)%*%dy)
    act_type <- "none"
    #(i,j) = (j,j')
    
    act_vec <- setdiff(1:K, act_group)
    v <- Inf
    for(istar in act_vec)
    { 
      #if ( istar %in% act_group) next
      sidx <- idx_gs[istar]:idx_ge[istar]
      for (i in sidx)
      {
        for (j in sidx)
        {
          if (i == j ) next
          cjj <- -corr_vec[i]+ corr_vec[j] - 2*rderiv[1+a+1]
          djj <- 2*lambda + grad_vec[i] - grad_vec[j]
          v0 <- djj/cjj
          if (cjj>=0 & djj < 0) break ("stop:: KKT violation error 1\n")
          if (cjj< 0 & djj < 0) break ("stop:: KKT violation error 2\n")
          if (cjj< 0 & djj > 0) v0 <- Inf
          
          if (v0 < v & v0>tol) 
          {
            v <- v0
            j1 <- i
            j2 <- j
          }
        }
      }
    }
    delta[2]<- v

    if (v<Inf)
    {
      g1 <- grad_vec[j1] + corr_vec[j1]*v
      g2 <- grad_vec[j2] + corr_vec[j2]*v
      # jstar1: positive activation ; jstar2: negative activation 
      if (g1 < g2)    
      {
        jstar1 <- j1
        jstar2 <- j2
      } else {
        jstar2 <- j1
        jstar1 <- j2
      }
      
      j1_range_max <- lambda + rderiv[1+a+1]*v - g1
      j1_range_min <-  -(lambda + rderiv[1+a+1]*v) - g1
      j2_range_max <- lambda + rderiv[1+a+1]*v - g2
      j2_range_min <- -(lambda + rderiv[1+a+1]*v) - g2
      range_diff <- abs(j1_range_min-j2_range_max)
      
      # j1_range_max = lambda + rderiv[1+a+1-1]*v - g1;
      # j1_range_min = -(lambda + rderiv[1+a+1-1]*v) - g1;
      # j2_range_max = lambda + rderiv[1+a+1-1]*v - g2;
      # j2_range_min = -(lambda + rderiv[1+a+1-1]*v) - g2;
      # tmp = j1_range_min- j2_range_max;
      # range_diff = abs_j(tmp);
      # 
      if ( range_diff  < tol )  
      {
        js1 <- j1
        js2 <- j2
      } else {
        js1 <- j2
        js2 <- j1
      }
      # mu_tmp is computed
      mu_tmp <- lambda + rderiv[1+a+1]*v -grad_vec[js2]-corr_vec[js2]*v 
      # mu_tmp should be equal to the following value
      # mu_tmp <- -lambda -rderiv[1+a+1]*v -grad_vec[js1] -corr_vec[js1]*v 
    }
    
    
    # C-code
    # d2_fun_C(act_vec, idx_gs, idx_ge,grad_vec, corr_vec,rderiv,lambda,a,tol)
    
    ### 2-2-b. activated individual variable in the active categorical variable
    v <- Inf
    q_istar <- 0
    q_istar1 <- 0 
    act_sign <- 0
    tmp_sign <- 0
    j1 <- NA
    # load.image("test2-2.R)  
    for(istar in act_group)
    { 
      q_istar = q_istar + 1
      sidx <- idx_gs[istar]:idx_ge[istar]
      #cat("q_istar:", q_istar,"\n")
      for (j in sidx)
      {
        if (beta_sign_vec[j]!=0) next
        
        cj1 <- corr_vec[j] - rderiv[1+a+1] + rderiv[1+a+1+q_istar]
        dj1 <- -grad_vec[j] + lambda - mu[istar]
        
        v1 <- dj1/cj1
        if (cj1>0 & dj1 < 0) 
        {
          if (v1 < -tol) 
          {
            stop ("stop:: KKT violation error 3\n") 
          } else {
            next
          }
          
        }
        if (cj1< 0 & dj1 < 0) 
        {
          if ( v1 > tol ) 
          {
            stop ("stop:: KKT violation error 4\n")
          } else {
            v1 <- Inf
          }
        }
        
        if (cj1< 0 & dj1 > 0) v1 <- Inf
        
        cj2 <- -corr_vec[j] -rderiv[1+a+1]- rderiv[1+a+1+q_istar]
        dj2 <- grad_vec[j] + lambda + mu[istar]
        v2 <- dj2/cj2
        if (cj2>=0 & dj2 < 0) 
        {
          if (v2 < -tol) 
          {
            stop ("stop:: KKT violation error 5\n")
          } else {
            next
          }
        }
          
        if (cj2< 0 & dj2 < 0) 
        {
          if (v2>tol) 
          {
            stop ("stop:: KKT violation error 6\n")
          } else {
            v2 <- Inf
          }
        }  
          
        if (cj2< 0 & dj2 > 0) v2 <- Inf
        
        v0 <- min(v1,v2)
        
        
        if(v0<v)
        {
          v <- v0
          j1 <- j
          q_istar1 <- q_istar
        }
      }  
    }

    if (v < 0) break
    
    delta[3]<- v
    if ( v< Inf)
    {
      tmp <-  grad_vec[j1] + corr_vec[j1]*v + mu[dict_idx_k[j1]] +  
        rderiv[1+a+1+q_istar1]*v
      if (tmp<0) act_sign = 1 else act_sign = -1
    }
    jstar3 <- j1
    
    #d3_fun_C(rderiv,corr_vec,grad_vec,mu,lambda,a,tol,beta_sign_vec,
    #         act_group, idx_gs, idx_ge, dict_idx_k)
    
    ### 2-3. Compute the distance needed for lambda to become 0
    if(rderiv[1+a+1]>0){
      delta[4] <- 0
      break
    }
    
    
    lam_final <- (lam_min-lambda)/rderiv[1+a+1]
    delta[4] = ifelse(lam_final < tol, Inf, lam_final)
    delta_f <- min(delta)
#if (iter == max_iter) break   
    #delta_f = 0.0001
    # Finding delta end --------------------------------------------------------->
    
    # update variable star ------------------------------------------------------>
    
    # update beta_vec
    beta_vec[beta_vec_A] <- beta_vec[beta_vec_A]+delta_f*rderiv[(1:a)+1] 
    beta0 <- beta0 + delta_f*rderiv[1]
    lambda <- lambda + delta_f*rderiv[1+a+1]
    mu[act_group]   <- mu[act_group] + delta_f*rderiv[(1+a+1+1):pA] 
    
    # residual update
    res <- y -beta0 - XA%*%beta_vec[beta_vec_A]
    grad_vec <- -drop(t(X)%*%res)
    #KKT_fun(beta0,beta_vec,mu)
    
    
    beta_mat[rec_idx, 1] <- beta0
    beta_mat[rec_idx, 1+which(beta_vec_A)] <- beta_vec[beta_vec_A]
    beta_mat[rec_idx, 1+p+1] <- lambda
    beta_mat[rec_idx, 1+p+1+act_group] <- mu[act_group]
    rec_idx <- rec_idx + 1
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
    if ((which.min(delta) == 4) | lambda < tol ) act_type <- "t_act"
    
    cat("Type of update:", act_type, "\n")
    
    if (act_type == "g_van")
    {
      num_g_A[k_A] <- 0
      sidx <- idx_gs[k_A]:idx_ge[k_A]
      beta_vec[sidx]<- 0
      beta_sign_vec[sidx] <- 0
      beta_vec_A[sidx] <- FALSE
      mu[k_A] <- NA
      # correct set_g_A !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      del_idx <- which(set_g_A==k_A)
      set_g_A[del_idx] <- 0;
      act_group <- sort(set_g_A[set_g_A!=0])
      i_g_A <- i_g_A - 1;
      set_g_A <- rep(0,K)
      set_g_A[1:i_g_A] <- act_group
    }
    
    if (act_type == "i_van")
    {
      beta_sign_vec[v_idx] <- 0
      beta_vec_A[v_idx] <- FALSE
      beta_vec[v_idx] <- 0
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
    
    ## Check KKT conditions:
    if (KKT_check)
    {
      check_vec <- KKT_fun(X,y, beta0,beta_vec,mu)
      cb1 <- abs(check_vec[beta_vec_A] + lambda*beta_sign_vec[beta_vec_A])
      if (max(cb1)>tol) 
      {
        cat("KKT  stationarity cond (active) violated!!\n")
        stop()
      }
      
      tmp_idx <- which ( (dict_idx_k%in%act_group) & beta_sign_vec == 0 ) 
      if (length(tmp_idx)>0)
      {
        if ( max(abs(check_vec[tmp_idx])) > (lambda + tol ) ) 
        {
          cat("KKT  stationarity cond (inactive) violated!!\n")
          stop()
        }
      }
      
      for (k in 1:K)
      {
        if (k %in% act_group) next
        sidx <- idx_gs[k]:idx_ge[k]
        Dvec  <- check_vec[sidx]
        Dmax <- Dvec + lambda
        Dmin <- Dvec - lambda
        if (max(Dmin) > min(Dmax)+ tol) 
        {
          cat("KKT  stationarity cond (inactive group) violated!!\n")
          stop()
        }
      }
      
      beta_gsum <-unname(unlist(by(beta_vec, dict_idx_k, sum)))
      if (any(abs(beta_gsum)>tol))
      {
        cat("KKT primal feasibility: violated!!\n")
        stop()
      }  
    }
    
    if (act_type == "t_act")
    {
      cat("lambda is zero!\n")
      break
    }
    # update variable end ------------------------------------------------------>
    
    # KKT condition check
    # stationary condition
    # LOOP procudure end ---------------------------------------------------------
  }  
  
  coefficients <- beta_mat[1:(rec_idx-1),1:(p+1)]
  colnames(coefficients) <- c("b0", paste("b", dict_idx_k,dict_idx_j,sep="_" ))
  lambda_vec <- beta_mat[1:(rec_idx-1),1+p+1]
  # mu
  lambda_vec
  
# KKT_fun(beta0,beta_vec,mu)
#   
# beta_mat[1:(rec_idx-1),]
# act_group
# sum(beta_vec_A)
# beta_mat[1:rec_idx,]
