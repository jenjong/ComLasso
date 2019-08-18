# Maintainer
# 

KKT_fun<- function(X,y,beta0, beta_vec, mu)
{
  res <- drop(y - (beta0 + X%*%beta_vec))
  mu_nonNA <- mu
  mu_nonNA[is.na(mu)] <- 0
  (-t(X)%*%res + rep(mu_nonNA, pk))
}

comLassoC <- function(X,y,pk,lam_min,max_iter=1e+5,tol=1e-8, KKT_check = TRUE)
{
  n = length(y)
  K = length(pk)
  p = sum(pk)
 # Setting index function #######################################################
  # starting and ending index in each group
  idx_gs <- cumsum(pk)-pk+1
  idx_ge <- cumsum(pk)
  # inverse search by matrix
  # single index ---> group index
  dict_idx_k <- rep(1:K, pk)
  # single index ---> within group index
  dict_idx_j <- unlist(lapply(pk, FUN = seq, from = 1, by = 1))
  
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
  init_fit <- init_fun_C(idx_gs, idx_ge, grad_vec)
  v <- init_fit$v
  j1 <- init_fit$jvec[1]
  j2 <- init_fit$jvec[2]
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
    if (iter == max_iter) break
    
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
    
    # Finding direction end ---------------------------------------------------->
    
    # Finding delta start ------------------------------------------------------>
    delta <- rep(Inf,4)
    ### 2-1. Derive when active dummy variables become inactive
    step_vec <- -beta_vec[beta_vec_A]/rderiv[(1:a)+1]
    step_vec[step_vec < tol] = Inf
    delta[1] <- min(step_vec)
    ### 2-2-a. Compute the distance needed to activate a new categorical variable
    dy <- rderiv[1] + XA%*%rderiv[(1:a)+1]
    corr_vec <- drop(t(X)%*%dy)
    act_type <- "none"
    act_vec <- setdiff(1:K, act_group)
    d2_fit <- d2_fun_C(act_vec, idx_gs, idx_ge,grad_vec, corr_vec,rderiv,lambda,a,tol)
    delta[2]<- d2_fit$v
    ### 2-2-b. activated individual variable in the active categorical variable
    d3_fit = d3_fun_C(rderiv,corr_vec,grad_vec,mu,lambda,a,tol,beta_sign_vec,
                      act_group, idx_gs, idx_ge, dict_idx_k)
    
    if (d3_fit$v < 0) 
    {
      break
    }
    delta[3] <- d3_fit$v
    ### 2-3. Compute the distance needed for lambda to become 0
    if(rderiv[1+a+1]>0){
      delta[4] <- 0
      break
    }
    lam_final <- (lam_min-lambda)/rderiv[1+a+1]
    delta[4] = ifelse(lam_final < tol, Inf, lam_final)
    delta_f <- min(delta)
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
    
    # check degree of freedom
    if ( n <= sum(beta_vec_A)-length(act_group)+1) break
    
    
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
    if ((which.min(delta) == 4) | lambda < tol) act_type <- "t_act"
    
     #cat("Type of update:", act_type, "\n")
    
    if (act_type == "g_van")
    {
      num_g_A[k_A] <- 0
      sidx <- idx_gs[k_A]:idx_ge[k_A]
      beta_vec[sidx]<- 0
      beta_sign_vec[sidx] <- 0
      beta_vec_A[sidx] <- FALSE
      mu[k_A] <- NA
      # correct set_g_A  for improving algorithm !!!
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
      jstar1 = d2_fit$jvec[1]
      jstar2 = d2_fit$jvec[2]
      beta_sign_vec[jstar1] <- 1
      beta_sign_vec[jstar2] <- -1
      beta_vec_A[jstar1] <- TRUE
      beta_vec_A[jstar2] <- TRUE
      k_A <- dict_idx_k[jstar1]
      mu[k_A] <- d2_fit$mu
      i_g_A <- i_g_A + 1
      set_g_A[i_g_A] <- k_A
      num_g_A[k_A] <- 2
      act_group <- sort(set_g_A[1:i_g_A])
    }
    
    if (act_type == "i_act")
    {
      jstar3 <- d3_fit$jstar3
      act_sign <- d3_fit$act_sign
      beta_sign_vec[jstar3] <- act_sign
      beta_vec_A[jstar3] <- TRUE
      k_A <- dict_idx_k[jstar3]
      num_g_A[k_A] <- num_g_A[k_A] + 1
    }

    if (KKT_check)
    {
      ## Check KKT conditions:
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
          cat("KKT  stationarity cond (inactive) violated!!\n")
          cat("max(Dmin)", max(Dmin),'\n')
          cat("min(Dmax)", min(Dmax),'\n')
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
      #cat("lambda is zero!\n")
      break
    }
    # update variable end ------------------------------------------------------>
    # LOOP procudure end ---------------------------------------------------------
  } 
#  cat("the df of the final model is", sum(beta_vec_A)-length(act_group)+1, "\n")
#  cat("the nobs is", n, "\n")
#  cat("p is", p,  "\n")  
  coefMat <- beta_mat[1:(rec_idx-1),1:(p+1)]
  colnames(coefMat) <- c("b0", paste("b", dict_idx_k,dict_idx_j,sep="_" ))
  lambda_vec <- beta_mat[1:(rec_idx-1),1+p+1]
  return(list(coefMat = coefMat, lambda_vec = lambda_vec))
}  
  
  

