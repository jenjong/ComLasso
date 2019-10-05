source('./library/src/constrsparsereg.R')
source('./library/src/find_rho_max.R')
if(!require("CVXR")) install.packages("CVXR")
library("CVXR")

zhou <- function(X, y, penwt, Aeq, beq, Aineq, bineq){
  # grhobla variable로 처리하면 좋을 것 같은 변수 *****
  # : Aeq, beq, Aineq, bineq, ρridge, penidx
  
  ### lsq_classopath 내부
  # X      :: AbstractMatrix{T},
  # y      :: AbstractArray{T};
  # Aeq    :: AbstractMatrix = zeros(T, 0, size(X, 2)),
  # beq    :: Union{AbstractArray, Number} = zeros(T, size(Aeq, 1)),
  # Aineq  :: AbstractMatrix = zeros(T, 0, size(X, 2)),
  # bineq  :: Union{AbstractArray, Number} = zeros(T, size(Aineq, 1)),
  # ρridge :: Number = zero(T),
  # penidx :: Array{Bool} = fill(true, size(X, 2)),
  # solver = ECOSSolver(maxit=10e8, verbose=0)
  n = dim(X)[1]
  p = dim(X)[2]
  
  rho_ridge = 0
  penidx = rep(TRUE, p)
  if(n < p) {
    print("Adding a small ridge penalty (default is 1e-4) since n < p")
    if(rho_ridge <= 0) {
      print("ridge must be positive, switching to default value (1e-4)")
      rho_ridge = 1e-4
    }
    # create augmented data
    y = matrix(c(y, rep(0, p)), nrow = n+p, ncol=1)
    X = rbind(X, sqrt(rho_ridge) * diag(rep(1, p)))
    # record original number of observations
    # n_orig = n
  } else {
    # make sure X is full column rank
    qrfact_X = qr(X)
    R = qr.R(qrfact_X)
    rankX = sum(abs(diag(R)) > (abs(R[1,1]) * max(n,p) * (.Machine$double.eps) ^ 2))  # 4.930380657631324e-32 *****
    
    if(rankX != p) {
      print("Adding a small ridge penalty (default is 1e-4) since X is rank deficient")
      if(rho_ridge <= 0) {
        print("ridge must be positive, switching to default value (1e-4)")
        rho_ridge = 1e-4
      }
      # create augmented data
      y = rbind(y, matrix(rep(0, p), nrow = p))
      X = rbind(X, sqrt(rho_ridge) * diag(rep(1, p)))
    }
  }
  
  # alrhocate variables arhong path
  neq = dim(Aeq)[1]
  nineq = dim(Aineq)[1]
  maxiters = 5 * (p + nineq) # max number of path segments to consider
  beta_path = matrix(rep(0, p * maxiters), nrow = p) 
  lambda_patheq = matrix(rep(0, neq * maxiters), nrow = neq) # dual variables for equality
  mu_pathineq = matrix(rep(0, nineq * maxiters), nrow = nineq) # dual variables for inequality
  rho_path = rep(0, maxiters) # tuning parameter
  df_path = rep(Inf, maxiters) # degree of freedom
  objval_path = rep(0, maxiters) # objective value
  objval_path2 = rep(0, maxiters) # objective value
  violation_path = rep(Inf, maxiters) 
  
  ### initialization
  # use LP to find ρmax
  H = t(X) %*% X 
  #sense = [repmat(['='], neq); repmat(['<'], nineq)]
  #β, _, problem = lsq_constrsparsereg(X, y, Inf;      # why necessary?
  #          [Aeq; Aineq], sense, [beq; bineq], lb, ub)#
  
  # find the maximum ρ (starting value)
  
  l = find_rho_max(X, Aeq, beq, Aineq, neq, nineq, penidx)
  
  rho_path[1] = l$rho_max
  #rho_path[1] = 144.0489670
  idx = l$ind_rho_max
  #cat("initial idx=", idx, "\n")
  
  # calculate at ρmax
  l = constrsparsereg(rho = rho_path[1], penwt = penidx, neq=neq, nineq=nineq)
  beta_path[, 1] = l$beta_value
  objval_path[1] = l$problem.optval
  problem = l$problem
  result = l$result
  
  for(i in 1:min(2, length(problem@constraints))) {
    if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinEqConstr") {
      # lambda_patheq[, 1] = result$getDualValue(problem@constraints[[i]]) # 이거 값 다름 - 왜 0 나오지?*****
      lambda_patheq[, 1] = result$getDualValue(problem@constraints[[i]]) # 이거 값 다름 - 왜 0 나오지?***** (부호변경)
    } else if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinLeqConstr") {
      mu_pathineq[, 1] = result$getDualValue(problem@constraints[[i]])
    }
  }
  
  #lambda_patheq[, 1] = -7.1592458 # 강제로 값을 지정*****
  
  mu_pathineq[mu_pathineq < 0] = 0
  setActive = (abs(beta_path[, 1]) > 1e-4) | (!penidx)
  beta_path[!setActive, 1] = 0
  
  residIneq = Aineq %*% beta_path[, 1] - bineq
  setIneqBorder = residIneq == 0
  nIneqBorder = sum(setIneqBorder != 0)
  
  # initialize subgradient vector
  resid = y - X %*% beta_path[, 1]
  
  if(neq > 0 & nineq > 0) {
    subgrad = t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1] - t(Aineq) %*% mu_pathineq[, 1]
  } else if(neq > 0 & nineq == 0) {
    subgrad = t(X) %*% resid - t(Aeq) %*% lambda_patheq[ ,1]
  } else if(neq == 0 & nieq > 0) {
    subgrad = t(X) %*% resid - t(Aineq) %*% mu_pathineq[, 1]
  }
  
  subgrad[setActive] = sign(beta_path[setActive, 1])
  subgrad[!setActive] = subgrad[!setActive] / rho_path[1]
  setActive[idx] = T
  nActive = sum(setActive != 0)
  
  ###
  #cat("idx=", idx, "\n")
  ###
  # calculate degrees of freedom
  rankAeq = Matrix::rankMatrix(Aeq) ### Brian's comment: need to make it more efficient
  df_path[1] = nActive - rankAeq - nIneqBorder
  
  # set initial violations counter to 0
  violation_path[1] = 0
  
  # sign for path direction (originally went both ways, but increasing was retired)
  dirsgn = -1
  
  ####################################
  ### main loop for path following ###
  ####################################
  
  k = 0
  
  for(k in 2:maxiters) {
    #cat("=====================\n")
    #cat(k, "th iteration", "\n")
    #cat(which(setActive), "\n")
    #print(beta_path[, k-1])
    #cat("=====================\n")
    # classo_for(k)_part
    
    # threshold near-zero rhos to zero and stop algorithm
    if(rho_path[k-1] <= 1e-4) {
      rho_path[k-1] = 0
      break
    }
    
    # calculate derivative for coefficients and multipliers
    # construct matrix
    activeCoeffs = which(setActive)
    inactiveCoeffs = which(!setActive)
    idxIneqBorder = which(setIneqBorder)
    
    # 여기 계산 불안정 - 제약조건이 1개만(등호 or 부등호) 있을 때*****
    M = cbind(H[activeCoeffs, activeCoeffs, drop=F], t(Aeq[, activeCoeffs, drop=F]), 
              t(Aineq[setIneqBorder, activeCoeffs]))
    M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
    M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs, drop=F], Aineq[idxIneqBorder, activeCoeffs])
    
    # calculate derivative
    tryCatch(
      expr = {
        dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
      }, 
      error = {
        dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
        # julia는 pinv 함수 씀
      }
    )
    
    if(any(is.nan(dir))) {
      dir = -(MASS::ginv(M) %*% rbind(subgrad[setActive], matrix(rep(0, neq + nIneqBorder), ncol = 1)))
    }
    
    # calculate the derivative for rho * subgradient
    # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
    temp = cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]))
    if(nineq == 0) {
      dirSubgrad = - temp %*% dir
    } else {
      dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
    }
    
    ### check additional events related to potential subgraient violations
    
    ## inactive coefficients moving too slowly
    # negative subgradient
    inactSlowNegIdx = which((1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                              subgrad[!setActive] <= (1 * dirsgn + 1e-8) &
                              1 * dirsgn < dirSubgrad)
    
    # positive subgradient
    inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                              subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                              -1 * dirsgn > dirSubgrad)
    
    ## "active" coeficients estimated as 0 with potential sign mismatch #%
    # positive subgrad but negative derivative
    signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                 subgrad[setActive] <= (1 + 1e-8) &
                                 dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                                 beta_path[activeCoeffs, k-1] == 0)
    
    # Negative subgradient but positive derivative
    signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                 subgrad[setActive] <= (0 + 1e-8) &
                                 dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                 beta_path[activeCoeffs, k-1] == 0)
    
    # reset violation counter (to avoid infinite loops)
    violateCounter = 0
    
    #print(sum(setActive))
    
    # outer while loop for checking all conditions together
    while(length(inactSlowNegIdx) > 0 | length(inactSlowPosIdx) > 0 |
          length(signMismatchPosIdx) > 0 | length(signMismatchNegIdx) > 0) {
      
      # 총 4개의 while문이 들어갈 예정
      
      # monitor and fix condition 1 violations
      while(length(inactSlowNegIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        inactiveCoeffs = which(!setActive)
        # identify problem coefficient
        viol_coeff = inactiveCoeffs[inactSlowNegIdx]
        setActive[viol_coeff] = T
        
        nActive = sum(setActive != 0)
        nIneqBorder = sum(setIneqBorder != 0)
        
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop = F], t(Aeq[, activeCoeffs, drop = F]), 
                  t(Aineq[setIneqBorder, activeCoeffs]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs], Aineq[idxIneqBorder, activeCoeffs])
        
        # calculate derivative
        tryCatch(
          expr = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
          }, 
          error = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
            # julia는 pinv 함수 씀
          }
        )
        
        # calculate the derivative for rho * subgradient
        # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
        temp = cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]))
        if(nineq == 0) {
          dirSubgrad = - temp %*% dir
        } else {
          dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
        }
        
        ## Misc. housekeeping #%
        # check for violations again
        
        ## inactive coefficients moving too slowly
        # negative subgradient
        inactSlowNegIdx = which((1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                                  subgrad[!setActive] <= (1 * dirsgn + 1e-8) &
                                  1 * dirsgn < dirSubgrad)
        
        # positive subgradient
        inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                                  subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                                  -1 * dirsgn > dirSubgrad)
        
        ## "active" coeficients estimated as 0 with potential sign mismatch #%
        # positive subgrad but negative derivative
        signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (1 + 1e-8) &
                                     dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # Negative subgradient but positive derivative
        signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (0 + 1e-8) &
                                     dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      # Monitor & fix subgradient condition 2 violations
      while(length(inactSlowPosIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        inactiveCoeffs = which(!setActive)
        # identify problem coefficient
        viol_coeff = inactiveCoeffs[inactSlowPosIdx]
        setActive[viol_coeff] = T
        
        nActive = sum(setActive != 0)
        nIneqBorder = sum(setIneqBorder != 0)
        
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop = F], t(Aeq[, activeCoeffs, drop = F]), 
                  t(Aineq[setIneqBorder, activeCoeffs]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs], Aineq[idxIneqBorder, activeCoeffs])
        
        # calculate derivative
        tryCatch(
          expr = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
          }, 
          error = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
            # julia는 pinv 함수 씀
          }
        )
        
        # calculate the derivative for rho * subgradient
        # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
        temp = cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]))
        if(nineq == 0) {
          dirSubgrad = - temp %*% dir
        } else {
          dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
        }
        
        ## Misc. housekeeping #%
        # check for violations again
        
        ## inactive coefficients moving too slowly
        # positive subgradient
        inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                                  subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                                  -1 * dirsgn > dirSubgrad)
        
        ## "active" coeficients estimated as 0 with potential sign mismatch #%
        # positive subgrad but negative derivative
        signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (1 + 1e-8) &
                                     dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # Negative subgradient but positive derivative
        signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (0 + 1e-8) &
                                     dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      # Monitor & fix condition 3 violations
      while(length(signMismatchPosIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        activeCoeffs = which(setActive)
        # identify problem coefficient
        viol_coeff = activeCoeffs[signMismatchPosIdx]
        # put problem coefficient back into inactive set;
        setActive[viol_coeff] = F
        # determine new number of active coefficients
        nActive = sum(setActive != 0)
        # determine number of active/binding inequality constraints
        nIneqBorder = sum(setIneqBorder != 0)
        
        ## Recalculate derivative for coefficients & multipliers #%
        # construct matrix
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop = F], t(Aeq[, activeCoeffs, drop = F]), 
                  t(Aineq[setIneqBorder, activeCoeffs]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs], Aineq[idxIneqBorder, activeCoeffs])
        
        # calculate derivative
        tryCatch(
          expr = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
          }, 
          error = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
            # julia는 pinv 함수 씀
          }
        )
        
        # calculate the derivative for rho * subgradient
        # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
        temp = cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop = F]))
        if(nineq == 0) {
          dirSubgrad = - temp %*% dir
        } else {
          dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
        }
        
        ## Misc. housekeeping #%
        # check for violations again
        ## "active" coeficients estimated as 0 with potential sign mismatch #%
        # positive subgrad but negative derivative
        signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (1 + 1e-8) &
                                     dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # Negative subgradient but positive derivative
        signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (0 + 1e-8) &
                                     dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      # Monitor & fix condition 3 violations
      while(length(signMismatchNegIdx) > 0) {
        ## Identify and move problem coefficient
        # indices corresponding to inactive coefficients
        activeCoeffs = which(setActive)
        # identify problem coefficient
        viol_coeff = activeCoeffs[signMismatchNegIdx]
        # put problem coefficient back into inactive set;
        setActive[viol_coeff] = F
        # determine new number of active coefficients
        nActive = sum(setActive != 0)
        # determine number of active/binding inequality constraints
        nIneqBorder = sum(setIneqBorder != 0)
        
        ## Recalculate derivative for coefficients & multipliers #%
        # construct matrix
        activeCoeffs = which(setActive)
        inactiveCoeffs = which(!setActive)
        idxIneqBorder = which(setIneqBorder)
        
        M = cbind(H[activeCoeffs, activeCoeffs, drop = F], t(Aeq[, activeCoeffs, drop = F]), 
                  t(Aineq[setIneqBorder, activeCoeffs]))
        M = rbind(M, matrix(rep(0, (neq + nIneqBorder) * dim(M)[2]), nrow = (neq + nIneqBorder)))
        M[(nrow(M) - neq - nIneqBorder + 1):nrow(M), 1:nActive] = rbind(Aeq[, activeCoeffs], Aineq[idxIneqBorder, activeCoeffs])
        
        # calculate derivative
        tryCatch(
          expr = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
          }, 
          error = {
            dir = -(MASS::ginv(M) %*% rbind(matrix(subgrad[setActive], ncol = 1), matrix(rep(0, neq + nIneqBorder), ncol = 1)))
            # julia는 pinv 함수 씀
          }
        )
        
        # calculate the derivative for rho * subgradient
        # 여기 계산 불안정 - 제약조건이 1개만 있을 때*****
        temp = cbind(H[inactiveCoeffs, activeCoeffs, drop=F], t(Aeq[, inactiveCoeffs, drop=F]))
        if(nineq == 0) {
          dirSubgrad = - temp %*% dir
        } else {
          dirSubgrad = - cbind(temp, Aineq[idxIneqBorder, inactiveCoeffs]) %*% dir  
        }
        
        ## Misc. housekeeping #%
        # check for violations again
        # Negative subgradient but positive derivative
        signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                     subgrad[setActive] <= (0 + 1e-8) &
                                     dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                     beta_path[activeCoeffs, k-1] == 0)
        
        # update violation counter
        violateCounter = violateCounter + 1
        # break loop if needed
        if(violateCounter >= maxiters) {
          break
        }
      }
      
      ## update violation trackers to see if any issues persist ##%
      activeCoeffs = which(setActive)
      inactiveCoeffs = which(!setActive)
      idxIneqBorder = which(setIneqBorder)
      
      ## inactive coefficients moving too slowly
      # negative subgradient
      inactSlowNegIdx = which((1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                                subgrad[!setActive] <= (1 * dirsgn + 1e-8) &
                                1 * dirsgn < dirSubgrad)
      
      # positive subgradient
      inactSlowPosIdx = which((-1 * dirsgn - 1e-8) <= subgrad[!setActive] &
                                subgrad[!setActive] <= (-1 * dirsgn + 1e-8) &
                                -1 * dirsgn > dirSubgrad)
      
      ## "active" coeficients estimated as 0 with potential sign mismatch #%
      # positive subgrad but negative derivative
      signMismatchPosIdx = which((0 - 1e-8) <= subgrad[setActive] &
                                   subgrad[setActive] <= (1 + 1e-8) &
                                   dirsgn * dir[1:nActive] <= (0 - 1e-8) &
                                   beta_path[activeCoeffs, k-1] == 0)
      
      # Negative subgradient but positive derivative
      signMismatchNegIdx = which((-1 - 1e-8) <= subgrad[setActive] &
                                   subgrad[setActive] <= (0 + 1e-8) &
                                   dirsgn * dir[1:nActive] >= (0 + 1e-8) &
                                   beta_path[activeCoeffs, k-1] == 0)
      
      # update violation counter
      violateCounter = violateCounter + 1
      # break loop if needed
      if(violateCounter >= maxiters) {
        break
      }
    }
    
    # after all while loop
    ###### after all while loop
    
    # store number of violations
    violation_path[k] = violateCounter
    
    # calculate derivative for residual inequality
    activeCoeffs = which(setActive)
    inactiveCoeffs = which(!setActive)
    idxIneqBorder = which(setIneqBorder)
    
    dirResidIneq = Aineq[which(!setIneqBorder), activeCoeffs, drop=F] %*% dir[1:nActive] # 계산 불안정*****
    
    ### Determine rho for next event (via delta rho) ###%
    ## Events based on coefficients changing activation status ##%
    next_rho_beta = matrix(rep(Inf, p), ncol = 1)
    
    # Active coefficient going inactive #%
    next_rho_beta[setActive] = -dirsgn * beta_path[activeCoeffs, k-1] / dir[1:nActive]
    
    # Inactive coefficient becoming positive #%
    t1 = dirsgn * rho_path[k-1] * (1 - subgrad[inactiveCoeffs]) / (dirSubgrad - 1)
    # threshold values hitting ceiling
    t1[t1 < 1e-8] = Inf
    
    # Inactive coefficient becoming negative #%
    t2 = -dirsgn * rho_path[k-1] * (1 + subgrad[!setActive]) / (dirSubgrad + 1)
    # threshold values hitting ceiling
    t2[t2 < 1e-8] = Inf
    
    # choose smaller delta rho out of t1 and t2
    next_rho_beta[!setActive] = pmin(t1, t2) # elementwise min
    next_rho_beta[(next_rho_beta <= 1e-8) | !penidx] = Inf
    
    ## Events based inequality constraints ##%
    # clear previous values
    next_rho_Ineq = matrix(rep(Inf, nineq), ncol = 1)
    
    # Inactive inequality constraint becoming active #%
    next_rho_Ineq[!setIneqBorder] = matrix(-dirsgn * residIneq[!setIneqBorder], sum(!setIneqBorder != 0), 1) /
      matrix(dirResidIneq, sum(!setIneqBorder != 0), 1)
    
    # Active inequality constraint becoming deactive #%
    if(length(mu_pathineq) > 0) {
      next_rho_Ineq[setIneqBorder] = -dirsgn * mu_pathineq[idxIneqBorder, k-1] / 
        matrix(dir[nActive + neq + 1, ], nIneqBorder, 1)
    }
    
    next_rho_Ineq[next_rho_Ineq <= 1e-8] = Inf
    
    ## determine next rho ##
    # find smallest rho
    chg_rho = min(rbind(next_rho_beta, next_rho_Ineq), na.rm = T)
    # find all indices corresponding to this chgρ
    idx = which(rbind(next_rho_beta, next_rho_Ineq) - chg_rho <= 1e-8)
    
    # terminate path following if no new event found
    if(chg_rho == Inf) {
      chg_rho = rho_path[k-1]
    }
    
    ## Update values at new rho ##%
    # move to next rho #%
    # make sure next rho isn't negative
    if(rho_path[k-1] + dirsgn * chg_rho < 0) {
      chg_rho = rho_path[k-1]
    }
    
    # calculate new value of rho
    rho_path[k] = rho_path[k-1] + dirsgn * chg_rho
    
    ## Update parameter and subgradient values #%
    # new coefficient estimates
    activeCoeffs = which(setActive)
    
    beta_path[activeCoeffs, k] = beta_path[activeCoeffs, k-1] + dirsgn * chg_rho * dir[1:nActive]
    # force near-zero coefficients to be zero (helps with numerical issues)
    beta_path[abs(beta_path[, k]) < 1e-12, k] = 0
    
    # new subgradient estimates
    subgrad[!setActive] = (rho_path[k-1] * subgrad[!setActive] + 
                             dirsgn * chg_rho * dirSubgrad) / rho_path[k]
    
    # Update dual variables #%
    # update lambda (lagrange multipliers for equality constraints)
    if(length(lambda_patheq) > 0) {
      lambda_patheq[, k] = lambda_patheq[, k-1] + 
        dirsgn * chg_rho * matrix(dir[(nActive + 1):(nActive + neq), ], neq, 1)
    }
    # update mu (lagrange multipliers for inequality constraints)
    if(length(mu_pathineq) > 0) {
      mu_pathineq[idxIneqBorder, k] = mu_pathineq[idxIneqBorder, k-1] + 
        dirsgn * chg_rho * matrix(dir[(nActive + neq + 1) : ncol(dir)], nIneqBorder, 1)
    }
    # update residual inequality
    residIneq = Aineq * beta_path[, k] - bineq
    
    ## update sets ##%
    for(j in 1:length(idx)) {
      curidx = idx[j]
      #cat("curidx", curidx, "\n")
      if(curidx <= p & setActive[curidx]) {
        # an active coefficient hits 0, or
        setActive[curidx] = F
      } else if(curidx <= p & !setActive[curidx]) {
        # a zero coefficient becomes nonzero
        setActive[curidx] = T
      } else if(curidx > p) {
        # an ineq on boundary becomes strict, or
        # a strict ineq hits boundary
        setIneqBorder[curidx - p] = !setIneqBorder[curidx - p]
      }
    }
    
    # determine new number of active coefficients
    nActive = sum(setActive != 0)
    # determine number of active/binding inequality constraints
    nIneqBorder = sum(nIneqBorder != 0)
    
    ## Calcuate and store values of interest along the path #%
    # calculate value of objective function
    objval_path[k] = norm(y - X %*% beta_path[, k], type = "F")^2 / 2 + 
      rho_path[k] * sum(abs(beta_path[, k]))
    
    objval_path[k] = mean((y - X %*% beta_path[, k])^2)/2 + rho_path[k] * sum(abs(beta_path[, k]))
    objval_path2[k] = mean((y - X %*% beta_path[, k])^2)/2
    # calculate degrees of freedom
    df_path[k] = nActive - rankAeq - nIneqBorder
    # break algorithm when df are exhausted
    if(df_path[k] >= n) {
      break
    }
    
    # end of big for loop
  }
  rho_path = rho_path[2:(k-1)]
  beta_path=beta_path[,2:(k-1)]
  objval_path = objval_path[2:(k-1)]
  objval_path2 = objval_path2[2:(k-1)]
  
  list(rho_path=rho_path, beta_path=beta_path, 
       objval_path=objval_path, objval_path2=objval_path2)
}
predict.zhou <- function (object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
                                                                                     "fraction", "norm", "lambda"), ...) 
{
  lams <- sort(object$rho_path, decreasing = TRUE)
  bb <- c()
  bb <- rbind(bb, object$beta_path)
  betas <- t(rbind(rep(0,ncol(bb)), bb))
  mode <- match.arg(mode)
  type <- match.arg(type)
  if (missing(newx) & type == "fit") {
    warning("Type=fit with no newx argument; type switched to coefficients")
    type <- "coefficients"
  }
  dimnames(betas) = list(NULL, dimnames(betas)[[2]])
  sbetas <- betas
  kp <- dim(sbetas)
  k <- kp[1]
  p <- kp[2]
  steps <- seq(k)
  if (missing(s)) {
    s <- steps
    mode <- "step"
  }
  sbeta <- switch(mode, step = {
    if (any(s < 0) | any(s > k)) {
      stop("Argument s out of range")
    }
    steps
  }, fraction = {
    if (any(s > 1) | any(s < 0)) {
      print(c(s, k))
      stop("Argument s out of range")
    }
    nbeta <- drop(abs(sbetas[, -1]) %*% rep(1, p - 1))
    nbeta/(nbeta[k])
  }, norm = {
    nbeta <- drop(abs(sbetas[, -1]) %*% rep(1, p - 1))
    if (any(s > nbeta[k]) | any(s < 0)) stop("Argument s out of range")
    nbeta
  }, lambda = {
    lambdas = lams
    s[s > max(lambdas)] = max(lambdas)
    s[s < 0] = 0
    lambdas
  })
  sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
  sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
  usbeta <- unique(sbeta)
  useq <- match(usbeta, sbeta)
  sbeta <- sbeta[useq]
  betas <- betas[useq, , drop = FALSE]
  coord <- approx(sbeta, seq(sbeta), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] + 
                 (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] - 
                                                                          sbeta[left])
  newbetas[left == right, ] <- betas[left[left == right], ]
  robject <- switch(type, coefficients = list(s = s, fraction = sfrac, 
                                              mode = mode, coefficients = drop(newbetas)), 
                    fit = list(s = s,  fraction = sfrac, mode = mode, fit = drop(cbind(1, newx) %*% 
                                                                                                                                                  t(newbetas))))
  robject
}



















