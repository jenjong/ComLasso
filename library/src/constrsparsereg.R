# Fit constrained lasso at fixed tuning parameter value(s) by minimizing
# `0.5sumabs2(√obswt .* (y - X * β)) + ρ * sumabs(penwt .* β)`
# subject to linear constraints, using `Convex.jl`.

################################################################
constrsparsereg = function(rho,
                           obswt = matrix(rep(1, n), nrow = n),
                           penwt, 
                           warmstart = F,
                           neq = 0, nineq=0) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  # function 인자값 setting
  # penwt = penidx
  beta_hat = matrix(rep(0, p*length(rho)), nrow = p, ncol = length(rho))
  optval_vec = rep(0, length(rho))
  prob_vec = list()
  result_vec = list()
  
  beta = Variable(p)
  loss = (1 / 2) * sum_squares(sqrt(obswt) * (y - X %*% beta)) # loss term
  pen = penwt %*% abs(beta)
  
  # get constraints
  if(neq > 0 & nineq > 0) {
    constraints = list(Aeq %*% beta == beq, Aineq %*% beta <= bineq)
  } else if(neq > 0 & nineq == 0) {
    constraints = list(Aeq %*% beta == beq)
  } else if(neq == 0 & nieq > 0) {
    constraints = list(Aineq %*% beta <= bineq)
  }
  
  for(i in 1:length(rho)) {
    rho_i = rho[i]
    if(rho_i == Inf) {
      objective = Minimize(pen)
      problem = Problem(objective, constraints)
    } else if(rho_i <= 0) {
      objective = Minimize(loss)
      problem = Problem(objective, constraints)
    } else {
      objective = Minimize(loss + rho_i * pen)
      problem = Problem(objective, constraints)
    }
    
    if(warmstart) {
      result = solve(problem, warm_start = i) # warm_start 잘 모르겠음...
      prob_vec[[i]] = problem
      result_vec[[i]] = result
      optval_vec[i] = result$value
    } else {
      result = solve(problem) 
      prob_vec[[i]] = problem
      result_vec[[i]] = result
      optval_vec[i] = result$value
    }
    
    if(length(rho) == 1) {
      return(list(beta_value = result$getValue(beta),
                  problem.optval = result$value,
                  problem = problem,
                  result = result))
    }
    
    beta_hat[, i] = result$getValue(beta)
  }
  
  return(list(beta_hat = beta_hat,
              optval_vec = optval_vec,
              prob_vec = prob_vec))
}

