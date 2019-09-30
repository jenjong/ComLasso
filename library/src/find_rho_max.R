### Find the maximum tuning parameter value `ρmax` to kick-start the solution path.
# X::AbstractMatrix,
# y::AbstractArray;
# Aeq::AbstractMatrix = zeros(eltype(X), 0, size(X, 2)),
# beq::Union{AbstractArray, Number} = zeros(eltype(X), size(Aeq, 1)),
# Aineq::AbstractMatrix = zeros(eltype(X), 0, size(X, 2)),
# bineq::Union{AbstractArray, Number} = zeros(eltype(X), size(Aineq, 1)),
# penidx::Array{Bool} = fill(true, size(X, 2)),
# solver = ECOSSolver(maxit=10e8, verbose=0)

find_rho_max = function(X, Aeq, beq, Aineq, neq, nineq, penidx) { # global variable로 처리한다면
  p = dim(X)[2]

  x = Variable(p)
  objective = Minimize(sum_entries(abs(x))) # sum_entries!
  if(neq > 0 & nineq > 0) {
    constraints = list(Aeq %*% x == beq, Aineq %*% x <= bineq)
  } else if(neq > 0 & nineq == 0) {
    constraints = list(Aeq %*% x == beq)
  } else if(neq == 0 & nieq > 0) {
    constraints = list(Aineq %*% x <= bineq)
  }
   
  problem = Problem(objective, constraints)
  result = solve(problem)
  
  beta = result$getValue(x)
  
  lambda_eq = matrix(0, 0, 1)
  mu_ineq = matrix(0, 0, 1)
  
  for(i in 1:min(2, length(problem@constraints))) {
    if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinEqConstr") {
      lambda_eq = result$getDualValue(problem@constraints[[i]])
    } else if(canonicalize(problem@constraints[[i]])[[2]][[1]]$class == "LinLeqConstr") {
      mu_ineq = result$getDualValue(problem@constraints[[i]])
    }
  }
  
  if(dim(mu_ineq)[1] != 0) {
    mu_ineq[mu_ineq < 0] = 0
  }
  
  setActive = abs(beta) > 1e-4 | !penidx
  beta[!setActive] = 0
  
  resid = y - X %*% beta
  #cat("lambda_eq", lambda_eq, "\n")
  #stop()
  subgrad = t(X) %*% resid - t(Aeq) %*% lambda_eq - t(Aineq) %*% mu_ineq
  rho_max = max(abs(subgrad))
  ind_rho_max = which.max(abs(subgrad))
  
  return(list(rho_max = rho_max,
              ind_rho_max = ind_rho_max,
              optval = result$value,
              lambda_eq = lambda_eq,
              mu_ineq = mu_ineq)) # optval :: problem.optval
} 

# l = find_rho_max()



















