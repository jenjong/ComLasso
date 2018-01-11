rm(list = ls())
library(Matrix)
library(Rcpp)
sourceCpp('D:/Jeon/rcode/ComLasso/library/inner.cpp')
source("D:/Jeon/rcode/ComLasso/library/ComLassoC.R")
source("D:/Jeon/rcode/ComLasso/library/ComLassoR.R")
n = 100
pk = rep(c(3,3,4),5000)
p = sum(pk)
set.seed(3)
X = matrix(rnorm(n*p), n, p)
y = rnorm(n)
lam_min=0
tol=1e-08
system.time({fit <- comLassoC(X,y,pk,lam_min,tol=1e-8, KKT_check = FALSE)})
system.time({fit <- comLasso(X,y,pk,lam_min,tol=1e-8, KKT_check = FALSE)})
#fit$coefMat[,1]
#fit$lambda_vec

