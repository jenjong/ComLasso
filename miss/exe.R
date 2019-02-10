# comparison of comlass with help of Rcpp
rm(list = ls())
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
#setwd("C:/Users/Jeon/Documents/GitHub/ComLasso")
sourceCpp('./library/inner.cpp')
source("./library/ComLassoC.R")
source("./library/ComLassoR.R")
n = 100
pk = rep(c(300,300,400),2)
p = sum(pk)
sim_iter = 20
r1 = r2 = c()
for (i in 1:sim_iter)
{
  cat("iter:", i,'\n')
  set.seed(i)
  
  X = matrix(rnorm(n*p), n, p)
  y = rnorm(n)
  lam_min=0
  tol=1e-08
  s1 = system.time({fit <- comLassoC(X,y,pk,lam_min,tol=1e-8, KKT_check = FALSE)})
  s2 = system.time({fit <- comLasso(X,y,pk,lam_min,tol=1e-8, KKT_check = FALSE)})
  r1[i]= s1[3]
  r2[i]= s2[3]
}
boxplot(r1,r2)

