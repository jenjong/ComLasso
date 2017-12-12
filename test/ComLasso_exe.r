rm(list = ls()) ; gc()
# Rtools is required and Path for rtools should be set.
library(Rcpp)
library(inline)
library(RcppArmadillo)
setwd("D:/Jeon/rcode/ComLasso/test")
sourceCpp('inner.cpp')

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
comlasso_l2(n=n,p=8,K=c(10,6,4,3,7,2,8,p-40),X.raw=z,y=y, 
                    weights=NULL,max.steps=8*min(n,40)+1,lam.min=0,tol=1e-08,trace=FALSE)
#### input
n = 10
p = 3
K = c(3,3,4)
y = rnorm(n)
X.raw= z = matrix(rnorm(n*length(K)), n, length(K))
weights=NULL
max.steps = 30
lam.min=0
tol=1e-08
trace=FALSE

################ classo

