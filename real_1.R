rm(list=ls())
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(devtools)
library(Rglpk)
library(MethylCapSig)
install_github("glmgen/genlasso")
library(genlasso)
setwd("~/github/ComLasso")
sourceCpp('./library/inner.cpp')
source("./library/ComLassoC.R")
# income
rX<- read.csv("./data/capital_season.csv")
ry<- read.csv("./data/income.csv")[,-1]
y = unlist(ry[5,]/ry[1,])
y = y[-length(y)] # 2003-1~2019-1
X = as.matrix(rX[,130:197]) # 2002-1~2018-4
X = t(X)
x = sweep(X,1,rowSums(X),"/")
idx = 1:nrow(x)
y
x1 = x[tail(idx, 65),]
x2 = x[tail(idx, 65)-1,]
x3 = x[tail(idx, 65)-2,]
x4 = x[tail(idx, 65)-3,]
xx = cbind(x1,x2,x3,x4)
xx = log(xx)
pk = rep(8,4)
fit = comLassoC(xx,y,pk=pk,lam_min=0,tol=1e-08,KKT_check=FALSE)
dim(fit$coefMat)
coefMat = fit$coefMat
est.var = sum((y-cbind(1,xx)%*%coefMat[nrow(coefMat),])^2)/(length(y)-1)

aic.vec = c()
i = 1
for (i in 1:nrow(coefMat))
{
  coefvec = coefMat[i,-1]
  df = 0
  for (j in 1:4)
  {
    subvec = coefvec[((j-1)*8+1):(j*8)]
    if (sum(subvec !=0)>0) df <- df + (sum(subvec !=0))-1
  }
  aic.vec[i]<-
  sum((y-cbind(1,xx)%*%coefMat[i,])^2)/est.var + log(length(y)) + 2*df
}
plot(aic.vec)
which.min(aic.vec)

## 17
k = 4
for ( i in 1:8)
{
  if (i == 1)
    plot(fit$coefMat[,i+8*(k-1)+1], ylim = c(-10,10), col = i, 
         type = 'l')
  if (i>1)
    lines(fit$coefMat[,i+8*(k-1)+1], ylim = c(-10,10),
          col = i)
  abline(v=17)
}  

plot(1:8,rep(0,8), col = 1:8)
xx = melt(x1)
colnames(xx)  = c("time", "index", "prop")
xx$index = as.factor(xx$index)
xx$time = as.integer(xx$time)
ggplot(data = xx) + geom_area(aes(x = time, y= prop, fill = index, group = index))


