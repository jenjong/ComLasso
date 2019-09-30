rm(list = ls())
gc()
#setwd("C:/Users/uos_stat/Dropbox/A composite lasso/prog/library")
setwd("/home/jeon/Dropbox/A composite lasso/prog/library")

tmp = as.matrix(read.csv('income.csv', header = T, stringsAsFactors = F)[,-1])
rtmp1 = as.matrix( tmp[2*(0:9)+1  , ] )
rtmp2 = as.matrix( tmp[2*(0:9)+2  , ] )
y =   colSums(rtmp2[10,,drop = F])/colSums(rtmp2[1,, drop = F])
plot(y)


tmp = read.csv('capital.csv', header = T, stringsAsFactors = F)
cname<-c('y',tmp[c(4,5,6,8,9,11,12),1])
tmp = as.matrix(tmp[c(4,5,6,8,9,11,12),-1])

cmean<-colSums(tmp)
for ( i in 1:ncol(tmp))
{
  tmp[,i] = tmp[,i]/cmean[i]
}

for ( i in 1:nrow(tmp))
{
  if (i == 1) plot(tmp[i,], type = 'l', col = i, ylim = c(0,0.5))  
  if (i != 1) lines(tmp[i,], col = i)
}

ctmp = tmp
x = t(ctmp)

x <- log(x)
z = x[,-ncol(x)]
for ( i in 1:ncol(z))
{
  z[,i] = z[,i]-x[,7]
}

ndata = data.frame(y,z)
names(ndata) = cname[-8]
fit = lm(y~., data = ndata)
plot(fit)
summary(fit)



##############
x = t(ctmp)
x <- log(x)
n <- nrow(x)
p <- ncol(x)
colnames(x) = cname[-1]
ltype="regression"

source("CLASSO_basic.R")
cfit1 <- HCLASSO(ltype="regression", n=104, p=1, K=7, X.raw=x, y=y, gam=100, weights=NULL, trace=FALSE)
matplot(log(cfit1$lambda), t(cfit1$beta.rec[[1]]), type="b")  

cfit2 <- HCLASSO(ltype="regression", n=104, p=1, K=7, X.raw=x, y=y, gam=2, weights=NULL, trace=FALSE)
matplot(log(cfit2$lambda), t(cfit2$beta.rec[[1]]), type="b")  

par(mfrow=c(1,2),mai=c(1.02,0.82,0.82,0.42)+c(0,0.15,-0.1,-0.2))
matplot(cfit1$lambda, t(cfit1$beta.rec[[1]]), 
       type="b", xlab=expression(lambda),ylab=expression(hat(beta)[j](lambda)))  
matplot(cfit2$lambda, t(cfit2$beta.rec[[1]]), 
       type="b", xlab=expression(lambda),ylab=expression(hat(beta)[j](lambda)))  

## adaptive
source("CLASSO_basic.R")
cfit3 <- HCLASSO(ltype="regression", n=104, p=1, K=7, X.raw=x, y=y, gam=100, 
                 weights=1/(abs(cfit1$beta[[1]])/sum(abs(cfit1$beta[[1]]))), trace=FALSE)

cfit4 <- HCLASSO(ltype="regression", n=104, p=1, K=7, X.raw=x, y=y, gam=2, 
                 weights=1/(abs(cfit1$beta[[1]])/sum(abs(cfit1$beta[[1]]))), trace=FALSE)

par(mfrow=c(2,2),mai=c(1.02,0.82,0.82,0.42)+c(0,0.15,-0.1,-0.2), 
mar=c(6.144578, 5.843373, 4.337349, 1.325301)+c(-2,0,-3,0))
graph_path(cfit1); graph_path(cfit2)
graph_path(cfit3); graph_path(cfit4)

