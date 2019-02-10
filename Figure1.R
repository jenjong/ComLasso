rm(list = ls())
library(comlasso)
data(sand)
mdat <- sand[,-1]
x <- as.matrix(mdat[,1:3])/100
x <- log(x); y <- log(mdat[,4])
ols <- lm((x[,1]-x[,3])~y)
ols <- lm((x[,1]-x[,2])~y)
n <- nrow(x); p <- ncol(x)
ltype="regression"

set.seed(1)
f1 <- comlasso("regression", n=n, p=1, K=p, X.raw=x, y=y, gam=1e5, weights=c(1,1,1))
w <- 1/(abs(f1$beta[[1]])/sum(abs(f1$beta[[1]])))
f2 <- comlasso("regression", n=n, p=1, K=p, X.raw=x, y=y, gam=1e5, weights=w)

op <- par(mfrow=c(1,2),mai=c(1.02,0.82-0.3,0.82,0.42),mar=c(5.1,4.1+0.7,4.1,2.1))
s = rowSums(abs(t(f1$beta.rec[[1]])))
matplot(s/s[3], t(f1$beta.rec[[1]]),type="b",xaxt="n",
  xlab=expression(paste(group("|", beta(lambda), "|")[1], "/" ,
    group("|", beta(0), "|")[1])),
  ylab=expression(hat(beta)[j](lambda)),ylim=c(-1,1)/2,main="lasso")#c(-16,16))
abline(h=0,lty=2,lwd=0.7)
abline(v=s[2]/s[3],lty=2,lwd=0.7)
legend(1,15,legend=c("comlasso","adaptive_comlasso"))
s = rowSums(abs(t(f2$beta.rec[[1]])))
matplot(c(s[1],s[3],s[2])/s[2],t(f2$beta.rec[[1]]),type="b",xaxt="n",
  xlab=expression(paste(group("|", beta(lambda), "||")[1], "/" ,
    group("||", beta(0), "|")[1])),
  ylab=expression(hat(beta)[j](lambda)),ylim=c(-1,1)/2,main="adaptive lasso")#c(-16,16))
abline(h=0,lty=2,lwd=0.7)
abline(v=s[3]/s[2],lty=2,lwd=0.7)
par(op)



