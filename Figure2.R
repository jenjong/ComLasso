rm(list = ls())
# load data
load("./bsi_new.RData"); 
library(comlasso)
n <- 28
gs <- matrix(0,107,4)
gname <- colnames(x)
wt <- function(x,y){wilcox.test(x[y==1], x[y==-1])$p.value}
pval <- apply(x[,gname],2,wt,y)
tmp <- data.frame(ugenera,pval)
aa <- tmp[sort(tmp[,2],decreasing=F,index=T)$ix,][1:10,]
names(ugenera) <- ugenera
selvar <- match(as.character(aa[,1]), ugenera)


###### Figure 2-(a) ######
fit2 <- comlasso(ltype="classification",n,p=1,K=107,max.steps=100,
                 x,y,gam=0,weights=NULL)
par(mar=c(5.1, 4.1, 4.1, 2.1)-c(0,0,2,-1.5))
can.lam <- exp(seq(log(max(fit2$lambda)),1e-5,length=50))
can.gam <- c(-0.2,-0.1,0)

mgs2 <- vector(mode="list",length=length(can.gam))
for(j in 1:length(can.gam))
  mgs2[[j]] <- matrix(0,28,length(can.lam))

for(j in 1:length(can.gam)){
  for(i in 1:28){
    loo <- comlasso(ltype="classification",n-1,p=1,K=107,max.steps=1e2,x[-i,],y[-i],gam=can.gam[j], weights=NULL)
    est <- predict.comlasso(loo,newx=x[i,,drop=F],s=can.lam,type="fit",mode="lambda")$fit
    mgs2[[j]][i,] <- est*y[i]
  }
}
par(mar=c(5,4,4,3)-c(0.8,-0.2,3,1.8))
plot(can.lam,apply(mgs2[[1]],2,function(x){mean(x<0)}),
  col=1,lty=1,type="s",lwd=2,
  ylab="leave-one-out error rate",
  xlab=expression(lambda),ylim=c(0.10,0.55),
  cex=1.5,cex.axis=1.2,cex.lab=1.5)
points(can.lam, apply(mgs2[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(can.lam, apply(mgs2[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(3,0.2,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)

###### Figure 2-(b) ######
fit3 <- comlasso(ltype="classification",n,p=length(K),K=K,
                 max.steps=100,x2,y,gam=0,weights=NULL)
coeffs <- predict.comlasso(fit3,newx=x2,s=fit3$lambda,type="coefficients",mode="lambda")$coefficients
b <- sum(abs(coeffs[71,]))
a <- rowSums(abs(coeffs))
can.lam <- exp(seq(log(max(fit3$lambda)),1e-5,length=50))
can.gam <- c(-0.2,-0.1,0)

mgs3 <- vector(mode="list",length=length(can.gam))
for(j in 1:length(can.gam))
  mgs3[[j]] <- matrix(0,28,length(can.lam))
for(j in 1:length(can.gam)){
  for(i in 1:28){
    loo <- comlasso(ltype="classification",n-1,p=length(K),K=K,
             max.steps=1e2,x2[-i,],y[-i],gam=can.gam[j], weights=NULL)
    print(i)
    est <- predict.comlasso(loo,newx=x2[i,,drop=F],s=can.lam,type="fit",mode="lambda")$fit
    mgs3[[j]][i,] <- est*y[i]
  }
}
par(mar=c(5,4,4,3)-c(0.8,-0.2,3,1.8))
plot(can.lam, apply(mgs3[[1]],2,function(x){mean(x<0)}),
     col=1,lty=1,type="s",lwd=2,
     ylab="leave-one-out error rate",
     xlab=expression(lambda),ylim=c(0.10,0.55),
     cex=1.5,cex.axis=1.2,cex.lab=1.5)
points(can.lam,apply(mgs3[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(can.lam,apply(mgs3[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(3,0.2,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)


###### Figure 2-(c) ######
fit1 <- comlasso(ltype="classification", n, p=1, K=length(selvar), x[,selvar], y, gam=0, weights=NULL)
can.lam <- exp(seq(log(max(fit1$lambda)),1e-5,length=50))
can.gam <- c(-0.2,-0.1,0)
mgs1 <- vector(mode="list",length=length(can.gam))
for(j in 1:length(can.gam))
  mgs1[[j]] <- matrix(0,28,length(can.lam))
for(j in 1:length(can.gam)){
  for(i in 1:28){
    loo <- comlasso(ltype="classification", n-1, p=1, K=length(selvar), x[-i,selvar], y[-i], gam=can.gam[j], weights=NULL)
    est <- predict.comlasso(loo,newx=x[i,selvar,drop=F],s=can.lam,type="fit",mode="lambda")$fit
    mgs1[[j]][i,] <- est*y[i]
  }
}

par(mar=c(5,4,4,3)-c(0.8,-0.2,3,1.8))
plot(can.lam,apply(mgs1[[1]],2,function(x){mean(x<0)}),
     col=1,lty=1,type="s",lwd=2,
     ylab="leave-one-out error rate",
     xlab=expression(lambda),ylim=c(0.10,0.55),
     cex=1.5,cex.axis=1.2,cex.lab=1.5)
points(can.lam, apply(mgs1[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(can.lam, apply(mgs1[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(3,0.2,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)








