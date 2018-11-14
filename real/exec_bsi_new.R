#############################################
#############################################
#setwd("~/comlasso/")
# load data
#load("./bsi_new.RData"); 
#n <- 28
#gs <- matrix(0,107,4)
#set.seed(1)
#gs[,1] <- stable_class_comlasso(x,y,gam=0,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
#set.seed(1)
#gs[,2] <- stable_class_comlasso(x,y,gam=-0.2,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
#set.seed(1)
#gs[,3] <- stable_class_comlasso(x2,y,gam=0,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
#set.seed(1)
#gs[,4] <- stable_class_comlasso(x2,y,gam=-0.2,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
#save.image(file="bsi_new_result.RData")


load("bsi_new_result.RData")
library(comlasso)
lab.fx <- vector(mode="list", 4)
lab.fx2 <- vector(mode="list", 4)
gname <- colnames(x)
lvs <- rep(0, length(gname))
for(j in 1:4){
  sfx <- sort(gs[,j], decreasing=T, index.return=T)
  sgs <- gs[sfx$ix,j]
  sgname <- gname[sfx$ix]
  num <- match(sgname,gname)
  pos <- 1:10
  lab.fx[[j]] <- rep(0, length(pos))

  for(i in 1:length(pos)){
    lab.fx[[j]][i] <- paste(i, ". ",sgname[pos[i]]," (",round(sgs[pos[i]],3),")",sep="")
    if(num[i]>=1 & num[i] <=20){
      lvs[i] <- 1  
    }else if(num[i]>=21 & num[i] <=30){
      lvs[i] <- 2
    }else if(num[i]>=31 & num[i] <=92){
      lvs[i] <- 3
    }else if(num[i]>=93 & num[i] <=107){
      lvs[i] <- 4
    }    
    lab.fx2[[j]][i] <- paste(uphylum[lvs[i]], ". ",sgname[pos[i]]," (",round(sgs[pos[i]],3),")",sep="")
  }
}


# barplot
op <- par(mar=c(5.1,4.1,4.1,2.1)-c(3.5,1.1,2.1,0.8))
plot(sgs[1:30], ylim=c(0,1), type="h", lwd=4, xaxt="n")#axis(1, at=pos)
abline(h=0.21, lty=2, lwd=2)
legend(horiz=F,2,0.95,legend=lab.fx, ncol=2, cex=1.20, bty="n")

# 1. classification with marginal association
wt <- function(x,y){wilcox.test(x[y==1], x[y==-1])$p.value}
pval <- apply(x[,gname],2,wt,y)
  sfx <- sort(pval, decreasing=F, index.return=T)
  sgs <- pval[sfx$ix]
  sgname <- gname[sfx$ix]
  num <- match(sgname, gname)
  pos <- 1:10
  lab.fx.ma <- rep(0, length(pos))
  lab.fx.ma2 <- rep(0, length(pos))
  for(i in 1:length(pos)){
    lab.fx.ma[i] <- paste(i, ". ",sgname[pos[i]]," (",round(sgs[pos[i]],3),")",sep="")
    if(num[i]>=1 & num[i] <=20){
      lvs[i] <- 1  
    }else if(num[i]>=21 & num[i] <=30){
      lvs[i] <- 2
    }else if(num[i]>=31 & num[i] <=92){
      lvs[i] <- 3
    }else if(num[i]>=93 & num[i] <=107){
      lvs[i] <- 4
    }    
    lab.fx.ma2[i] <- paste(uphylum[lvs[i]], ". ",sgname[pos[i]]," (",round(sgs[pos[i]],3),")",sep="")
  }
  library(xtable); xtable(cbind(lab.fx2[[2]][1:10],lab.fx2[[4]][1:10],lab.fx.ma2[1:10]))
  library(xtable); xtable(cbind(lab.fx[[2]][1:13],lab.fx[[4]][1:13],lab.fx.ma[1:13]))

###
tmp <- data.frame(ugenera,pval)
aa <- tmp[sort(tmp[,2],decreasing=F,index=T)$ix,][1:10,]
names(ugenera) <- ugenera
selvar <- match(as.character(aa[,1]), ugenera)
xtable(aa,digits=3)

####################################################
##1. gamma=0 & marginal association으로 선택된 변수 10개   #
####################################################
fit1 <- comlasso(ltype="classification", n, p=1, K=length(selvar), x[,selvar], y, gam=0, weights=NULL)
#matplot(t(fit$beta.rec[[1]]),type="b",lty=2)
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
#par(mar=c(5,4,4,3)-c(1,-0.2,1.5,0.5),mai=c(0.8,0.8,0.5,0.5)-c(0.01,0.3,0.3,0.3))
par(mar=c(5,4,4,3)-c(3.8,-0.2,3,1.8))
plot(apply(mgs1[[1]],2,function(x){mean(x<0)}),
     col=1,lty=1,type="s",lwd=2,xaxt="n",
     ylab="leave-one-out error rate",
     xlab=expression(lambda),ylim=c(0.10,0.55),
     cex=1.5,cex.axis=1.2)
points(apply(mgs1[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(apply(mgs1[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(20,0.55,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)

####################################
# 2. classification with complasso #
####################################
# comlassoA
# fit2는 107개 통으로 적합
fit2 <- comlasso(ltype="classification",n,p=1,K=107,max.steps=100,
                 x,y,gam=0,weights=NULL)
par(mar=c(5.1, 4.1, 4.1, 2.1)-c(0,0,2,-1.5))
matplot.comlasso1(fit2$lambda,t(fit2$beta.rec[[1]][sort(gs[,1],decreasing=T,index=T)$ix,]), 
  1:10, ylab="coeffients",cex.axis=0.8, cex.lab=0.7,
  xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, breaks=FALSE)

###
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
par(mar=c(5,4,4,3)-c(3.8,-0.2,3,1.8))
plot(apply(mgs2[[1]],2,function(x){mean(x<0)}),col=1,lty=1,type="s",lwd=2,xaxt="n",
     ylab="leave-one-out error rate",
     xlab=expression(lambda),ylim=c(0.10,0.55),cex=1.5,cex.axis=1.2)
points(apply(mgs2[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(apply(mgs2[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(20,0.55,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)

#######################################################
### comlassoB
### fit3는 4개의 대 그룹으로 묶어서 적합
#######################################################
fit3 <- comlasso(ltype="classification",n,p=length(K),K=K,max.steps=100,
  x2,y,gam=0,weights=NULL)
  
# 해의 solution path  
coeffs <- predict.comlasso(fit3,newx=x2,s=fit3$lambda,type="coefficients",mode="lambda")$coefficients
b <- sum(abs(coeffs[71,]))
a <- rowSums(abs(coeffs))

par(mfrow=c(2,2),mar=c(5.1, 4.1, 4.1, 2.1)-c(0,2,2,1.5))
matplot(a/b,t(fit3$beta.rec[[1]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[1])
matplot(a/b,t(fit3$beta.rec[[2]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[2])
matplot(a/b,t(fit3$beta.rec[[3]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[3])
matplot(a/b,t(fit3$beta.rec[[4]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[4])

can.lam <- exp(seq(log(max(fit3$lambda)),1e-5,length=50))
can.gam <- c(-0.2,-0.1,0)

mgs3 <- vector(mode="list",length=length(can.gam))
for(j in 1:length(can.gam))
  mgs3[[j]] <- matrix(0,28,length(can.lam))
for(j in 1:length(can.gam)){
  for(i in 1:28){
loo <- comlasso(ltype="classification",n-1,p=length(K),K=K,max.steps=1e2,x2[-i,],y[-i],gam=can.gam[j], weights=NULL)
    est <- predict.comlasso(loo,newx=x2[i,,drop=F],s=can.lam,type="fit",mode="lambda")$fit
    mgs3[[j]][i,] <- est*y[i]
  }
}
#par(mar=c(5,4,4,3)-c(1,-0.2,1.5,0.5),
#    mai=c(0.8,0.8,0.5,0.5)-c(0.01,0.3,0.3,0.3))
par(mar=c(5,4,4,3)-c(3.8,-0.2,3,1.8))
plot(apply(mgs3[[1]],2,function(x){mean(x<0)}),col=1,lty=1,type="s",lwd=2,xaxt="n",
     ylab="leave-one-out error rate",
     xlab=expression(lambda),ylim=c(0.10,0.55),cex=1.5,cex.axis=1.2)
points(apply(mgs3[[2]],2,function(x){mean(x<0)}),col=2,lty=2,type="s",lwd=2,cex=2)
points(apply(mgs3[[3]],2,function(x){mean(x<0)}),col=3,lty=3,type="s",lwd=2,cex=2)
legend(20,0.55,legend=c(expression(gamma==-0.2),expression(gamma==-0.1),expression(gamma==0)),col=c(1,2,3),lty=c(1,2,3),bty="n",ncol=1,cex=1.5)

#
#save.image(file="bsi_result_full.RData")








