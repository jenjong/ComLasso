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

###### Figure 3 (left panel) ######
fit2 <- comlasso(ltype="classification",n,p=1,K=107,max.steps=100,
                 x,y,gam=0,weights=NULL)
coeffs <- predict.comlasso(fit2,newx=x2,s=fit2$lambda,type="coefficients",mode="lambda")$coefficients
b <- sum(abs(coeffs[nrow(coeffs),]))
a <- rowSums(abs(coeffs))
matplot(a/b,t(fit2$beta.rec[[1]][1:20,]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[1])
matplot(a/b,t(fit2$beta.rec[[1]][21:30,]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[2])
matplot(a/b,t(fit2$beta.rec[[1]][31:92,]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[3])
matplot(a/b,t(fit2$beta.rec[[1]][93:107,]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[4])



###### Figure 3 (right panel) ######
fit3 <- comlasso(ltype="classification",n,p=length(K),K=K,
                 max.steps=100,x2,y,gam=0,weights=NULL)
# solution path  
coeffs <- predict.comlasso(fit3,newx=x2,s=fit3$lambda,type="coefficients",mode="lambda")$coefficients
b <- sum(abs(coeffs[71,]))
a <- rowSums(abs(coeffs))
matplot(a/b,t(fit3$beta.rec[[1]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[1])
matplot(a/b,t(fit3$beta.rec[[2]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[2])
matplot(a/b,t(fit3$beta.rec[[3]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[3])
matplot(a/b,t(fit3$beta.rec[[4]]),type="l",ylim=c(-0.5,0.5),ylab="coeffients",cex.axis=0.8, cex.lab=0.9,
        xlab=expression(sum(hat(beta)[j](lambda),j=1,107)/sum(abs(hat(beta)[j](0)),j==1,107)), lwd=3, main=uphylum[4])





