#############################################
#############################################
rm(list = ls())
setwd("./GitHub/ComLasso/real")
load("bsi_result_full.RData")
library(comlasso)
# 
fit = fit2
coeffs <- predict.comlasso(fit,newx=x2,s=fit$lambda,type="coefficients",mode="lambda")$coefficients
tmp = fit$beta.rec[[1]];dim(tmp)
cumsum(K)
idx1 = 1:20;idx2 = 21:30;idx3 = 31:92;idx4 = 93:107
tmplist = list()
tmplist[[1]] <- tmp[idx1,]
tmplist[[2]] <- tmp[idx2,]
tmplist[[3]] <- tmp[idx3,]
tmplist[[4]] <- tmp[idx4,]
fit$beta.rec <- tmplist
b <- sum(abs(coeffs[nrow(coeffs),]))
a <- rowSums(abs(coeffs))
frac_lambda =  a/b
for (knum in 1:4)
{
  png(file = paste0("comlassoA_",knum,'.png'), width = 720, height = 480)
  bmat = t(fit$beta.rec[[knum]])
  matplot(frac_lambda, bmat ,type="l", ylim=c(-0.5,0.5), ylab="coefficient",
          cex.main = 1.8, cex.axis=1.8, cex.lab=1.8,
          xlab= expression(paste(group("|", beta(lambda), "|")[1], "/" ,
                                 group("|", beta(lambda[max]), "|")[1])),
          lwd=3, main=uphylum[knum])
  
  dev.off()
}

# comlassoB
fit = fit3
coeffs <- predict.comlasso(fit,newx=x2,s=fit$lambda,type="coefficients",mode="lambda")$coefficients
b <- sum(abs(coeffs[nrow(coeffs),]))
a <- rowSums(abs(coeffs))
frac_lambda =  a/b
for (knum in 1:4)
{
  bmat = t(fit$beta.rec[[knum]])
  png(file = paste0("comlassoB_",knum,'.png'), width = 720, height = 480)
  bmat = t(fit$beta.rec[[knum]])
  matplot(frac_lambda, bmat ,type="l", ylim=c(-0.5,0.5), ylab="coefficient",
          cex.main = 1.8, cex.axis=1.8, cex.lab=1.8,
          xlab= expression(paste(group("|", beta(lambda), "|")[1], "/" ,
                                 group("|", beta(lambda[max]), "|")[1])),
          lwd=3, main=uphylum[knum])
  
  dev.off()
  
}

