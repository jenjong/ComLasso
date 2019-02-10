rm(list = ls())
# load data
load("./bsi_new.RData"); 
# Implement stability Selection by fitting comlasso 
# It takes a few minutes to obtain the result
library(comlasso)
n <- 28
gs <- matrix(0,107,4)
set.seed(1)
gs[,1] <- stable_class_comlasso(x,y,gam=0,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
set.seed(1)
gs[,2] <- stable_class_comlasso(x,y,gam=-0.2,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
set.seed(1)
gs[,3] <- stable_class_comlasso(x2,y,gam=0,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 
set.seed(1)
gs[,4] <- stable_class_comlasso(x2,y,gam=-0.2,nbootstrap=1e3,nsteps=n/2,alpha=0.5,plotme=F,trace=F) 

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

# classification with marginal association
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
### Table 3  is generated
library(xtable)
xtable(cbind(lab.fx2[[2]][1:10],lab.fx2[[4]][1:10],lab.fx.ma2[1:10]))

