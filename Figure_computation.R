rm(list=ls())
library(ggplot2)
setwd("~/github/comlasso")
load("table1.rdata")
v = matrix(0,3,6)
for(i in 1:length(runtime.list))
{
  for (j in 1:3)  v[j,i]<- median(runtime.list[[i]][,j], na.rm = T)
}
v
fv = factor(paste0("p=", c(20,50,200,500,1000,2000)), 
            levels = paste0("p=", c(20,50,200,500,1000,2000)))
       
       
vv = data.frame(method = rep(c("comlasso", "genlasso", "constrained lasso"), 6),
p = rep(fv, each = 3), time = c(v))
ggplot(data = vv) + geom_bar(stat='identity', position = 'dodge',
                             aes(x = p, y = time, fill = method))

rm(list=ls())
setwd("~/github/comlasso")
load("table1-2.rdata")
v = matrix(0,3,6)
for(i in 1:length(runtime.list))
{
  for (j in 1:3)  v[j,i]<- median(runtime.list[[i]][,j], na.rm = T)
}
v
fv = factor(paste0("n=", c(20,50,100,200,500,1000)), 
            levels = paste0("n=", c(20,50,100,200,500,1000)))


vv = data.frame(method = rep(c("comlasso", "genlasso", "constrained lasso"), 6),
                n = rep(fv, each = 3), time = c(v))
ggplot(data = vv) + geom_bar(stat='identity', position = 'dodge',
                             aes(x = n, y = time, fill = method))




# K = 5
rm(list=ls())
setwd("~/github/comlasso")
load("table2.rdata")
v = matrix(0,3,6)
for(i in 1:length(runtime.list))
{
  for (j in 1:3)  v[j,i]<- median(runtime.list[[i]][,j], na.rm = T)
}
v = v[,-2]

fv = factor(paste0("K=", c(1,5,10,25,50)), 
            levels = paste0("K=", c(1,5,10,25,50)))


vv = data.frame(method = rep(c("comlasso", "genlasso", "constrained lasso"), 5),
                K = rep(fv, each = 3), time = c(v))
ggplot(data = vv) + geom_bar(stat='identity', position = 'dodge',
                             aes(x = K, y = time, fill = method))


