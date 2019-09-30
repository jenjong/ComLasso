#
slump_test <- read.csv('slump_flow.csv')

mdat <- slump_test[,-1]
head(mdat)
x <- mdat[,1:7]
x[x==0] = 0.1
for ( i in 1:nrow(x))
{
  x[i,] = x[i,]/sum(x[i,])
}
x <- log(x)
y <- mdat[,8]

z = x[,-7]
for ( i in 1:ncol(z))
{
  z[,i] = z[,i]-x[,7]
}
ndata = data.frame(y,z)
fit = lm(y~., data = ndata)
plot(fit)

#
slump_test <- read.csv('slump_flow.csv')

mdat <- slump_test[,-1]
head(mdat)
x <- mdat[,1:7]
x <- x[,-c(1,3)]
x[x==0] = 1
for ( i in 1:nrow(x))
{
  x[i,] = x[i,]/sum(x[i,])
}
x <- log(x)
y <- mdat[,8]

z = x[,-ncol(x)]
for ( i in 1:ncol(z))
{
  z[,i] = z[,i]-x[,ncol(x)]
}
ndata = data.frame(y,z)
fit = lm(y~., data = ndata)
plot(fit)
summary(fit)
