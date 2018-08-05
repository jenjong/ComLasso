n = 100
M1 = matrix(1, n, 1)
M2 = cbind( c(rep(1,30), rep(0,70)),
            c(rep(0,30), rep(1,70)) )
M3 = cbind( c(rep(1,10), rep(0,20), rep(0,70)),
            c(rep(0,10), rep(1,20), rep(0,70)),
            c(rep(0,30), rep(1,30), rep(0,40)),
            c(rep(0,30), rep(0,30), rep(1,40)))
proj = function(x)
{
  x%*%solve(t(x)%*%x)%*%t(x)
}
P1 = proj(M1)
P2 = proj(M2)
P3 = proj(M3)
I = diag(1,n)
sum(abs((I-P1)))
sum(abs(P1-P2))
sum((abs(P2-P3)))
image(P3)

# simulation model

