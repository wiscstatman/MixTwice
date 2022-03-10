l = 3000 ## number of testing units

m = 4 ## number of groups

neach = 20 ## number of subjects in each group

group = rep(factor(1:m), each = neach)

n = m*neach ## total number

pi0 = 0.8

signal2 = rep(0, l)

signal2[1:round((1-pi0)*l)] = rnorm(round((1-pi0)*l), mean = 0, sd = 2)

signal = cbind(rep(0, l), signal2)

## I will calculate the lambda parameter (true signal, H0: lambda = 0)

lambda = apply(signal, 1, function(x){return(sum((x - mean(x))^2)*neach)})

## I will generate the sigma^2 parameter

sigma2 = rep(0.5,l)

sigma = sqrt(sigma2)

## Then I can generate data, and calculate SSE and SSB

data = matrix(NA, nrow = l, ncol = n)

for (i in 1:l) {
  
  data[i,] = rnorm(n, mean = rep(signal[i,], each = neach), sd = sigma[i])
  
}

SSB = apply(data, 1, 
            function(y){
              fit = aov(y~group)
              fit.summary = summary(fit)[[1]]
              return(fit.summary$`Sum Sq`[1])
            })

SSE = apply(data, 1, 
            function(y){
              fit = aov(y~group)
              fit.summary = summary(fit)[[1]]
              return(fit.summary$`Sum Sq`[2])
            })

################################################################################
################################################################################

Blambda = 20

Bsigma2 = 6

df1 = m-1
df2 = n-m

kk = SSB*(n-m-2)/SSE - m + 1
kk2 = SSB - (m-1)*SSE/(n-m)

grid.lambda = seq(0, sqrt(max(lambda)), length = Blambda)^2

grid.sigma2 = seq(min(SSE/(n-m))/1.1, 1.1*max(SSE/(n-m)), length = Bsigma2)

## conditional likelihood

likelihood.SSE = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
likelihood.SSB = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)

grid.lambda2 = rep(grid.lambda, each = Bsigma2)
grid.sigma22 = rep(grid.sigma2, Blambda)

t1 = rep(c(1:length(grid.lambda)), each = Bsigma2)
t2 = rep(c(1:length(grid.sigma2)), Blambda)

for(i in 1:l){
  
  likelihood.SSE[i,] = dchisq(SSE[i]/grid.sigma22, df = n-m)
  
  likelihood.SSB[i,] = dchisq(SSB[i]/grid.sigma22, df = m-1, ncp = grid.lambda2/grid.sigma22)
  
}

likelihood.conditional = likelihood.SSB * likelihood.SSE

loglik = log(likelihood.conditional)

L = function(x) {
  xlambda = x[1:Blambda]
  xsigma2 = x[(Blambda + 1):(Blambda + Bsigma2)]
  yy = array(x[(Blambda + 1):(Blambda + Bsigma2)] %o% x[1:Blambda])
  return(-sum(log(yy %*% t(likelihood.conditional))))
}

## point mass estimation of lambda

est.lambda = rep(1, length(grid.lambda)); est.lambda = est.lambda/sum(est.lambda)
est.lambda = abs(rnorm(length(grid.lambda))); est.lambda = est.lambda/sum(est.lambda)

## point mass estimation of sigma2

est.sigma2 = rep(1, length(grid.sigma2)); est.sigma2 = est.sigma2/sum(est.sigma2)

count = 0
maxit = 100
notdone = TRUE

LL = NULL

while (notdone) {
  
  ll = L(c(est.lambda, est.sigma2))
  LL = c(LL, ll)
  
  est = as.numeric(array(est.sigma2 %o% est.lambda))
  vv = t(t(loglik) + log(est))
  uu = exp(vv)
  vv = uu/rowSums(uu)
  est = colMeans(vv)
  
  ## de-compose est into est.lambda and est.sigma2
  
  est.lambda = as.numeric(tapply(est, t1, sum))
  est.sigma2 = as.numeric(tapply(est, t2, sum))
  
  ## use Pava to get est.lambda monotone decreasing
  
  tmp = pava(est.lambda, decreasing = T)
  
  est.lambda = tmp/sum(tmp)
  
  notdone = (count < maxit )
  count = count + 1
  
  
}

### Optimization goal

## first order derivative of objective function

G = function(x) {
  g = h = NULL
  
  xlambda = x[1:Blambda]
  xsigma2 = x[(Blambda + 1):(Blambda + Bsigma2)]
  
  yy = array(x[(Blambda + 1):(Blambda + Bsigma2)] %o% x[1:Blambda])
  
  d = yy %*% t(likelihood.conditional)
  
  for (i in (1:Blambda)) {
    g[i] = -sum((xsigma2 %*% t(likelihood.conditional[, c((Bsigma2 * (i - 
                                                                        1) + 1):(Bsigma2 * i))]))/d)
  }
  
  for (j in (1:Bsigma2)) {
    h[j] = -sum((xlambda %*% t(likelihood.conditional[, seq(j, (j + (Blambda - 
                                                                       1) * (Bsigma2)), by = Bsigma2)]))/d)
  }
  return(c(g, h))
}

## equality constraint, summation is 1

heq = function(x) {
  h = NULL
  h[1] = sum(x[1:Blambda]) - 1
  h[2] = sum(x[(Blambda + 1):(Blambda + Bsigma2)]) - 1
  return(h)
}

heq.jac.fun = function(x) {
  hh1 = c(rep(1, Blambda), rep(0, Bsigma2))
  hh2 = c(rep(0, Blambda), rep(1, Bsigma2))
  heq.jac = rbind(hh1, hh2)
  return(heq.jac)
}

## inequality constraint, all between 0 and 1 (since we already have sum = 1, it is equivalent only saying all larger than 0)

hin = function(x) {
  h = NULL
  for (i in 1:((Blambda) + (Bsigma2))) {
    h[i] = x[i]
  }
  
  h2 = NULL
  for(i in 1:(Blambda - 1)){
    h2[i] = x[i] - x[i+1]
  }
  return(c(h,h2))
}

hin.jac1 = diag(1, nrow = (Blambda + Bsigma2), ncol = (Blambda + Bsigma2))

hin.jac2 = matrix(0, nrow = Blambda - 1, ncol = Blambda)

for (i in 1:(Blambda - 1)) {
  hin.jac2[i, i] = 1
  hin.jac2[i, i + 1] = -1
}
hin.jac3 = matrix(0, nrow = Blambda - 1, ncol = Bsigma2)

hin.jac.fun = function(x) {
  
  return(rbind(hin.jac1, cbind(hin.jac2, hin.jac3)))
}

## initial guess

a1 = seq(10, 1, length = Blambda)
a1 = a1/sum(a1)
a2 = rep(1, Bsigma2)
a2 = a2/sum(a2)
a = c(a1, a2)

try1 = suppressWarnings(alabama::auglag(par = a, fn = L, 
                                        gr = G, heq = heq, hin = hin, heq.jac = heq.jac.fun, 
                                        hin.jac = hin.jac.fun,control.outer = list(trace = F)))


est.lambda.opt = try1$par[1:Blambda]

est.lambda.opt[est.lambda.opt<0] = 0

est.sigma2.opt = try1$par[(Blambda + 1):(Blambda + Bsigma2)]

est.sigma2.opt[est.sigma2.opt<0] = 0

save.image(file = "./check for loglik convergence.RData")

