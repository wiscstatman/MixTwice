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

## point mass estimation of lambda

est.lambda = rep(1, length(grid.lambda)); est.lambda = est.lambda/sum(est.lambda)
est.lambda = abs(rnorm(length(grid.lambda))); est.lambda = est.lambda/sum(est.lambda)

## point mass estimation of sigma2

est.sigma2 = rep(1, length(grid.sigma2)); est.sigma2 = est.sigma2/sum(est.sigma2)

count = 0
maxit = 50
notdone = TRUE

DD = NULL

while (notdone) {
  
  ee.old = c(est.lambda, est.sigma2)
  
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
  
  ee.new = c(est.lambda, est.sigma2)
  dd = sum((ee.new - ee.old)^2)
  DD = c(DD, dd)
  
  notdone = (count < maxit )
  count = count + 1
  
  
}

save.image(file = "./check for error convergence.RData")

