l = 3000 ## number of testing units

m = 5 ## number of groups

neach = 4 ## number of subjects in each group

group = rep(factor(1:m), each = neach)

n = m*neach ## total number

simu = function(pi0, ## number of testings on null
                t, ## signal
                s ## error
){
  
  signal = matrix(0, nrow = l, ncol = m) ## mean parameters of each group for each unit
  
  ## generate a covariate
  
  covar = covar = 1 / (1+exp(-runif(l, 0, 10) + 5))
  
  covar_non = runif(l)
  
  ss = sample(l, l - round(l*pi0), prob = covar)
  
  for (i in ss) {
    
    signal[i,] = seq(-t, t, length = m)
    
  }
  
  truth = rep(0, l)
  truth[ss] = 1
  
  ## I will calculate the lambda parameter (true signal, H0: lambda = 0)
  
  lambda = apply(signal, 1, function(x){return(sum((x - mean(x))^2)*neach)})
  
  ## I will generate the sigma^2 parameter
  
  sigma2 = rep(s,l)
  
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
  
  Fstat = (SSB/(m-1))/(SSE/(n-m))
  
  pval = 1-pf(Fstat, df1 = m-1, df2 = n-m)
  
  ## bonferonni
  
  p_bonferroni = p.adjust(pval, method = "bonferroni")
  
  ## BH
  
  p_BH = p.adjust(pval, method = "BH")
  
  ## qvalue
  
  p_qvalue = qvalue(pval)$qvalue
  
  ## IHW
  
  library(IHW)
  
  p_ihw_s = adj_pvalues(ihw(pval, covariates = covar, alpha = 0.05))
  
  p_ihw_w = adj_pvalues(ihw(pval, covariates = covar_non, alpha = 0.05))
  
  ##  BL
  
  p_BL_s = lm_pi0(pval, X = covar)$pi0 * p_ihw_s
  
  p_BL_w = lm_pi0(pval, X = covar_non)$pi0 * p_ihw_w
  
  ## AdaPT (GLM)
  
  xx = adapt_glm(x = data.frame(icov = covar), 
                 pvals = pval, 
                 pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
                 mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)
  
  p_AdaPT_s = xx$qvals
  
  xx = adapt_glm(x = data.frame(icov = covar_non), 
                 pvals = pval, 
                 pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
                 mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)
  
  p_AdaPT_w = xx$qvals
  
  ## LFDR
  
  LFDR_s = clfdr_hickswrapper(unadj_p = pval,
                     groups = IHW::groups_by_filter(covar, 20), 
                     lfdr_estimation="fdrtool")
  
  LFDR_w = clfdr_hickswrapper(unadj_p = pval,
                              groups = IHW::groups_by_filter(covar_non, 20), 
                              lfdr_estimation="fdrtool")
  
  ## let's try mixtwice approach
  
  Blambda = 20
  
  Bsigma2 = 8
  
  kk = SSB*(n-m-2)/SSE - m + 1
  kk2 = SSB - (m-1)*SSE/(n-m)
  
  grid.lambda = seq(0, max(lambda), length = Blambda)
  
  grid.sigma2 = seq(min(SSE/(n-m))/1.1, 1.1*max(SSE/(n-m)), length = Bsigma2)
  
  ## conditional likelihood
  
  likelihood.SSE = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  likelihood.SSB = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  
  grid.lambda2 = rep(grid.lambda, each = Bsigma2)
  grid.sigma22 = rep(grid.sigma2, Blambda)
  
  for(i in 1:l){
    
    likelihood.SSE[i,] = dchisq(SSE[i]/grid.sigma22, df = n-m)
    
    likelihood.SSB[i,] = dchisq(SSB[i]/grid.sigma22, df = m-1, ncp = grid.lambda2/grid.sigma22)
    
  }
  
  likelihood.conditional = likelihood.SSB * likelihood.SSE
  
  ## objective function
  
  L = function(x) {
    xlambda = x[1:Blambda]
    xsigma2 = x[(Blambda + 1):(Blambda + Bsigma2)]
    yy = array(x[(Blambda + 1):(Blambda + Bsigma2)] %o% x[1:Blambda])
    return(-sum(log(yy %*% t(likelihood.conditional))))
  }
  
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
    return(h)
  }
  hin.jac.fun = function(x) {
    hin.jac = diag(1, nrow = (Blambda + Bsigma2), 
                   ncol = (Blambda + Bsigma2))
    return(hin.jac)
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
  
  est.lambda = try1$par[1:Blambda]
  est.lambda[est.lambda < 0] = 0
  
  est.sigma2 = try1$par[(Blambda + 1):(Blambda + Bsigma2)]
  est.sigma2[est.sigma2 < 0] = 0
  
  #plot(grid.lambda, cumsum(est.lambda), ylim = c(0, 1), type = "s")
  #lines(ecdf(lambda))
  
  #plot(grid.sigma2, cumsum(est.sigma2), ylim = c(0, 1), type = "s")
  #lines(ecdf(sigma2))
  
  est.matrix = outer(est.lambda, est.sigma2)
  
  est.array = NULL
  for (i in 1:Blambda) {
    est.array = c(est.array, est.matrix[i, ])
  }
  
  LFDR = matrix(NA, ncol = Blambda, nrow = l)
  
  for (i in 1:l) {
    
    ddd = likelihood.conditional[i, ] * est.array
    ddd = ddd/sum(ddd)
    UUU = NULL
    
    for (j in 1:Blambda) {
      begin = (j - 1) * Bsigma2 + 1
      end = j * Bsigma2
      uuu = sum(ddd[begin:end])
      UUU = c(UUU, uuu)
    }
    
    LFDR[i, ] = UUU
  }
  
  
  lfdr = LFDR[, 1]
  
  summary.table = data.frame(p_bonferroni, p_BH, p_qvalue,
                             p_ihw_s, p_ihw_w,
                             p_BL_s, p_BL_w,
                             p_AdaPT_s, p_AdaPT_w,
                             LFDR_s, LFDR_w,
                             mixtwice = lfdr)
  
  return(list(summary.table = summary.table,
              truth = truth))
  
}


pi0_true = seq(0.1, 0.9, length = 30)

tt = c(1, 2, 3)

ss = c(1, 2, 3)

i = 1

out = NULL

for(t in tt){
  
  for (s in ss) {
    
    for (pi0 in pi0_true) {
      
      xx = simu(pi0, t, s)
      
      out[[i]] = xx
      
      i = i + 1
      
      print(i)
      
    }
    
  }
}
