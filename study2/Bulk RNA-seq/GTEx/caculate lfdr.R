load("E:/github_local/MixTwice/study2/Bulk RNA-seq/GTEx/GTEx_processed data.RData")

get_mixing = function(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df, prop = 1){
  
  theta0 = thetaHat
  s20 = s2
  
  ok.sample = sample(length(theta0), length(s20) * prop)
  theta = theta0[ok.sample]
  s2 = s20[ok.sample]
  p0 = length(theta0)
  p = length(theta)
  cc = max(abs(theta)) * 1.1
  gridtheta = (cc/Btheta) * seq(-Btheta, Btheta, by = 1)
  gridsigma = seq(sqrt(min(s2)), sqrt(max(s2)), by = (sqrt(max(s2)) - 
                                                        sqrt(min(s2)))/Bsigma2)
  ltheta = length(gridtheta)
  lsigma = length(gridsigma)
  grid1 = rep(gridtheta, each = length(gridsigma))
  grid2 = rep(gridsigma, length(gridtheta))
  rbind(grid1, grid2)
  lik1 = t(exp(-0.5 * (t((outer(theta, grid1, "-"))^2)/(grid2)^2))/(grid2 * 
                                                                      sqrt(2 * pi)))
  y = outer(df * s2, (1/gridsigma^2), "*")
  m = df/2
  lik2 = y^(m - 1) * exp(-0.5 * y)/((2^m) * (gamma(m)))
  lik22 = matrix(rep(lik2, ltheta), nrow = nrow(lik2))
  lik = lik1 * lik22
  L = function(x) {
    xtheta = x[1:ltheta]
    xsigma = x[(ltheta + 1):(ltheta + lsigma)]
    yy = array(x[(ltheta + 1):(ltheta + lsigma)] %o% x[1:ltheta])
    return(-sum(log(yy %*% t(lik))))
  }
  G = function(x) {
    g = h = NULL
    xtheta = x[1:ltheta]
    xsigma = x[(ltheta + 1):(ltheta + lsigma)]
    yy = array(x[(ltheta + 1):(ltheta + lsigma)] %o% x[1:ltheta])
    d = yy %*% t(lik)
    for (i in (1:ltheta)) {
      g[i] = -sum((xsigma %*% t(lik[, c((lsigma * (i - 
                                                     1) + 1):(lsigma * i))]))/d)
    }
    for (j in (1:lsigma)) {
      h[j] = -sum((xtheta %*% t(lik[, seq(j, (j + (ltheta - 
                                                     1) * (lsigma)), by = lsigma)]))/d)
    }
    return(c(g, h))
  }
  heq <- function(x) {
    h = NULL
    h[1] = sum(x[1:ltheta]) - 1
    h[2] = sum(x[(ltheta + 1):(lsigma + ltheta)]) - 1
    return(h)
  }
  hh1 = c(rep(1, ltheta), rep(0, lsigma))
  hh2 = c(rep(0, ltheta), rep(1, lsigma))
  heq.jac = rbind(hh1, hh2)
  heq.jac.fun = function(x) {
    j = heq.jac
    return(j)
  }
  hin <- function(x) {
    h1 = NULL
    for (i in 1:((ltheta) + (lsigma))) {
      h1[i] = x[i]
    }
    h2 = NULL
    for (i in 1:(Btheta)) {
      h2[i] = x[i + 1] - x[i]
    }
    for (i in (Btheta + 1):(ltheta - 1)) {
      h2[i] = x[i] - x[i + 1]
    }
    h = c(h1, h2)
    return(h)
  }
  hin.jac1 = diag(1, nrow = (ltheta + lsigma), ncol = (ltheta + 
                                                         lsigma))
  hin.jac2 = matrix(0, ncol = ltheta, nrow = ltheta - 1)
  for (i in 1:(Btheta)) {
    hin.jac2[i, i] = -1
    hin.jac2[i, i + 1] = 1
  }
  for (i in (Btheta + 1):(ltheta - 1)) {
    hin.jac2[i, i] = 1
    hin.jac2[i, i + 1] = -1
  }
  hin.jac3 = matrix(0, nrow = ltheta - 1, ncol = lsigma)
  hin.jac = rbind(hin.jac1, cbind(hin.jac2, hin.jac3))
  hin.jac.fun = function(x) {
    j = hin.jac
    return(j)
  }
  a1 = rep(1, ltheta)
  a1 = a1/sum(a1)
  a2 = rep(1, lsigma)
  a2 = a2/sum(a2)
  a = c(a1, a2)
  try1 = suppressWarnings(alabama::auglag(par = a, fn = L, 
                                          gr = G, heq = heq, hin = hin, heq.jac = heq.jac.fun, 
                                          hin.jac = hin.jac.fun, control.outer = list(trace = F)))
  
  est.theta = try1$par[1:ltheta]
  est.theta[est.theta < 0] = 0
  est.sigma = try1$par[(ltheta + 1):(ltheta + lsigma)]
  est.sigma[est.sigma < 0] = 0
  
  return(list(grid.theta = gridtheta, grid.sigma2 = gridsigma^2, 
              mix.theta = est.theta, mix.sigma2 = est.sigma))
  
}

mix.mixtwice = get_mixing(thetaHat = res$effect_size,
                          s2 = res$SE^2,
                          Btheta = 30,
                          Bsigma = 10,
                          df = 9, 
                          prop = 0.01)

plot(mix.mixtwice$grid.theta, mix.mixtwice$mix.theta)

mixtwice_lfdr = function(thetaHat, s2, df, mix.mixtwice){
  
  gridtheta = mix.mixtwice$grid.theta
  est.theta = mix.mixtwice$mix.theta
  
  gridsigma = sqrt(mix.mixtwice$grid.sigma2)
  est.sigma = mix.mixtwice$mix.sigma2
  
  est.matrix = outer(est.theta, est.sigma)
  
  ltheta = length(gridtheta)
  lsigma = length(gridsigma)
  
  est.array = NULL
  
  for (i in 1:ltheta) {
    est.array = c(est.array, est.matrix[i, ])
  }
  
  p0 = length(thetaHat)
  
  LFDR = matrix(NA, ncol = ltheta, nrow = p0)
  
  theta = thetaHat
  s2 = s2
  
  grid1 = rep(gridtheta, each = length(gridsigma))
  grid2 = rep(gridsigma, length(gridtheta))
  
  lik1 = t(exp(-0.5 * (t((outer(theta, grid1, "-"))^2)/(grid2)^2))/(grid2 * 
                                                                      sqrt(2 * pi)))
  y = outer(df * s2, (1/gridsigma^2), "*")
  m = df/2
  lik2 = y^(m - 1) * exp(-0.5 * y)/((2^m) * (gamma(m)))
  lik22 = matrix(rep(lik2, ltheta), nrow = nrow(lik2))
  lik1 = lik1 * lik22
  
  for (i in 1:p0) {
    ddd = lik1[i, ] * est.array
    UUU = NULL
    for (j in 1:ltheta) {
      begin = (j - 1) * lsigma + 1
      end = j * lsigma
      uuu = sum(ddd[begin:end])
      UUU = c(UUU, uuu)
    }
    UUU = UUU/sum(UUU)
    LFDR[i, ] = UUU
  }
  
  lfdr = LFDR[, ((length(gridtheta)-1)/2 + 1)]
  
  return(lfdr = lfdr)
  
}

lfdr.mt = mixtwice_lfdr(thetaHat = res$effect_size,
                        s2 = res$SE^2,
                        df = 9,
                        mix.mixtwice = mix.mixtwice)

save(mix.mixtwice, lfdr.mt, file = "F://Harddrive-Jul-17-2021//MixTwice follow up//RNAseq//mixtwice_result.RData")

## consider the rest of the methods.

## unadjust

unadjust = res$pval

## bonferonni

bonferroni = p.adjust(res$pval, method = "bonferroni")

## BH

BH = p.adjust(res$pval, method = "BH")

## storey's qvalue
library(qvalue)

q.value = qvalue(res$pval)$qvalue

## IHW
library(IHW)

ihw = adj_pvalues(ihw(res$pval, covariates = res$ind_covariate, alpha = 0.05))

## ash
library(ashr)

ash = get_lfsr(ash(betahat = res$effect_size, sebetahat = res$SE, df = 9))

## BL

library(swfdr)

BL = lm_pi0(res$pval, X = res$ind_covariate)$pi0 * BH

## AdaPT (GLM)

library(adaptMT)

xx = adapt_glm(x = data.frame(icov = res$ind_covariate), 
               pvals = res$pval, 
               pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
               mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)

AdaPT = xx$qvals

AdaPT[AdaPT>=1] = 1


## LFDR

clfdr_hickswrapper <- function(unadj_p, groups, lfdr_estimation="fdrtool") {
  
  # Exclude this method if there are fewer than 200 tests within a grouping 
  # (fdrtool is applied separately to each group, and throws a warning in such case)
  #if(min(table(groups)) < 200)
  #  stop("Not enough tests to apply this method. Require at least 200.")
  
  ## estimate local fdr within each stratum first
  lfdr_res <- lfdr_fit(unadj_p, groups, lfdr_estimation=lfdr_estimation)
  lfdrs <- lfdr_res$lfdr
  
  ## now use the rejection rule described in Cai's paper
  
  ## Remark:
  ## When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  ## we get monotonic adjusted p-values as a function of the p-values
  ## This is mainly needed for grenander based lfdrs, with most other
  ## lfdr estimation methods lfdr ties are not a problem usually
  
  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]
  return(adj_p)
}

lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){
  
  pvals_list <- split(unadj_p, group)
  
  if (lfdr_estimation == "fdrtool"){
    lfdr_fun <- function(pv) fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
    
  } else if (lfdr_estimation == "locfdr"){
    if (!requireNamespace("locfdr", quietly=TRUE)){
      stop("locfdr package required for this function to work.")
    }
    lfdr_fun <- function(pv) locfdr::locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
  } else {
    stop("This lfdr estimation method is not available.")
  }
  
  lfdr_list <- lapply(pvals_list, lfdr_fun)
  lfdrs <- unsplit(lfdr_list, group)
  
  fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
  fit_obj
}

LFDR = clfdr_hickswrapper(unadj_p = res$pval,
                          groups = IHW::groups_by_filter(res$ind_covariate, 20), 
                          lfdr_estimation="fdrtool")

lfdr.summary = data.frame(unadjust, bonferroni, BH, q.value, ihw, ash,
                          BL, LFDR, AdaPT)

save(lfdr.summary,
     file = "F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/lfdr_all.RData")


