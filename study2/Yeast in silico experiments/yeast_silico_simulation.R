setwd("F:/Harddrive-Jul-17-2021/MixTwice follow up")

X = 1

group = 'Snf2'

rseed = 1

nDE = 500

sampleSize = 10

signal = 3 ## the fold value that I will add to one conditions, for generating non-null cases

uninformativeCovariate = FALSE
pvalHists = FALSE
strongCovariate = TRUE

load("F:/Harddrive-Jul-17-2021/MixTwice follow up/yeast_fulldata_DESeqdata.RData")

dds_full <- dds_full[,colData(dds_full)$condition == group] ## make a subset group

# set random seed

set.seed(as.numeric(X)*as.numeric(rseed))

covar = 1 / (1+exp(-runif(nrow(dds_full), 0, 10) + 5))

# select a random subset of samples of size sampleSize each, with in WT or Snf2 group

dds_test = dds_full[,sample(1:ncol(dds_full), sampleSize*2)]

# add a fake condition column to coldat
colData(dds_test)$fake =  factor(c(rep("A", sampleSize), 
                                   rep("B", sampleSize))[sample(1:(sampleSize*2), 
                                                                sampleSize*2)])
design(dds_test) = ~fake

pzero = rowSums(counts(dds_test)==0)/ncol(counts(dds_test))
dds_test = dds_test[pzero < 0.5,]
covar = covar[pzero < 0.5]

truth <- rep(FALSE, nrow(dds_test))

DE <- sample(1:nrow(dds_test), nDE, prob = covar)
truth[DE] <- TRUE

### I will have a signal vector of length length(truth)

multiple = rep(1, length(truth))

multiple[DE] = signal

counts_new <- counts(dds_test)

# randomize which condition is shifted up or down
ran <- runif(nrow(dds_test)) 
refcond <- ifelse(ran < 0.5, "A", "B")
down <- which(ran < 0.5)

counts_new[down,colData(dds_test)$fake==unique(refcond[down])] <- 
  counts(dds_test)[down, colData(dds_test)$fake==unique(refcond[down])] *
  multiple[down]
counts_new[-down,colData(dds_test)$fake==unique(refcond[-down])] <- 
  counts(dds_test)[-down, colData(dds_test)$fake==unique(refcond[-down])] *
  multiple[-down]

counts_new <- apply(counts_new, 2, as.integer)


counts(dds_test) <- counts_new

# replace existing size factors 
dds_test <- estimateSizeFactors(dds_test)

dds_test <- DESeq(dds_test, parallel = FALSE)
resTEST <- results(dds_test, name="fake_B_vs_A", independentFiltering = FALSE)

geneExp <- tbl_df(data.frame(pval=resTEST$pvalue, 
                             SE=resTEST$lfcSE,                 
                             ind_covariate = covar,
                             effect_size = resTEST$log2FoldChange, 
                             test_statistic = resTEST$stat,
                             qvalue = truth))

if (uninformativeCovariate){
  geneExp <- mutate(geneExp, ind_covariate = runif(length(covar)))
}else if(!strongCovariate){
  geneExp <- mutate(geneExp, ind_covariate = pmin(1, abs(covar + rnorm(length(covar), 0, 0.25))))
}

geneExp <-  geneExp %>% dplyr::filter(!is.na(pval))

### after that, I will call every methods for testing comparison

mm = function(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df, prop = 1){
  
  if (length(thetaHat) == 0) {
    stop("the input is empty")
  }
  if (!length(thetaHat) == length(s2)) {
    stop("the length of estimated effect sizes must be equal to the estimated squared standard errors")
  }
  if (sum(s2 < 0) > 0) {
    stop("estimated squared standard errors can not be negative")
  }
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
  est.matrix = outer(est.theta, est.sigma)
  est.array = NULL
  for (i in 1:ltheta) {
    est.array = c(est.array, est.matrix[i, ])
  }
  LFDR = matrix(NA, ncol = ltheta, nrow = p0)
  theta = theta0
  s2 = s20
  lik1 = t(exp(-0.5 * (t((outer(theta, grid1, "-"))^2)/(grid2)^2))/(grid2 * 
                                                                      sqrt(2 * pi)))
  y = outer(df * s2, (1/gridsigma^2), "*")
  m = df/2
  lik2 = y^(m - 1) * exp(-0.5 * y)/((2^m) * (gamma(m)))
  lik22 = matrix(rep(lik2, ltheta), nrow = nrow(lik2))
  lik = lik1 * lik22
  
  for (i in 1:p0) {
    ddd = lik[i, ] * est.array
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
  
  lfdr = LFDR[, (Btheta + 1)]
  lfsr = NULL
  for (i in 1:p0) {
    xi = theta[i]
    if (xi > 0) {
      lfsr[i] = sum(LFDR[i, 1:(Btheta + 1)])
    }
    if (xi < 0) {
      lfsr[i] = sum(LFDR[i, (Btheta + 1):ltheta])
    }
  }
  
  return(list(grid.theta = gridtheta, grid.sigma2 = gridsigma^2, 
              mix.theta = est.theta, mix.sigma2 = est.sigma, lfdr = lfdr, 
              lfsr = lfsr))
  
}

mt = mm(geneExp$effect_size, geneExp$SE^2, df = sampleSize*2-2 , prop = 0.01, Btheta = 30)

lfdr.mixtwice = mt$lfdr

## try different methods

covar = geneExp$ind_covariate
covar2 = runif(length(covar))

## unadjust

unadjust = geneExp$pval

## bonferonni

bonferroni = p.adjust(geneExp$pval, method = "bonferroni")

## BH

BH = p.adjust(geneExp$pval, method = "BH")

## storey's qvalue

q.value = qvalue(geneExp$pval)$qvalue

## IHW

ihw_strong = adj_pvalues(ihw(geneExp$pval, covariates = covar, alpha = 0.05))

ihw_weak = adj_pvalues(ihw(geneExp$pval, covariates = covar2, alpha = 0.05))

## ash

ash = get_lfsr(ash(betahat = geneExp$effect_size, sebetahat = geneExp$SE, df = 2*sampleSize - 2))

## BL

BL_strong = lm_pi0(geneExp$pval, X = covar)$pi0 * BH

BL_weak = lm_pi0(geneExp$pval, X = covar2)$pi0 * BH

## AdaPT (GLM)

xx = adapt_glm(x = data.frame(icov = covar), 
               pvals = geneExp$pval, 
               pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
               mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)

AdaPT_strong = xx$qvals

AdaPT_strong[AdaPT_strong>=1] = 1

xx = adapt_glm(x = data.frame(icov = covar2), 
               pvals = geneExp$pval, 
               pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
               mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)

AdaPT_weak = xx$qvals

AdaPT_weak[AdaPT_weak>=1] = 1

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

LFDR_strong = clfdr_hickswrapper(unadj_p = geneExp$pval,
                                  groups = IHW::groups_by_filter(covar, 20), 
                                  lfdr_estimation="fdrtool")

LFDR_weak = clfdr_hickswrapper(unadj_p = geneExp$pval,
                                groups = IHW::groups_by_filter(covar2, 20), 
                                lfdr_estimation="fdrtool")

lfdr.summary = data.frame(unadjust, bonferroni, BH, q.value, ihw_strong, ihw_weak,
                          ash, BL_strong, BL_weak, AdaPT_strong, AdaPT_weak, LFDR_strong, LFDR_weak, lfdr.mixtwice)

save(nDE, sampleSize, signal, geneExp, lfdr.summary,
     file = "./yeast_silico_simulation.RData")

## analysis the data

analysis_FDR=function(lfdr, FDR, truth){
  
  C=FDR
  C_a=TP_a=FP_a=NULL
  a=lfdr
  t=truth
  
  for (cutoff in C){
    
    ## FDR
    
    c_a = sum(t == 0 & a <= cutoff)/sum(a<=cutoff)
    
    ## ture positive
    
    tp_a = sum(t == 1 & a <= cutoff)/p1
    
    C_a=c(C_a, c_a)
    TP_a=c(TP_a, tp_a)
  }
  
  return(list(ControlledFDR=FDR, localFDR=lfdr, 
              EmpiricalFDR=C_a, TruePositive=TP_a))
}

analysis_FDR2=function(lfdr, FDR, truth){
  
  C=FDR
  C_a=TP_a=FP_a=NULL
  a=lfdr
  t=truth
  
  for (cutoff in C){
    
    aa=a[order(a, decreasing = FALSE)]
    
    aaa=cumsum(aa)/c(1:length(aa))
    
    ok.a=aaa<=cutoff
    
    p_a=t[order(a, decreasing = FALSE)][ok.a]
    
    c_a=sum(p_a==0)/length(p_a) ## empirical FDR
    tp_a=sum(p_a==1)/p1 ## true positive rate
    fp_a=sum(p_a==0)/p2 ## false positive rate
    
    C_a=c(C_a, c_a)
    
    TP_a=c(TP_a, tp_a)
    FP_a=c(FP_a, fp_a)
  }
  
  return(list(ControlledFDR=FDR, localFDR=lfdr, 
              EmpiricalFDR=C_a, TruePositive=TP_a, FalsePositive=FP_a))
}

FDR = seq(10^-3, 0.2, by = 0.01)

p1 = sum(geneExp$qvalue)
p2 = length(geneExp$qvalue) - p1

EFDR = matrix(NA, nrow = length(FDR), ncol = dim(lfdr.summary)[2])

colnames(EFDR) = colnames(lfdr.summary)

TPR = EFDR

for(j in 1:dim(EFDR)[2]){
  
  yy = analysis_FDR(lfdr.summary[,j], FDR, geneExp$qvalue)
  
  EFDR[,j] = yy$EmpiricalFDR
  
  TPR[,j] = yy$TruePositive
  
}

## some analysis based on lfdr, not qvalue..

for(j in c(7, 13, 14)){
  
  yy = analysis_FDR2(lfdr.summary[,j], FDR, geneExp$qvalue)
  
  EFDR[,j] = yy$EmpiricalFDR
  
  TPR[,j] = yy$TruePositive
  
}

data = data.frame(efdr = as.numeric(EFDR),
                  etpr = as.numeric(TPR),
                  FDR = rep(FDR, dim(EFDR)[2]),
                  method = rep(colnames(EFDR), each = length(FDR)))

name0 = colnames(EFDR)[order(colnames(EFDR))]

name1 = name0

labels=c("LFDR-s", "LFDR-w", "mixtwice", "AdaPT-s", "AdaPT-w", "BH", "BL-s", "BL-w", "Bonferroni", "ihw-s", "ihw-w", "qvalue")
values = c("tomato2", "springgreen2", 
           "black", 
           "yellow2", "cyan",
           "thistle3", 
           "red", "steelblue2", 
           "purple", 
           "orange","slateblue", 
           "gold")

name1[1] = "AdaPT-s"
name1[2] = "AdaPT-w"
name1[5] = "BL-s"
name1[6] = "BL-w"
name1[7] = "Bonferroni"
name1[8] = "ihw-s"
name1[9] = "ihw-w"
name1[10] = "mixtwice"
name1[11] = "LFDR-s"
name1[12] = "LFDR-w"
name1[13] = 'qvalue'
name1[14] = "two step ash"

col = rep(NA, length(name1))

mm = match(labels, name1)

col[mm] = values

col[3] = "coral"

col[14] = "coral4"

col[15] = "darkgoldenrod"

p1 = ggplot(data = data, aes(x = FDR, y = efdr, col = method))+
  geom_point()+
  geom_smooth(lwd = 1.5)+
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1)+
  scale_color_manual("method",
                     breaks=name0,
                     labels=name1,
                     values = col)+
  labs(title = "FDR control",
       x = "Controlled FDR",
       y = "Empirical FDR")+
  theme(plot.title = element_text(color = "black", size = 16, hjust = 0.5),
        strip.text.x = element_text(size = 12, colour = "black", angle = 0),
        legend.position = "top",
        plot.caption = element_text(color = "black", size = 16, face = "italic", hjust = 1),
        axis.text = element_text(size = 12),
        axis.title.x  = element_text(size = 16, angle = 0),
        axis.title.y  = element_text(size = 16, angle = 90),
        legend.title = element_text(size = 12, angle = 0),
        legend.text = element_text(size = 12, angle = 0),
        axis.line = element_line(linetype = "solid"),
        panel.border = element_rect(linetype = "solid", size = 1.5, fill = NA))

## I calculate the positives, using upset figure

## under the level of 0.01

alpha = 0.05

uu = lfdr.summary

for(j in 1:dim(uu)[2]){
  
  uu[,j] = as.numeric(lfdr.summary[,j]<=alpha)
}

for(j in c(7,13,14)){
  
  x = lfdr.summary[,j]
  
  oo = order(x)
  
  x = x[oo]
  
  pos = c(1:length(x))[oo][(cumsum(x)/c(1:length(x))) <= alpha]
  
  uu[,j] = 0
  
  uu[pos,j] = 1
  
}

colnames(uu) = name1[match(colnames(uu), name0)]

uu$true = as.numeric(geneExp$qvalue)

library(UpSetR)

color2 = rep("black", dim(uu)[2])

name3 = colnames(uu)

name3 = name3[order(colSums(uu), decreasing = T)]

color2[name3 == "true"] = "gold"

color2[name3 == "mixtwice"] = "red"

kk = upset(uu, nsets = dim(uu)[2],order.by = "freq",sets.bar.color = color2)

pdf(file = "./yeast_silico_simulation_figure1.pdf", height = 10, width = 12)

p1

dev.off()

pdf(file = "./yeast_silico_simulation_figure2.pdf", height = 10, width = 14)

kk

dev.off()




