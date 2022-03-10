load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjexample.RData")

## calculate some useful statistics

d1 = dat[,group == "SSA-"]
d2 = dat[,group == "SSA+"]
d3 = dat[,group == "Control"]

mu1 = rowMeans(d1)
mu2 = rowMeans(d2)
mu3 = rowMeans(d3)

mu = rowMeans(dat)

SSB = 8*((mu1 - mu)^2 + (mu2 - mu)^2 + (mu3 - mu)^2)

SST = apply(as.matrix(dat), 1, sd)^2*(23)

SSE = SST - SSB

n = length(group)
m = nlevels(as.factor(group))

df1 = m - 1
df2 = n - m

Fstat = (SSB/df1)/(SSE/df2)

pval = 1-pf(Fstat, df1 = df1, df2 = df2)

covar = sqrt(SST/23)

a = list(SSB = SSB, SSE = SSE, 
         df1 = df1, df2 = df2, 
         Fstat = Fstat, pval = pval, 
         covar = covar)

save(a, file = "C:/Users/appli/Desktop/mixtwice anova/SSjStat.RData")


## try MixTwice-ANOVA..

ok.noPA27<-info$CONTAINER=="CTRL_PLACEHOLDER"|info$CONTAINER=="EXPT_HUMAN_PROTEOME"|info$CONTAINER=="EXPT_PROTEOMES_K"|info$CONTAINER=="EXPT_PROTEOMES_NO_RK"|info$CONTAINER=="EXPT_PROTEOMES_R"|info$CONTAINER=="EXPT_PROTS_K"|info$CONTAINER=="EXPT_PROTS_NO_RK"|info$CONTAINER=="EXPT_PROTS_R"
ok=!ok.noPA27

mixtwice_ANOVA = function(SSB, SSE, df1, df2, Blambda = 20, Bsigma2 = 8, prop = 0.1){
  
  SSB0 = SSB
  SSE0 = SSE
  
  l = length(SSB)
  
  ok.sample = sample(l, l * prop)
  
  SSB = SSB0[ok.sample]
  SSE = SSE0[ok.sample]
  
  prop = 0.01
  
  l = length(SSB)
  
  Blambda = 20
  
  Bsigma2 = 8
  
  kk1 = SSB - (df1)*SSE/(df2)
  kk2 = (SSB/SSE)*(df2 - 1)*SSE/(df2) - df1*SSE/df2
  
  grid.lambda = seq(0, max(kk1), length = Blambda)
  
  grid.sigma2 = seq(min(SSE/(df2))/1.1, 1.1*max(SSE/(df2)), length = Bsigma2)
  
  ## conditional likelihood
  
  likelihood.SSE = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  likelihood.SSB = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  
  grid.lambda2 = rep(grid.lambda, each = Bsigma2)
  grid.sigma22 = rep(grid.sigma2, Blambda)
  
  for(i in 1:l){
    
    likelihood.SSE[i,] = dchisq(SSE[i]/grid.sigma22, df = df2)/grid.sigma22
    
    likelihood.SSB[i,] = dchisq(SSB[i]/grid.sigma22, df = df1, ncp = grid.lambda2/grid.sigma22)/grid.sigma22
    
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
  
  #try1 = suppressWarnings(alabama::auglag(par = a, fn = L, 
  #                                        gr = G, heq = heq, hin = hin, heq.jac = heq.jac.fun, 
  #                                        hin.jac = hin.jac.fun,control.outer = list(trace = F)))
  
  try1 = suppressWarnings(alabama::auglag(par = a, fn = L, 
                                          gr = G, heq = heq, hin = hin, heq.jac = heq.jac.fun, 
                                          hin.jac = hin.jac.fun))
  
  est.lambda = try1$par[1:Blambda]
  est.lambda[est.lambda < 0] = 0
  
  est.sigma2 = try1$par[(Blambda + 1):(Blambda + Bsigma2)]
  est.sigma2[est.sigma2 < 0] = 0
  
  est.matrix = outer(est.lambda, est.sigma2)
  
  est.array = NULL
  for (i in 1:Blambda) {
    est.array = c(est.array, est.matrix[i, ])
  }
  
  ## conditional likelihood on the whole data
  
  SSB = SSB0
  SSE = SSE0
  
  l = length(SSB)
  
  likelihood.SSE = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  likelihood.SSB = matrix(NA, nrow = l, ncol = Blambda*Bsigma2)
  
  for(i in 1:l){
    
    likelihood.SSE[i,] = dchisq(SSE[i]/grid.sigma22, df = df2)/grid.sigma22
    
    likelihood.SSB[i,] = dchisq(SSB[i]/grid.sigma22, df = df1, ncp = grid.lambda2/grid.sigma22)/grid.sigma22
    
  }
  
  likelihood.conditional = likelihood.SSB * likelihood.SSE
  
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
  
  return(list(grid.lambda = grid.lambda, grid.sigma2 = grid.sigma2, 
              mix.lambda = est.lambda, mix.sigma2 = est.sigma2, lfdr = lfdr))
  
}

att2 = mixtwice_ANOVA(SSB = a$SSB[ok], SSE = a$SSE[ok],
                      df1 = a$df1, df2 = a$df2,
                      Blambda = 20,
                      Bsigma2 = 8,
                      prop = 0.1)

save(att, file = "./SSjMixTwice.RData")

pval = a$pval[ok]
covar = a$covar[ok]

## bonferonni

p_bonferroni = p.adjust(pval, method = "bonferroni")

## BH

p_BH = p.adjust(pval, method = "BH")

## qvalue

p_qvalue = qvalue(pval)$qvalue

## IHW

library(IHW)

p_ihw = adj_pvalues(ihw(pval, covariates = covar, alpha = 0.05))

##  BL

p_BL = lm_pi0(pval, X = covar)$pi0 * p_ihw

## AdaPT (GLM)

## this does not work

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

LFDR = clfdr_hickswrapper(unadj_p = pval,
                           groups = IHW::groups_by_filter(covar, 20), 
                           lfdr_estimation="fdrtool")

out = list(pval = a$pval[ok],
           bonferroni = p_bonferroni,
           BH = p_BH,
           qvalue = p_qvalue,
           ihw = p_ihw,
           BL = p_BL,
           LFDR = LFDR,
           MixTwice = att$lfdr)

save(out, file = "./SSjall.RData")

ok.noPA27<-info$CONTAINER=="CTRL_PLACEHOLDER"|info$CONTAINER=="EXPT_HUMAN_PROTEOME"|info$CONTAINER=="EXPT_PROTEOMES_K"|info$CONTAINER=="EXPT_PROTEOMES_NO_RK"|info$CONTAINER=="EXPT_PROTEOMES_R"|info$CONTAINER=="EXPT_PROTS_K"|info$CONTAINER=="EXPT_PROTS_NO_RK"|info$CONTAINER=="EXPT_PROTS_R"
ok=!ok.noPA27

okok = out$MixTwice <=0.05

dat = dat[ok,]
info = info[ok,]

dat.sig = dat[out$MixTwice<=0.05,]

mu1 = rowMeans(dat.sig[,1:8])
mu2 = rowMeans(dat.sig[,9:16])
mu3 = rowMeans(dat.sig[,17:24])

okok2 = mu2>=mu1 & mu2>=mu3

dat.sig = dat.sig[okok2,]

info.sig = info[,c(2, 5, 6, 14)][okok,][okok2,]

output = cbind(info.sig, dat.sig)

## write.csv(output, row.names = F, file = "./SSj_sig.csv")

rr = rank(a$Fstat[ok])/sum(ok)

col = rep("100%", sum(ok))

col[rr<=0.8] = "80%"
col[rr<=0.6] = "60%"
col[rr<=0.4] = "40%"
col[rr<=0.2] = "20%"

d = data.frame(SSB = a$SSB[ok],
               SSE = a$SSE[ok],
               quantitle = col)

A = ggplot()+
  geom_point(data = d[sample(dim(d)[1], 0.1*dim(d)[1]),], 
             aes(x = SSB, y = SSE, col = quantitle),
             cex = 1)+
  labs(title = "A: Scatter plot between SSB and SSE",
    x = "SSB",
    y = "SSE")+
  scale_y_log10()+
  scale_x_log10()+
  geom_point(data =d[okok,], aes(x = SSB, y = SSE),
             cex = 1, pch = 15)+
  theme(plot.title = element_text(color = "red", size = 16, hjust = 0),
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


## let's see the boxplot representation...

d1 = data.frame(signal = as.numeric(as.matrix(dat)),
                group = rep(c("SSA-", "SSA+", "Control"), each = dim(dat)[1]*8))

d1 = d1[sample(dim(d1)[1], 0.01*dim(d1)[1]),]

d2 = data.frame(signal = as.numeric(as.matrix(dat.sig)),
                group = rep(c("SSA-", "SSA+", "Control"), each = dim(dat.sig)[1]*8),
                peptide = rep(as.character(info.sig$PROBE_SEQUENCE)), 24)

B = ggplot(data = d2, aes(x = group, y = signal, fill = group))+
  geom_boxplot()+
  facet_wrap(~peptide)+
  labs(title = "B: Boxplot for significant peptides",
       x = "group",
       y = "signal in double log scale")+
  theme(plot.title = element_text(color = "red", size = 16, hjust = 0),
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

# d = rbind(d1,d2)
# 
# d$type = c(rep("all peptides", dim(d1)[1]), rep("significant peptides", dim(d2)[1]))
# 
# C = ggplot(data = d, aes(x = group, y = signal, fill = group))+
#   #geom_boxplot(outlier.colour = NA)+
#   geom_boxplot()+
#   #scale_y_continuous(limits=c(1.1,1.7))+
#   facet_grid(~type)+
#   labs(title = "C: Boxplot for all peptides and significant peptides",
#        x = "group",
#        y = "signal in double log scale")+
#   theme(plot.title = element_text(color = "red", size = 16, hjust = 0),
#         strip.text.x = element_text(size = 12, colour = "black", angle = 0),
#         legend.position = "top",
#         plot.caption = element_text(color = "black", size = 16, face = "italic", hjust = 1),
#         axis.text = element_text(size = 12),
#         axis.title.x  = element_text(size = 16, angle = 0),
#         axis.title.y  = element_text(size = 16, angle = 90),
#         legend.title = element_text(size = 12, angle = 0),
#         legend.text = element_text(size = 12, angle = 0),
#         axis.line = element_line(linetype = "solid"),
#         panel.border = element_rect(linetype = "solid", size = 1.5, fill = NA))


pdf(file = "C://Users//appli//Desktop//mixtwice anova//SSj example//f.pdf", 
    height = 12, width = 12)

ggarrange(A, B, nrow = 2,
          heights = c(1,2))

dev.off()
