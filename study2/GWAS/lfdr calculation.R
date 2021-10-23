### bmi data

load("E:/github_local/MixTwice/study2/GWAS/bmi.RData")

bmi <- bmi %>% 
  dplyr::rename(pval = p_value,
                SE = stderr,
                effect_size = effect,
                ind_covar_N = N,
                ind_covar_AF = Freq_Allele1_HapMapCEU) %>%
  mutate(test_statistic = effect_size / SE) %>%
  select(pval, SE, effect_size, ind_covar_N, ind_covar_AF, test_statistic)

set.seed(39580)
bmi <- bmi %>% mutate(ind_covar_uninf = rnorm(nrow(bmi)))

p_effect = ggplot(data=bmi) +
  geom_histogram(bins=100, aes(x = effect_size, y = ..density..))+
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


p_test = ggplot(data=bmi) +
  geom_histogram(bins=100, aes(x = test_statistic, y = ..density..))+
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

p_pval = ggplot(data=bmi) +
  geom_histogram(bins=100, aes(x = pval, y = ..density..))+
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

p_N = ggplot(data=bmi) +
  geom_histogram(bins=100, aes(x = ind_covar_N, y = ..density..))+
  labs(x = "N")+
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

pdf(file = "F://Harddrive-Jul-17-2021//MixTwice follow up//GWAS_visualization.pdf", height = 8, width = 14)

figure

dev.off()

## unadjust

unadjust = bmi$pval

## bonferonni

bonferroni = p.adjust(bmi$pval, method = "bonferroni")

## BH

BH = p.adjust(bmi$pval, method = "BH")

## storey's qvalue

library(qvalue)

q.value = qvalue(bmi$pval)$qvalue

## IHW

library(IHW)

ihw = adj_pvalues(ihw(bmi$pval, covariates = bmi$ind_covar_N, alpha = 0.05))

## BL

library(swfdr)

BocaLeek = lm_pi0(bmi$pval, X = bmi$ind_covar_N)$pi0 * BH


### try ash

ash.fit = ash(bmi$effect_size, bmi$SE)

ash.lfdr = get_lfdr(ash.fit)

fit.g = get_fitted_g(ash.fit)

ash = ash.lfdr

### MixTwice

theta0 = bmi$effect_size
s20 = bmi$SE^2

Btheta = 30
Bsigma = 10

prop = 10^-3

df = 100

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

lfdr = rep(NA, p0)

N = 10

oo = sample(1:N, dim(bmi)[1], replace = T)

for(nn in 1:10){
  
  LFDR = matrix(NA, ncol = ltheta, nrow = sum(oo==nn))
  
  theta = theta0[oo == nn]
  s2 = s20[oo == nn]
  lik1 = t(exp(-0.5 * (t((outer(theta, grid1, "-"))^2)/(grid2)^2))/(grid2 * 
                                                                      sqrt(2 * pi)))
  for (i in 1:sum(oo == nn)) {
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
  
  lfdr[oo==nn] = LFDR[, (Btheta + 1)]
  
}

mt = lfdr.mixtwice


lfdr.summary = data.frame(unadjust,
                          bonferroni,
                          BH,
                          q.value,
                          ihw,
                          BocaLeek,
                          ash,
                          lfdr.mixtwice)

save(lfdr.summary, file = "E://github_local//MixTwice//study2//GWAS//lfdr.summary.RData")

lfdr2 = lfdr.summary

for (j in 1:dim(lfdr2)[2]) {
  
  lfdr2[,j] = as.numeric(lfdr2[,j] <= 0.05)
  
}

color = rep("black", 7)
color[6] = "red"

pdf(file = "F://Harddrive-Jul-17-2021//MixTwice follow up//GWAS_upset.pdf", height = 8, width = 14)

upset(lfdr2[,-1], nsets = dim(lfdr2)[2], order.by = 'freq', sets.bar.color = color)

dev.off()
