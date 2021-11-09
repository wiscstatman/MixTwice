load("C:/Users/appli/Desktop/compare MixTwice-ANOVA with all others.RData")

l = 3000 ## number of testing units

pi0_true = seq(0.1, 0.9, length = 30)

tt = c(1, 2, 3)

ss = c(1, 2, 3)

input = matrix(NA, nrow = length(out), ncol = 3)

colnames(input) = c("s", "t", "pi0")

i = 1

for(t in tt){
  
  for (s in ss) {
    
    for (pi0 in pi0_true) {
      
      input[i,] = c(t, s, pi0)
      
      i = i + 1
      
      print(i)
      
    }
    
  }
}

positive = matrix(NA, nrow = dim(input)[1], ncol = dim(out[[1]]$summary.table)[2])

colnames(positive) = colnames(out[[1]]$summary.table)

for (i in 1:dim(positive)[1]) {
  
  tt = out[[i]]$summary.table
  
  rr = out[[i]]$truth
  
  ## others are based on qvalue, so we just cut by 0.01 FDR
  
  for(j in (1:(dim(positive)[2]-1))){
    
    positive[i,j] = sum(tt[rr==1,j] <= 0.1)/l
  }
  
  ## for mixtwice, we ouput the lfdr
  
  tt = tt$mixtwice
  
  oo = order(tt)
  
  tt = tt[oo]
  
  tt2 = cumsum(tt)/c(1:l)
  
  rr = rr[oo]
  
  positive[i,dim(positive)[2]] = sum(tt2[rr == 1] <= 0.1)/l
}

dat = NULL

for (j in 1:dim(positive)[2]) {
  
  dat = rbind(dat, input)
  
}

dat = data.frame(dat)

dat$method = rep(colnames(positive), each = dim(positive)[1])

dat$truepositive = as.numeric(positive)

ok = dat$s <= 2 & dat$t %in% c(2,3)

dat = dat[ok,]

dat$s[dat$s == 1] = "small signal"
dat$s[dat$s == 2] = "large signal"

dat$t[dat$t == 2] = "small variance"
dat$t[dat$t == 3] = "large variance"

ok = dat$method %in% colnames(positive)[c(10,11,12)]

sensitivity = ggplot(data = dat, aes(x = pi0, y = truepositive))+
  ## geom_point(aes(color = method))+
  geom_smooth(aes(color = method), lwd = 1.5, method = "lm", level = 0.5)+
  facet_grid(s~t)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  geom_abline(slope = -1, intercept = 1, lty = 2, lwd = 1)+
  scale_color_manual("method",
                     breaks=c("LFDR_s", "LFDR_w", "mixtwice", "p_AdaPT_s", "p_AdaPT_w", "p_BH", "p_BL_s", "p_BL_w", "p_bonferroni", "p_ihw_s", "p_ihw_w", "p_qvalue"),
                     labels=c("LFDR-s", "LFDR-w", "mixtwice-ANOVA", "AdaPT-s", "AdaPT-w", "BH", "BL-s", "BL-w", "Bonferroni", "ihw-s", "ihw-w", "qvalue"),
                     values = c("tomato2", "springgreen2", 
                                "black", 
                                "yellow2", "cyan",
                                "thistle3", 
                                "red", "steelblue2", 
                                "purple", 
                                "orange","slateblue", 
                                "gold"))+
  labs(title = "Sensitivity analysis",
       x = "Simulated null propotion (pi0)",
       y = "Proportion of true positives")+
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

p_mixtwice = positive[,12]
p_no = rowMeans(positive[,c(1,2,3)])
p_strong = rowMeans(positive[,c(4,6,8,10)])
p_weak = rowMeans(positive[,c(5,7,9,11)])

positive2 = cbind(p_mixtwice, p_no, p_strong, p_weak)

dat2 = NULL

for (j in 1:dim(positive2)[2]) {
  
  dat2 = rbind(dat2, input)
  
}

dat2 = data.frame(dat2)

dat2$method = rep(colnames(positive2), each = dim(positive2)[1])

dat2$truepositive = as.numeric(positive2)

ok2 = dat2$s <= 2 & dat2$t %in% c(2,3)

dat2 = dat2[ok2,]

dat2$s[dat2$s == 1] = "small signal"
dat2$s[dat2$s == 2] = "large signal"

dat2$t[dat2$t == 2] = "small variance"
dat2$t[dat2$t == 3] = "large variance"

sensitivity2 = ggplot(data = dat2, aes(x = pi0, y = truepositive))+
  ## geom_point(aes(color = method))+
  geom_smooth(aes(color = method), lwd = 1.5, method = "lm", level = 0.5)+
  facet_grid(s~t)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  geom_abline(slope = -1, intercept = 1, lty = 2, lwd = 1)+
  scale_color_manual("method",
                     breaks=c("p_mixtwice", "p_no", "p_strong", "p_weak"),
                     labels=c("mixtwice-ANOVA", "methods without covariate", "methods with strong covariate", "methods with mis-specified covariate"),
                     values = c("black", "springgreen2", 
                                "tomato2", 
                                "yellow2"))+
  labs(title = "Sensitivity analysis",
       x = "Simulated null propotion (pi0)",
       y = "Proportion of true positives")+
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

FDR = seq(0.001, 0.2, by = 0.01)

eFDR2 = NULL

for(kk in 1:30){
  
  okok = input[,1] == 2 & input[,2] == 1 & input[,3] == seq(0.1, 0.9, length = 30)[kk]
  
  oo = out[okok][[1]]
  
  pi0 = input[okok,3]
  
  truth = oo$truth
  
  p1 = sum(truth)
  p2 = length(truth) - p1
  
  a = NULL
  
  for (j in 1:dim(positive)[2]) {
    
    a[[j]] = analysis_FDR(oo$summary.table[,j], FDR, truth)
    
  }
  
  a[[12]] = analysis_FDR2(oo$summary.table[,12], FDR, truth)
  
  eFDR = NULL
  
  for(j in 1:dim(positive)[2]){
    
    eFDR = c(eFDR, a[[j]]$EmpiricalFDR)
  }
  
  eFDR2 = rbind(eFDR2, eFDR)
  
}

d = data.frame(cFDR = rep(FDR, dim(positive)[2]),
               eFDR = colMeans(eFDR2, na.rm = T),
               method = rep(colnames(positive), each = length(FDR)))

p_FDR = ggplot(data = d, aes(x = cFDR, y = eFDR))+
  geom_smooth(lwd = 1.5, aes(color = method), method = "lm")+
  ##geom_point(aes(color = method))+
  #xlim(c(0,0.1))+
  #ylim(c(0,0.1))+
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1)+
  scale_color_manual("method",
                     breaks=c("LFDR_s", "LFDR_w", "mixtwice", "p_AdaPT_s", "p_AdaPT_w", "p_BH", "p_BL_s", "p_BL_w", "p_bonferroni", "p_ihw_s", "p_ihw_w", "p_qvalue"),
                     labels=c("LFDR-s", "LFDR-w", "mixtwice-ANOVA", "AdaPT-s", "AdaPT-w", "BH", "BL-s", "BL-w", "Bonferroni", "ihw-s", "ihw-w", "qvalue"),
                     values = c("tomato2", "springgreen2", 
                                "black", 
                                "yellow2", "cyan",
                                "thistle3", 
                                "red", "steelblue2", 
                                "purple", 
                                "orange","slateblue", 
                                "gold"))+
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

eTPR2 = NULL

for(kk in 1:30){
  
  okok = input[,1] == 2 & input[,2] == 1 & input[,3] == seq(0.1, 0.9, length = 30)[kk]
  
  oo = out[okok][[1]]
  
  pi0 = input[okok,3]
  
  truth = oo$truth
  
  p1 = sum(truth)
  p2 = length(truth) - p1
  
  a = NULL
  
  for (j in 1:dim(positive)[2]) {
    
    a[[j]] = analysis_FDR(oo$summary.table[,j], FDR, truth)
    
  }
  
  a[[12]] = analysis_FDR2(oo$summary.table[,12], FDR, truth)
  
  eTPR = NULL
  
  for(j in 1:dim(positive)[2]){
    
    eTPR = c(eTPR, a[[j]]$TruePositive)
  }
  
  eTPR2 = rbind(eTPR2, eTPR)
  
}
  
d = data.frame(cFDR = rep(FDR, dim(positive)[2]),
               eTPR = colMeans(eTPR2, na.rm = T),
               method = rep(colnames(positive), each = length(FDR)))

p_power = ggplot(data = d, aes(x = cFDR, y = eTPR))+
  geom_smooth(lwd = 1.5, aes(color = method), level = 0.5)+
  ## geom_point(aes(color = method))+
  scale_color_manual("method",
                     breaks=c("LFDR_s", "LFDR_w", "mixtwice", "p_AdaPT_s", "p_AdaPT_w", "p_BH", "p_BL_s", "p_BL_w", "p_bonferroni", "p_ihw_s", "p_ihw_w", "p_qvalue"),
                     labels=c("LFDR-s", "LFDR-w", "mixtwice-ANOVA", "AdaPT-s", "AdaPT-w", "BH", "BL-s", "BL-w", "Bonferroni", "ihw-s", "ihw-w", "qvalue"),
                     values = c("tomato2", "springgreen2", 
                                "black", 
                                "yellow2", "cyan",
                                "thistle3", 
                                "red", "steelblue2", 
                                "purple", 
                                "orange","slateblue", 
                                "gold"))+
  labs(title = "Power",
       x = "Controlled FDR",
       y = "Empirical True Positive Rate")+
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

load("C:/Users/appli/Desktop/compare MixTwice-ANOVA with BH.RData")

rate1 = NULL
rate2 = NULL

pi0est = NULL

for (i in 1:length(out)) {
  
  x = out[[i]]
  
  x1 = x$padj
  x2 = x$lfdr
  
  pi0i = input[i,3]
  
  rate1[i] = sum(x1[((round(pi0i*l) + 1):l)] <= 0.01)/l
  rate2[i] = sum(x2[((round(pi0i*l) + 1):l)] <= 0.01)/l
  
  pi0est[i] = x$pi0est
  
  
}

dat = data.frame(rbind(input, input))

dat$method = rep(c("BH", "MixTwice-ANOVA"), each = dim(input)[1])

dat$truepositive = c(rate1, rate2)

ok = dat$s >= 2 & dat$t %in% c(1,3)

dat = dat[ok,]

dat$s[dat$s == 2] = "small signal"
dat$s[dat$s == 3] = "large signal"

dat$t[dat$t == 1] = "small variance"
dat$t[dat$t == 3] = "large variance"

dat = data.frame(input, pi0est = pi0est)

ok = dat$s >= 2 & dat$t %in% c(1,3)

dat = dat[ok,]

dat$s[dat$s == 2] = "small signal"
dat$s[dat$s == 3] = "large signal"

dat$t[dat$t == 1] = "small variance"
dat$t[dat$t == 3] = "large variance"

p_pi0 = ggplot(data = dat, aes(x = pi0, y = pi0est))+
  geom_point()+
  geom_smooth(lwd = 1.5)+
  facet_grid(s~t)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1)+
  labs(title = "Consistency of non-null propotion estimation",
       x = "Simulated null proportion (pi0)",
       y = "Estimated null proportion (pi0)")+
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

pdf(file = "C://Users//appli//Desktop//f1.pdf", height = 12, width = 10)

p_pi0

ggarrange(sensitivity, sensitivity2, nrow = 2)

ggarrange(p_FDR, p_power, nrow = 2)

dev.off()

