group = 'Snf2'

rseed = 1

nDE = 500

sampleSize = 10

signal = 10

load("F:/Harddrive-Jul-17-2021/MixTwice follow up/yeast_fulldata_DESeqdata.RData")

EFDR2 = NULL
ETPR2 = NULL

N = 200

for(X in 1:N){
  
  x = simulateSplit(X, group, rseed, nDE, sampleSize, signal)
  
  ## I am trying to get the t.test
  
  dA = x$data[,x$group == "A"]
  dB = x$data[,x$group == "B"]
  
  muA = rowMeans(dA)
  muB = rowMeans(dB)
  
  sdA = apply(dA, 1, sd)
  sdB = apply(dB, 1, sd)
  
  effectsize = muA - muB
  
  se = sqrt(sdA^2/sampleSize + sdB^2/sampleSize)
  
  statvalue = effectsize/se
  
  pvalue = 2*(1-pt(abs(statvalue), df = 2*sampleSize - 2))
  
  out = calculate_lfdr(effectsize,
                       se,
                       statvalue,
                       pvalue,
                       covar = x$covar,
                       sampleSize)
  
  out[is.na(out)] = 1
  
  out = out[,-9] ## AdaPT_GLM never worked...
  
  truth = x$truth
  
  p1 = sum(x$truth)
  p2 = length(x$truth) - p1
  
  FDR = seq(0.01, 0.99, by = 0.01)
  
  EFDR = matrix(NA, ncol = dim(out)[2], nrow = length(FDR))
  
  ETPR = EFDR
  
  colnames(EFDR) = colnames(out)
  
  for (j in 1:dim(out)[2]) {
    
    lfdr = out[,j]
    
    x.j = analysis_FDR(lfdr, FDR, truth)
    
    efdr = x.j$EmpiricalFDR
    etpr = x.j$TruePositive
    
    EFDR[,j] = efdr
    ETPR[,j] = etpr
    
  }
  
  EFDR2 = rbind(EFDR2, EFDR)
  ETPR2 = rbind(ETPR2, ETPR)
  
  print(X)
  
}

FF = matrix(0, nrow = length(FDR), ncol = dim(EFDR2)[2])
TT = FF

for(i in 1:N){
  
  FF = FF + EFDR2[(((length(FDR))*(i-1)+1):(length(FDR)*i)),]
  
  TT = TT + ETPR2[(((length(FDR))*(i-1)+1):(length(FDR)*i)),]
}

FF = FF/N
TT = TT/N

dd = data.frame(empirical_fdr = as.numeric(FF),
                controlled_fdr = rep(FDR, dim(out)[2]),
                method = rep(colnames(out), each = length(FDR)))

library(ggplot2)

ggplot(data = dd,aes(x = controlled_fdr, y = empirical_fdr, col = method))+
  geom_line(lwd = 1)+
  geom_point(cex = 0.5)+
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.5)+
  labs(title = "Yeast silico 10 vs 10 in Snf2 condition",
       caption = "Only 500/6000 non-null",
       x = "Controlled False Discovery Rate",
       y = "Empirical False Discovery Rate") + 
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

