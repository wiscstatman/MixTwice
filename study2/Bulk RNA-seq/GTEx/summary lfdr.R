load("F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/lfdr_all.RData")
load("F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/mixtwice_result.RData")

lfdr.mt[is.na(lfdr.mt)] = 1 ## there is one NA, its s.e. is 0

lfdr.summary$mixtwice = lfdr.mt

lfdr.summary = lfdr.summary[,-1] ## I dont want the unadjust one

lfdr.summary = lfdr.summary[,order(colnames(lfdr.summary))]

color_code = c("yellow2", "green", "thistle3", "red", "purple",
               "orange", "tomato2", "black", "gold", "blue")

cut = seq(0, 0.05, length = 30)

out = matrix(NA, nrow = length(cut), ncol = dim(lfdr.summary)[2])

colnames(out) = colnames(lfdr.summary)

for (j in 1:dim(lfdr.summary)[2]) {
  
  ll = lfdr.summary[,j]
  
  for (i in 1:length(cut)) {
    
    out[i,j] = mean(ll<=cut[i])
    
  }
  
}

for (j in c(2,8,10)) { ## for those who output lfdr rather than qvalue
  
  ll = lfdr.summary[,j]
  
  for(i in 1:length(cut)){
    
    ll = ll[order(ll,decreasing = F)]
    
    out[i,j] = mean(cumsum(ll)/c(1:length(ll)) <= cut[i])
  }
  
}

reject_prop = as.numeric(out)

dd = data.frame(reject_prop,
                method = rep(colnames(lfdr.summary),each = length(cut)),
                FDR_cutoff = rep(cut, dim(lfdr.summary)[2]))

rejection = ggplot(data = dd,aes(x = FDR_cutoff, y = reject_prop, col = method))+
  geom_smooth(lwd = 1.5, se = F)+
  ## geom_point(cex = 0.5)+
  labs(title = "RNA-seq GTEx",
       x = "FDR cutoff",
       y = "Proportion of significance")+
  scale_color_manual("method",
                     breaks=colnames(lfdr.summary),
                     labels=colnames(lfdr.summary),
                     values = color_code)+
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

cut = 0.05

upset = lfdr.summary

for(j in 1:dim(lfdr.summary)[2]){
  
  upset[,j] = as.numeric(lfdr.summary[,j] <= cut)
  
}

for(j in c(3,8,10)){
  
  ll = lfdr.summary[,j]
  
  ref = c(1:length(ll))
  
  oo = order(ll)
  
  ll = ll[oo]
  
  ok = cumsum(ll)/c(1:length(ll)) <= cut
  
  ref.ok = ref[oo][ok]
  
  upset[,j] = 0
  
  upset[ref.ok,j] = 1
  
  
}

library(UpSetR)

cc = rep("black", dim(upset)[2])

cc[dim(upset)[2] - rank(colSums(upset))[colnames(upset) == "mixtwice"] + 1] = "red"

upsetfigure = upset(upset, nsets = dim(upset)[2], nintersects = 20,
                    order.by = "freq", sets.bar.color = cc)

## visualization of mixtwice mixing distribution

grid = c(mix.mixtwice$grid.theta, mix.mixtwice$grid.sigma2)

mix = c(cumsum(mix.mixtwice$mix.theta), cumsum(mix.mixtwice$mix.sigma2))

type = c(rep("effect size (g)", length(mix.mixtwice$grid.theta)), 
         rep("squared standard error (h)", length(mix.mixtwice$grid.sigma2)))

d3 = data.frame(grid, mix, type)

mtfit = ggplot(data = d3, aes(x = grid, y = mix))+
  geom_step(lwd = 2)+
  facet_wrap(~type, scales="free")+
  ylim(c(0,1.0001))+
  labs(x = "grid", y = "cummulative distribution function (cdf)",
       title = "RNA-seq GTEx: Estimated mixing distribution")+
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

## visualize p-value

load("E:/github_local/MixTwice/study2/Bulk RNA-seq/GTEx/GTEx_processed data.RData")

hist_p = ggplot(res, aes(x = pval))+
  geom_histogram(bins = 100, aes(y = ..density..))+
  geom_abline(slope = 0, intercept = 1, lwd = 2, col = "red", lty = 2)+
  labs(x = "p-value", y = "density",
       title = "RNA-seq GTEx: histogram of p-value")+
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

save(rejection, upsetfigure, mtfit, hist_p,
     file = "F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/graphic summary.RData")

pdf("F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/figure_draft.pdf",
    height = 10, width = 14)

hist_p

mtfit

rejection

upsetfigure

dev.off()