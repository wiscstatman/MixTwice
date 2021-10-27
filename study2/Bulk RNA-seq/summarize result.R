load("F:/Harddrive-Jul-17-2021/MixTwice follow up/RNAseq/lfdr.summary.RData")

lfdr.summary$AdaPT_GLM[lfdr.summary$AdaPT_GLM>1] = 1

ok = !is.na(lfdr.summary$lfdr.mixtwice)

lfdr.summary = lfdr.summary[ok,]

## cutoff by 0.05

dat = lfdr.summary

for(j in 1:dim(lfdr.summary)[2]){
  
  dat[,j] = as.numeric(lfdr.summary[,j] <= 0.05)
}

color = rep("black", dim(lfdr.summary)[2])

color[6] = "red"

p.reject = reject_figure(lfdr.summary, title = "rejection figure for Bulk RNA-seq data")

load("E:/github_local/MixTwice/study2/Bulk RNA-seq/GTEx_processed data.RData")

p_effect = ggplot(data=res) +
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


p_test = ggplot(data=res) +
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

p_pval = ggplot(data=res) +
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

library(ggpubr)

figure = ggarrange(p_effect, p_test, p_pval)

pdf(file = "F://Harddrive-Jul-17-2021//MixTwice follow up//RNAseq//RNAseq_visualization.pdf", height = 8, width = 14)

figure

dev.off()

pdf(file = "F://Harddrive-Jul-17-2021//MixTwice follow up///RNAseq///RNAseq_upset&rejection.pdf", height = 8, width = 14)

upset(dat, nsets = dim(lfdr.summary)[2], order.by = 'freq', sets.bar.color = color)

p.reject

dev.off()


