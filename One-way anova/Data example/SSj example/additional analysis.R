## some initial figure and statistic

load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjexample.RData")
load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjStat.RData")

ok.noPA27<-info$CONTAINER=="CTRL_PLACEHOLDER"|info$CONTAINER=="EXPT_HUMAN_PROTEOME"|info$CONTAINER=="EXPT_PROTEOMES_K"|info$CONTAINER=="EXPT_PROTEOMES_NO_RK"|info$CONTAINER=="EXPT_PROTEOMES_R"|info$CONTAINER=="EXPT_PROTS_K"|info$CONTAINER=="EXPT_PROTS_NO_RK"|info$CONTAINER=="EXPT_PROTS_R"
ok=!ok.noPA27

p = a$pval[ok]

q = qvalue(p)

q$pi0

est.pi0 = lm_pi0(p, X = a$covar[ok])$pi0

## I made some figures of pvalue hist and distribution of pi0 estimation

d = data.frame(p = p, pi0 = est.pi0)

a = ggplot(data = d, aes(x = p))+
  geom_histogram(aes(y = ..density..), bins = 100)+
  geom_hline(yintercept = 1, col = "red", lwd = 2, lty = 2)+
  labs(##title = "B",
    x = "p-value",
    y = "Density")+
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

b = ggplot(data = d, aes(x = pi0))+
  geom_histogram(aes(y = ..density..), bins = 100)+
  geom_vline(xintercept = q$pi0, col = "red", lwd = 2, lty = 2)+
  labs(##title = "B",
    x = "estimation of unit-specific null proportion",
    y = "Density")+
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

### something else on the boxplot for those two good peptides

### array data

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

rownames(dat.sig) = info.sig$PROBE_SEQUENCE

dat.sig2 = dat.sig[info.sig$SEQ_ID %in% c("P08603", "P11940"),][-1,]

dd.array = data.frame(signal = as.numeric(as.matrix(dat.sig2)),
                      peptide = rep(rownames(dat.sig2), dim(dat.sig2)[2]),
                      group = rep(c("SSA-", "SSA+", "Control"), each = 16))

c.array = ggplot(data = dd.array, aes(x = group, y = signal, fill = group))+
  geom_boxplot()+
  facet_grid(~peptide)+
  labs(##title = "B: Boxplot for significant peptides",
       x = "group",
       y = "array signal intensity in double log scale")+
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

### elisa data

dat = read.delim("C:\\Users\\appli\\Desktop\\mixtwice anova\\SSj example\\ELISA data.txt", header = T)

l = levels(as.factor(dat$peptide))

i = 4

d.4 = dat[dat$peptide == l[i],]

i = 5

d.5 = dat[dat$peptide == l[i],]

d.elisa = rbind(d.4, d.5)

d.elisa$Subject.Designation[d.elisa$Subject.Designation == "Healthy Control"] = "Control"

d.elisa$peptide[d.elisa$peptide == "P08603"] = "PSQIAQLRPSPR"
d.elisa$peptide[d.elisa$peptide == "P11940"] = "SFTMIGHRSITC"

colnames(d.elisa)[2] = "group"

c.elisa = ggplot(data = d.elisa, aes(x = group, y = ELISA.OD.Value, fill = group))+
  geom_boxplot()+
  facet_grid(~peptide)+
  labs(##title = "B: Boxplot for significant peptides",
    x = "group",
    y = "ELISA signal intensity")+
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

pdf(file = "C://Users//appli//Desktop//mixtwice anova//SSj example//Sjogren data.pdf",
    height = 10, width = 12)

ggarrange(a,b, c.array, c.elisa,
          labels = c("A", "B", "C", "D"), font.label = list(size = 16, color = "red"))

dev.off()