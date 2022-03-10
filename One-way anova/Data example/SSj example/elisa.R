dat = read.delim("C:\\Users\\appli\\Desktop\\mixtwice anova\\SSj example\\ELISA data.txt", header = T)

#dat$ELISA.OD.Value = round(dat$ELISA.OD.Value,2)

l = levels(as.factor(dat$peptide))

i = 2

d = dat[dat$peptide == l[i],]

boxplot(d$ELISA.OD.Value~d$Subject.Designation)

fit = aov(d$ELISA.OD.Value~d$Subject.Designation)

summary(fit)

kruskal.test(d$ELISA.OD.Value~d$Subject.Designation)