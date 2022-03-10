load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjall.RData")
load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjexample.RData")
load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjMixTwice.RData")
load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjStat.RData")

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

F.sig = a$Fstat[ok][okok][okok2]

dat.sig = exp(exp(dat.sig))

mu1 = rowMeans(dat.sig[,1:8])
mu2 = rowMeans(dat.sig[,9:16])
mu3 = rowMeans(dat.sig[,17:24])

fold.SSApositive = mu2/mu3
fold.SSAnegative = mu1/mu3

se1 = apply(dat.sig[,1:8], 1, sd)
se2 = apply(dat.sig[,9:16], 1, sd)
se3 = apply(dat.sig[,17:24], 1, sd)

t.SSApositive = (mu2 - mu3)/(sqrt(se2^2/8 + se3^2/8))
t.SSAnegative = (mu1 - mu3)/(sqrt(se1^2/8 + se3^2/8))

z.SSApositive = qnorm(pt(t.SSApositive, df = 14))
z.SSAnegative = qnorm(pt(t.SSAnegative, df = 14))

p.SSApositive = 2*(1-pnorm(abs(z.SSApositive)))
p.SSAnegative = 2*(1-pnorm(abs(z.SSAnegative)))

output = data.frame(info.sig,
                    Fstat = F.sig,
                    fold.SSApositive,
                    fold.SSAnegative,
                    z.SSApositive,
                    z.SSAnegative,
                    p.SSApositive,
                    p.SSAnegative,
                    dat.sig)

output[,5:35] = round(output[,5:35],3)

write.csv(output, row.names = F, file = "C:/Users/appli/Desktop/mixtwice anova/SSj example/summary_SSj_sig.csv")
