### limma and sva check for sjogren disease (8 vs 8 vs 8)

load("C:/Users/appli/Desktop/mixtwice anova/SSj example/SSjexample.RData")

ok.noPA27<-info$CONTAINER=="CTRL_PLACEHOLDER"|info$CONTAINER=="EXPT_HUMAN_PROTEOME"|info$CONTAINER=="EXPT_PROTEOMES_K"|info$CONTAINER=="EXPT_PROTEOMES_NO_RK"|info$CONTAINER=="EXPT_PROTEOMES_R"|info$CONTAINER=="EXPT_PROTS_K"|info$CONTAINER=="EXPT_PROTS_NO_RK"|info$CONTAINER=="EXPT_PROTS_R"
ok=!ok.noPA27

dat = dat[ok,]

info = info[ok,]

### try limma first

library(edgeR)
library(limma)

dgecounts = calcNormFactors(DGEList(counts=dat, group=group))

design = model.matrix(~group)
v = voom(dgecounts,design,plot=FALSE)

lim = lmFit(v)

## use Limma topTable

lim = eBayes(lim)
ans = topTable(lim, number = 100,
               adjust.method = "BH", sort.by = "F")

summary(ans$adj.P.Val)

### Then I try sva

library(sva)

mod = model.matrix(~as.factor(group), data = data.frame(group))
mod0 = model.matrix(~1, data = data.frame(group))

n.sv.be = num.sv(as.matrix(dat),mod,method="be")

n.sv.leek = num.sv(as.matrix(dat),mod,method="leek")

sv.be = sva(as.matrix(dat),mod, mod0 = mod0, n.sv = 2)

summary(sv.be$pprob.b)

## sv.leek = sva(as.matrix(dat),mod, mod0 = mod0, n.sv = n.sv.leek) ## since n.sv.leek = 0, not useful