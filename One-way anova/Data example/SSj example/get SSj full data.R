load("F:/Harddrive-Jul-17-2021/R02_02/try2 and clustering_merged data/try3_merge_Ops.RData")

info=info3

dat = cbind(tpdat.disease[,25:40],
            tpdat.control2[,25:32])

dat = round(dat, 3)

group = rep(c("SSA-", "SSA+", "Control"), each = 8)

colnames(dat) = group

save(info, dat, group,
     file = "C://Users//appli//Desktop//mixtwice anova//SSjexample.RData")