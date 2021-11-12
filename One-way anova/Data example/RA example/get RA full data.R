### array1, full data

load("F:/Harddrive-Jul-17-2021/R01/data&code/Rdata/try2-G.RData")
ok = c(1:12, 37:48, 55:66)
group1 = rep(c("CCP+RF+", "CCP-RF-", "Control"), each = 12)
dat.full1 = tpdat2[,ok]

load("F:/Harddrive-Jul-17-2021/R01/data&code/Rdata/info.RData")
info1 = info

### array2, full data

load("F:/Harddrive-Jul-17-2021/R02_02/try2 and clustering_merged data/try3_merge_Ops.RData")

info=info3
ok.noPA27<-info$CONTAINER=="CTRL_PLACEHOLDER"|info$CONTAINER=="EXPT_HUMAN_PROTEOME"|info$CONTAINER=="EXPT_PROTEOMES_K"|info$CONTAINER=="EXPT_PROTEOMES_NO_RK"|info$CONTAINER=="EXPT_PROTEOMES_R"|info$CONTAINER=="EXPT_PROTS_K"|info$CONTAINER=="EXPT_PROTS_NO_RK"|info$CONTAINER=="EXPT_PROTS_R"
ok.PA27=!ok.noPA27

info2=info[ok.PA27,]

dat.full2 = cbind(tpdat.disease[ok.PA27, 41:48], 
                  tpdat.disease[ok.PA27, 49:56],
                  tpdat.control2[ok.PA27, 41:48])

group2 = rep(c("CCP+RF+", "CCP-RF-", "Control"), each = 8)

#################################################################################
#################################################################################

### delete duplicated and match

du1 = duplicated(info1$PROBE_SEQUENCE)

info1 = info1[!du1,]

dat.full1 = dat.full1[!du1,]

du2 = duplicated(info2$PROBE_SEQUENCE)

info2 = info2[!du2,]

dat.full2 = dat.full2[!du2,]

mm = match(info1$PROBE_SEQUENCE, info2$PROBE_SEQUENCE)

info2 = info2[mm,]

dat.full2 = dat.full2[mm,]

info = info1

save(info, dat1, dat2, group1, group2,
     file = "C://Users//appli//Desktop//mixtwice anova//RAexample.RData")
