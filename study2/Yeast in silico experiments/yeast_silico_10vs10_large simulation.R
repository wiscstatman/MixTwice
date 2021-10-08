FDR = seq(0.01, 0.99, by = 0.01)

EFDR = matrix(NA, ncol = dim(fdr_summary)[2], nrow = length(FDR))

ETPR = EFDR

colnames(EFDR) = colnames(fdr_summary)

for (j in 1:dim(fdr_summary)[2]) {
  
  lfdr = fdr_summary[,j]
  
  x.j = analysis_FDR(lfdr, FDR, truth)
  
  efdr = x.j$EmpiricalFDR
  etpr = x.j$TruePositive
  
  EFDR[,j] = efdr
  ETPR[,j] = etpr
  
}

dd = data.frame(empirical_fdr = as.numeric(EFDR),
                controlled_fdr = rep(FDR, dim(fdr_summary)[2]),
                method = rep(colnames(fdr_summary), each = length(FDR)))

group = 'Snf2'

rseed = 1

nDE = 500

sampleSize = 10

signal = 10 ## the fold value that I will add to one conditions, for generating non-null cases

uninformativeCovariate = FALSE
pvalHists = FALSE
strongCovariate = TRUE

load("F:/Harddrive-Jul-17-2021/MixTwice follow up/yeast_fulldata_DESeqdata.RData")

dds_full <- dds_full[,colData(dds_full)$condition == group] ## make a subset group

# set random seed

EFDR2 = NULL
ETPR2 = NULL

for(X in 1:100){
  
  set.seed(as.numeric(X)*as.numeric(rseed))
  
  covar = 1 / (1+exp(-runif(nrow(dds_full), 0, 10) + 5))
  
  # select a random subset of samples of size sampleSize each, with in WT or Snf2 group
  
  dds_test = dds_full[,sample(1:ncol(dds_full), sampleSize*2)]
  
  # add a fake condition column to coldat
  colData(dds_test)$fake =  factor(c(rep("A", sampleSize), 
                                     rep("B", sampleSize))[sample(1:(sampleSize*2), 
                                                                  sampleSize*2)])
  design(dds_test) = ~fake
  
  pzero = rowSums(counts(dds_test)==0)/ncol(counts(dds_test))
  dds_test = dds_test[pzero < 0.5,]
  covar = covar[pzero < 0.5]
  
  truth <- rep(FALSE, nrow(dds_test))
  
  DE <- sample(1:nrow(dds_test), nDE, prob = covar)
  truth[DE] <- TRUE
  
  ### I will have a signal vector of length length(truth)
  
  multiple = rep(1, length(truth))
  
  multiple[DE] = signal
  
  counts_new <- counts(dds_test)
  
  # randomize which condition is shifted up or down
  ran <- runif(nrow(dds_test)) 
  refcond <- ifelse(ran < 0.5, "A", "B")
  down <- which(ran < 0.5)
  
  counts_new[down,colData(dds_test)$fake==unique(refcond[down])] <- 
    counts(dds_test)[down, colData(dds_test)$fake==unique(refcond[down])] *
    multiple[down]
  counts_new[-down,colData(dds_test)$fake==unique(refcond[-down])] <- 
    counts(dds_test)[-down, colData(dds_test)$fake==unique(refcond[-down])] *
    multiple[-down]
  
  counts_new <- apply(counts_new, 2, as.integer)
  
  
  counts(dds_test) <- counts_new
  
  # replace existing size factors 
  dds_test <- estimateSizeFactors(dds_test)
  
  dds_test <- DESeq(dds_test, parallel = FALSE)
  resTEST <- results(dds_test, name="fake_B_vs_A", independentFiltering = FALSE)
  
  geneExp <- tbl_df(data.frame(pval=resTEST$pvalue, 
                               SE=resTEST$lfcSE,                 
                               ind_covariate = covar,
                               effect_size = resTEST$log2FoldChange, 
                               test_statistic = resTEST$stat,
                               qvalue = truth))
  
  if (uninformativeCovariate){
    geneExp <- mutate(geneExp, ind_covariate = runif(length(covar)))
  }else if(!strongCovariate){
    geneExp <- mutate(geneExp, ind_covariate = pmin(1, abs(covar + rnorm(length(covar), 0, 0.25))))
  }
  
  geneExp <-  geneExp %>% dplyr::filter(!is.na(pval))
  
  ### after that, I will call every methods for testing comparison
  
  mt = mixtwice(geneExp$effect_size, geneExp$SE^2, df = 100, prop = 0.01, Btheta = 30)
  
  bd <- initializeBenchDesign()
  
  bd <- addBMethod(bd, "fdrreg-t",
                   FDRreg::FDRreg,
                   function(x) { x$FDR },
                   z = test_statistic,
                   features = model.matrix( ~  splines::bs(ind_covariate, df = 3) - 1),
                   nulltype = 'theoretical',
                   control = list(lambda = 0.01))
  bd <- addBMethod(bd, "fdrreg-e",
                   FDRreg::FDRreg,
                   function(x) { x$FDR },
                   z = test_statistic,
                   features = model.matrix( ~  splines::bs(ind_covariate, df = 3) - 1),
                   nulltype = 'empirical',
                   control = list(lambda = 0.01))
  
  sb <- bd %>% buildBench(data=geneExp, parallel = FALSE)
  
  ll = assay(sb)
  
  ll = ll[,c(0:15, 20, 21)] ## there are several testing procedure not available here
  
  ll[ll[,17] >= 1,17] = 1 ## some lfdr calculation above 1, curious!
  
  ll = ll[,-c(6:14)] ## different parameters for ihw, but almost identical to each other
  
  ll = data.frame(ll, MixTwice = mt$lfdr)
  
  colnames(ll)[5] = "ihw"
  
  fdr_summary = ll
  
  p1 = sum(truth)
  p2 = length(truth) - p1
  
  FDR = seq(0.01, 0.99, by = 0.01)
  
  EFDR = matrix(NA, ncol = dim(fdr_summary)[2], nrow = length(FDR))
  
  ETPR = EFDR
  
  colnames(EFDR) = colnames(fdr_summary)
  
  for (j in 1:dim(fdr_summary)[2]) {
    
    lfdr = fdr_summary[,j]
    
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

for(i in 1:100){
  
  FF = FF + EFDR2[(((length(FDR))*(i-1)+1):(length(FDR)*i)),]
  
  TT = TT + ETPR2[(((length(FDR))*(i-1)+1):(length(FDR)*i)),]
}

FF = FF/100
TT = TT/100

FF = FF[,-7] ## some error with lfdr, since discovery list is most of time, empty

dd = data.frame(empirical_fdr = as.numeric(FF),
                controlled_fdr = rep(FDR, dim(fdr_summary)[2]),
                method = rep(colnames(fdr_summary), each = length(FDR)))

dd = dd[dd$method != "lfdr",] ## some error with lfdr, since discovery list is most of time, empty

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


