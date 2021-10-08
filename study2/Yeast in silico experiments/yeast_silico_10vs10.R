X = 1

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

save(fdr_summary, truth, file = "F:/Harddrive-Jul-17-2021/MixTwice follow up/yeast_silico_10vs10.RData")
