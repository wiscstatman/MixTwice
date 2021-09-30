simulateSplit = function(X, 
                         group = 'Snf2',
                         rseed, 
                         nDE, 
                         sampleSize,
                         uninformativeCovariate = FALSE, 
                         pvalHists = FALSE,
                         strongCovariate = TRUE){
  
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
  
  # make sure null comparison is truly null
  
  # if PC1 or PC2 sig different, or if test of join means of PCs 1:4 is significant, reshuffle sample labels
  dds_test <- estimateSizeFactors(dds_test)
  x <- t(counts(dds_test, normalize=TRUE))
  pc <- prcomp(log(x + 0.5), scale.=TRUE)
  
  a1 <- pc$x[colData(dds_test)$fake=="A",1]
  b1 <- pc$x[colData(dds_test)$fake=="B",1] 
  a2 <- pc$x[colData(dds_test)$fake=="A",2]
  b2 <- pc$x[colData(dds_test)$fake=="B",2] 
  p1 <- t.test(a1, b1)$p.value
  p2 <- t.test(a2, b2)$p.value
  tries <- 0
  
  while(p1 < 0.10 || p2 < 0.10 && tries < 10){
    colData(dds_test)$fake <- sample(colData(dds_test)$fake, ncol(dds_test))
    x <- t(counts(dds_test, normalize=TRUE))
    pc <- prcomp(log(x + 0.5), scale.=TRUE)
    
    a1 <- pc$x[colData(dds_test)$fake=="A",1]
    b1 <- pc$x[colData(dds_test)$fake=="B",1] 
    a2 <- pc$x[colData(dds_test)$fake=="A",2]
    b2 <- pc$x[colData(dds_test)$fake=="B",2] 
    p1 <- t.test(a1, b1)$p.value
    p2 <- t.test(a2, b2)$p.value
    
    tries <- tries + 1
  }
  
  if(nDE > 0){
    # pick random set of nDE genes to add signal to
    DE <- sample(1:nrow(dds_test), nDE, prob = covar)
    truth[DE] <- TRUE
    
    # randomly sample a log2FC from original FCs (without regard to DE)
    counts_new <- counts(dds_test)
    log2FC <- rep(0, nrow(dds_test))
    
    if(sampleSize == 5){
      log2FC[DE] <- res_5$log2FoldChange[pzero < 0.5][DE]
    }else if (sampleSize == 10){
      log2FC[DE] <- res_10$log2FoldChange[pzero < 0.5][DE]
    }else{
      stop("Only sample sizes 5 and 10 are currently supported with pre-",
           "computed fold changes for sampling from.")
    }
    
    # randomize which condition is shifted up or down
    ran <- runif(nrow(dds_test)) 
    refcond <- ifelse(ran < 0.5, "A", "B")
    down <- which(ran < 0.5)
    
    counts_new[down,colData(dds_test)$fake==unique(refcond[down])] <- 
      counts(dds_test)[down, colData(dds_test)$fake==unique(refcond[down])] *
      2^log2FC[down]
    counts_new[-down,colData(dds_test)$fake==unique(refcond[-down])] <- 
      counts(dds_test)[-down, colData(dds_test)$fake==unique(refcond[-down])] *
      2^log2FC[-down]
    
    counts_new <- apply(counts_new, 2, as.integer)
    
    
    counts(dds_test) <- counts_new
  }
  
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
  
  return(geneExp)
  
}



