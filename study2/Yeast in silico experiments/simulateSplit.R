simulateSplit = function(X, 
                         group = 'Snf2',
                         rseed, 
                         nDE, 
                         sampleSize,
                         signal ## fold value with DE
){
  
  dds_full <- dds_full[,colData(dds_full)$condition == group] ## make a subset group
  
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
  
  return(list(data = counts_new,
              truth = truth,
              group =colData(dds_test)$fake,
              covar = covar))
}