clfdr_hickswrapper <- function(unadj_p, groups, lfdr_estimation="fdrtool") {
  
  # Exclude this method if there are fewer than 200 tests within a grouping 
  # (fdrtool is applied separately to each group, and throws a warning in such case)
  #if(min(table(groups)) < 200)
  #  stop("Not enough tests to apply this method. Require at least 200.")
  
  ## estimate local fdr within each stratum first
  lfdr_res <- lfdr_fit(unadj_p, groups, lfdr_estimation=lfdr_estimation)
  lfdrs <- lfdr_res$lfdr
  
  ## now use the rejection rule described in Cai's paper
  
  ## Remark:
  ## When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  ## we get monotonic adjusted p-values as a function of the p-values
  ## This is mainly needed for grenander based lfdrs, with most other
  ## lfdr estimation methods lfdr ties are not a problem usually
  
  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]
  return(adj_p)
}

lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){
  
  pvals_list <- split(unadj_p, group)
  
  if (lfdr_estimation == "fdrtool"){
    lfdr_fun <- function(pv) fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
    
  } else if (lfdr_estimation == "locfdr"){
    if (!requireNamespace("locfdr", quietly=TRUE)){
      stop("locfdr package required for this function to work.")
    }
    lfdr_fun <- function(pv) locfdr::locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
  } else {
    stop("This lfdr estimation method is not available.")
  }
  
  lfdr_list <- lapply(pvals_list, lfdr_fun)
  lfdrs <- unsplit(lfdr_list, group)
  
  fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
  fit_obj
}