##### I want to write a function to help me analysis the FDR control

### The input should be localLDR and FDR as a controlled FDR and the truth

analysis_FDR=function(lfdr, FDR, truth){
  
  C=FDR
  C_a=TP_a=FP_a=NULL
  a=lfdr
  t=truth
  
  for (cutoff in C){
    
    aa=a[order(a, decreasing = FALSE)]
    
    aaa=cumsum(aa)/c(1:length(aa))
    
    ok.a=aaa<=cutoff
    
    p_a=t[order(a, decreasing = FALSE)][ok.a]
    
    c_a=sum(p_a==0)/length(p_a) ## empirical FDR
    tp_a=sum(p_a==1)/p1 ## true positive rate
    fp_a=sum(p_a==0)/p2 ## false positive rate
    
    C_a=c(C_a, c_a)
    
    TP_a=c(TP_a, tp_a)
    FP_a=c(FP_a, fp_a)
  }
  
  return(list(ControlledFDR=FDR, localFDR=lfdr, 
              EmpiricalFDR=C_a, TruePositive=TP_a, FalsePositive=FP_a))
}
