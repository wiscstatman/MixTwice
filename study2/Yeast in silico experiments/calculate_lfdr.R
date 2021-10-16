calculate_lfdr = function(effectsize,
                          se,
                          statvalue,
                          pvalue,
                          covar,
                          sampleSize){
  
  ## unadjust
  
  unadjust = pvalue
  
  ## bonferonni
  
  bonferroni = p.adjust(pvalue, method = "bonferroni")
  
  ## BH
  
  BH = p.adjust(pvalue, method = "BH")
  
  ## storey's qvalue
  
  library(qvalue)
  
  q.value = qvalue(pvalue)$qvalue
  
  ## IHW
  
  library(IHW)
  
  ihw = adj_pvalues(ihw(pvalue, covariates = covar, alpha = 0.05))
  
  ## ash
  
  library(ashr)
  
  ash = get_lfsr(ash(betahat = effectsize, sebetahat = se, df = 2*sampleSize - 2))
  
  ## MixTwice
  
  library(MixTwice)
  
  mt = mixtwice(effectsize, se^2, df = 2*sampleSize - 2, prop = 0.05, Btheta = 30)$lfdr
  
  ## BL
  
  library(swfdr)
  
  BocaLeek = lm_pi0(pvalue, X = covar)$pi0 * BH
  
  ## AdaPT (GLM)
  
  library(adaptMT)
  
  xx = adapt_glm(x = data.frame(icov = covar), 
                 pvals = pvalue, 
                 pi_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), 
                 mu_formulas = paste0("splines::ns(icov, df = ", seq(2, 10, 2), ")"), alphas = 0)
  
  AdaPT_GLM = xx$qvals
  
  return(data.frame(unadjust,
                    bonferroni,
                    BH,
                    q.value,
                    ihw,
                    ash,
                    MixTwice = mt,
                    BocaLeek,
                    AdaPT_GLM))

}

