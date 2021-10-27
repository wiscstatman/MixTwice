reject_figure = function(lfdr.summary, cut = seq(0, 1, by = 0.01),
                         title = "rejection figure for Bulk RNA-seq data"){
  
  ## the input is a lfdr.summary matrix
  
  out = matrix(NA, nrow = length(cut), ncol = dim(lfdr.summary)[2])
  
  colnames(out) = colnames(lfdr.summary)
  
  for (j in 1:dim(lfdr.summary)[2]) {
    
    ll = lfdr.summary[,j]
    
    for (i in 1:length(cut)) {
      
      out[i,j] = mean(ll<=cut[i])
      
    }
    
  }
  
  reject_prop = as.numeric(out)
  
  dd = data.frame(reject_prop,
                  method = rep(colnames(lfdr.summary),each = length(cut)),
                  FDR_cutoff = rep(cut, dim(lfdr.summary)[2]))
  ggplot(data = dd,aes(x = FDR_cutoff, y = reject_prop, col = method))+
    geom_line(lwd = 1.5)+
    geom_point(cex = 0.5)+
    labs(title = title,
         x = "FDR cutoff",
         y = "Discovery rate") + 
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
  
}