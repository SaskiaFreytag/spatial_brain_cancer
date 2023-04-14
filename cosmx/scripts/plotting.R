
make_log_plot <- function(sce, feature, df, order=FALSE) {
  
  df$feature <- as.numeric(logcounts(sce[feature,]))
  
  if(order) {
    
    ggplot(df[order(df$feature),], aes(x=UMAP1, y=UMAP2, color=feature)) + geom_point(size=1, alpha=0.7) +
      theme_classic() + scale_color_viridis(feature) + ggtitle(paste(feature)) 
    
  } else {
    
    ggplot(df[sample(1:nrow(df)),], aes(x=UMAP1, y=UMAP2, color=feature)) + geom_point(size=1, alpha=0.7) +
      theme_classic() + scale_color_viridis(feature) + ggtitle(paste(feature)) 
    
  }
  
}

make_ind_plot <- function(sce, feature, df) {
  
  df$feature <- as.numeric(logcounts(sce[feature,]))
  df$feature <- df$feature>0
  
  ggplot(df[order(df$feature),], aes(x=UMAP1, y=UMAP2, color=feature)) + geom_point(size=1, alpha=0.7) +
    theme_classic() + scale_color_manual(values=c("TRUE"="red", "FALSE"="lightgrey"), name=paste0(feature)) 
  
}


find_ranks <- function(markers_tmp, gene_sets, go_term){
  
  gene_set <- gene_sets[[go_term]]
  gene_set <- intersect(rownames(markers_tmp), gene_set)
  signs <- sign(markers_tmp$summary.logFC) 
  markers_tmp$Rank <- log10(markers_tmp$FDR)*-1*signs
  markers_tmp$Rank <- rank(markers_tmp$Rank)
  tmp <- markers_tmp[gene_set,]
  data.frame(GO=go_term, gene=rownames(tmp), Rank=tmp$Rank)
  
}