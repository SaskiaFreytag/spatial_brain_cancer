l_r_interest_df <- function(cell_index, ligand_cells, receptor_cells, 
        l_r_interest, res_means, res_pvalues, res_gene_comb){
  
  res_cell_comb <- cell_index
  res_cell_comb_pa <- apply(res_cell_comb, 1, function(x) paste0(x[1], "|", x[2]))
  
  cell_interest <- res_cell_comb[,1] %in% ligand_cells & res_cell_comb[,2] %in% receptor_cells 
  
  ligand_index <- which(res_gene_comb[,1] %in% l_r_interest)
  receptor_index <- which(res_gene_comb[,2] %in% l_r_interest)
  
  res_means_sub <- res_means[c(ligand_index, receptor_index), cell_interest]
  res_pvalues_sub <- res_pvalues[c(ligand_index, receptor_index), cell_interest]
  l_r_sub <- apply(res_gene_comb[c(ligand_index, receptor_index),], 1, function(x) 
    paste0(x[1], "|", x[2]))
  rownames(res_means_sub) <- l_r_sub
  rownames(res_pvalues_sub) <- l_r_sub
  colnames(res_means_sub) <- res_cell_comb_pa[cell_interest]
  colnames(res_pvalues_sub) <- res_cell_comb_pa[cell_interest]
  
  sub <- apply(res_pvalues_sub, 1, function(x) sum(x<0.001, na.rm = TRUE))!=0
  
  if(all(!sub)) return("None found") else{
    
    res_pvalues_sub <- res_pvalues_sub[sub,, drop=FALSE]
    res_means_sub <- res_means_sub[sub,, drop=FALSE]
    
    df <- reshape2::melt(res_pvalues_sub)
    df$mean <- reshape2::melt(res_means_sub)[,3]
    colnames(df)[1:3] <- c("L-R", "Cells", "pvalue")
    df <- df[!is.na(df$pvalue),]
    df <- df[df$pvalue<0.001,]
    
    return(df)
    
  }
  
}


find_specific_comb <- function(res_cell_comb, ligand_interest, receptor_interest,
   population_interest, type_interest, res_gene_comb, res_pvalues, res_means){
  
  all_combs_recp <- (res_cell_comb[,1] %in% ligand_interest & 
                       res_cell_comb[,2] %in% receptor_interest)
  
  pop_interest <- res_cell_comb[all_combs_recp, type_interest] == population_interest
  
  true_stuff <- apply(res_pvalues[, all_combs_recp] < 0.05, 1, function(x) 
    identical(x, pop_interest))
  
  df_l_r <- l_r_interest_df(res_cell_comb, 
       ligand_interest,  receptor_interest,
       unlist(res_gene_comb[true_stuff,]), 
       res_means, res_pvalues, res_gene_comb)
  
  sub <- paste0(res_gene_comb[true_stuff,1], "|", 
                res_gene_comb[true_stuff,2])
  
  df_l_r_sub <- df_l_r[df_l_r$`L-R` %in% sub, ]
  return(df_l_r_sub)
  
}
