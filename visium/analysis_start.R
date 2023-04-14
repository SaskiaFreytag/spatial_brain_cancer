library(Seurat)
library(ggplot2)
library(spatstat)
library(paletteer)
library(wesanderson)

setwd("gbm_spatial/doi_10.5061_dryad.h70rxwdmj__v11/10XVisium_2/")

features <- c("APOE", "A2M", "VEGFA", "EGR1", "GLIS3", "CHIC2", "C4B")

meta <- readRDS("Features_all_new.RDS")
names(meta) <- sapply(meta, function(x) x$sample[1])

sc <- list()

for(i in 1:length(meta)) {

  sc[[i]] <- Load10X_Spatial(
    paste0("#UKF", names(meta)[i], "_ST/outs"),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1"
  )

  sc[[i]]@meta.data <- cbind(sc[[i]]@meta.data, meta[[i]][
    match(colnames(sc[[i]]), meta[[i]]$barcodes),])
  sc[[i]]$MES_like <- sc[[i]]$MES_like_hypoxia_independent +
    sc[[i]]$MES_like_hypoxia_MHC
  sc[[i]]$NPC_like <- sc[[i]]$NPC_like_neural + sc[[i]]$NPC_like_OPC +
    sc[[i]]$NPC_like_Prolif
  sc[[i]]$TAM <- sc[[i]]$TAM_BDM_anti_infl + sc[[i]]$TAM_BDM_hypoxia_MES + 
    sc[[i]]$TAM_BDM_INF + sc[[i]]$TAM_BDM_MHC 
  
  sc[[i]] <- NormalizeData(sc[[i]], assay = "Spatial", verbose = FALSE)
  
}


find_range <- function(sc, feature){
  
  tmp <- sapply(sc, function(x) as.numeric(x@assays$Spatial[feature]))
  range(unlist(tmp))
  
}

vegfa_range <- find_range(sc, "VEGFA")
apoe_range <- find_range(sc, "APOE")

wesanderson::wes_palette("Zissou1", 10, type = "continuous")

for(i in 1:length(sc)){

  gg1 <- SpatialFeaturePlot(sc[[i]], features = c("VEGFA"), pt.size.factor = 1.6, 
      ncol = 1, crop = TRUE) & 
    scale_fill_gradientn(colours=rev(brewer.rdylbu(15)), limits=vegfa_range)
  gg2 <- SpatialFeaturePlot(sc[[i]], features = c("APOE"), pt.size.factor = 1.6, 
                            ncol = 1, crop = TRUE) & 
    scale_fill_gradientn(colours=rev(brewer.rdylbu(15)),
                      limits=vegfa_range)
  gg3 <- SpatialPlot(sc[[i]], features=c("MES_like"), pt.size.factor = 1.6, 
                     ncol = 1, crop = TRUE) &
    paletteer::scale_fill_paletteer_c("viridis::plasma", limits=c(0,1))
    
  ggsave(gg1, file=paste0("plots/", names(meta)[i], "_VEGFA.pdf"), 
         width=5, height=5)
  ggsave(gg2, file=paste0("plots/", names(meta)[i], "_APOE.pdf"), 
         width=5, height=5)
  ggsave(gg3, file=paste0("plots/", names(meta)[i], "_MES.pdf"), 
         width=5, height=5)
  
}
 
}

 

