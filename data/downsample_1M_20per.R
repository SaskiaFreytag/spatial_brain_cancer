library(reticulate)
reticulate::use_virtualenv("r-reticulate")
library(scater)

setwd("/stornext/Bioinf/data/lab_brain_cancer/public_datasets/sc_gbm_big/")

sc <- import("scanpy")
ad <- import("anndata", convert = FALSE)
np <- import("numpy") 
pd <- import("pandas") 
geosketch <- import("geosketch", convert=FALSE)

# malignant cell - M                 macrophage - I                    
# mural cell  - M                   dendritic cell - I                
# microglial cell - I                 monocyte - I                    
# oligodendrocyte                endothelial cell - M              
# mature T cell - I                  oligodendrocyte precursor cell - M
# mast cell - I                     B cell - I                        
# plasma cell  - I                   natural killer cell -I           
# astrocyte - M                     radial glial cell - M             
# neuron - M                         neutrophil - I     

adata <- sc$read("local.h5ad", backed="r+")
X_umap <- adata$obsm["X_umap"]
sketch.size <- as.integer(ceiling(adata$n_obs*0.2))
sketch.indices <- geosketch$gs(X_umap, sketch.size, one_indexed = FALSE)
sketch.indices <- py_to_r(sketch.indices)
sketch.indices <- unlist(sketch.indices)

adata <- adata[sketch.indices]

sce <- SingleCellExperiment(
  assays      = list(logcounts = t(adata$X)),
  colData     = adata$obs,
  rowData     = adata$var
)

reducedDim(sce, "UMAP") <- X_umap[sketch.indices,]

saveRDS(sce, file="sce_downsampled.rds")

all_mal <- c('malignant cell', 'mural cell', 'oligodendrocyte', 'endothelial cell', 
             'oligodendrocyte precursor cell', 'astrocyte', 'radial glial cell', 'neuron')
all_imm <- c('microglial cell', 'macrophage', 'dendritic cell', 'monocyte', 'mature T cell', 
             'mast cell', 'B cell', 'plasma cell', 'natural killer cell', 'neutrophil')

sce_mal <- sce[, sce$cell_type %in% all_mal]
saveRDS(sce_mal, file="sce_downsampled_malignant.rds")
rm(sce_mal)
gc()

sce_imm <- sce[, sce$cell_type %in% all_imm]
saveRDS(sce_imm, file="sce_downsampled_immune.rds")