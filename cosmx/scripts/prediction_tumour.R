
library(Seurat)
library(SeuratObject)
library(SingleR)
library(ggplot2)
library(SingleCellExperiment)
library(pals)
library(intrinsicDimension)
library(reticulate)
library(readxl)
library(viridis)
library(dplyr)
library(magrittr)
library(scater)
library(scuttle)
library(scran)

sce_tumour <- readRDS("data/sce_tumour_jj98.rds")

sc_couturier <- readRDS("../data/Couturier_dataset_tumor.RDS")
rownames(sc_couturier) <- rowData(sc_couturier)$Symbol
sc_couturier <- logNormCounts(sc_couturier)

available <- intersect(rownames(sc_couturier), rownames(sce_tumour))

pred_labels <- SingleR(test=sce_tumour[available,], ref=sc_couturier[available,], 
                       labels=sc_couturier$cluster, de.method="wilcox")
table(pred_labels$labels)

saveRDS(pred_labels, file="data/prediction_couturier.rds")