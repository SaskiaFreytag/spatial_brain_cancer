library(SpatialExperiment)
library(SPOTlight)
library(ggplot2)
library(scran)
library(scater)
library(spatstat)
library(paletteer)
library(wesanderson)

setwd("/stornext/Bioinf/data/lab_brain_cancer/public_datasets/gbm_spatial/doi_10.5061_dryad.h70rxwdmj__v11/10XVisium_2/")

sce <- read10xVisium(
  samples = "#UKF243_T_ST/outs",
  sample_id = paste0("sample", sprintf("%02d", seq_along(samples))),
  type = c("HDF5", "sparse"),
  data = c("filtered"),
  images = "lowres",
  load = FALSE
)
meta <- readRDS("Features_all_new.RDS")
names(meta) <- sapply(meta, function(x) x$sample[1])

colData(sce) <- cbind(colData(sce) , meta[["243_T"]][
  match(colnames(sce), meta[["243_T"]]$barcodes),])

sce$AC_like <- sce$AC_like_Prolif + sce$AC_like 
sce$AC_like <- sce$AC_like*sce$Nr_of_cells

sce$OPC_like <- sce$OPC_like + sce$OPC_like_Prolif
sce$OPC_like <- sce$OPC_like*sce$Nr_of_cells


sce$NPC_like <- sce$NPC_like_neural + sce$NPC_like_OPC + sce$NPC_like_Prolif
sce$NPC_like <- sce$NPC_like*sce$Nr_of_cells


sce$MES_like <- sce$MES_like_hypoxia_independent + sce$MES_like_hypoxia_MHC
sce$MES_like <- sce$MES_like*sce$Nr_of_cells


xy <- spatialCoords(sce)
sce$x <- xy[, 1]
sce$y <- xy[, 2]

gg_mes <- ggcells(sce, aes(x, y)) +
  geom_point(shape = 1, size = 1,colour = "black") +
  geom_point(aes(color = MES_like), size=0.9) +
  scale_colour_gradient( low = "white", high = "#E2726E") +
  coord_fixed() +
  theme_void()

ggsave(gg_mes, file="plots1/MES_UKF243_T.jpg", height=3, width=3)

gg_ac <- ggcells(sce, aes(x, y)) +
  geom_point(shape = 1, size = 1,colour = "black")  +
  geom_point(aes(color = AC_like), size=0.9) +
  scale_colour_gradient( low = "white", high = "#8FBAEE") +
  coord_fixed() +
  theme_void()

ggsave(gg_ac, file="plots1/AC_UKF243_T.jpg", height=3, width=3)

gg_opc <- ggcells(sce, aes(x, y)) +
  geom_point(shape = 1, size = 1,colour = "black")  +
  geom_point(aes(color = OPC_like), size=0.9) +
  scale_colour_gradient( low = "white", high = "#B4F0D5") +
  coord_fixed() +
  theme_void()

ggsave(gg_opc, file="plots1/OPC_UKF243_T.jpg", height=3, width=3)

gg_npc <- ggcells(sce, aes(x, y)) +
  geom_point(shape = 1, size = 1,colour = "black")  +
  geom_point(aes(color = NPC_like), size=0.9) +
  scale_colour_gradient( low = "white", high = "#C6A8F6") +
  coord_fixed() +
  theme_void()

ggsave(gg_npc, file="plots1/NPC_UKF243_T.jpg", height=3, width=3)

