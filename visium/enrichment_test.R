library(igraph)
library(ggplot2)
library(spatstat)
library(SpatialExperiment)
library(parallel)


setwd("gbm_spatial/doi_10.5061_dryad.h70rxwdmj__v11/10XVisium_2/")


find_observed <- function(x, df, a="TAM", b="MES_like"){
  
  barcodes <- names(x)[x==1]
  sum((df[barcodes, a]+df[barcodes, b])==2)
  
}

permute <- function(adj, df, a="TAM", b="MES_like"){
  
  df <- df[, c(a,b)]
  df <- apply(df, 2, function(x) sample(x))
  colnames(df) <- c(a,b)
  rownames(df) <- colnames(adj)
  
  exp <- apply(adj, 1, function(x) find_observed(x, df, a=a, b=b))
  sum(exp)
  
}

enrichment <- function(spe10x, a="TAM", b="MES_like", n_perm=1000, max_dist=100, cores=10){
  
  dist_spe10x <- dist(spatialCoords(spe10x), method = "euclidean", diag = FALSE, upper = FALSE)
  adj <- as.matrix(dist_spe10x)
  adj <-  apply(adj, 2, function(x) as.numeric(x<=max_dist))
  
  df <- colData(spe10x)
  
  obs <- apply(adj, 1, function(x) find_observed(x, colData(spe10x), a=a, b=b))
  obs <- sum(obs)

  exp <- mclapply(1:n_perm, function(x) permute(adj, df, a=a, b=b), mc.cores=cores)

  zscore <- (obs-mean(unlist(exp)))/(sd(unlist(exp)))
  zscore
  
}

all_zscores <- list()
all_sums <- list()

meta <- readRDS("Features_all_new.RDS")
names(meta) <- sapply(meta, function(x) x$sample[1])

all <- length(meta)
#all <- 2

for(i in 1:all){

  spe10x <- read10xVisium( paste0("#UKF", names(meta)[i], "_ST/outs"), 
    "slice1", type = "sparse", data = "filtered",
     images = "lowres", load = FALSE)

  spe10x <- spe10x[, colnames(spe10x) %in% meta[[i]]$barcodes]
  colData(spe10x) <- cbind(colData(spe10x),  
    meta[[i]][match(colnames(spe10x), meta[[i]]$barcodes),])

  spe10x$MES_like <- spe10x$MES_like_hypoxia_independent + spe10x$MES_like_hypoxia_MHC
  spe10x$MES_like <- spe10x$MES_like*spe10x$Nr_of_cells
  spe10x$MES_like <- as.numeric(spe10x$MES_like>1)
  
  spe10x$Mono <- spe10x$Mono_anti_infl + spe10x$Mono_hypoxia + 
    spe10x$Mono_naive 
  spe10x$Mono <- spe10x$Mono*spe10x$Nr_of_cells
  spe10x$Mono<- as.numeric(spe10x$Mono>1)
  
  spe10x$Macro <- spe10x$TAM_BDM_anti_infl +  spe10x$TAM_BDM_hypoxia_MES +        
    spe10x$TAM_BDM_INF + spe10x$TAM_BDM_MHC
  spe10x$Macro <- spe10x$Macro*spe10x$Nr_of_cells
  spe10x$Macro<- as.numeric(spe10x$Macro>1)
  
  spe10x$Micro <- spe10x$TAM_MG_aging_sig + spe10x$TAM_MG_pro_infl_I +
    spe10x$TAM_MG_pro_infl_II + spe10x$TAM_MG_prolif
  spe10x$Micro <- spe10x$Micro*spe10x$Nr_of_cells
  spe10x$Micro<- as.numeric(spe10x$Micro>1)
  
  spe10x$AC_like <- spe10x$AC_like_Prolif + spe10x$AC_like 
  spe10x$AC_like <- spe10x$AC_like*spe10x$Nr_of_cells
  spe10x$AC_like<- as.numeric(spe10x$AC_like>1)
  
  spe10x$OPC_like <- spe10x$OPC_like + spe10x$OPC_like_Prolif
  spe10x$OPC_like <- spe10x$OPC_like*spe10x$Nr_of_cells
  spe10x$OPC_like <- as.numeric(spe10x$OPC_like>1)
  
  spe10x$NPC_like <- spe10x$NPC_like_neural + spe10x$NPC_like_OPC + spe10x$NPC_like_Prolif
  spe10x$NPC_like <- spe10x$NPC_like*spe10x$Nr_of_cells
  spe10x$NPC_like <- as.numeric(spe10x$NPC_like>1)
  
  spe10x$VEGFA <- counts(spe10x)["ENSG00000112715",]>1
  
  all_sums[[names(meta)[i]]] <- c(MES_like=sum(spe10x$MES_like), 
    Mono=sum(spe10x$Mono), Macro=sum(spe10x$Macro),
    NPC_like=sum(spe10x$NPC_like),
        VEGFA=sum(spe10x$VEGFA))
  
  all_tests <- list(
    Mono_MES = enrichment(spe10x, a="Mono", b="MES_like", n_perm=1000, cores=21),
    NPC_Micro = enrichment(spe10x, a="NPC_like", b="Micro", n_perm=1000, cores=21),
    OPC_Micro = enrichment(spe10x, a="OPC_like", b="Micro", n_perm=1000, cores=21),
    AC_Micro = enrichment(spe10x, a="AC_like", b="Micro", n_perm=1000, cores=21),
    MES_Micro = enrichment(spe10x, a="MES_like", b="Micro", n_perm=1000, cores=21),
    VEGFA_MES = enrichment(spe10x, a="VEGFA", b="MES_like", n_perm=1000, cores=21),
    NPC_Macro = enrichment(spe10x, a="NPC_like", b="Macro", n_perm=1000, cores=21),
    OPC_Macro = enrichment(spe10x, a="OPC_like", b="Macro", n_perm=1000, cores=21),
    AC_Macro = enrichment(spe10x, a="AC_like", b="Macro", n_perm=1000, cores=21),
    MES_Macro = enrichment(spe10x, a="MES_like", b="Macro", n_perm=1000, cores=21)
    )

  all_zscores[[names(meta)[i]]] <- all_tests
  
  rm(spe10x)
  gc()
  
}


saveRDS(all_zscores, file="all_zscores.rds")
saveRDS(all_sums, file="all_sums.rds")

all_sums <- readRDS("all_sums.rds")  
all_zscores <- readRDS("all_zscores.rds")


all_zscores <- lapply(all_zscores, function(x) unlist(x))
all_pvalues <- lapply(all_zscores, function(x) pnorm(abs(x), lower.tail=FALSE)*2)
all_pvalues <- lapply(all_pvalues, function(x) p.adjust(x, method="fdr"))

all_res <- list()

for(i in 1:length(all_zscores)){
  
  all_res[[i]] <- data.frame(test=names(all_zscores[[i]]),
    sample=names(all_zscores)[i], zscore=all_zscores[[i]],
                             FDR=all_pvalues[[i]]) 
  
}

all_res <- do.call(rbind, all_res)
rownames(all_res)<-NULL

write.csv(all_res, file="Test_results.csv", quote=FALSE)



