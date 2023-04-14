setwd("/stornext/Bioinf/data/lab_brain_cancer/public_datasets/gbm_spatial/doi_10.5061_dryad.h70rxwdmj__v11/10XVisium_2/")

library(tidyverse)

all_res <- read.csv("Test_results.csv", header=TRUE)
all_res <- all_res %>% mutate(log10_PValue=-log10(FDR)) %>% mutate(direction=sign(zscore))
all_res <- all_res %>% mutate(Result = case_when(direction < 0 & FDR < 0.05 ~ "down", 
                                                 direction > 0 & FDR < 0.05 ~ "up",
                                                 FDR > 0.05 ~ "not significant"))
all_res$Result <- factor(all_res$Result, levels=c("up", "down", "not significant"))
all_res$test <- factor(all_res$test, levels=c(paste0(rep(c("MES", "AC", "OPC", "NPC"), 2), 
        "_", rep(c("Micro", "Macro"), each=4)), "Mono_MES"))
all_res <- all_res %>% group_by(test, Result,.drop = FALSE) %>% count() %>% filter(test != "VEGFA_MES")

gg1 <- ggplot(all_res, aes(x=Result, y=test)) + geom_tile(aes(fill=n)) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_distiller(palette = "RdPu", direction = 1)

ggsave(gg1, file="Results_associations_Ravi.pdf", height=5, width=3)




