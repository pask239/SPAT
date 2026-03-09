# script to perform spatial variable genes between spatial domains
library(DESpace)
library(SpatialExperiment)

spe <- readRDS("data/spe_Banksy.rds")

multi_results <- svg_test(spe = spe,
                          cluster_col = "Banksy_smooth",
                          sample_col = "sample_id",
                          replicates = TRUE)











