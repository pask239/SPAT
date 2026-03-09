suppressMessages({
  library(DESpace)
  library(ggplot2)
  library(SpatialExperiment)
  library(reshape2)
  library(tidyverse)
  library(patchwork)
  library(splines)
  library(edgeR)
})
set.seed(123)

spe <- readRDS("data/spe_Banksy.rds")

condition <- factor(condition,levels = c("Wt","Tg"))
spe$condition <- condition

results <- dsp_test(spe = spe,
                    cluster_col = "Banksy_smooth",
                    sample_col = "sample_id",
                    condition_col = "condition",
                    test = 'QLF',
                    verbose = TRUE)

cluster_results <- individual_dsp(spe,
                                  cluster_col = "Banksy_smooth",
                                  sample_col = "sample_id",
                                  condition_col = "condition")

res <- cluster_results$`9` #|> head(n=4)

sample_ids <- c("9-1","9-3","9-5","9-4")

plots <- lapply(sample_ids, function(sample_id) {
  
  # Subset spe for each sample
  spe_j <- spe[, colData(spe)$sample_id == sample_id]
  
  # Create FeaturePlot for the sample
  plot <- FeaturePlot(spe_j,feature = "Layn", 
                      cluster_col = "Banksy_smooth",
                      cluster = '9',
                      diverging = TRUE,
                      point_size = 0.1,
                      linewidth = 0.6
                      ,assay.type = "counts"
                      ,platform = "OpenST") +
    theme(legend.position = "right",
          legend.key.size = unit(0.5, 'cm')) +
    labs(color = "")
  
  return(plot)
  
})

combined_plot <- wrap_plots(plots, ncol = 3) + 
  # common legend
  plot_layout(guides = 'collect')  

