library(ggspavis)

sample_ids <- c("9-1","9-3","9-5","9-4")

plots <- lapply(sample_ids, function(sample_id) {
  
  # Subset spe for each sample
  spe_j <- spe[, colData(spe)$sample_id == sample_id]
  
  plotCoords(spe = spe_j,in_tissue = NULL,annotate = "Aldh1a1")
  
  
  
})  

plots_2 <- lapply(sample_ids, function(sample_id) {
  
  # Subset spe for each sample
  spe_j <- spe[, colData(spe)$sample_id == sample_id]
  
  plotCoords(spe = spe_j,in_tissue = NULL,annotate = "Calb1")
  
  
  
})

plots_3 <- lapply(sample_ids, function(sample_id) {
  
  # Subset spe for each sample
  spe_j <- spe[, colData(spe)$sample_id == sample_id]
  
  plotCoords(spe = spe_j,in_tissue = NULL,annotate = "Calb2")
  
})  
