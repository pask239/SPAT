library(SingleCellExperiment)
library(SpatialExperiment)
library(S4Vectors)

reticulate::use_condaenv("single_cell", required = TRUE)

sample_id <- list("9-1" = "9-1"
                  ,"9-3" = "9-3"
                  ,"9-5" = "9-5"
                  ,"9-4" = "9-4")

spe_list <- lapply(sample_id,function(sample){
  
  sce <- anndataR::read_h5ad(path = file.path("data/processed_data"
                                              ,sample
                                              ,"multimodal/stitched_segmented.h5ad")
                             ,as = "SingleCellExperiment")
  
  sce@int_colData$reducedDims@listData$spatial <- NULL 
  
  sce <- sce[,sce$cell_ID_mask != 0]
  
  colnames(sce) <- paste(sample,colnames(sce),sep = "_")
  
  spe <- SpatialExperiment(
    assays = assays(sce)
    , rowRanges = sce@rowRanges
    , colData = colData(sce)
    , spatialCoordsNames = c("x_pos","y_pos")
    ,sample_id = sample
  )
  
  names(assays(spe)) <- "counts"
  
  return(spe)
  
})

spe_list <- lapply(spe_list, function(spe){
  
  mcols(spe) <- NULL
  return(spe)
  
})

common_gene <- Reduce(intersect, lapply(spe_list, rownames))
spe_all <- do.call(cbind, lapply(spe_list, function(spe) spe[common_gene,]))
