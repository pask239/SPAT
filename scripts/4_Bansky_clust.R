library(SpatialExperiment)
library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(harmony)
library(ggplot2)
library(dplyr)

spe <- readRDS("data/spe_staggered_cordinates.rds")-
sample_ids <- unique(spe$sample_id)
spe_list <- lapply(sample_ids, function(x) spe[, spe$sample_id == x])

# Normalize data ---------------------------------------------------------------
seu_list <- lapply(spe_list, function(x) {
  x <- as.Seurat(x, data = NULL)
  scale_factor <-  median(colSums(x@assays[["originalexp"]]@counts))
  NormalizeData(x, scale.factor = scale_factor, normalization.method = 'RC')
})

# Compute HVGs -----------------------------------------------------------------
hvgs <- lapply(seu_list, function(x) {
  VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)

aname <- "normcounts"
spe_list <- Map(function(spe, seu) {
  assay(spe, aname, withDimnames=FALSE) <- GetAssayData(seu)
  spe[hvgs,]
}, spe_list, seu_list)

invisible(gc())

# Banksy cluster ---------------------------------------------------------------
compute_agf <- FALSE
k_geom <- 18
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, 
                   compute_agf = compute_agf, k_geom = k_geom)

spe_joint <- do.call(cbind, spe_list)
seu_joint <- as.Seurat(spe_joint, data = NULL)

group_name = c("sample_id")
dimx = "sdimx"; dimy = "sdimy"
col_data <- cbind(colData(spe_joint), spatialCoords(spe_joint))
seu_joint <- AddMetaData(seu_joint, metadata = as.data.frame(col_data))

# Fit Banksy
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)

seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 20)

seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='pca', reduction.save='harmony')

seu_joint = FindNeighbors(seu_joint, dims = 1:20, reduction = 'harmony')

seu_joint = FindClusters(seu_joint, resolution = 0.3, graph.name = 'BANKSY_snn')

# Banksy smooth
# spe_joint = spe
# spe_joint$BANKSY_snn_res0.1 <- seu_joint@meta.data$BANKSY_snn_res.0.1
# spe_list <- lapply(sample_ids, function(x) spe_joint[, spe_joint$sample_id == x])
# spe_list <- lapply(spe_list, smoothLabels, cluster_names = "BANKSY_snn_res0.1"
#                    , k = 6L, verbose = FALSE)
# names(spe_list) <- paste0("sample_", sample_ids)
# spe_Banksy <- do.call(cbind, spe_list)

spe_joint = spe
spe_joint$BANKSY_snn_res0.3 <- seu_joint@meta.data$BANKSY_snn_res.0.3
spe_list <- lapply(sample_ids, function(x) spe_joint[, spe_joint$sample_id == x])
spe_list <- lapply(spe_list, smoothLabels, cluster_names = "BANKSY_snn_res0.3"
                   , k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_ids)
spe_Banksy <- do.call(cbind, spe_list)


# create final spe object -------------------------------------------------

colData(spe_Banksy) <- cbind(colData(spe_Banksy), 
                             spatialCoords(spe_Banksy))

colData(spe_Banksy) <- subset(colData(spe_Banksy), select = c(sample_id
                                                              ,condition
                                                              ,time
                                                              ,BANKSY_snn_res0.3
                                                              ,BANKSY_snn_res0.3_smooth
                                                              ,sdimx
                                                              ,sdimy))

colnames(colData(spe_Banksy))[4:5] <- c(
  "Banksy", "Banksy_smooth" 
)

spe <- spe_Banksy; rm(spe_Banksy)

saveRDS(spe, file = "data/spe_Banksy.rds")

CD <- colData(spe) |> as.data.frame()
ggplot(CD, aes(x = sdimx, y = sdimy, color = factor(Banksy_smooth))) +
  geom_point(size = 0.25) +
  facet_wrap(~sample_id,scales = "free") +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(color = NULL, title = "Banksy Spatial Clusters")

