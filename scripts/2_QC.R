library(scuttle)
library(scater)
library(ggspavis)
library(SpotSweeper)
library(patchwork)

source("scripts/functions/detect_global_outliers.R")

# -------------------------------------------------------------------------
spe_list <- readRDS("data/spe_list.rds")

# ------------------------- ADD per cell QC ------------------------------------
spe_list <- lapply(spe_list, function(spe){
  
  is_mito <- grepl("(^MT-)|(^mt-)", rownames(spe))
  table(is_mito)
  spe <- addPerCellQC(spe, subsets=list(mito = is_mito))
  return(spe)
  
})

# Global outlier ---------------------------------------------------------------
# select QC thresholds for library size, 
# detected features & mito. proportion

spe_list <- lapply(spe_list, detect_global_outlier,perc_lib_size = 0.2
                   ,perc_detected = 0.2,perc_mito_percent = 0.99)

lapply(spe_list, function(spe){
  
  par(mfrow=c(1, 3))
  thresh_qc_libsize <- quantile(spe$sum,probs = 0.2)
  hist(spe$sum, xlab="sum", main="UMIs per cells")
  abline(v = thresh_qc_libsize, col = "red", lwd = 2)
  
  thresh_qc_detected <- quantile(spe$detected,probs = 0.2,na.rm= T)
  hist(spe$detected, xlab="detected genes", main="Genes per cells")
  abline(v = thresh_qc_detected , col = "red", lwd = 2)
  
  thresh_qc_mito_percent <- quantile(spe$subsets_mito_percent,probs = 0.99,na.rm= T)
  hist(spe$subsets_mito_percent, xlab="pct mito", main="Percent mito UMIs")
  abline(v = thresh_qc_mito_percent, col = "red", lwd = 2)
  
  
})


lapply(spe_list,function(spe) {
  
  p1 <- plotObsQC(spe, 
                  plot_type="spot", annotate="qc_lib_size") + 
    ggtitle("Library size")
  
  p2 <- plotObsQC(spe, 
                  plot_type="spot", annotate="qc_detected") + 
    ggtitle("Detected genes")
  
  p3 <- plotObsQC(spe, 
                  plot_type="spot", annotate="qc_mito_prop") + 
    ggtitle("Mitocondrial proportion")
  
  p1 | p2 | p3
  
})

# Local Outliers ---------------------------------------------------------------
spe_list <- lapply(spe_list, function(spe){

  spe <- localOutliers(spe, metric="sum", direction="lower", log=TRUE
                       ,workers = 24)
  spe <- localOutliers(spe, metric="detected", direction="lower", log=TRUE
                       ,workers = 24)
  spe <- localOutliers(spe, metric="subsets_mito_percent", direction="higher"
                       , log=FALSE,workers = 24)

  return(spe)

})

lapply(spe_list, function(spe){

  p1 <- plotCoords(spe,
                   annotate="sum_log",in_tissue = NULL) +
    ggtitle("log2(Library Size)")

  p2 <- plotObsQC(spe,
                  plot_type="spot", in_tissue= NULL,
                  annotate="sum_outliers", point_size=0.2) +
    ggtitle("Local Outliers (Library Size)")

  # spot plot of log-transformed detected genes
  p3 <- plotCoords(spe,
                   annotate="detected_log",in_tissue = NULL) +
    ggtitle("log2(Detected)")

  p4 <- plotObsQC(spe,
                  plot_type="spot", in_tissue=NULL,
                  annotate="detected_outliers", point_size=0.2) +
    ggtitle("Local Outliers (Detected)")

  # spot plot of mitochondrial proportion
  p5 <- plotCoords(spe,
                   annotate="subsets_mito_percent",in_tissue = NULL) +
    ggtitle("Mito Proportion")

  p6 <- plotObsQC(spe,
                  plot_type="spot", in_tissue=NULL,
                  annotate="subsets_mito_percent_outliers", point_size=0.2) +
    ggtitle("Local Outliers (Mito Prop)")

  # plot using patchwork
  (p1 / p2) | (p3 / p4) | (p5 / p6)

})

spe_list <- lapply(spe_list, function(spe){
  
  spe$global_outliers <- spe$qc_lib_size | spe$qc_detected | spe$qc_mito_prop 
  
  spe$local_outliers <-
    spe$sum_outliers |
    spe$detected_outliers |
    spe$subsets_mito_percent_outliers

  print(rbind(global=table(spe$global_outliers)
              ,local=table(spe$local_outliers
              )
  )
  )
  
  return(spe)
  
})

spe_list <- lapply(spe_list, function(spe){
  
  spe$discard <- 
    spe$global_outliers |
    spe$local_outliers
  
  spe <- spe[, !spe$discard]
  
  return(spe)
  
})

# Discard lowly abundant genes, which were detected in less than 20 spots.
qc_low_genes <- sapply(spe_list, function(spe){
  
  qc_low_gene <- rowSums(assays(spe)$counts) > 0
  # # Select QC threshold for lowly expressed genes: at least 20 non-zero cells:
  # qc_low_gene <- rowSums(assays(spe)$counts > 0) >= 20
  return(rownames(spe)[qc_low_gene])
  
})


# Find common genes across all samples
common_genes <- Reduce(intersect, qc_low_genes)

spe_list <- lapply(spe_list, function(spe){
  
  spe <- spe[common_genes, ]
  
  return(spe)
  
})

lapply(spe_list, function(spe){
  
  plotCoords(spe,in_tissue = NULL)
  
})


# Hangnails artifact -----------------------------------------------------------
# spe_list <- lapply(spe_list, function(spe){
#   
#   spe <- localVariance(spe, 
#                        n_neighbors=36, 
#                        metric="subsets_mito_percent", 
#                        name="local_mito_variance_k36")
#   return(spe)
#   
# })
# 
# lapply(spe_list, function(spe){
#   
#   plotCoords(spe, annotate="subsets_mito_percent",in_tissue = NULL)
#   
# })
# 
# lapply(spe_list, function(spe){
#   
#   plotCoords(spe, annotate="local_mito_variance_k36", point_size=1
#              ,in_tissue = NULL)
#   
# })
# 
# lapply(spe_list, function(spe){
#   
#   hangnail_sample <- data.frame(x=spe$local_mito_variance_k36)
#   
#   ggplot(hangnail_sample, aes(x)) + ggtitle("Hangnail Sample") +
#     geom_density(fill="gray", alpha=0.5)
#   
# })

spe <- do.call(cbind,spe_list)

spe$condition <- ifelse(stringr::str_detect(spe$sample_id,"9-1|9-3"),"Wt","Tg")
spe$time <- ifelse(stringr::str_detect(spe$sample_id,"9-3|9-4"),"3m","1m")

saveRDS(spe,file = "data/spe_qc.rds")

