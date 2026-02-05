# function to detect global outlier

detect_global_outlier <- function(spe
                                  ,perc_lib_size
                                  ,perc_detected
                                  ,perc_mito_percent) {
  
  thresh_qc_libsize <- quantile(spe$sum,probs = perc_lib_size)
  cat("Threshold lib size:", thresh_qc_libsize,"\n")
  
  thresh_qc_mito_percent <- quantile(spe$subsets_mito_percent,probs = perc_mito_percent
                                     ,na.rm = T)
  cat("Threshold qc mito percent: ", thresh_qc_mito_percent,"\n")
  
  thresh_qc_detected <- quantile(spe$detected, probs = perc_detected)
  cat("Threshold detected genes :", thresh_qc_detected,"\n")
  
  
  spe$qc_lib_size <- spe$sum < thresh_qc_libsize
  spe$qc_detected <- spe$detected < thresh_qc_detected
  spe$qc_mito_prop <- spe$subsets_mito_percent > thresh_qc_mito_percent
  
  qc <- grep("^qc", names(colData(spe)))
  
  print(sapply(colData(spe)[qc], table))
  
  return(spe)
  
  
}
















