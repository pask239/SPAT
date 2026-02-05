# Create Spe object with staggered coordinates

library(data.table)
library(SpatialExperiment)

spe_list <- readRDS("data/spe_list_qc.rds")

spe <- do.call(cbind,spe_list)
spe$condition <- ifelse(stringr::str_detect(spe$sample_id,"9-1|9-3"),"Wt","Tg")
spe$time <- ifelse(stringr::str_detect(spe$sample_id,"9-5|9-4"),"3m","1m")

# Stagger spatial coordinates
# stagger the spatial coordinates across the samples so that spots from different samples do not overlap.
locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = as.numeric(factor(spe$sample_id)))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs

saveRDS(spe,"data/spe_staggered_cordinates.rds")















