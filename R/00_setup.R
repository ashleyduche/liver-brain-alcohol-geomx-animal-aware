suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(pheatmap)
  library(readxl)
  library(yaml)
  library(tidyverse)   # loaded last so dplyr verbs (slice, select, filter) are not masked by Bioconductor
})

set.seed(1)
cfg <- read_yaml(Sys.getenv("GEOMX_CONFIG", "config.yml"))   # config_real.yml for the real data

for (d in c("processed", "qc", "dge_tables", "volcano", "heatmaps", "ranked_lists", "pseudobulk")) {
  dir.create(file.path(cfg$paths$results_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# Subset metadata + counts to one ROI type; columns of counts follow coldata rows.
# Genes kept if above each ROI's LOQ in >= gene_detection_rate of the ROI type's ROIs
# (GeoMx standard), and optionally if their mean count is >= min_rowmean.
subset_roi <- function(meta, counts, roi_def,
                       min_rowmean = cfg$thresholds$min_rowmean,
                       gene_detection_rate = cfg$geomx$gene_detection_rate) {
  keep <- meta$tissue == roi_def$tissue & meta$region == roi_def$region
  if (!is.null(roi_def$pathology)) keep <- keep & meta$pathology == roi_def$pathology
  coldata <- as.data.frame(meta[keep, ])
  rownames(coldata) <- coldata$sample_id
  roi_counts <- counts[, coldata$sample_id]
  above_loq <- sweep(roi_counts, 2, coldata$loq, ">")
  keep <- rowMeans(above_loq) >= gene_detection_rate & rowMeans(roi_counts) >= min_rowmean
  roi_counts <- roi_counts[keep, ]
  list(coldata = coldata, counts = roi_counts)
}
