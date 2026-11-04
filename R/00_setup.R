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
# Genes are kept if their mean count across the ROI type's ROIs is >= min_rowmean.
subset_roi <- function(meta, counts, roi_def, min_rowmean = cfg$thresholds$min_rowmean) {
  keep <- meta$tissue == roi_def$tissue & meta$region == roi_def$region
  if (!is.null(roi_def$pathology)) keep <- keep & meta$pathology == roi_def$pathology
  coldata <- as.data.frame(meta[keep, ])
  rownames(coldata) <- coldata$sample_id
  roi_counts <- counts[, coldata$sample_id]
  roi_counts <- roi_counts[rowMeans(roi_counts) >= min_rowmean, ]
  list(coldata = coldata, counts = roi_counts)
}
