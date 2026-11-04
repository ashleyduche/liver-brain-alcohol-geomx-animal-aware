suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(pheatmap)
  library(readxl)
  library(yaml)
  library(tidyverse)   # loaded last so dplyr verbs like slice, select and filter are not masked by Bioconductor
})

set.seed(1)

# Settings file, which defaults to the example data (set GEOMX_CONFIG to use the real data)
cfg <- read_yaml(Sys.getenv("GEOMX_CONFIG", "config.yml"))

# Create the output folders
for (d in c("processed", "qc", "dge_tables", "volcano", "heatmaps", "ranked_lists", "pseudobulk")) {
  dir.create(file.path(cfg$paths$results_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# Subset metadata and counts to one ROI type, with count columns in the same order as coldata rows
# Genes are kept when they are above each ROI's LOQ in enough of the ROI type's ROIs (GeoMx standard)
# An optional mean count filter can also be applied with min_rowmean
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
