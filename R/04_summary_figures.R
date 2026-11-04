# Summary barplot of significant genes per model and animal level heatmaps of the top genes
summary_all <- read_csv(file.path(cfg$paths$results_dir, "dge_summary.csv"), show_col_types = FALSE)

# One row per model and one column per ROI type
bar_mat <- summary_all %>%
  select(roi, method, n_fdr) %>%
  pivot_wider(names_from = roi, values_from = n_fdr) %>%
  column_to_rownames("method") %>%
  as.matrix()
bar_mat <- bar_mat[, names(cfg$rois)]

png(file.path(cfg$paths$results_dir, "fdr_hits_by_method.png"), width = 2000, height = 1000, res = 200)
par(mar = c(10, 5, 3, 1))
bp <- barplot(bar_mat, beside = TRUE, las = 2, col = c("grey60", "skyblue3", "darkorange"),
              ylab = paste0("Genes with FDR < ", cfg$thresholds$fdr),
              main = "Genes with FDR < 0.05 by model and ROI type")
text(bp, bar_mat, labels = bar_mat, pos = 3, cex = 0.7, xpd = TRUE)
legend("topleft", legend = rownames(bar_mat), fill = c("grey60", "skyblue3", "darkorange"), bty = "n", cex = 0.8)
dev.off()

# Heatmap of the top genes for each ROI type
for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]

  tab <- read_csv(file.path(cfg$paths$results_dir, "dge_tables", paste0(roi_name, "_pseudobulk_edgeR_animals.csv")),
                  show_col_types = FALSE)
  logcpm <- read_csv(file.path(cfg$paths$results_dir, "pseudobulk", paste0(roi_name, "_logCPM.csv")),
                     show_col_types = FALSE)

  top_genes <- tab %>% arrange(pvalue) %>% slice_head(n = cfg$thresholds$top_genes_heatmap) %>% pull(gene)
  mat <- logcpm %>% filter(gene %in% top_genes) %>% column_to_rownames("gene") %>% as.matrix()

  # Look up each animal's group and put the controls first
  animal_group <- meta %>% distinct(animal, Group) %>% arrange(Group, animal)
  mat <- mat[, animal_group$animal]
  col_anno <- data.frame(Group = animal_group$Group, row.names = animal_group$animal)

  pheatmap(mat, scale = "row", cluster_cols = FALSE, annotation_col = col_anno,
           main = paste0(roi_def$label, ": top ", nrow(mat), " genes by nominal p (animal-level)"),
           filename = file.path(cfg$paths$results_dir, "heatmaps", paste0(roi_name, "_top_genes.png")),
           width = 6, height = 8)
}
