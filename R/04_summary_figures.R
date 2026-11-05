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

  # Z-score each gene across animals so the color shows up or down relative to that gene's mean
  z_mat <- t(scale(t(mat)))

  # Group bar on top in the same colors as the PCA plots
  top_anno <- HeatmapAnnotation(Group = animal_group$Group,
                                col = list(Group = c("Control" = "forestgreen", "Alcohol-Treated" = "firebrick")))

  # Leaving col unset keeps the ComplexHeatmap default of blue, white and red
  ht <- Heatmap(z_mat, name = "Row z-score",
                top_annotation = top_anno,
                cluster_columns = FALSE,
                column_split = animal_group$Group,
                row_names_gp = gpar(fontsize = 8),
                column_names_rot = 0,
                column_names_centered = TRUE,
                column_title = paste0(roi_def$label, "\ntop ", nrow(z_mat), " genes by pseudobulk p value"),
                column_title_gp = gpar(fontsize = 10, fontface = "bold"))

  png(file.path(cfg$paths$results_dir, "heatmaps", paste0(roi_name, "_top_genes.png")),
      width = 1200, height = 1800, res = 200)
  draw(ht, merge_legend = TRUE)
  dev.off()
}
