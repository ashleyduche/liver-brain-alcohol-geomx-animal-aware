# Three models per ROI type, each comparing Alcohol-Treated with Control
#   Model 1  DESeq2_ROIs               every ROI is a sample, as in the published analysis
#   Model 2  voom_dupCor_ROIs          ROI level limma voom with animal as a correlated block
#   Model 3  pseudobulk_edgeR_animals  ROIs summed per animal and tested with the edgeR QL F test (n = 2 vs 2), the primary result
dge_dir <- file.path(cfg$paths$results_dir, "dge_tables")
summary_all <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def    <- cfg$rois[[roi_name]]
  roi_data   <- subset_roi(meta, counts, roi_def)
  coldata    <- roi_data$coldata
  roi_counts <- roi_data$counts

  cat(sprintf("\n=== %s: %d ROIs from %d animals ===\n", roi_name, nrow(coldata), n_distinct(coldata$animal)))

  # ---- Model 1: DESeq2 with each ROI as a sample (published approach) ----
  dds <- DESeqDataSetFromMatrix(countData = roi_counts, colData = coldata, design = ~ Group)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds)

  deseq2_rois <- as_tibble(as.data.frame(res), rownames = "gene") %>%
    transmute(gene, logFC = log2FoldChange, pvalue, FDR = padj)

  # ---- Model 2: ROI level voom with animal as a correlated block ----
  y <- DGEList(roi_counts)
  y <- y[filterByExpr(y, group = coldata$Group), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "upperquartile")
  design <- model.matrix(~ Group, data = coldata)

  # voom and duplicateCorrelation are run twice, as recommended by limma
  v <- voom(y, design)
  dupcor <- duplicateCorrelation(v, design, block = coldata$animal)
  v <- voom(y, design, block = coldata$animal, correlation = dupcor$consensus.correlation)
  dupcor <- duplicateCorrelation(v, design, block = coldata$animal)

  fit <- lmFit(v, design, block = coldata$animal, correlation = dupcor$consensus.correlation)
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  voom_dupcor <- tibble(gene = rownames(tt), logFC = tt$logFC, pvalue = tt$P.Value, FDR = tt$adj.P.Val)
  cat(sprintf("  within-animal ROI correlation: %.2f\n", dupcor$consensus.correlation))

  # ---- Model 3: pseudobulk per animal with edgeR quasi likelihood ----
  # Sum the ROI counts for each animal
  animal_f <- factor(coldata$animal)
  pb_counts <- roi_counts %*% model.matrix(~ 0 + animal_f)
  colnames(pb_counts) <- levels(animal_f)

  pb_group <- coldata %>% distinct(animal, Group) %>% arrange(match(animal, levels(animal_f)))
  stopifnot(identical(pb_group$animal, colnames(pb_counts)))

  yb <- DGEList(pb_counts, group = pb_group$Group)
  yb <- yb[filterByExpr(yb, group = pb_group$Group), , keep.lib.sizes = FALSE]
  yb <- calcNormFactors(yb, method = "upperquartile")
  design_pb <- model.matrix(~ Group, data = pb_group)

  yb <- estimateDisp(yb, design_pb, robust = TRUE)
  fit_pb <- glmQLFit(yb, design_pb, robust = TRUE)
  qlf <- glmQLFTest(fit_pb, coef = 2)
  tt_pb <- topTags(qlf, n = Inf, sort.by = "none")$table

  pseudobulk <- tibble(gene = rownames(tt_pb), logFC = tt_pb$logFC, pvalue = tt_pb$PValue, FDR = tt_pb$FDR)

  # Animal level logCPM, used for the heatmaps
  write_csv(as_tibble(cpm(yb, log = TRUE, prior.count = 2), rownames = "gene"),
            file.path(cfg$paths$results_dir, "pseudobulk", paste0(roi_name, "_logCPM.csv")))

  # Ranked gene list (signed -log10 p) for optional GSEA
  pseudobulk %>%
    mutate(rank_stat = sign(logFC) * -log10(pvalue)) %>%
    arrange(dplyr::desc(rank_stat)) %>%
    select(gene, rank_stat) %>%
    write_tsv(file.path(cfg$paths$results_dir, "ranked_lists", paste0(roi_name, "_pseudobulk.rnk")))

  # ---- Save tables, summary counts and volcano panels ----
  method_tables <- list(DESeq2_ROIs = deseq2_rois,
                        voom_dupCor_ROIs = voom_dupcor,
                        pseudobulk_edgeR_animals = pseudobulk)

  png(file.path(cfg$paths$results_dir, "volcano", paste0(roi_name, "_Volcano.png")), width = 2400, height = 850, res = 200)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 4.5, 1))

  for (method in names(method_tables)) {
    tab <- method_tables[[method]]
    write_csv(tab, file.path(dge_dir, paste0(roi_name, "_", method, ".csv")))

    summary_all <- bind_rows(summary_all, tibble(
      roi            = roi_name,
      method         = method,
      genes_tested   = sum(!is.na(tab$pvalue)),
      n_nominal_p    = sum(tab$pvalue < cfg$thresholds$nominal_p, na.rm = TRUE),
      n_fdr          = sum(tab$FDR <= cfg$thresholds$fdr, na.rm = TRUE),
      n_deg          = sum(tab$FDR <= cfg$thresholds$fdr & abs(tab$logFC) > cfg$thresholds$abs_log2fc, na.rm = TRUE),
      min_fdr        = min(tab$FDR, na.rm = TRUE)
    ))

    # Red for up and blue for down when a gene passes FDR and the fold change cutoff
    is_up   <- tab$FDR <= cfg$thresholds$fdr & tab$logFC >  cfg$thresholds$abs_log2fc
    is_down <- tab$FDR <= cfg$thresholds$fdr & tab$logFC < -cfg$thresholds$abs_log2fc
    pt_col <- case_when(is_up ~ "red", is_down ~ "blue", TRUE ~ "grey70")

    # Draw the grey genes first so the colored genes sit on top
    draw_order <- order(pt_col != "grey70")

    plot(tab$logFC[draw_order], -log10(tab$pvalue[draw_order]), col = pt_col[draw_order], pch = 20,
         xlab = "log2 fold change (Alcohol vs Control)", ylab = "-log10(p value)",
         main = paste0(method, "\n", sum(tab$FDR <= cfg$thresholds$fdr, na.rm = TRUE), " genes FDR <= ", cfg$thresholds$fdr,
                       "\n", sum(is_up, na.rm = TRUE), " up (red) and ", sum(is_down, na.rm = TRUE),
                       " down (blue) past |log2FC| ", cfg$thresholds$abs_log2fc),
         cex.main = 0.85)
    abline(v = c(-cfg$thresholds$abs_log2fc, cfg$thresholds$abs_log2fc), lty = 2, col = "grey40")

    # The horizontal line marks the largest p value that still passes FDR, when any gene passes
    if (any(tab$FDR <= cfg$thresholds$fdr, na.rm = TRUE)) {
      abline(h = -log10(max(tab$pvalue[which(tab$FDR <= cfg$thresholds$fdr)])), lty = 2, col = "grey40")
    }
  }
  dev.off()
}

write_csv(summary_all, file.path(cfg$paths$results_dir, "dge_summary.csv"))
print(as.data.frame(summary_all))
