# Three models per ROI type, same genes/contrast (Alcohol-Treated vs Control):
#   1. naive_DESeq2_ROIs        every ROI is a sample (the original analysis; pseudoreplicated)
#   2. voom_dupCor_ROIs         ROI-level limma-voom, animal as a correlated block
#   3. pseudobulk_edgeR_animals ROIs summed per animal, edgeR QL F-test (n = 2 vs 2); PRIMARY result
dge_dir <- file.path(cfg$paths$results_dir, "dge_tables")
summary_all <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def  <- cfg$rois[[roi_name]]
  roi_data <- subset_roi(meta, counts, roi_def)
  coldata  <- roi_data$coldata
  roi_counts <- roi_data$counts
  cat(sprintf("\n=== %s: %d ROIs from %d animals ===\n", roi_name, nrow(coldata), n_distinct(coldata$animal)))

  # ---- 1. naive DESeq2, ROIs as independent samples ----
  dds <- DESeqDataSetFromMatrix(countData = roi_counts, colData = coldata, design = ~ Group)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds)
  naive <- as_tibble(as.data.frame(res), rownames = "gene") %>%
    transmute(gene, logFC = log2FoldChange, pvalue, FDR = padj)

  # ---- 2. ROI-level voom with animal as random-like block ----
  y <- DGEList(roi_counts)
  y <- y[filterByExpr(y, group = coldata$Group), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "upperquartile")
  design <- model.matrix(~ Group, data = coldata)
  v <- voom(y, design)
  dupcor <- duplicateCorrelation(v, design, block = coldata$animal)
  v <- voom(y, design, block = coldata$animal, correlation = dupcor$consensus.correlation)
  dupcor <- duplicateCorrelation(v, design, block = coldata$animal)
  fit <- lmFit(v, design, block = coldata$animal, correlation = dupcor$consensus.correlation)
  fit <- eBayes(fit, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  voom_dupcor <- tibble(gene = rownames(tt), logFC = tt$logFC, pvalue = tt$P.Value, FDR = tt$adj.P.Val)
  cat(sprintf("  within-animal ROI correlation: %.2f\n", dupcor$consensus.correlation))

  # ---- 3. pseudobulk per animal, edgeR quasi-likelihood ----
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

  write_csv(as_tibble(cpm(yb, log = TRUE, prior.count = 2), rownames = "gene"),
            file.path(cfg$paths$results_dir, "pseudobulk", paste0(roi_name, "_logCPM.csv")))

  # Ranked list (signed -log10 p) for optional GSEA; ties/zeros are fine for pre-ranked tools
  pseudobulk %>%
    mutate(rank_stat = sign(logFC) * -log10(pvalue)) %>%
    arrange(desc(rank_stat)) %>%
    select(gene, rank_stat) %>%
    write_tsv(file.path(cfg$paths$results_dir, "ranked_lists", paste0(roi_name, "_pseudobulk.rnk")))

  # ---- save tables, summary counts, volcano panels ----
  method_tables <- list(naive_DESeq2_ROIs = naive,
                        voom_dupCor_ROIs = voom_dupcor,
                        pseudobulk_edgeR_animals = pseudobulk)

  png(file.path(cfg$paths$results_dir, "volcano", paste0(roi_name, "_Volcano.png")), width = 2400, height = 850, res = 200)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

  for (method in names(method_tables)) {
    tab <- method_tables[[method]]
    write_csv(tab, file.path(dge_dir, paste0(roi_name, "_", method, ".csv")))

    summary_all <- bind_rows(summary_all, tibble(
      roi            = roi_name,
      method         = method,
      genes_tested   = sum(!is.na(tab$pvalue)),
      n_nominal_p    = sum(tab$pvalue < cfg$thresholds$nominal_p, na.rm = TRUE),
      n_fdr          = sum(tab$FDR < cfg$thresholds$fdr, na.rm = TRUE),
      n_fdr_original = sum(tab$FDR <= cfg$thresholds$fdr_original & abs(tab$logFC) > cfg$thresholds$abs_log2fc, na.rm = TRUE),
      min_fdr        = min(tab$FDR, na.rm = TRUE)
    ))

    pt_col <- case_when(
      tab$pvalue < cfg$thresholds$nominal_p & tab$logFC >  cfg$thresholds$abs_log2fc ~ "red",
      tab$pvalue < cfg$thresholds$nominal_p & tab$logFC < -cfg$thresholds$abs_log2fc ~ "blue",
      TRUE ~ "grey70")
    plot(tab$logFC, -log10(tab$pvalue), col = pt_col, pch = 20,
         xlab = "log2 fold change (Alcohol vs Control)", ylab = "-log10(nominal p)",
         main = paste0(method, "\n", sum(tab$FDR < cfg$thresholds$fdr, na.rm = TRUE), " genes FDR < ", cfg$thresholds$fdr),
         cex.main = 0.85)
    abline(v = c(-cfg$thresholds$abs_log2fc, cfg$thresholds$abs_log2fc), lty = 2, col = "grey40")
    abline(h = -log10(cfg$thresholds$nominal_p), lty = 2, col = "grey40")
    points(tab$logFC[which(tab$FDR < cfg$thresholds$fdr)], -log10(tab$pvalue[which(tab$FDR < cfg$thresholds$fdr)]), pch = 1, cex = 1.3)
  }
  dev.off()
}

write_csv(summary_all, file.path(cfg$paths$results_dir, "dge_summary.csv"))
print(as.data.frame(summary_all))
