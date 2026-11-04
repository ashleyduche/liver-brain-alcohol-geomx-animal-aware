# How many genes survive each rowMeans cutoff, and how many pass FDR in the pseudobulk model.
cutoffs <- c(0, 10, 20, 30, 50, 100)
filter_sens <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]

  for (cutoff in cutoffs) {
    roi_data <- subset_roi(meta, counts, roi_def, min_rowmean = cutoff)
    coldata  <- roi_data$coldata
    animal_f <- factor(coldata$animal)
    pb_counts <- roi_data$counts %*% model.matrix(~ 0 + animal_f)
    colnames(pb_counts) <- levels(animal_f)
    pb_group <- coldata %>% distinct(animal, Group) %>% arrange(match(animal, levels(animal_f)))

    yb <- DGEList(pb_counts, group = pb_group$Group)
    yb <- yb[filterByExpr(yb, group = pb_group$Group), , keep.lib.sizes = FALSE]
    yb <- calcNormFactors(yb, method = "upperquartile")
    design_pb <- model.matrix(~ Group, data = pb_group)
    yb <- estimateDisp(yb, design_pb, robust = TRUE)
    fit_pb <- glmQLFit(yb, design_pb, robust = TRUE)
    tt <- topTags(glmQLFTest(fit_pb, coef = 2), n = Inf, sort.by = "none")$table

    filter_sens <- bind_rows(filter_sens, tibble(
      roi = roi_name, min_rowmean = cutoff,
      genes_kept  = nrow(roi_data$counts),
      genes_tested = nrow(tt),
      n_nominal_p = sum(tt$PValue < cfg$thresholds$nominal_p),
      n_fdr       = sum(tt$FDR < cfg$thresholds$fdr)))
  }
}

write_csv(filter_sens, file.path(cfg$paths$results_dir, "filter_sensitivity.csv"))
print(as.data.frame(filter_sens %>% select(roi, min_rowmean, genes_kept, n_fdr) %>%
                      pivot_wider(names_from = min_rowmean, values_from = c(genes_kept, n_fdr))))
