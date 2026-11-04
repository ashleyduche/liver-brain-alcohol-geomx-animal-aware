# With 2 vs 2 animals there are only 3 distinct ways to split them into two pairs.
# Rerun the pseudobulk model on the 2 WRONG splits: any "significant" gene is a false positive,
# so this shows how much to trust FDR calls from the true split.
animal_ids <- c("6868", "6870", "6876", "6880")
true_split   <- c("g1", "g1", "g2", "g2")   # 6868 + 6870 (alcohol) vs 6876 + 6880 (control)
decoy_split1 <- c("g1", "g2", "g1", "g2")
decoy_split2 <- c("g1", "g2", "g2", "g1")
splits <- list(true = true_split, decoy1 = decoy_split1, decoy2 = decoy_split2)

swap_summary <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def  <- cfg$rois[[roi_name]]
  roi_data <- subset_roi(meta, counts, roi_def)
  animal_f <- factor(roi_data$coldata$animal)
  stopifnot(identical(levels(animal_f), animal_ids))

  pb_counts <- roi_data$counts %*% model.matrix(~ 0 + animal_f)
  colnames(pb_counts) <- levels(animal_f)

  for (split_name in names(splits)) {
    split_group <- factor(splits[[split_name]])
    yb <- DGEList(pb_counts, group = split_group)
    yb <- yb[filterByExpr(yb, group = split_group), , keep.lib.sizes = FALSE]
    yb <- calcNormFactors(yb, method = "upperquartile")
    design_pb <- model.matrix(~ split_group)
    yb <- estimateDisp(yb, design_pb, robust = TRUE)
    fit_pb <- glmQLFit(yb, design_pb, robust = TRUE)
    tt <- topTags(glmQLFTest(fit_pb, coef = 2), n = Inf, sort.by = "none")$table

    swap_summary <- bind_rows(swap_summary, tibble(
      roi = roi_name, split = split_name,
      n_nominal_p = sum(tt$PValue < cfg$thresholds$nominal_p),
      n_fdr = sum(tt$FDR < cfg$thresholds$fdr)))
  }
}

write_csv(swap_summary, file.path(cfg$paths$results_dir, "label_swap_check.csv"))
print(as.data.frame(swap_summary %>% pivot_wider(names_from = split, values_from = c(n_nominal_p, n_fdr))))
