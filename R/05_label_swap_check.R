# With 2 vs 2 animals there are only 3 ways to split the mice into two pairs
# The two incorrect splits put one alcohol mouse and one control in each group
# Any gene they call significant is a false positive, which shows how much to trust the true split
treated_animals <- sort(unique(meta$animal[meta$Group == "Alcohol-Treated"]))
control_animals <- sort(unique(meta$animal[meta$Group == "Control"]))
stopifnot(length(treated_animals) == 2, length(control_animals) == 2)
animal_ids <- c(treated_animals, control_animals)

# Group labels in the same order as animal_ids
true_split   <- c("g1", "g1", "g2", "g2")   # the 2 alcohol mice vs the 2 controls
decoy_split1 <- c("g1", "g2", "g1", "g2")
decoy_split2 <- c("g1", "g2", "g2", "g1")
splits <- list(true = true_split, decoy1 = decoy_split1, decoy2 = decoy_split2)

swap_summary <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def  <- cfg$rois[[roi_name]]
  roi_data <- subset_roi(meta, counts, roi_def)

  # Sum the ROI counts for each animal
  animal_f <- factor(roi_data$coldata$animal)
  pb_counts <- roi_data$counts %*% model.matrix(~ 0 + animal_f)
  colnames(pb_counts) <- levels(animal_f)
  pb_counts <- pb_counts[, animal_ids]

  # Rerun the pseudobulk model for each split
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
