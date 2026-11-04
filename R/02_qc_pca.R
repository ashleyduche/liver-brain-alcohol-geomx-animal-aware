# ROI level QC and PCA to see whether the animal explains more variance than the treatment

# Library size per ROI
qc <- meta %>%
  mutate(lib_size = colSums(counts))
write_csv(qc, file.path(cfg$paths$results_dir, "qc", "roi_qc_metrics.csv"))

# Q3 against negative probe background, the standard GeoMx QC plot
# Points close to the y = x line are close to background
tissue_col <- ifelse(qc$tissue == "brain", "purple", "darkorange")

png(file.path(cfg$paths$results_dir, "qc", "q3_vs_negprobe.png"), width = 1000, height = 900, res = 150)
plot(qc$neg_probe, qc$q3, log = "xy", col = tissue_col, pch = 16,
     xlab = "NegProbe-WTX count (background)", ylab = "Q3 count (signal)",
     main = "ROI signal vs background")
abline(0, 1, lty = 2, col = "grey40")
legend("topleft", legend = c("brain", "liver"), col = c("purple", "darkorange"), pch = 16, bty = "n")
dev.off()

# Library size for each animal
png(file.path(cfg$paths$results_dir, "qc", "library_size_by_animal.png"), width = 1600, height = 900, res = 200)
par(mar = c(9, 5, 2, 1))
boxplot(lib_size ~ interaction(animal, Group), data = qc, las = 2, xlab = "",
        ylab = "Library size (counts)", col = "grey85")
dev.off()

# PCA for each ROI type
pca_variance <- tibble()

for (roi_name in names(cfg$rois)) {
  roi_def  <- cfg$rois[[roi_name]]
  roi_data <- subset_roi(meta, counts, roi_def)
  coldata  <- roi_data$coldata

  y <- DGEList(roi_data$counts)
  y <- y[filterByExpr(y, group = coldata$Group), ]
  y <- calcNormFactors(y, method = "upperquartile")   # Q3, the NanoString GeoMx convention
  logcpm <- cpm(y, log = TRUE, prior.count = 2)

  pca <- prcomp(t(logcpm), scale. = TRUE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  # Color by group and shape by animal
  animals    <- factor(coldata$animal)
  group_col  <- ifelse(coldata$Group == "Control", "forestgreen", "firebrick")
  animal_pch <- c(15, 16, 17, 18)[as.integer(animals)]

  png(file.path(cfg$paths$results_dir, "qc", paste0(roi_name, "_PCA.png")), width = 1000, height = 900, res = 150)
  plot(pca$x[, 1], pca$x[, 2], col = group_col, pch = animal_pch, cex = 1.6,
       xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
       main = paste0(roi_def$label, " (ROI-level, Q3 logCPM)"))
  legend("topright", legend = c(levels(coldata$Group), paste("animal", levels(animals))),
         col = c("forestgreen", "firebrick", rep("black", 4)),
         pch = c(16, 16, 15, 16, 17, 18)[1:(2 + nlevels(animals))], bty = "n", cex = 0.8)
  dev.off()

  # Variance of PC1 and PC2 explained by animal compared with treatment (animal is nested in treatment)
  for (pc in c(1, 2)) {
    r2_animal <- summary(lm(pca$x[, pc] ~ animals))$r.squared
    r2_group  <- summary(lm(pca$x[, pc] ~ coldata$Group))$r.squared

    pca_variance <- bind_rows(pca_variance,
      tibble(roi = roi_name, PC = paste0("PC", pc), pct_var = pct[pc],
             r2_animal = round(r2_animal, 3), r2_treatment = round(r2_group, 3)))
  }
}

write_csv(pca_variance, file.path(cfg$paths$results_dir, "qc", "pca_variance_animal_vs_treatment.csv"))
print(as.data.frame(pca_variance))
