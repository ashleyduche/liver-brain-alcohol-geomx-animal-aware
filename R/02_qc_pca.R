# ROI level QC and PCA to see whether the animal explains more variance than the treatment

# Library size per ROI
qc <- meta %>%
  mutate(lib_size = colSums(counts))
write_csv(qc, file.path(cfg$paths$results_dir, "qc", "roi_qc_metrics.csv"))

# Q3 against negative probe background, the standard GeoMx QC plot
# Points close to the y = x line are close to background
tissue_col <- ifelse(qc$tissue == "brain", "purple", "darkorange")

# Both axes share limits that start at 1 so the y = x line is always inside the plot
# A NegProbe count of 0 is drawn at 1 since a log axis cannot show 0
neg_plot  <- pmax(qc$neg_probe, 1)
axis_lims <- c(1, max(c(neg_plot, qc$q3)) * 1.5)

png(file.path(cfg$paths$results_dir, "qc", "q3_vs_negprobe.png"), width = 1000, height = 900, res = 150)
plot(neg_plot, qc$q3, log = "xy", col = tissue_col, pch = 16,
     xlim = axis_lims, ylim = axis_lims,
     xlab = "NegProbe-WTX count (background)", ylab = "Q3 count (signal)",
     main = "ROI signal vs background")
abline(0, 1, lty = 2, col = "grey40")
legend("bottomright", legend = c("brain", "liver", "y = x (signal equals background)"),
       col = c("purple", "darkorange", "grey40"), pch = c(16, 16, NA), lty = c(NA, NA, 2), bty = "n")
dev.off()

# Library size for each animal, ordered and colored by group
animal_order <- qc %>% distinct(animal, Group) %>% arrange(Group, animal)
qc_box <- qc %>% mutate(animal = factor(animal, levels = animal_order$animal))
box_col <- ifelse(animal_order$Group == "Control", "forestgreen", "firebrick")

png(file.path(cfg$paths$results_dir, "qc", "library_size_by_animal.png"), width = 1600, height = 900, res = 200)
par(mar = c(4.5, 6, 2, 1))
boxplot(lib_size / 1e6 ~ animal, data = qc_box, col = adjustcolor(box_col, alpha.f = 0.5),
        xlab = "Animal", ylab = "Library size (million counts)",
        main = "Library size per ROI, by animal")
legend("topright", legend = levels(qc$Group), fill = adjustcolor(c("forestgreen", "firebrick"), alpha.f = 0.5),
       bty = "n", cex = 0.8)
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

  # Each animal's shape in the legend is drawn in its group color
  legend_group <- coldata %>% distinct(animal, Group) %>% arrange(animal)
  legend_col   <- ifelse(legend_group$Group == "Control", "forestgreen", "firebrick")

  # Legend sits in the right margin so it never covers points
  png(file.path(cfg$paths$results_dir, "qc", paste0(roi_name, "_PCA.png")), width = 1200, height = 900, res = 150)
  par(mar = c(5, 5, 4, 11))
  plot(pca$x[, 1], pca$x[, 2], col = group_col, pch = animal_pch, cex = 1.6,
       xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
       main = paste0(roi_def$label, " (ROI-level, Q3 logCPM)"))
  legend("topleft", inset = c(1.02, 0), xpd = TRUE, title = "Animal",
         legend = paste0(legend_group$animal, " (", legend_group$Group, ")"),
         col = legend_col, pch = c(15, 16, 17, 18)[seq_len(nrow(legend_group))],
         pt.cex = 1.4, bty = "n", cex = 0.8)
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
