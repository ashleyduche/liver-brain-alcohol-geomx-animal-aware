# Read the GEO metadata sheet, starting from its "*library name" header row
raw_meta <- read_excel(cfg$paths$metadata_xlsx, sheet = cfg$paths$metadata_sheet,
                       col_names = FALSE, .name_repair = "minimal")
header_row <- which(raw_meta[[1]] == "*library name")
stopifnot(length(header_row) == 1)
colnames(raw_meta) <- as.character(raw_meta[header_row, ])

# Keep one row per ROI and pull the animal ID out of each sample name
meta <- raw_meta %>%
  select(`*library name`, `**tissue`, treatment, region, pathology) %>%   # the sheet has blank and duplicate headers
  slice((header_row + 1):n()) %>%
  filter(!is.na(`**tissue`)) %>%
  transmute(
    sample_id = `*library name`,
    tissue    = `**tissue`,
    region,
    pathology,
    animal    = str_extract(sample_id, "[0-9]{4}"),
    Group     = factor(if_else(tolower(treatment) == "ethanol", "Alcohol-Treated", "Control"),
                       levels = c("Control", "Alcohol-Treated"))
  )

# Read the counts and set the negative probe row aside
counts <- read.csv(cfg$paths$raw_counts_csv, row.names = 1, check.names = FALSE)
neg_probe <- counts[cfg$paths$negative_probe, ]
counts <- counts[rownames(counts) != cfg$paths$negative_probe, ]
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# Match the count columns to the metadata rows
stopifnot(setequal(meta$sample_id, colnames(counts)))
counts <- counts[, meta$sample_id]

# GeoMx ROI QC with the LOQ from the negative probe, Q3 signal against background and gene detection rate
meta <- meta %>%
  mutate(neg_probe       = as.numeric(neg_probe[sample_id]),
         q3              = apply(counts, 2, quantile, probs = 0.75),
         q3_to_neg       = q3 / neg_probe,
         loq             = pmax(cfg$geomx$min_loq, neg_probe),
         genes_above_loq = colMeans(sweep(counts, 2, loq, ">")))

# Remove ROIs that detect too few genes above their LOQ
low_quality <- meta$genes_above_loq < cfg$geomx$roi_min_gene_detection
cat(sprintf("ROI QC: %d of %d ROIs removed (< %.0f%% of genes above LOQ)\n",
            sum(low_quality), nrow(meta), 100 * cfg$geomx$roi_min_gene_detection))
meta   <- meta[!low_quality, ]
counts <- counts[, meta$sample_id]

# Design check, since the animal and not the ROI is the experimental unit
animal_design <- meta %>% distinct(animal, Group)
cat("Animals per group:\n")
print(table(animal_design$Group))

cat("\nROIs per animal x ROI type:\n")
print(table(meta$animal, paste(meta$region, meta$pathology, sep = "_")))

saveRDS(list(meta = meta, counts = counts, neg_probe = neg_probe),
        file.path(cfg$paths$results_dir, "processed", "prepared_data.rds"))
