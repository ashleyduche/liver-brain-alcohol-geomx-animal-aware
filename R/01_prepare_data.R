# Parse GEO-style metadata and attach animal ID (embedded in every sample name).
raw_meta <- read_excel(cfg$paths$metadata_xlsx, sheet = cfg$paths$metadata_sheet,
                       col_names = FALSE, .name_repair = "minimal")
header_row <- which(raw_meta[[1]] == "*library name")
stopifnot(length(header_row) == 1)
colnames(raw_meta) <- as.character(raw_meta[header_row, ])

meta <- raw_meta %>%
  select(`*library name`, `**tissue`, treatment, region, pathology) %>%   # sheet has blank/duplicate headers
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

counts <- read.csv(cfg$paths$raw_counts_csv, row.names = 1, check.names = FALSE)
neg_probe <- counts[cfg$paths$negative_probe, ]
counts <- counts[rownames(counts) != cfg$paths$negative_probe, ]
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

stopifnot(setequal(meta$sample_id, colnames(counts)))
counts <- counts[, meta$sample_id]

# Design check: the animal, not the ROI, is the experimental unit.
animal_design <- meta %>% distinct(animal, Group)
cat("Animals per group:\n"); print(table(animal_design$Group))
cat("\nROIs per animal x ROI type:\n")
print(table(meta$animal, paste(meta$region, meta$pathology, sep = "_")))

saveRDS(list(meta = meta, counts = counts, neg_probe = neg_probe),
        file.path(cfg$paths$results_dir, "processed", "prepared_data.rds"))
