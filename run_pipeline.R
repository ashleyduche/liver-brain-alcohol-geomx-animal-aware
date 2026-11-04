# Run from the project root with: Rscript run_pipeline.R
source("R/00_setup.R")
source("R/01_prepare_data.R")
source("R/02_qc_pca.R")
source("R/03_run_dge.R")
source("R/04_summary_figures.R")
source("R/05_label_swap_check.R")
source("R/06_filter_sensitivity.R")

# Save package versions so the run can be reproduced
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$results_dir, "sessionInfo.txt"))

cat("\nDone. See results/.\n")
