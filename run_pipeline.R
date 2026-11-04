# Run from the project root: Rscript run_pipeline.R
source("R/00_setup.R")
source("R/01_prepare_data.R")
source("R/02_qc_pca.R")
source("R/03_run_dge.R")
source("R/04_summary_figures.R")
source("R/05_label_swap_check.R")

writeLines(capture.output(sessionInfo()), file.path(cfg$paths$results_dir, "sessionInfo.txt"))
cat("\nDone. See results/.\n")
