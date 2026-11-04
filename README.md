# APP/PS1 Alcohol Liver and Brain Spatial Transcriptomics: Animal Level Analysis

This repository contains a reanalysis of the GeoMx spatial transcriptomics data published on GEO under [GSE324193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324193). The study profiled the liver and brain of APP/PS1 mice that were fed alcohol intragastrically for 5 weeks. In the brain, regions of interest (ROIs) were collected from the cortex and hippocampus, both with and without amyloid plaques. In the liver, ROIs were collected from the periportal and perivenous zones.

The original analysis in [liver-brain-alcohol-dge-pipeline](https://github.com/ashleyduche/liver-brain-alcohol-dge-pipeline) treated every ROI as an independent sample. However, ROIs collected from the same mouse are more similar to each other than to ROIs from other mice, so they are not independent replicates. Treating them as independent is called pseudoreplication, and it makes p values look stronger than the data support. To avoid this, this pipeline treats each animal as the replicate. It does so by adding together the counts of all ROIs from the same animal and ROI type, which is an approach called pseudobulk.

> [!WARNING]
> **This is an exploratory analysis.** The study includes four mice, with two controls and two mice fed alcohol. With only two animals per group, there are just six ways to assign the four mice to two groups. As a result, a permutation test can never reach p < 0.05, because the smallest possible p value is 0.33. The p values reported here depend heavily on the statistical model, so they are best used to rank candidate genes rather than as proof of differential expression.

## Analysis Approach

Each ROI type is analyzed with three models so that their results can be compared directly. The first model repeats the original analysis and is included only for comparison. The second model keeps each ROI as a sample but tells limma that ROIs from the same animal are correlated. The third model is the primary analysis, and it uses pseudobulk counts so that each animal contributes one sample.

| Model | What counts as a replicate | Role |
|---|---|---|
| `naive_DESeq2_ROIs` | Each ROI | Repeats the original analysis for comparison |
| `voom_dupCor_ROIs` | Each ROI, with ROIs from the same animal modeled as correlated (limma `duplicateCorrelation`) | Sensitivity analysis |
| `pseudobulk_edgeR_animals` | Each animal, with ROI counts summed per animal and tested with the edgeR quasi likelihood F test | **Primary analysis** |

## Processing

According to the GEO record, the sequencing reads were processed with the GeoMx NGS Pipeline to produce DCC files. These files were then converted into a table of raw counts for every ROI, which is the `RawCounts_AllSamples.csv` file in `data/raw/`. That table also includes `NegProbe-WTX`, which measures the background signal from probes that do not target any gene.

Starting from those counts, this pipeline follows the standard GeoMx workflow from NanoString's GeomxTools package.

1. **ROI quality control.** Each ROI is given a limit of quantification (LOQ), which is the larger of 2 and its `NegProbe-WTX` count. Any ROI that detects fewer than 10% of genes above its LOQ is removed. Every ROI in this study passed this check.
2. **Gene filtering.** A gene is kept for an ROI type when it is above the LOQ in at least 10% of that type's ROIs. The remaining genes then pass through the edgeR `filterByExpr` function.
3. **Normalization.** Counts are normalized with the upper quartile (Q3) method, which is the GeoMx standard.
4. **Differential expression.** Each model compares the mice fed alcohol with the controls. Every gene is reported with its log2 fold change and p value, and p values are adjusted with the Benjamini Hochberg false discovery rate (FDR). The second model, which treats ROIs from the same animal as correlated, matches the linear mixed model that NanoString recommends, where each animal is given its own random intercept.

The published `NegProbe-WTX` count is already combined across the individual negative probes. Because of this, the LOQ uses that count directly and cannot include the spread between probes that GeomxTools normally adds. All of these thresholds can be changed in `config.yml` under `geomx:`.

The quality control plot in `results/qc/q3_vs_negprobe.png` shows that brain ROIs sit much closer to background than liver ROIs. In the brain, the Q3 signal is only about 1.2 to 1.5 times the negative probe, while in the liver it reaches about 6 times. Gene level results from the brain are therefore noisier.

## Checking the Results

Two additional checks help judge how much to trust the primary results. First, a PCA for each ROI type measures how much of the variation is explained by the animal compared with the treatment. Second, a label swap check reruns the pseudobulk model using the two incorrect ways of pairing the mice. Because those groupings have no biological meaning, any gene they call significant is a false positive. Comparing these numbers with the true grouping shows whether the true result stands out from chance.

## Results

The table below shows the results from the primary pseudobulk model.

| ROI type | Genes tested | Genes with FDR < 0.05 | Genes with FDR < 0.05 using the two incorrect pairings |
|---|---|---|---|
| plaqueCortex | 15570 | 1 | 3 / 1 |
| nonplaqueCortex | 18549 | 1 | 19 / 8 |
| plaqueHippo | 16624 | 1993 | 0 / 0 |
| nonplaqueHippo | 17024 | 0 | 0 / 0 |
| periportal | 19303 | 9 | 823 / 8 |
| perivenous | 18235 | 78 | 12 / 23 |

The plaque bearing hippocampus is the only ROI type where the true grouping produces far more significant genes than the incorrect pairings. In the periportal liver, one incorrect pairing produced 823 significant genes compared with 9 for the true grouping, so those results cannot be trusted. Even in the hippocampus, two animals per group is not enough to confirm differential expression, so these genes should be treated as candidates for validation.

The script `R/06_filter_sensitivity.R` also tests an optional filter that keeps genes based on their mean count. This filter is set by `thresholds: min_rowmean` and is turned off by default. The liver results change with this filter. For example, the number of perivenous genes with FDR < 0.05 drops from 78 to 31 when genes need a mean count of at least 50.

## Running the Pipeline

The pipeline runs on a synthetic example dataset by default. To run it on the published data instead, use `config_real.yml`.

```bash
Rscript run_pipeline.R                                      # synthetic example data
GEOMX_CONFIG=config_real.yml Rscript run_pipeline.R         # published data from GSE324193
```

The `data/raw/` folder holds the published counts and metadata. The `results/` folder holds the figures, quality control files, summaries, and the primary pseudobulk tables. The larger intermediate tables are not stored in the repository, but they are created again each time the pipeline runs. The pipeline needs R with the `tidyverse`, `readxl`, `yaml`, `edgeR`, `limma`, `DESeq2`, and `pheatmap` packages.

```
R/00_setup.R               loads packages, settings, and the ROI subsetting helper
R/01_prepare_data.R        reads metadata and counts, then adds animal ID, LOQ, and ROI quality control
R/02_qc_pca.R              plots Q3 against background and runs PCA for each ROI type
R/03_run_dge.R             runs the three models and saves tables, volcano plots, and ranked gene lists
R/04_summary_figures.R     plots significant genes per model and heatmaps of the top genes
R/05_label_swap_check.R    reruns the primary model with the incorrect animal pairings
R/06_filter_sensitivity.R  tests the optional mean count filter
```

Each run saves its outputs to `results/`. The `dge_tables/` folder holds one table per ROI type and model, and `dge_summary.csv` summarizes them. The `volcano/`, `heatmaps/`, and `qc/` folders hold the figures. The `ranked_lists/` folder holds ranked gene lists that can be used for GSEA. The label swap and filter checks are saved as `label_swap_check.csv` and `filter_sensitivity.csv`, and the R session details are saved in `sessionInfo.txt`.

## Interpreting the Results

Because the study has only two animals per group, the genes reported here should be described as candidates rather than findings. The most useful outputs are the effect sizes and the top ranked genes at the animal level. Before relying on any result, check `label_swap_check.csv`. If the incorrect pairings produce as many significant genes as the true grouping, then the FDR from the true grouping is not informative. Confirming these candidates will require more animals or an independent method such as RNAscope or qPCR.

## References

- Lazic SE, Clarke-Williams CJ, Munafò MR (2018). What exactly is "N" in cell culture and animal experiments? *PLoS Biol* 16:e2005282.
- Zimmerman KD, Espeland MA, Langefeld CD (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738.
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- NanoString. GeomxTools: NanoString GeoMx tools (Bioconductor package and workflow vignette).
- Smyth GK, Michaud J, Scott HS (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21:2067-75.
