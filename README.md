# Animal-aware GeoMx DGE: alcohol exposure in APP/PS1 liver and brain

Data: GEO [GSE324193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324193). Re-analysis of a GeoMx whole-transcriptome spatial dataset (APP/PS1 mice, intragastric alcohol vs. control, 5 weeks; brain plaque-bearing/plaque-free cortex and hippocampus, liver periportal/perivenous) that treats the **animal**, not the ROI, as the experimental unit.

> [!WARNING]
> **Exploratory analysis.** The study has **4 animals (2 control, 2 alcohol-treated)**. With n = 2 vs. 2 there are only 6 possible group labelings, so no permutation-based animal-level test can reach p < 0.05 (the smallest attainable two-sided p is 2/6 = 0.33). Parametric p-values at this n lean heavily on modeled variance and should be read as a ranking of candidate genes, not as evidence. See [Interpreting the results](#interpreting-the-results).

Companion repo: [`liver-brain-alcohol-dge-pipeline`](https://github.com/ashleyduche/liver-brain-alcohol-dge-pipeline) (original analysis, every ROI treated as an independent sample).

## Why a second analysis

ROIs from the same animal are correlated, so 16 periportal ROIs from 4 mice are not 16 independent replicates (pseudoreplication; Lazic et al. 2018, *PLoS Biol*; Zimmerman et al. 2021, *Nat Commun*; Squair et al. 2021, *Nat Commun*). This pipeline runs three models side by side, per ROI type:

| Model | Replicate unit | Role |
|---|---|---|
| `naive_DESeq2_ROIs` | each ROI | reproduces the original analysis, for comparison only |
| `voom_dupCor_ROIs` | ROI, with animal as a correlated block (limma `duplicateCorrelation`) | sensitivity analysis |
| `pseudobulk_edgeR_animals` | animal (ROI counts summed per animal; edgeR quasi-likelihood F-test) | **primary** |

Shared choices: upper-quartile (Q3) normalization, the GeoMx standard; a mean-count filter (`thresholds: min_rowmean` in the config, default 50: keep genes whose mean count across an ROI type's ROIs is >= 50) followed by `filterByExpr`; Alcohol-Treated vs. Control contrast; effect size (log2FC) reported with every p-value; Benjamini-Hochberg FDR.

Additional checks: PCA with variance explained by animal vs. treatment (`results/qc/`), and a **label-swap check** (`R/05_label_swap_check.R`) that reruns the pseudobulk model on the two *wrong* animal pairings. Genes called significant there are false positives, which calibrates how far to trust FDR from the true split.

## Gene filtering

`R/06_filter_sensitivity.R` reruns the primary model at mean-count cutoffs of 0, 10, 20, 30, 50 and 100 (`results/filter_sensitivity.csv`). The cutoff matters most in the liver:

| ROI type | genes kept at 30 / 50 / 100 | pseudobulk FDR < 0.05 at 30 / 50 / 100 |
|---|---|---|
| periportal | 19724 / 13982 / 7117 | 10 / 0 / 0 |
| perivenous | 19929 / 17423 / 8580 | 166 / 34 / 0 |
| plaqueHippo | 19961 / 19951 / 17814 | 2295 / 2298 / 2142 |

For the original ROI-as-sample DESeq2 model, `rowMeans >= 30` gives 189 periportal DEGs and `rowMeans >= 50` gives 97 plaque-bearing hippocampus DEGs, matching the original README. Perivenous (802 DEGs, 17,325 genes) could not be matched exactly: 802 genes was not reached at any cutoff tried (804-805 at 10-30), and `rowMeans >= 50` keeps 17,423 genes.

## Run

```bash
Rscript run_pipeline.R                                      # synthetic example data (default)
GEOMX_CONFIG=config_real.yml Rscript run_pipeline.R         # real data, once placed in data/raw/
```

`data/raw/` holds the published study counts and metadata. `results/` tracks figures, QC, summaries and the primary pseudobulk tables; the larger intermediate tables regenerate when you run the pipeline. Needs R with `tidyverse`, `readxl`, `yaml`, `edgeR`, `limma`, `DESeq2`, `pheatmap`.

```
R/00_setup.R               packages, config, ROI subsetting helper
R/01_prepare_data.R        metadata + counts, animal ID from sample name, design table
R/02_qc_pca.R              library size, per-ROI-type PCA, animal vs treatment variance
R/03_run_dge.R             the three models, tables, volcano panels, ranked lists
R/04_summary_figures.R     hits-per-method barplot, animal-level heatmaps of top genes
R/05_label_swap_check.R    wrong-grouping false-positive check
```

Outputs in `results/`: `dge_tables/` (per ROI per model), `dge_summary.csv`, `volcano/`, `heatmaps/`, `qc/`, `ranked_lists/` (pre-ranked `.rnk` for GSEA), `label_swap_check.csv`, `sessionInfo.txt`.

## Interpreting the results

- Report effect sizes and the top-ranked genes at the animal level; call them candidates, not findings.
- Check `label_swap_check.csv`: if wrong pairings produce as many "significant" genes as the true one, the FDR from the true split is not informative.
- Confirmation needs more animals (or an orthogonal assay such as RNAscope/qPCR on candidate genes).

## References

- Lazic SE, Clarke-Williams CJ, Munafò MR (2018). What exactly is "N" in cell culture and animal experiments? *PLoS Biol* 16:e2005282.
- Zimmerman KD, Espeland MA, Langefeld CD (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738.
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- Smyth GK, Michaud J, Scott HS (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21:2067-75.
