# Animal-aware GeoMx DGE: alcohol exposure in APP/PS1 liver and brain

Data: GEO [GSE324193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324193). Re-analysis of a GeoMx whole-transcriptome spatial dataset (APP/PS1 mice, intragastric alcohol vs. control, 5 weeks; brain plaque-bearing/plaque-free cortex and hippocampus, liver periportal/perivenous) that treats the **animal**, not the ROI, as the experimental unit. ROIs are aggregated per animal (**pseudobulk**), which avoids pseudoreplication.

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

## Processing

Upstream (as described in GSE324193): GeoMx NGS Pipeline → DCC files → per-ROI digital counts (`RawCounts_AllSamples.csv`, including the aggregated `NegProbe-WTX` background probe). Downstream, this repo follows the standard GeoMx workflow (NanoString GeomxTools):

1. **ROI QC** — limit of quantification (LOQ) per ROI = max(2, `NegProbe-WTX` count); ROIs with < 10% of genes above LOQ are removed (none were); Q3 signal vs. background is plotted (`results/qc/q3_vs_negprobe.png`).
2. **Gene filtering** — a gene is kept for an ROI type if it is above LOQ in ≥ 10% of that ROI type's ROIs, then `filterByExpr`.
3. **Normalization** — upper-quartile (Q3).
4. **Differential expression** — the three models above; Alcohol-Treated vs. Control; log2FC reported with every p-value; Benjamini-Hochberg FDR. The ROI-level model with animal as a correlated block corresponds to NanoString's recommended linear mixed model with a random intercept per animal/slide.

The `NegProbe-WTX` row in the published counts is already aggregated, so the LOQ uses it directly; the per-probe geometric SD that GeomxTools also uses is not available. All thresholds are in `config.yml` under `geomx:`.

QC note: brain ROIs sit closer to background than liver ROIs (Q3 is about 1.2–1.5× the negative probe in brain vs. up to ~6× in liver), so brain gene-level estimates are noisier.

Additional checks: PCA with variance explained by animal vs. treatment (`results/qc/`), and a **label-swap check** (`R/05_label_swap_check.R`) that reruns the pseudobulk model on the two *wrong* animal pairings. Genes called significant there are false positives, which calibrates how far to trust FDR from the true split.

## Results (primary pseudobulk model)

| ROI type | genes tested | FDR < 0.05 | FDR < 0.05 under the two wrong animal pairings |
|---|---|---|---|
| plaqueCortex | 15570 | 1 | 3 / 1 |
| nonplaqueCortex | 18549 | 1 | 19 / 8 |
| plaqueHippo | 16624 | 1993 | 0 / 0 |
| nonplaqueHippo | 17024 | 0 | 0 / 0 |
| periportal | 19303 | 9 | 823 / 8 |
| perivenous | 18235 | 78 | 12 / 23 |

Only plaque-bearing hippocampus gives more hits for the true grouping than for wrong groupings by a wide margin. Even there, with 2 animals per group, these are candidates for validation, not confirmed DEGs.

`R/06_filter_sensitivity.R` adds an optional mean-count filter (`thresholds: min_rowmean`, off by default) at 0–100 on top of the LOQ filter (`results/filter_sensitivity.csv`); liver results are sensitive to it (e.g. perivenous 78 → 31 FDR hits at `rowMeans >= 50`).

## Run

```bash
Rscript run_pipeline.R                                      # synthetic example data (default)
GEOMX_CONFIG=config_real.yml Rscript run_pipeline.R         # published data (GSE324193), in data/raw/
```

`data/raw/` holds the published study counts and metadata. `results/` tracks figures, QC, summaries and the primary pseudobulk tables; the larger intermediate tables regenerate when you run the pipeline. Needs R with `tidyverse`, `readxl`, `yaml`, `edgeR`, `limma`, `DESeq2`, `pheatmap`.

```
R/00_setup.R               packages, config, ROI subsetting helper
R/01_prepare_data.R        metadata + counts, animal ID, LOQ and ROI QC, design table
R/02_qc_pca.R              ROI QC (Q3 vs background), per-ROI-type PCA, animal vs treatment variance
R/03_run_dge.R             the three models, tables, volcano panels, ranked lists
R/04_summary_figures.R     hits-per-method barplot, animal-level heatmaps of top genes
R/05_label_swap_check.R    wrong-grouping false-positive check
R/06_filter_sensitivity.R  optional mean-count filter sensitivity
```

Outputs in `results/`: `dge_tables/` (per ROI per model), `dge_summary.csv`, `volcano/`, `heatmaps/`, `qc/`, `ranked_lists/` (pre-ranked `.rnk` for GSEA), `label_swap_check.csv`, `filter_sensitivity.csv`, `sessionInfo.txt`.

## Interpreting the results

- Report effect sizes and the top-ranked genes at the animal level; call them candidates, not findings.
- Check `label_swap_check.csv`: if wrong pairings produce as many "significant" genes as the true one, the FDR from the true split is not informative.
- Confirmation needs more animals (or an orthogonal assay such as RNAscope/qPCR on candidate genes).

## References

- Lazic SE, Clarke-Williams CJ, Munafò MR (2018). What exactly is "N" in cell culture and animal experiments? *PLoS Biol* 16:e2005282.
- Zimmerman KD, Espeland MA, Langefeld CD (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nat Commun* 12:738.
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- NanoString. GeomxTools: NanoString GeoMx tools (Bioconductor package and workflow vignette).
- Smyth GK, Michaud J, Scott HS (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21:2067-75.
