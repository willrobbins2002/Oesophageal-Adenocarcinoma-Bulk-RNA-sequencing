# Oesophageal Adenocarcinoma Bulk RNA-sequencing

An analysis pipeline for bulk RNA-seq data from oesophageal adenocarcinoma (OAC) patients. Analyses are performed on pre-produced count matrices, one per sample batch, combined with a shared sample metadata file.

## Contents

```
scripts/
  01_pca/                         — Principle Component Analysis, coloured by tissue type and batch, for all QC-pass samples
  02_dge/                         — Differential Gene Expression analysis for patients who had normal oesophagus AND OAC
  03a_top_genes/                  — Plotted the expression of the top 10 most downregulated and upregulated genes from DGE 
  03b_sense_check_genes/          — Plotted the expression of 10-well documented genes associated with BO and OAC
  04_gsea/                        — Gene Set Enrichment analysis, associating gene expression from DGE to hallmarks and pathways 
  05_treatment_modality_analysis/ — Naïve vs Chemo vs CRT DESeq2 and Hallmark GSEA across all QC-pass OAC samples
  06a_full_triad_volcano_plots/   — Paired Normal vs Barrett's vs OAC DESeq2 and volcano plots for the matched triad cohort
  06b_full_triad_heatmap/         — Heatmap of the top DEGs for the matched triad cohort
```

Each script folder has its own README describing that analysis's inputs, outputs and any analysis-specific settings.

## How the analyses relate

`01_pca` runs on all samples in the cohort. `02_dge`, `03a_top_genes` and `03b_sense_check_genes` work from a shared subset of 52 patients, all of whom have QC-pass normal oesophagus and OAC samples. `04_gsea` is run based on the results from DGE. `05_treatment_modality_analysis` uses 94 QC-pass OAC samples with a recorded treatment modality (Naïve, Chemo, CRT), independent of the DGE/GSEA sample set. `06a_full_triad_volcano_plots` and `06b_full_triad_heatmap` both use a separate, smaller cohort: 7 patients with a complete matched triad of Normal, Barrett's and OAC samples, analysed with a paired design.

## Count matrix files

The analysis scripts are written to work with one or more count matrix files.

By default, the scripts expect the following three files:

```r
count_files <- c(
  "Count Matrix Batch 1.xlsx",
  "Count Matrix Batch 2.xlsx",
  "Count Matrix Batch 3.xlsx"
)
```

If your data are provided as a different number of count matrix files, edit this section of the script so that it contains the filenames of your available count matrices.

For example, if you have one count matrix file:

```r
count_files <- "Count Matrix.xlsx"
```

If you have two count matrix files:

```r
count_files <- c(
  "Count Matrix 1.xlsx",
  "Count Matrix 2.xlsx"
)
```

If you have three count matrix files:

```r
count_files <- c(
  "Count Matrix Batch 1.xlsx",
  "Count Matrix Batch 2.xlsx",
  "Count Matrix Batch 3.xlsx"
)
```

The supplied count matrix file(s) must collectively contain the samples required for the analysis. Each script automatically combines the supplied matrices and checks that all samples identified from the metadata are present.

## Metadata file

All scripts read a shared sample metadata file, `Moremetadata.xlsx`, which must include:

| Column | Description |
|---|---|
| `SampleName` | Longer sample identifier which includes patient and tissue type, not used in analyses  |
| `SampleID` | Abbreviated sample identified, used to match rows against count matrices |
| `Batch` | Sequencing batch (1, 2 or 3), used to select the correct sample-matching rule |
| `TissueType` | Tissue type (e.g. Normal, Barretts, OAC) |
| `QCStatus` | QC outcome; scripts filter to `Pass` |
| `Treatment Modality` | Treatment group (Naïve, Chemo, CRT), used by the treatment modality analysis |

Batch 1 samples are matched to count-matrix columns using a short code derived from `SampleName`; Batches 2 and 3 are matched directly on `SampleName`. See individual script READMEs for the exact matching logic.

## Requirements

R version: 2026.09.0+174

Each script's own README lists the R packages it needs

## How to run

1. Place the required count matrix file(s) and `Moremetadata.xlsx` in your working directory
2. Install the R packages listed in the relevant script's README
3. Source or run the script for the analysis you want (e.g. `scripts/06_full_triad/06_full_triad_deseq2.R`)
4. Outputs are written to `results/` (tables, `.rds` objects, session info) and `figures/` (plots) inside the working directory

## Data availability

The data used in this analysis has been kept private due to patient confidentiality, institutional data-governance restrictions, and ethical approval

## Contact

Will Robbins
willrobbins2002@gmail.com
