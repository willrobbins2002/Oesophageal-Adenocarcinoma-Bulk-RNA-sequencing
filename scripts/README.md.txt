# Full Triad RNA-seq Analysis (DESeq2)

Differential gene expression analysis across the full matched triad cohort of **Normal Oesophagus, Barrett's Oesophagus and Oesophageal Adenocarcinoma (OAC)** using DESeq2.

The analysis identifies genes differentially expressed between each pair of tissue types while accounting for the matched/repeated-measures nature of the samples.

## What the script does

1. Reads the sample metadata and cleans the tissue-type and batch information
2. Identifies patients who have all three required tissue types: Normal Oesophagus, Barrett's Oesophagus and OAC
3. Filters samples to those that have passed QC
4. Automatically identifies the count-matrix sample identifier (`MatrixID`), accounting for the different naming convention used by Batch 1
5. Reads and combines one or more count matrices
6. Removes Ensembl version suffixes and collapses duplicated Ensembl gene IDs by summing counts
7. Checks that all samples identified from the metadata are present in the supplied count matrices
8. Filters genes to retain those with at least 10 counts in at least 3 samples
9. Fits a DESeq2 model using the design `~ PatientID + TissueType`
10. Accounts for patient matching while comparing tissue types
11. Extracts three pairwise contrasts:

    * Barrett's Oesophagus vs Normal Oesophagus
    * OAC vs Barrett's Oesophagus
    * OAC vs Normal Oesophagus
12. Maps Ensembl gene IDs to HGNC gene symbols
13. Classifies genes as upregulated, downregulated or not significant based on an adjusted p-value < 0.05
14. Generates volcano plots for each comparison, highlighting the five most significant upregulated and downregulated genes
15. Combines the three volcano plots into a single figure
16. Saves the DESeq2 results, processed metadata, count matrix, plots and fitted DESeq2 object
17. Saves R and package version information to support reproducibility

## Requirements

R packages:

`readxl`, `dplyr`, `stringr`, `DESeq2`, `ggplot2`, `ggrepel`, `tibble`, `patchwork`, `org.Hs.eg.db`, `AnnotationDbi`

The required packages should be installed before running the script. See the main repository README for installation instructions.

## Input files

The analysis requires a metadata file and one or more count matrix files.

### Metadata

`Moremetadata.xlsx`

The metadata file should contain the following columns:

* `SampleName`
* `SampleID`
* `Batch`
* `TissueType`
* `QCStatus`

The script uses:

* `SampleID` to identify patients
* `TissueType` to identify Normal Oesophagus, Barrett's Oesophagus and OAC samples
* `QCStatus` to retain QC-pass samples
* `Batch` and `SampleName` to construct the count-matrix sample identifier (`MatrixID`)

### Count matrices

By default, the script expects:

* `Count Matrix Batch 1.xlsx`
* `Count Matrix Batch 2.xlsx`
* `Count Matrix Batch 3.xlsx`

Each count matrix should contain **gene IDs in the first column** and sample identifiers in the remaining columns.

The script can also work with a different number of count matrix files. If your data contain one or more count matrices with different filenames, edit the `count_files` section of the script:

```r
count_files <- c(
  "Count Matrix Batch 1.xlsx",
  "Count Matrix Batch 2.xlsx",
  "Count Matrix Batch 3.xlsx"
)
```

For example, if all samples are contained in a single count matrix:

```r
count_files <- "Count Matrix.xlsx"
```

If samples are distributed across two matrices:

```r
count_files <- c(
  "Count Matrix 1.xlsx",
  "Count Matrix 2.xlsx"
)
```

The supplied count matrix file(s) must collectively contain all samples required for the analysis. The script automatically checks that all samples identified from the metadata are present.

## Outputs

All outputs are written to:

`Full_Triad_Analysis/`

The script produces:

* `full_triad_metadata.csv` — metadata for the samples included in the final triad cohort
* `full_triad_counts.csv` — processed count matrix for the triad samples
* `DESeq2_Barrett_vs_Normal.csv` — gene-level results for Barrett's Oesophagus vs Normal Oesophagus
* `DESeq2_OAC_vs_Barrett.csv` — gene-level results for OAC vs Barrett's Oesophagus
* `DESeq2_OAC_vs_Normal.csv` — gene-level results for OAC vs Normal Oesophagus
* `DESeq2_summary.csv` — summary of significant, upregulated and downregulated genes for each comparison
* `Volcano_Barrett_vs_Normal.png` — volcano plot for Barrett's Oesophagus vs Normal Oesophagus
* `Volcano_OAC_vs_Barrett.png` — volcano plot for OAC vs Barrett's Oesophagus
* `Volcano_OAC_vs_Normal.png` — volcano plot for OAC vs Normal Oesophagus
* `Full_Triad_Volcano_Plots.png` — combined figure containing all three volcano plots
* `dds_full_triad.rds` — fitted DESeq2 object
* `session_info.txt` — R version and package information used for the analysis

## Notes

* Only patients with **all three required tissue types** are included.
* A patient may also have a Normal Stomach sample; this does not prevent inclusion in the triad cohort.
* Only samples marked as `Pass` in `QCStatus` are considered.
* The analysis uses the matched patient structure through the design `~ PatientID + TissueType`.
* **Batch is not included as a covariate and no batch correction is performed.**
* **Treatment modality is not included in the analysis.**
* The three tissue types are ordered as Normal Oesophagus, Barrett's Oesophagus and OAC, with Normal Oesophagus used as the reference level.
* DESeq2 significance is defined as an adjusted p-value (`padj`) < 0.05 using the Benjamini-Hochberg method.
* No additional log2 fold-change threshold is applied when defining significant genes.
* Positive log2 fold changes indicate higher expression in the first-named group of each comparison.
* Negative log2 fold changes indicate higher expression in the second-named group.
* The volcano plots label the five most significant upregulated and five most significant downregulated genes with available gene symbols.
* The analysis is standalone and does **not** require any previous analysis script or intermediate `.rds` file.

## Usage

Place `Moremetadata.xlsx` and the required count matrix file(s) in the working directory.

Then run:

```r
source("scripts/XX_full_triad/XX_full_triad_deseq2.R")
```

The exact script path should be updated to match the location and filename used in this repository.

The script will automatically identify the QC-pass patients with all three required tissue types, match their samples to the supplied count matrices, perform the DESeq2 analysis and generate the output files and figures.
