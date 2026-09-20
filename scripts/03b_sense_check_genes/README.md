## OAC RNA-seq Sense Check: Marker Gene Expression

Expression of literature-reported marker genes across normal stomach, normal oesophagus, Barrett's oesophagus and OAC, used to check that the data behave as expected biologically

## What the script does
1. Loads and cleans the sample metadata
2. Merges three count matrices and aligns them to the metadata
3. Filters to samples that passed QC
4. Defines a panel of 11 marker genes (TFF3, ANXA10, LGALS4, ERBB2, GATA4, GATA6, VEGFA, CDX2, MUC2, KRT7, KRT20)(well documented markers in OAC development) and maps them to the Ensembl IDs in the count matrix
5. Normalises counts with DESeq2 size factors and converts to log2
6. Plots each marker as a boxplot with individual samples across the four tissue groups

## Requirements
R packages: DESeq2, org.Hs.eg.db, AnnotationDbi, readxl, dplyr, tidyr, ggplot2, stringr, tibble

## Input files (not included)
- `Moremetadata.xlsx` with columns `SampleName`, `TissueType`, `Batch`, `QCStatus`
- `Count Matrix Batch 1.xlsx`, `Count Matrix Batch 2.xlsx`, `Count Matrix Batch 3.xlsx`, each with a `GeneID` column

## Outputs (not included)
Written to a `figures/` folder:
- `sense_check_genes.png`

## Notes
- Expression is exploratory, no batch correction is applied, and no statistical tests are run
- The script warns if a gene is missing from the count matrix or matches to more than one Ensembl ID

## Usage
Place the input files in the working directory and run `scripts/04_sense_check_genes/04_sense_check_genes.R`. This script is independent of the DGE analysis, so it can be run without running `02_dge` first.
