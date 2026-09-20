## OAC RNA-seq PCA plot

PCA plot of QC-pass RNA-seq samples from normal stomach, normal oesophagus, Barrett's oesophagus, and oesophageal adenocarcinoma (OAC)

## What the script does
1. Loads and cleans the sample metadata
2. Merges three count matrices (Batches 1-3)
3. Filters to samples that passed QC
4. Runs a variance-stabilising transformation (DESeq2)
5. Plots a PCA of the top 1000 most variable genes, coloured by tissue type

## Requirements
R packages: DESeq2, readxl, dplyr, ggplot2, stringr

## Input files (not included)
- `Moremetadata.xlsx` with columns `SampleName`, `TissueType`, `Batch`, `QCStatus`
- `Count Matrix Batch 1.xlsx`, `Count Matrix Batch 2.xlsx`, `Count Matrix Batch 3.xlsx`, each with a `GeneID` column

## Usage
Place the input files in the working directory and run `scripts/01_pca.R`
