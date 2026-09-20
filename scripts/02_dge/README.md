## OAC RNA-seq Differential Gene Expression (DGE) Analysis

DGE of normal oesophagus versus OAC, with both tissues coming from the same patient in each test 

## What the script does
1. Loads and cleans the sample metadata
2. Merges three count matrices (Batches 1-3) and aligns them to the metadata
3. Subsets to QC-pass Normal Oesophagus and OAC samples with documented treatment modality, keeping only patients with both tissues
4. Builds a paired DESeq2 model (`~ PatientID + TissueType`, reference level Normal Oesophagus) and removes low-count genes (at least 10 counts in at least 3 samples)
5. Runs DESeq2 and extracts OAC vs Normal Oesophagus results (adjusted p-value threshold 0.05)
6. Applies ashr log2 fold-change shrinkage for ranking and plotting
7. Saves the fitted model and results for downstream scripts

## Requirements
R packages: DESeq2, ashr, readxl, dplyr, stringr, tibble

## Input files (not included)
- `Moremetadata.xlsx` with columns `SampleName`, `SampleID`, `TissueType`, `Batch`, `QCStatus`, `TreatmentModality`
- `Count Matrix Batch 1.xlsx`, `Count Matrix Batch 2.xlsx`, `Count Matrix Batch 3.xlsx`, each with a `GeneID` column

## Outputs (not included)
Written to a `results/` folder:
- `dds_de.rds`: fitted DESeq2 object
- `res_shrunk_df.rds`: shrunken results table
- `dge_oac_vs_normal_oesophagus_shrunk.csv` and `dge_oac_vs_normal_oesophagus_unshrunk.csv`: full results tables
- `session_info_dge.txt`: R and package versions

## Usage
Place the input files in the working directory and run `scripts/02_dge_paired/02_dge_paired_oac_vs_normal.R'. Run this script before the downstream scripts, which read the saved results.
