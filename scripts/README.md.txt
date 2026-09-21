## OAC RNA-seq Treatment Modality Analysis (DESeq2 and Hallmark GSEA)

Pairwise comparison of treatment modalities (Naïve, Chemo, CRT) across all QC-pass OAC samples, using DESeq2 for gene-level differential expression and pre-ranked GSEA on the Hallmark gene sets

## What the script does
1. Reads and combines the three count matrices
2. Reads the sample metadata and matches sample names to count-matrix columns (Batch 1 uses short codes; Batches 2 and 3 use `SampleName`)
3. Filters to QC-pass OAC samples with treatment modality Naïve, Chemo or CRT, and resolves duplicated Batch 1 entries 
4. Removes Ensembl version suffixes and filters low-count genes (at least 10 counts in as many samples as the smallest treatment group)
5. Fits a DESeq2 model with design `~ TreatmentModality` 
6. Extracts three pairwise contrasts: Naïve vs Chemo, Naïve vs CRT and Chemo vs CRT
7. Maps Ensembl IDs to gene symbols
8. Runs pre-ranked GSEA on the MSigDB Hallmark sets, ranking all mapped genes by signed Wald statistic without filtering on gene-level significance
9. Plots the top 10 Hallmark pathways for each comparison (five enriched in each group) and saves all tables and plots

## Requirements
R packages: readxl, dplyr, tibble, stringr, readr, ggplot2, ggrepel, DESeq2, fgsea, msigdbr, AnnotationDbi, org.Hs.eg.db

## Input files (not included)
- `Moremetadata.xlsx` with columns `SampleName`, `SampleID`, `Batch`, `TissueType`, `QCStatus`, `TreatmentModality`
- `Count Matrix Batch 1.xlsx`, `Count Matrix Batch 2.xlsx`, `Count Matrix Batch 3.xlsx`, each with gene IDs in the first column

## Outputs (not included)
Written to `results/` and `figures/`:
- `figures/*_Hallmark_GSEA_top10.png`: top 10 Hallmark pathways for each of the three comparisons
- `results/*_Hallmark_GSEA.csv`: full Hallmark GSEA results for each comparison
- `results/*_DESeq2_gene_results.csv`: gene-level DESeq2 results for each comparison
- `results/*.rds`: fitted DESeq2 object, DESeq2 and GSEA result lists, and the metadata used
- `results/session_info_treatment_modality.txt`: R and package versions

## Notes
- The first-named group is the numerator in each comparison, so positive log2 fold changes and positive normalised enrichment scores (NES) mean higher or enriched in that group (Naïve for Naïve vs Chemo and Naïve vs CRT, Chemo for Chemo vs CRT)
- DESeq2 significance uses FDR < 0.05. Pathways are also counted at the GSEA FDR < 0.25 threshold (more lenient)
- This analysis is independent of the paired normal oesophagus vs OAC analysis (`02_dge`) and does not need any earlier script to be run first

## Usage
Place the input files in the working directory and run `scripts/05_treatment_modality/05_treatment_modality_deseq2_gsea.R`.