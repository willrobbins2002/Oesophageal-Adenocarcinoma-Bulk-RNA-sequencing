## OAC RNA-seq Gene Set Analysis (CAMERA)

Gene set analysis of OAC versus normal oesophagus, using the same paired design as the DGE analysis (each patient's normal oesophagus and OAC samples are compared)

## What the script does
1. Loads the fitted paired DESeq2 object (`02_dge`) and extracts raw counts and sample metadata
2. Removes Ensembl version suffixes so gene IDs match the gene sets
3. Builds the paired design (`~ PatientID + TissueType`, reference level normal oesophagus), filters low-expression genes (edgeR `filterByExpr`) and applies TMM normalisation
4. Applies a voom transformation
5. Retrieves Hallmark, Reactome and Gene Ontology (BP, CC, MF) gene sets from MSigDB via `msigdbr`
6. Runs CAMERA (Wu & Smyth, 2012) on each collection for the OAC vs normal oesophagus contrast
7. Reports gene sets with FDR < 0.05, saves all results and plots the top Hallmark sets

## Requirements
R packages: DESeq2, edgeR, limma, msigdbr, dplyr, tibble, ggplot2

## Input files (not included)
- `results/dds_de.rds`, created by `scripts/02_dge/02_dge_paired_oac_vs_normal.R`

## Outputs (not included)
- `results/GSEA_CAMERA_all_results.csv`: CAMERA results for every gene set tested
- `results/GSEA_CAMERA_significant_results.csv`: gene sets with FDR < 0.05
- `results/session_info_gsea.txt`: R and package versions, including `msigdbr`
- `figures/voom_mean_variance.png`: voom mean-variance trend
- `figures/top_hallmark_sets.png`: top significant Hallmark gene sets

## Notes
- CAMERA is a test that accounts for inter-gene correlation
- Gene set contents depend on the MSigDB release. Newer `msigdbr` versions rename the `category` and `subcategory` arguments to `collection` and `subcollection`, so you will need to update these in the script if you see an error or warning. For me, the msigdbr version is: 26.1.0
- Gene sets are matched on Ensembl IDs. The script stops if none of the genes match

## Usage
Run `scripts/02_dge/02_dge_paired_oac_vs_normal.R` first, then run `scripts/05_gsea/05_gsea_camera.R` from the same working directory
