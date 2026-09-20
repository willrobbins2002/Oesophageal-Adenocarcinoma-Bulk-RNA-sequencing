## OAC RNA-seq Top Up- and Down-regulated Genes

Tables of the most up- and down-regulated genes in OAC compared with normal oesophagus, taken from the paired DGE analysis

## What the script does
1. Loads the shrunken DGE results from the paired analysis (`02_dge`)
2. Maps Ensembl IDs to gene symbols, with the version suffix removed 
3. Selects genes with adjusted p-value < 0.05 and absolute shrunken log2 fold change > 1
4. Ranks them by shrunken log2 fold change and reports the top 10 upregulated and top 10 downregulated genes
5. Saves the annotated results and both tables

## Requirements
R packages: org.Hs.eg.db, AnnotationDbi, dplyr, stringr, tibble

## Input files (not included)
- `results/res_shrunk_df.rds`, created by `scripts/02_dge/02_dge_paired_oac_vs_normal.R`

## Outputs (not included)
Written to a `results/` folder:
- `res_annotated.rds`: DGE results with gene symbols added
- `top_up_genes.csv` and `top_down_genes.csv`: top 10 genes per direction

## Usage
Run `scripts/02_dge/02_dge_paired_oac_vs_normal.R` first, then run `scripts/03_top_genes/03a_top_genes_table.R` from the same working directory.