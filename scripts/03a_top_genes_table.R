## Install the annotation package 
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}

library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)      
library(stringr)
library(tibble)

## Parameters 
results_dir <- "results"
padj_cutoff <- 0.05   # adjusted p-value threshold
lfc_cutoff  <- 1      # minimum |shrunken log2 fold change|
n_top       <- 10     # genes reported per direction

## Load DGE results 
res_path <- file.path(results_dir, "res_shrunk_df.rds")
if (!file.exists(res_path)) {
  stop("Can't find ", res_path, ". Run 02_dge_paired_oac_vs_normal.R first.")
}
res_shrunk_df <- readRDS(res_path)

## Map Ensembl IDs to gene symbols 
## Ensembl IDs include a version suffix (e.g. ".16") which org.Hs.eg.db doesn't recognise, so it is stripped before mapping
res_annotated <- res_shrunk_df %>%
  rownames_to_column("EnsemblID_full") %>%
  mutate(EnsemblID = str_remove(EnsemblID_full, "\\..*"))

res_annotated$GeneSymbol <- unname(mapIds(org.Hs.eg.db,
                                          keys      = res_annotated$EnsemblID,
                                          column    = "SYMBOL",
                                          keytype   = "ENSEMBL",
                                          multiVals = "first"))

## Some Ensembl IDs won't map (non-coding RNAs, pseudogenes, or IDs missing from this annotation build), so a number of NAs is expected
cat("Genes without a gene symbol:", sum(is.na(res_annotated$GeneSymbol)),
    "of", nrow(res_annotated), "\n")

## Select top genes 
format_table <- function(df) {
  df %>%
    dplyr::select(GeneSymbol, EnsemblID_full, baseMean, log2FoldChange, padj) %>%
    mutate(across(c(baseMean, log2FoldChange), ~ round(.x, 2)),
           padj = signif(padj, 3))
}

top_up <- res_annotated %>%
  filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = n_top) %>%
  format_table()

top_down <- res_annotated %>%
  filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = n_top) %>%
  format_table()

cat("Genes passing the rule - up:",
    sum(res_annotated$padj < padj_cutoff & res_annotated$log2FoldChange >  lfc_cutoff, na.rm = TRUE),
    "| down:",
    sum(res_annotated$padj < padj_cutoff & res_annotated$log2FoldChange < -lfc_cutoff, na.rm = TRUE), "\n")

print(top_up)
print(top_down)

## Save outputs 
dir.create(results_dir, showWarnings = FALSE)
saveRDS(res_annotated, file.path(results_dir, "res_annotated.rds"))   # reused by later scripts
write.csv(top_up,   file.path(results_dir, "top_up_genes.csv"),   row.names = FALSE)
write.csv(top_down, file.path(results_dir, "top_down_genes.csv"), row.names = FALSE)