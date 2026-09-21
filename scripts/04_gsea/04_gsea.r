library(DESeq2)
library(edgeR)
library(limma)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)

## Parameters 
results_dir    <- "results"
figures_dir    <- "figures"
contrast_name  <- "TissueTypeOAC"
fdr_cutoff     <- 0.05             ## Benjamini-Hochberg FDR threshold
inter_gene_cor <- 0.01             ## limma's default; set to NA to estimate per gene set
n_top_plot     <- 15               ## Hallmark sets shown in the barplot

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

## ---- 1. Load counts and metadata from the paired DESeq2 object ----
dds_path <- file.path(results_dir, "dds_de.rds")
if (!file.exists(dds_path)) {
  stop("Can't find ", dds_path, ". Run 02_dge_paired_oac_vs_normal.R first.")
}
dds_de <- readRDS(dds_path)

raw_counts <- counts(dds_de, normalized = FALSE)
coldata    <- as.data.frame(colData(dds_de))

## Strip Ensembl version suffixes (e.g. "ENSG00000141510.15" -> "ENSG00000141510"). Msigdbr's ensembl_gene column has no version numbers, so gene set matching fails without this
rownames(raw_counts) <- sub("\\..*$", "", rownames(raw_counts))

## If stripping versions made two IDs identical, sum their counts
if (any(duplicated(rownames(raw_counts)))) {
  cat("Duplicate Ensembl IDs found after stripping versions - collapsing by sum\n")
  raw_counts <- rowsum(raw_counts, group = rownames(raw_counts))
}

## Sample order must match exactly between counts and metadata
stopifnot(all(colnames(raw_counts) == rownames(coldata)))

cat("Samples in GSEA input:", ncol(raw_counts), "\n")
cat("Patients:", n_distinct(coldata$PatientID), "\n")
print(table(coldata$TissueType))

## Design matrix (paired: PatientID + TissueType)
## PatientID and TissueType are already factors, with Normal Oesophagus as the reference level, from dds_de
design <- model.matrix(~ PatientID + TissueType, data = coldata)

contrast_coef <- which(colnames(design) == contrast_name)
stopifnot(length(contrast_coef) == 1)   # contrast_name must match colnames(design)
cat("Testing design column:", colnames(design)[contrast_coef], "\n")

## Filter low-expression genes and normalise 
dge <- DGEList(counts = raw_counts)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes retained after filtering:", sum(keep), "out of", length(keep), "\n")

dge <- calcNormFactors(dge)   

## Voom transformation 
## CAMERA can be run on a voom EList, which incorporates the mean-variance relationship, rather than on plain logCPM
png(file.path(figures_dir, "voom_mean_variance.png"), width = 1600, height = 1200, res = 300)
v <- voom(dge, design, plot = TRUE)
dev.off()

## Retrieve gene sets from msigdbr 
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
go_bp_sets    <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_cc_sets    <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
go_mf_sets    <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")

## Convert msigdbr's long format to a named list of gene IDs per set. The count matrix uses Ensembl IDs, so the default is "ensembl_gene"
make_gene_set_list <- function(msigdbr_df, id_col = "ensembl_gene") {
  split(msigdbr_df[[id_col]], msigdbr_df$gs_name)
}

hallmark_list <- make_gene_set_list(hallmark_sets)
reactome_list <- make_gene_set_list(reactome_sets)
go_bp_list    <- make_gene_set_list(go_bp_sets)
go_cc_list    <- make_gene_set_list(go_cc_sets)
go_mf_list    <- make_gene_set_list(go_mf_sets)

## Map gene sets to row indices of the voom object 
## Check that gene IDs match before going any further: if none of the voom row names appear in the gene sets, the ID types don't match
n_matched <- sum(rownames(v) %in% unlist(hallmark_list))
cat("Genes in the voom object found in Hallmark sets:", n_matched, "of", nrow(v), "\n")
stopifnot(n_matched > 0)

idx_hallmark <- ids2indices(hallmark_list, identifiers = rownames(v))
idx_reactome <- ids2indices(reactome_list, identifiers = rownames(v))
idx_go_bp    <- ids2indices(go_bp_list,    identifiers = rownames(v))
idx_go_cc    <- ids2indices(go_cc_list,    identifiers = rownames(v))
idx_go_mf    <- ids2indices(go_mf_list,    identifiers = rownames(v))

## Run CAMERA for each gene set collection 
run_camera <- function(idx_list, label) {
  camera(v, idx_list, design,
         contrast = contrast_coef,
         inter.gene.cor = inter_gene_cor) %>%
    rownames_to_column("GeneSet") %>%
    mutate(Collection = label) %>%
    arrange(FDR)
}

camera_hallmark <- run_camera(idx_hallmark, "Hallmark")
camera_reactome <- run_camera(idx_reactome, "Reactome")
camera_go_bp    <- run_camera(idx_go_bp,    "GO_BP")
camera_go_cc    <- run_camera(idx_go_cc,    "GO_CC")
camera_go_mf    <- run_camera(idx_go_mf,    "GO_MF")

## Significant gene sets 
sig_hallmark <- camera_hallmark %>% filter(FDR < fdr_cutoff)
sig_reactome <- camera_reactome %>% filter(FDR < fdr_cutoff)
sig_go_bp    <- camera_go_bp    %>% filter(FDR < fdr_cutoff)
sig_go_cc    <- camera_go_cc    %>% filter(FDR < fdr_cutoff)
sig_go_mf    <- camera_go_mf    %>% filter(FDR < fdr_cutoff)

cat("Significant Hallmark sets:", nrow(sig_hallmark), "\n")
cat("Significant Reactome sets:", nrow(sig_reactome), "\n")
cat("Significant GO:BP sets:",    nrow(sig_go_bp), "\n")
cat("Significant GO:CC sets:",    nrow(sig_go_cc), "\n")
cat("Significant GO:MF sets:",    nrow(sig_go_mf), "\n")

## Combine and export results 
all_camera_results <- bind_rows(
  camera_hallmark, camera_reactome, camera_go_bp, camera_go_cc, camera_go_mf
)
write.csv(all_camera_results,
          file.path(results_dir, "GSEA_CAMERA_all_results.csv"), row.names = FALSE)

all_sig_results <- all_camera_results %>% filter(FDR < fdr_cutoff)
write.csv(all_sig_results,
          file.path(results_dir, "GSEA_CAMERA_significant_results.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "session_info_gsea.txt"))

## Barplot of the top significant Hallmark sets 
top_hallmark <- sig_hallmark %>%
  arrange(FDR) %>%
  slice_head(n = n_top_plot)

if (nrow(top_hallmark) > 0) {
  p_hallmark <- ggplot(top_hallmark,
                       aes(x = reorder(GeneSet, -log10(FDR)), y = -log10(FDR),
                           fill = Direction)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "-log10(FDR)", title = "Top significant Hallmark gene sets") +
    theme_minimal()

  print(p_hallmark)
  ggsave(file.path(figures_dir, "top_hallmark_sets.png"), p_hallmark,
         width = 9, height = 6, dpi = 300)
} else {
  cat("No significant Hallmark sets at FDR <", fdr_cutoff, "- barplot skipped\n")
}
