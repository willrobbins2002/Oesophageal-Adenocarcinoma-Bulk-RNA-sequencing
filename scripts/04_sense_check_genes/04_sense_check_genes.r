
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)
library(dplyr)     
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)

## Parameters 
data_dir    <- "."         # folder containing the Excel input files
figures_dir <- "figures"   # created if it doesn't exist

tissue_colours <- c("Normal Stomach"       = "#4E79A7",
                    "Normal Oesophagus"    = "#59A14F",
                    "Barrett's Oesophagus" = "#F28E2B",
                    "OAC"                  = "#E15759")

dir.create(figures_dir, showWarnings = FALSE)

## Load and clean metadata 
metadata <- read_excel(file.path(data_dir, "Moremetadata.xlsx"), sheet = "Sheet1") %>%
  mutate(
    TissueType = str_trim(TissueType),   # e.g. "Normal Oesophagus " -> "Normal Oesophagus"
    TissueType = factor(TissueType, levels = names(tissue_colours)),
    Batch = factor(Batch)
  )

stopifnot(all(c("SampleName", "TissueType", "Batch", "QCStatus") %in% colnames(metadata)))
table(metadata$TissueType, useNA = "always")

## Build a MatrixID column matching the count matrix columns.
## Batch 1 columns are short codes 
## Batch 2/3 columns match SampleName exactly
metadata <- metadata %>%
  mutate(MatrixID = if_else(Batch == "1",
                            str_extract(SampleName, "^[^_]+"),
                            SampleName))

## Batch 1 has duplicate lane entries sharing one MatrixID (already summed in the count matrix), so collapse to one row each
metadata <- metadata %>%
  distinct(MatrixID, .keep_all = TRUE)

## Load and merge the three count matrices 
counts1 <- read_excel(file.path(data_dir, "Count Matrix Batch 1.xlsx"), sheet = "batch_1_count_matrix 1")
counts2 <- read_excel(file.path(data_dir, "Count Matrix Batch 2.xlsx"), sheet = "batch_2_count_matrix 1")
counts3 <- read_excel(file.path(data_dir, "Count Matrix Batch 3.xlsx"), sheet = "batch_3_count_matrix 1")

stopifnot(identical(counts1$GeneID, counts2$GeneID))
stopifnot(identical(counts1$GeneID, counts3$GeneID))

all_sample_cols <- c(setdiff(colnames(counts1), "GeneID"),
                     setdiff(colnames(counts2), "GeneID"),
                     setdiff(colnames(counts3), "GeneID"))
stopifnot(!anyDuplicated(all_sample_cols))

counts_merged <- bind_cols(counts1,
                           counts2 %>% dplyr::select(-GeneID),
                           counts3 %>% dplyr::select(-GeneID))

gene_ids <- counts_merged$GeneID
counts_merged <- counts_merged %>% dplyr::select(-GeneID) %>% as.data.frame()
rownames(counts_merged) <- gene_ids

## Align metadata to counts and keep QC-pass samples 
metadata <- metadata %>% filter(MatrixID %in% colnames(counts_merged))
counts_merged <- counts_merged[, metadata$MatrixID]
stopifnot(all(colnames(counts_merged) == metadata$MatrixID))

metadata_pass <- metadata %>% filter(QCStatus == "Pass") %>% as.data.frame()
rownames(metadata_pass) <- metadata_pass$MatrixID
counts_pass <- counts_merged[, metadata_pass$MatrixID]
stopifnot(all(colnames(counts_pass) == metadata_pass$MatrixID))

cat("QC-pass samples:", ncol(counts_pass), "\n")
print(table(metadata_pass$TissueType))

## Define the marker panel 
marker_genes <- tibble::tibble(
  Symbol = c("TFF3", "ANXA10", "LGALS4", "ERBB2", "GATA4", "GATA6", "VEGFA",
             "CDX2", "MUC2", "KRT7", "KRT20"),
  Source = c("Upregulated in BO",
             "Upregulated in BO",
             "Upregulated in BO",
             "Amplified in OAC",
             "Amplified in OAC",
             "Amplified in OAC",
             "Amplified in OAC",
             "Master regulator of intestinal differentiation - canonical BO marker",
             "Intestinal goblet-cell mucin - hallmark of BO metaplasia",
             "Columnar/glandular cytokeratin - BO/gastric differential marker",
             "Gastric-type cytokeratin - BO/gastric differential marker")
)

## Map gene symbols to the Ensembl IDs in the count matrix 
ensembl_map <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys    = marker_genes$Symbol,
                                     keytype = "SYMBOL",
                                     columns = "ENSEMBL")

## Ensembl IDs in the count matrix carry a version suffix, which org.Hs.eg.db doesn't recognise, so strip it before matching
gene_lookup <- data.frame(EnsemblID_full = rownames(counts_pass),
                          EnsemblID      = str_remove(rownames(counts_pass), "\\..*")) %>%
  inner_join(ensembl_map, by = c("EnsemblID" = "ENSEMBL")) %>%
  inner_join(marker_genes, by = c("SYMBOL" = "Symbol"))

## Every marker should be found once, so investigate anything flagged here
missing_genes <- setdiff(marker_genes$Symbol, gene_lookup$SYMBOL)
if (length(missing_genes) > 0) {
  warning("Marker genes not found in the count matrix: ",
          paste(missing_genes, collapse = ", "), immediate. = TRUE)
}
dup_symbols <- unique(gene_lookup$SYMBOL[duplicated(gene_lookup$SYMBOL)])
if (length(dup_symbols) > 0) {
  warning("Symbols matching more than one Ensembl ID (they would share a panel): ",
          paste(dup_symbols, collapse = ", "), immediate. = TRUE)
}
print(gene_lookup)

## Normalise counts (all QC-pass samples, all tissues) 
## The design is only needed to build the object; size factors are estimated from all genes, with no low-count filter
dds_full <- DESeqDataSetFromMatrix(countData = round(counts_pass),
                                   colData   = metadata_pass,
                                   design    = ~ TissueType)
dds_full <- estimateSizeFactors(dds_full)
norm_counts <- counts(dds_full, normalized = TRUE)

## Extract marker expression and reshape for plotting 
marker_expr <- as.data.frame(norm_counts[gene_lookup$EnsemblID_full, , drop = FALSE])
marker_expr$EnsemblID_full <- rownames(marker_expr)
marker_expr <- marker_expr %>% left_join(gene_lookup, by = "EnsemblID_full")

marker_long <- marker_expr %>%
  pivot_longer(cols = all_of(metadata_pass$MatrixID),
               names_to = "MatrixID", values_to = "NormCount") %>%
  left_join(metadata_pass %>% dplyr::select(MatrixID, TissueType), by = "MatrixID") %>%
  mutate(log2Count = log2(NormCount + 1),
         SYMBOL = factor(SYMBOL, levels = intersect(marker_genes$Symbol, unique(SYMBOL))))

## Plot expression of each marker across tissue groups 
p <- ggplot(marker_long, aes(x = TissueType, y = log2Count, fill = TissueType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  scale_fill_manual(values = tissue_colours) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "log2(normalized count + 1)",
       title = "Expression of literature-reported marker genes across tissue groups")

print(p)
ggsave(file.path(figures_dir, "sense_check_genes.png"), p,
       width = 12, height = 9, dpi = 300)
