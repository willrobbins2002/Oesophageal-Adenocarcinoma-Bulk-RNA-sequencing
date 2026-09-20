## A PCA plot of QC-pass RNA-seq samples, including OAC, Barrett's, normal oesophagus, and normal stomach 
## Inputs: Moremetadata.xlsx, Count Matrix Batches 1-3.xlsx
## Outputs : PCA plot
## Requires: DESeq2, readxl, dplyr, ggplot2, stringr

library(DESeq2)
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

select <- dplyr::select

## Load in and clean metadata ----
metadata <- read_excel("Moremetadata.xlsx", sheet = "Sheet1") %>%
  mutate(
    TissueType = str_trim(TissueType),   # this fixes names, for example; "Normal Oesophagus " or " Normal Oesophagus" -> "Normal Oesophagus"
    TissueType = factor(TissueType, levels = c("Normal Stomach", "Normal Oesophagus",
                                               "Barrett's Oesophagus", "OAC")),
    Batch = factor(Batch)
  )

table(metadata$TissueType, useNA = "always")

## Building a MatrixID column that matches each count matrix's column names
## (For the count matrix data I had) Batch 1 columns are short codes (prefix before first "_"); Batch 2/3 columns match the SampleName exactly
metadata <- metadata %>%
  mutate(MatrixID = if_else(Batch == "1",
                            str_extract(SampleName, "^[^_]+"),
                            SampleName))

## Batch 1 has duplicate lane entries sharing one MatrixID 
## Confirmed TissueType and QCStatus are identical within each duplicate pair, so reduce to 1 row each 
metadata <- metadata %>%
  distinct(MatrixID, .keep_all = TRUE)

## Load and merge the three count matrices

counts1 <- read_excel("Count Matrix Batch 1.xlsx", sheet = "batch_1_count_matrix 1")
counts2 <- read_excel("Count Matrix Batch 2.xlsx", sheet = "batch_2_count_matrix 1")
counts3 <- read_excel("Count Matrix Batch 3.xlsx", sheet = "batch_3_count_matrix 1")

# GeneID order is identical across all three - this is confirming this before cbind
stopifnot(identical(counts1$GeneID, counts2$GeneID))
stopifnot(identical(counts1$GeneID, counts3$GeneID))

counts_merged <- bind_cols(counts1, counts2 %>% select(-GeneID), counts3 %>% select(-GeneID))

gene_ids <- counts_merged$GeneID
counts_merged <- counts_merged %>% select(-GeneID) %>% as.data.frame()
rownames(counts_merged) <- gene_ids

## Align metadata to count matrix columns
metadata <- metadata %>% filter(MatrixID %in% colnames(counts_merged))
counts_merged <- counts_merged[, metadata$MatrixID]  # reorder columns to match metadata row order

stopifnot(all(colnames(counts_merged) == metadata$MatrixID))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$MatrixID

## Filter for only samples which  passed QC (QC information was in my metadata)
metadata_pass <- metadata %>% filter(QCStatus == "Pass")
counts_pass <- counts_merged[, metadata_pass$MatrixID]

cat("Samples going into PCA:", ncol(counts_pass), "\n")
table(metadata_pass$TissueType)

## Build DESeq2 DataSet ----
dds <- DESeqDataSetFromMatrix(countData = round(counts_pass),
                              colData = metadata_pass,
                              design = ~ TissueType)

## Remove low-count genes before transformation
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

## Variance-stabilising transformation
vsd <- vst(dds, blind = TRUE)

## Producing a PCA plot on top 1000 variable genes
pcaData <- plotPCA(vsd, intgroup = c("TissueType", "Batch"),
                   ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p<- ggplot(pcaData, aes(x = PC1, y = PC2, color = TissueType)) +
  geom_point(size = 2.5, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("Normal Stomach" = "#4E79A7",
                                "Normal Oesophagus" = "#59A14F",
                                "Barrett's Oesophagus" = "#F28E2B",
                                "OAC" = "#E15759")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("PCA of QC-pass samples (top 1000 variable genes) by tissue type")

ggsave("pca_tissue_type.png", p, width = 7, height = 5, dpi = 300)