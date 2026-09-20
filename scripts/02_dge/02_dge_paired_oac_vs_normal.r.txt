## Paired differential gene expression of OAC vs normal oesophagus (OAC and normal tissue coming from the same patients)
## DESeq2, paired by patient (design = ~ PatientID + TissueType)
## Inputs:   Moremetadata.xlsx
##           Count Matrix Batch 1.xlsx, Count Matrix Batch 2.xlsx, Count Matrix Batch 3.xlsx
##           Moremetadata.xlsx must contain: SampleName, SampleID, TissueType, Batch, QCStatus, TreatmentModality
## Outputs:  results/dds_de.rds
##           results/res_shrunk_df.rds (ashr-shrunken results)
##           results/dge_oac_vs_normal_oesophagus_shrunk.csv
##           results/dge_oac_vs_normal_oesophagus_unshrunk.csv
##           results/session_info_dge.txt
## Requires: DESeq2, ashr, readxl, dplyr, stringr, tibble
## Downstream scripts (top-gene tables, per-patient plots, plots of well-established genes in OAC) 
## Start from the saved .rds files, so run this script first

library(DESeq2)
library(ashr)
library(readxl)
library(dplyr)
library(stringr)
library(tibble)

## Set the parameters
data_dir    <- "."        ## folder containing the Excel input files
results_dir <- "results"  ## created if it doesn't exist
min_count   <- 10         ## minimum count per gene, set as 10 here
min_samples <- 3          ## ...in at least this many samples. set as 3 here
alpha       <- 0.05       ## adjusted p-value threshold for results, set as <0.05 here

dir.create(results_dir, showWarnings = FALSE)

## Load and clean metadata 
metadata <- read_excel(file.path(data_dir, "Moremetadata.xlsx"), sheet = "Sheet1") %>%
  mutate(
    TissueType = str_trim(TissueType),   # e.g. "Normal Oesophagus " -> "Normal Oesophagus"
    TissueType = factor(TissueType, levels = c("Normal Stomach", "Normal Oesophagus",
                                               "Barrett's Oesophagus", "OAC")),
    Batch = factor(Batch)
  )

required_cols <- c("SampleName", "SampleID", "TissueType", "Batch",
                   "QCStatus", "TreatmentModality")
stopifnot(all(required_cols %in% colnames(metadata)))

table(metadata$TissueType, useNA = "always")

## Build a MatrixID column matching the count matrix columns
## Batch 1 columns are short codes 
## Batch 2/3 columns match SampleName exactly
metadata <- metadata %>%
  mutate(MatrixID = if_else(Batch == "1",
                            str_extract(SampleName, "^[^_]+"),
                            SampleName))

## Batch 1 has duplicate lane entries sharing one MatrixID. TissueType and QCStatus are identical within each pair, so collapse to one row each.
metadata <- metadata %>%
  distinct(MatrixID, .keep_all = TRUE)

## Load and merge the three count matrices 
counts1 <- read_excel(file.path(data_dir, "Count Matrix Batch 1.xlsx"), sheet = "batch_1_count_matrix 1")
counts2 <- read_excel(file.path(data_dir, "Count Matrix Batch 2.xlsx"), sheet = "batch_2_count_matrix 1")
counts3 <- read_excel(file.path(data_dir, "Count Matrix Batch 3.xlsx"), sheet = "batch_3_count_matrix 1")

## GeneID order must be identical across all three before binding columns
stopifnot(identical(counts1$GeneID, counts2$GeneID))
stopifnot(identical(counts1$GeneID, counts3$GeneID))

## bind_cols() silently renames duplicated column names, so check first
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

## Align metadata to count matrix columns 
metadata <- metadata %>% filter(MatrixID %in% colnames(counts_merged))
counts_merged <- counts_merged[, metadata$MatrixID]   ## reorder the columns to match the metadata
stopifnot(all(colnames(counts_merged) == metadata$MatrixID))

## Subset to normal oesophagus + OAC, QC pass, treatment documented
## TreatmentModality, "Unknown" excludes samples whose treatment status was never recorded, and keeps chemoradiotherapy, chemotherapy, and treatment-naïve 
de_meta <- metadata %>%
  filter(QCStatus == "Pass") %>%
  filter(TissueType %in% c("Normal Oesophagus", "OAC")) %>%
  filter(TreatmentModality != "Unknown") %>%
  mutate(PatientID = SampleID) %>%
  mutate(TissueType = droplevels(TissueType))

## Keep only patients with both normal oesophagus AND OAC present
paired_patients <- de_meta %>%
  dplyr::count(PatientID, TissueType) %>%
  dplyr::count(PatientID) %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(PatientID)

de_meta <- de_meta %>% filter(PatientID %in% paired_patients)

cat("Matched patients:", length(paired_patients), "\n")
cat("Samples going into paired DE:", nrow(de_meta), "\n")
print(table(de_meta$TissueType))

## Duplicate patient/tissue combinations (e.g. two OAC biopsies for one patient)
dup_check <- de_meta %>% dplyr::count(PatientID, TissueType) %>% dplyr::filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Some patients have more than one sample of the same tissue type; decide how to handle these.")
  print(dup_check)
}

## Align counts to this subset
de_meta <- as.data.frame(de_meta)
rownames(de_meta) <- de_meta$MatrixID
counts_de <- counts_merged[, de_meta$MatrixID]
stopifnot(all(colnames(counts_de) == de_meta$MatrixID))

## Build paired DESeqDataSet 
de_meta$PatientID  <- factor(de_meta$PatientID)
de_meta$TissueType <- relevel(factor(de_meta$TissueType), ref = "Normal Oesophagus")

dds_de <- DESeqDataSetFromMatrix(countData = round(counts_de),
                                 colData   = de_meta,
                                 design    = ~ PatientID + TissueType)

## Drop low-count genes
keep <- rowSums(counts(dds_de) >= min_count) >= min_samples
dds_de <- dds_de[keep, ]
cat("Genes retained after low-count filter:", nrow(dds_de), "\n")

## Run DESeq2 
dds_de <- DESeq(dds_de)

## Results: OAC vs Normal Oesophagus 
res <- results(dds_de,
               contrast = c("TissueType", "OAC", "Normal Oesophagus"),
               alpha = alpha)
summary(res)

## Shrink log2FC (this is more reliable for ranking or plotting) 
res_shrunk <- lfcShrink(dds_de,
                        contrast = c("TissueType", "OAC", "Normal Oesophagus"),
                        res = res, type = "ashr")

res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df <- res_shrunk_df[order(res_shrunk_df$padj), ]
head(res_shrunk_df, 20)

## Save outputs for downstream scripts 
saveRDS(dds_de,        file.path(results_dir, "dds_de.rds"))
saveRDS(res_shrunk_df, file.path(results_dir, "res_shrunk_df.rds"))

## Both versions are saved: shrunken (ashr) for ranking or plotting, and the unshrunk MLE results from results()
write.csv(rownames_to_column(res_shrunk_df, "EnsemblID_full"),
          file.path(results_dir, "dge_oac_vs_normal_oesophagus_shrunk.csv"),
          row.names = FALSE)

res_unshrunk_df <- as.data.frame(res)
res_unshrunk_df <- res_unshrunk_df[order(res_unshrunk_df$padj), ]
write.csv(rownames_to_column(res_unshrunk_df, "EnsemblID_full"),
          file.path(results_dir, "dge_oac_vs_normal_oesophagus_unshrunk.csv"),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "session_info_dge.txt"))
