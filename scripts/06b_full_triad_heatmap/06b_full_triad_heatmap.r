## LOAD PACKAGES

library(DESeq2)
library(dplyr)
library(pheatmap)

## FILE LOCATIONS

input_dir <- "Full_Triad_Analysis"

dds_file <- file.path(
  input_dir,
  "dds_full_triad.rds"
)

results_file <- file.path(
  input_dir,
  "DESeq2_OAC_vs_Normal.csv"
)

metadata_file <- file.path(
  input_dir,
  "full_triad_metadata.csv"
)

output_file <- file.path(
  input_dir,
  "Heatmap_Top20_DEGs_Full_Triad.png"
)

## CHECK INPUT FILES

required_files <- c(
  dds_file,
  results_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Missing files:\n\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )

}

## LOAD DATA

dds <- readRDS(dds_file)

deseq_results <- read.csv(
  results_file,
  stringsAsFactors = FALSE
)

triad_metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

## CHECK SAMPLE MATCHING

cat(
  "\nSamples in DESeq2 object:",
  ncol(dds),
  "\n"
)

cat(
  "Metadata rows:",
  nrow(triad_metadata),
  "\n"
)

missing_metadata <- setdiff(
  colnames(dds),
  triad_metadata$MatrixID
)

if (length(missing_metadata) > 0) {

  stop(
    "Samples missing from metadata:\n\n",
    paste(
      missing_metadata,
      collapse = "\n"
    )
  )

}

## VST TRANSFORMATION

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(vsd)

## SELECT TOP 20 DEGs

top_genes <- deseq_results %>%

  filter(
    !is.na(padj),
    !is.na(log2FoldChange)
  ) %>%

  arrange(
    padj
  ) %>%

  slice_head(
    n = 20
  )

cat(
  "\nTop genes selected:\n"
)

print(
  top_genes %>%
    select(
      EnsemblID,
      Gene,
      log2FoldChange,
      padj
    )
)

## EXTRACT EXPRESSION MATRIX

heatmap_ensembl_ids <- top_genes$EnsemblID[
  top_genes$EnsemblID %in%
    rownames(vst_matrix)
]

if (length(heatmap_ensembl_ids) == 0) {

  stop(
    "No selected genes found in VST matrix."
  )

}

heatmap_mat <- vst_matrix[
  heatmap_ensembl_ids,
  ,
  drop = FALSE
]

## CREATE GENE LABELS

gene_labels <- top_genes %>%

  filter(
    EnsemblID %in%
      heatmap_ensembl_ids
  ) %>%

  arrange(
    match(
      EnsemblID,
      heatmap_ensembl_ids
    )
  ) %>%

  mutate(

    display_label = ifelse(
      !is.na(Gene) &
        Gene != "",
      Gene,
      EnsemblID
    )

  ) %>%

  pull(
    display_label
  )

rownames(
  heatmap_mat
) <- gene_labels

## PREPARE METADATA

annotation_col <- triad_metadata %>%

  select(
    MatrixID,
    SampleID,
    TissueType
  ) %>%

  mutate(

    TissueType = factor(
      TissueType,
      levels = c(
        "Normal Oesophagus",
        "Barrett's Oesophagus",
        "OAC"
      )
    )

  )

rownames(
  annotation_col
) <- annotation_col$MatrixID

annotation_col <- annotation_col[
  colnames(heatmap_mat),
  ,
  drop = FALSE
]

## ORDER SAMPLES BY PATIENT

sample_order <- annotation_col %>%

  arrange(
    SampleID,
    TissueType
  ) %>%

  rownames()

heatmap_mat <- heatmap_mat[
  ,
  sample_order,
  drop = FALSE
]

annotation_col <- annotation_col[
  sample_order,
  ,
  drop = FALSE
]

## CREATE PATIENT GAPS

patient_sizes <- table(
  annotation_col$SampleID
)

gaps_col <- cumsum(
  patient_sizes
)

gaps_col <- gaps_col[
  -length(gaps_col)
]

## SIMPLIFY ANNOTATIONS

annotation_display <- annotation_col %>%

  select(
    SampleID,
    TissueType
  )

annotation_colours <- list(

  TissueType = c(
    "Normal Oesophagus" = "#1F78B4",
    "Barrett's Oesophagus" = "#9467BD",
    "OAC" = "#E377C2"
  )

)

## COLOUR PALETTE

heatmap_colours <- colorRampPalette(

  c(
    "#2166AC",
    "white",
    "#B2182B"
  )

)(100)

## GENERATE HEATMAP

pheatmap(

  mat = heatmap_mat,

  scale = "row",

  cluster_rows = TRUE,

  cluster_cols = FALSE,

  annotation_col = annotation_display,

  annotation_colors = annotation_colours,

  gaps_col = gaps_col,

  show_colnames = FALSE,

  border_color = NA,

  color = heatmap_colours,

  fontsize_row = 10,

  fontsize_col = 8,

  main =
    "Top 20 DEGs Across Matched Normal, Barrett's and OAC Samples",

  filename = output_file,

  width = 12,

  height = 8

)

## SAVE SESSION INFO

writeLines(

  capture.output(
    sessionInfo()
  ),

  file.path(
    input_dir,
    "session_info_heatmap.txt"
  )

)

## SUMMARY

cat(
  "\n====================================\n"
)

cat(
  "FULL TRIAD HEATMAP COMPLETE\n"
)

cat(
  "====================================\n\n"
)

cat(
  "Patients:",
  length(
    unique(
      triad_metadata$SampleID
    )
  ),
  "\n"
)

cat(
  "Samples:",
  ncol(heatmap_mat),
  "\n"
)

cat(
  "Genes displayed:",
  nrow(heatmap_mat),
  "\n\n"
)

cat(
  "Top genes selected from:\n"
)

cat(
  "OAC vs Normal Oesophagus DESeq2 comparison\n\n"
)

cat(
  "Output:\n",
  output_file,
  "\n"
)
