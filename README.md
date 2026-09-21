# Oesophageal-Adenocarcinoma-Bulk-RNA-sequencing
An analysis pipeline of bulk RNA-seq data from oesophageal adenocarcinoma patients. Performed using three pre-produced count matrices, one for each sample batch in data


### Count matrix files

The analysis scripts are written to work with **one or more count matrix files**.

By default, the scripts expect the following three files:

```r
count_files <- c(
  "Count Matrix Batch 1.xlsx",
  "Count Matrix Batch 2.xlsx",
  "Count Matrix Batch 3.xlsx"
)
```

If your data are provided as a different number of count matrix files, edit this section of the script so that it contains the filenames of your available count matrices.

For example, if you have **one count matrix file**:

```r
count_files <- "Count Matrix.xlsx"
```

If you have **two count matrix files**:

```r
count_files <- c(
  "Count Matrix 1.xlsx",
  "Count Matrix 2.xlsx"
)
```

If you have **three count matrix files**:

```r
count_files <- c(
  "Count Matrix Batch 1.xlsx",
  "Count Matrix Batch 2.xlsx",
  "Count Matrix Batch 3.xlsx"
)
```

The supplied count matrix file(s) must collectively contain the samples required for the analysis. The script will automatically combine the supplied matrices and check that all samples identified from the metadata are present.
