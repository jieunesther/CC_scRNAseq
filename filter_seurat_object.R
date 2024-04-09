library(rlang)
library(vctrs)
library(backports)
library(cli)
library(Seurat)
library(rstudioapi)
library(withr)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(farver)
library(usethis)
library(labeling)
library(Matrix)
library(stringdist)

CC_combined = readRDS("/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined_filtered_220325.rds")
CC_combined <- PercentageFeatureSet(object = CC_combined, pattern = "^mt-", col.name = "percent.mt")
combined.filtered <- subset(CC_combined, subset = nCount_RNA > 1000 & nCount_RNA < 100000 & nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
combined.filtered

# Gene-level filtering (remove genes that have zero expression in all cells, then keep only genes which are expressed in 30 or more cells)
# remove genes not expressed in at least 30 cells, if less than 3 transcripts (have not done this)
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = combined.filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 30 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 30

# Only keeping those genes expressed in more than 30 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
combined.filtered <- CreateSeuratObject(filtered_counts, meta.data = combined.filtered@meta.data)
