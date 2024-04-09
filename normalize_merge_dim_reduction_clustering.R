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

#options(future.globals.maxSize=5242880000)
all.list <- SplitObject(combined.filtered, split.by = "sampleNames")

#for (i in 1:length(x = all.list)) {
  all.list[[i]] <- PercentageFeatureSet(object = all.list[[i]], pattern = "^mt-", col.name = "percent.mt")
  all.list[[i]] <- SCTransform(all.list[[i]], vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes=FALSE)
}

saveRDS(all.list, "/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined_filtered_220325_SCTransform.rds")
rm(list=ls())

all.list = readRDS("/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined_filtered_220325_SCTransform.rds")
combined = merge(all.list[[1]], y = c(all.list[[2]], all.list[[3]], all.list[[4]], all.list[[5]], 
                                     all.list[[6]], all.list [[7]], all.list[[8]], all.list[[9]], all.list[[10]]))
combined = merge(combined, y = c(all.list[[11]], all.list[[12]], all.list[[13]], all.list[[14]], all.list[[15]], 
                                 all.list[[16]], all.list[[17]], all.list[[18]], all.list[[19]], all.list[[20]]))
combined = merge(combined, y = c(all.list[[21]], all.list[[22]], all.list[[23]], all.list[[24]], all.list[[25]], 
                                 all.list[[26]], all.list[[27]], all.list[[28]], all.list[[29]], all.list[[30]]))
combined = merge(combined, y = c(all.list[[31]], all.list[[32]], all.list[[33]], all.list[[34]], all.list[[35]], 
                                 all.list[[36]], all.list[[37]], all.list[[38]], all.list[[39]], all.list[[40]]))
combined = merge(combined, y = c(all.list[[41]], all.list[[42]], all.list[[43]], all.list[[44]], all.list[[45]], 
                                 all.list[[46]], all.list[[47]], all.list[[48]], all.list[[49]], all.list[[50]]))
combined = merge(combined, y = c(all.list[[51]], all.list[[52]], all.list[[53]], all.list[[54]], all.list[[55]], 
                                 all.list[[56]], all.list[[57]], all.list[[58]], all.list[[59]], all.list[[60]]))
combined = merge(combined, y = c(all.list[[61]], all.list[[62]], all.list[[63]], all.list[[64]], all.list[[65]], 
                                 all.list[[66]], all.list[[67]], all.list[[68]], all.list[[69]], all.list[[70]]))
combined = merge(combined, y = c(all.list[[71]], all.list[[72]], all.list[[73]], all.list[[74]], all.list[[75]], 
                                 all.list[[76]], all.list[[77]], all.list[[78]], all.list[[79]], all.list[[80]]))
combined = merge(combined, y = c(all.list[[81]], all.list[[82]], all.list[[83]], all.list[[84]], all.list[[85]], 
                                 all.list[[86]], all.list[[87]], all.list[[88]], all.list[[89]], all.list[[90]]))
combined = merge(combined, y = c(all.list[[91]], all.list[[92]], all.list[[93]], all.list[[94]]))

saveRDS(combined, "/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined_filtered_220325_SCTransform_merged.rds")

CC.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 5000)
VariableFeatures(combined) <- CC.features

combined <- RunPCA(combined, verbose = TRUE, npcs = 100)
combined <- RunUMAP(combined, dims = 1:100)
combined <- FindNeighbors(combined, dims = 1:100, verbose = TRUE)
combined <- FindClusters(combined, verbose = TRUE, resolution = c(0.2, 0.3, 0.4, 0.5, 0.8, 1), algorithm=2)

saveRDS(combined,"/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined_filtered_220325_SCTransform_merged_clustered.rds")
