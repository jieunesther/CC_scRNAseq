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
library(fishpond)
library(Seurat)
library(SingleCellExperiment)

# Read in zUMIs data
S1 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S1_counts", outputFormat = "snRNA")
S2 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S2_counts", outputFormat = "snRNA")
S3 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S3_counts", outputFormat = "snRNA")
S4 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S4_counts", outputFormat = "snRNA")
S5 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S5_counts", outputFormat = "snRNA")
S6 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S6_counts", outputFormat = "snRNA")
S7 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S7_counts", outputFormat = "snRNA")
S8 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S8_counts", outputFormat = "snRNA")

# Rename cells to include sample number, and divide barcodes
S1.cnames = gsub(colnames(S1),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S1_\\1-\\2-\\3",perl=T)
colnames(S1) = S1.cnames
S2.cnames = gsub(colnames(S2),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S2_\\1-\\2-\\3",perl=T)
colnames(S2) = S2.cnames
S3.cnames = gsub(colnames(S3),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S3_\\1-\\2-\\3",perl=T)
colnames(S3) = S3.cnames
S4.cnames = gsub(colnames(S4),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S4_\\1-\\2-\\3",perl=T)
colnames(S4) = S4.cnames
S5.cnames = gsub(colnames(S5),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S5_\\1-\\2-\\3",perl=T)
colnames(S5) = S5.cnames
S6.cnames = gsub(colnames(S6),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S6_\\1-\\2-\\3",perl=T)
colnames(S6) = S6.cnames
S7.cnames = gsub(colnames(S7),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S7_\\1-\\2-\\3",perl=T)
colnames(S7) = S7.cnames
S8.cnames = gsub(colnames(S8),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="S8_\\1-\\2-\\3",perl=T)
colnames(S8) = S8.cnames


head(S1.cnames)
head(S2.cnames)
head(S3.cnames)
head(S4.cnames)
head(S5.cnames)
head(S6.cnames)
head(S7.cnames)
head(S8.cnames)

# Collapse transcripts to genes
tx2gene = read.table("/proj/zylkalab/Esther/genomeAnnotation/gencode.vM25.annotation.lookup.txt",header=F,sep="\t",col.names=c("tx","gene"))

S1.txId = rownames(S1)
S1.geneId <- as.vector(tx2gene$gene[match(S1.txId, tx2gene$tx)])
S1.tx.grp = t(sparse.model.matrix(~ 0 + S1.geneId))
S1.summarized <- S1.tx.grp %*% assay(S1)
rownames(S1.summarized) <- rownames(S1.summarized) %>% str_replace_all(".+.geneId","")

S2.txId = rownames(S2)
S2.geneId <- as.vector(tx2gene$gene[match(S2.txId, tx2gene$tx)])
S2.tx.grp = t(sparse.model.matrix(~ 0 + S2.geneId))
S2.summarized <- S2.tx.grp %*% assay(S2)
rownames(S2.summarized) <- rownames(S2.summarized) %>% str_replace_all(".+.geneId","")

S3.txId = rownames(S3)
S3.geneId <- as.vector(tx2gene$gene[match(S3.txId, tx2gene$tx)])
S3.tx.grp = t(sparse.model.matrix(~ 0 + S3.geneId))
S3.summarized <- S3.tx.grp %*% assay(S3)  
rownames(S3.summarized) <- rownames(S3.summarized) %>% str_replace_all(".+.geneId","")

S4.txId = rownames(S4)
S4.geneId <- as.vector(tx2gene$gene[match(S4.txId, tx2gene$tx)])
S4.tx.grp = t(sparse.model.matrix(~ 0 + S4.geneId))
S4.summarized <- S4.tx.grp %*% assay(S4)
rownames(S4.summarized) <- rownames(S4.summarized) %>% str_replace_all(".+.geneId","")

S5.txId = rownames(S5)
S5.geneId <- as.vector(tx2gene$gene[match(S5.txId, tx2gene$tx)])
S5.tx.grp = t(sparse.model.matrix(~ 0 + S5.geneId))
S5.summarized <- S5.tx.grp %*% assay(S5)
rownames(S5.summarized) <- rownames(S5.summarized) %>% str_replace_all(".+.geneId","")

S6.txId = rownames(S6)
S6.geneId <- as.vector(tx2gene$gene[match(S6.txId, tx2gene$tx)])
S6.tx.grp = t(sparse.model.matrix(~ 0 + S6.geneId))
S6.summarized <- S6.tx.grp %*% assay(S6)  
rownames(S6.summarized) <- rownames(S6.summarized) %>% str_replace_all(".+.geneId","")

S7.txId = rownames(S7)
S7.geneId <- as.vector(tx2gene$gene[match(S7.txId, tx2gene$tx)])
S7.tx.grp = t(sparse.model.matrix(~ 0 + S7.geneId))
S7.summarized <- S7.tx.grp %*% assay(S7)
rownames(S7.summarized) <- rownames(S7.summarized) %>% str_replace_all(".+.geneId","")

S8.txId = rownames(S8)
S8.geneId <- as.vector(tx2gene$gene[match(S8.txId, tx2gene$tx)])
S8.tx.grp = t(sparse.model.matrix(~ 0 + S8.geneId))
S8.summarized <- S8.tx.grp %*% assay(S8)  
rownames(S8.summarized) <- rownames(S8.summarized) %>% str_replace_all(".+.geneId","")

# Import data into Seurat objects
S1.seurat = CreateSeuratObject(S1.summarized,project="S1")
S2.seurat = CreateSeuratObject(S2.summarized,project="S2")
S3.seurat = CreateSeuratObject(S3.summarized,project="S3")
S4.seurat = CreateSeuratObject(S4.summarized,project="S4")
S5.seurat = CreateSeuratObject(S5.summarized,project="S5")
S6.seurat = CreateSeuratObject(S6.summarized,project="S6")
S7.seurat = CreateSeuratObject(S7.summarized,project="S7")
S8.seurat = CreateSeuratObject(S8.summarized,project="S8")

# Check raw dimensions
dim(S1.seurat)
dim(S2.seurat)
dim(S3.seurat)
dim(S4.seurat)
dim(S5.seurat)
dim(S6.seurat)
dim(S7.seurat)
dim(S8.seurat)

combined_1 = merge(x=S1.seurat, y=c(S2.seurat,S3.seurat,S4.seurat,S5.seurat,S6.seurat,S7.seurat,S8.seurat))

# Look up barcodes (last 8, this orientation) to get sample name
# Searches for a hamming distance within 1, but for this dataset, it seems to be all exact matches
# If no match, cell name will be prepended with "TRASH" for subsequent filtering
# Matches will have cell name prepended with GTrep-ampBC

bclookup=read.table("/proj/zylkalab/Esther/splitseq/SSSC081621/barcode_files/CC_day1_BC1.txt",header=F,sep="\t")

ids = colnames(combined_1@assays$RNA)
ids_l8 = substr(ids,22,29)
head(ids_l8)


newNames = ids
for(i in 1:length(ids_l8)) {
  ham = stringdist(ids_l8[i],bclookup$V1,method="hamming")
  ct = length(ham[ham<2])
  if(ct==0){
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else if(ct > 1){
    print(paste0("Barcode ",ids_l8[i]," had ",ct," matches within a Hamming distance of 1"))
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else{
    bestMatch = which.min(ham)
    matchingBC = as.character(bclookup[bestMatch,1])
    newNames[i] = paste0(as.character(bclookup[bestMatch,2]),"_",newNames[i])
  }
}
combined_1 = RenameCells(combined_1,new.names = newNames)

dim(combined_1)
saveRDS(combined_1,"/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined1_220316.rds")
combined_1 = readRDS("/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined1_220316.rds") # 55291 genes, 3280203 nuclei
combined_1 = subset(combined_1, subset = nCount_RNA > 500 & nFeature_RNA > 250)
dim(combined_1)# [1]  55291 238428
combined_1_1 = subset(combined_1, subset = nCount_RNA > 1000 & nFeature_RNA > 500)
combined_1_2 = subset(combined_1, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)
dim(combined_1_1) # [1]  55291 193720
dim(combined_1_2) # [1]  55291 155314
# Parse out sample metadata, drop TRASH cells
cnames = colnames(combined_1)
genotype = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\1",cnames,perl=T))
replicate = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\2",cnames,perl=T))
sex = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\3",cnames,perl=T))
wellID = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\4",cnames,perl=T))
ampBC = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\5",cnames,perl=T))
sublibrary = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\6",cnames,perl=T))
sampleNames = as.character(gsub("(.+)_(.+)_(.+)_(.+)_(.+)_(.+)_(.+)-(.+)-(.+)","\\1_\\2",cnames,perl=T)) 

combined_1 = AddMetaData(combined_1,metadata=genotype,col.name="genotype")
combined_1 = AddMetaData(combined_1,metadata=replicate,col.name="replicate")
combined_1 = AddMetaData(combined_1,metadata=sex,col.name="sex")
combined_1 = AddMetaData(combined_1,metadata=ampBC,col.name="ampBC")
combined_1 = AddMetaData(combined_1,metadata=sampleNames,col.name="sampleNames")
combined_1 = AddMetaData(combined_1,metadata=wellID,col.name="wellID")
combined_1 = AddMetaData(combined_1,metadata=sublibrary,col.name="sublibrary")
table(as.data.frame(ampBC_wellID))

# Read in zUMIs data
T1 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S9_counts", outputFormat = "snRNA")
T2 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S10_counts", outputFormat = "snRNA")
T3 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S11_counts", outputFormat = "snRNA")
T4 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S12_counts", outputFormat = "snRNA")
T5 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S13_counts", outputFormat = "snRNA")
T6 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S14_counts", outputFormat = "snRNA")
T7 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S15_counts", outputFormat = "snRNA")
T8 = loadFry(fryDir = "/pine/scr/j/i/jieunp/alevin-fry/CC/S16_counts", outputFormat = "snRNA")

# Rename cells to include sample number, and divide barcodes
T1.cnames = gsub(colnames(T1),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T1_\\1-\\2-\\3",perl=T)
colnames(T1) = T1.cnames
T2.cnames = gsub(colnames(T2),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T2_\\1-\\2-\\3",perl=T)
colnames(T2) = T2.cnames
T3.cnames = gsub(colnames(T3),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T3_\\1-\\2-\\3",perl=T)
colnames(T3) = T3.cnames
T4.cnames = gsub(colnames(T4),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T4_\\1-\\2-\\3",perl=T)
colnames(T4) = T4.cnames
T5.cnames = gsub(colnames(T5),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T5_\\1-\\2-\\3",perl=T)
colnames(T5) = T5.cnames
T6.cnames = gsub(colnames(T6),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T6_\\1-\\2-\\3",perl=T)
colnames(T6) = T6.cnames
T7.cnames = gsub(colnames(T7),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T7_\\1-\\2-\\3",perl=T)
colnames(T7) = T7.cnames
T8.cnames = gsub(colnames(T8),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="T8_\\1-\\2-\\3",perl=T)
colnames(T8) = T8.cnames


head(T1.cnames)
head(T2.cnames)
head(T3.cnames)
head(T4.cnames)
head(T5.cnames)
head(T6.cnames)
head(T7.cnames)
head(T8.cnames)

# Collapse transcripts to genes
tx2gene = read.table("/proj/zylkalab/Esther/genomeAnnotation/gencode.vM25.annotation.lookup.txt",header=F,sep="\t",col.names=c("tx","gene"))

T1.txId = rownames(T1)
T1.geneId <- as.vector(tx2gene$gene[match(T1.txId, tx2gene$tx)])
T1.tx.grp = t(sparse.model.matrix(~ 0 + T1.geneId))
T1.summarized <- T1.tx.grp %*% assay(T1)
rownames(T1.summarized) <- rownames(T1.summarized) %>% str_replace_all(".+.geneId","")

T2.txId = rownames(T2)
T2.geneId <- as.vector(tx2gene$gene[match(T2.txId, tx2gene$tx)])
T2.tx.grp = t(sparse.model.matrix(~ 0 + T2.geneId))
T2.summarized <- T2.tx.grp %*% assay(T2)
rownames(T2.summarized) <- rownames(T2.summarized) %>% str_replace_all(".+.geneId","")

T3.txId = rownames(T3)
T3.geneId <- as.vector(tx2gene$gene[match(T3.txId, tx2gene$tx)])
T3.tx.grp = t(sparse.model.matrix(~ 0 + T3.geneId))
T3.summarized <- T3.tx.grp %*% assay(T3)  
rownames(T3.summarized) <- rownames(T3.summarized) %>% str_replace_all(".+.geneId","")

T4.txId = rownames(T4)
T4.geneId <- as.vector(tx2gene$gene[match(T4.txId, tx2gene$tx)])
T4.tx.grp = t(sparse.model.matrix(~ 0 + T4.geneId))
T4.summarized <- T4.tx.grp %*% assay(T4)
rownames(T4.summarized) <- rownames(T4.summarized) %>% str_replace_all(".+.geneId","")

T5.txId = rownames(T5)
T5.geneId <- as.vector(tx2gene$gene[match(T5.txId, tx2gene$tx)])
T5.tx.grp = t(sparse.model.matrix(~ 0 + T5.geneId))
T5.summarized <- T5.tx.grp %*% assay(T5)
rownames(T5.summarized) <- rownames(T5.summarized) %>% str_replace_all(".+.geneId","")

T6.txId = rownames(T6)
T6.geneId <- as.vector(tx2gene$gene[match(T6.txId, tx2gene$tx)])
T6.tx.grp = t(sparse.model.matrix(~ 0 + T6.geneId))
T6.summarized <- T6.tx.grp %*% assay(T6)  
rownames(T6.summarized) <- rownames(T6.summarized) %>% str_replace_all(".+.geneId","")

T7.txId = rownames(T7)
T7.geneId <- as.vector(tx2gene$gene[match(T7.txId, tx2gene$tx)])
T7.tx.grp = t(sparse.model.matrix(~ 0 + T7.geneId))
T7.summarized <- T7.tx.grp %*% assay(T7)
rownames(T7.summarized) <- rownames(T7.summarized) %>% str_replace_all(".+.geneId","")

T8.txId = rownames(T8)
T8.geneId <- as.vector(tx2gene$gene[match(T8.txId, tx2gene$tx)])
T8.tx.grp = t(sparse.model.matrix(~ 0 + T8.geneId))
T8.summarized <- T8.tx.grp %*% assay(T8)  
rownames(T8.summarized) <- rownames(T8.summarized) %>% str_replace_all(".+.geneId","")

# Import data into Seurat objects
T1.seurat = CreateSeuratObject(T1.summarized,project="T1")
T2.seurat = CreateSeuratObject(T2.summarized,project="T2")
T3.seurat = CreateSeuratObject(T3.summarized,project="T3")
T4.seurat = CreateSeuratObject(T4.summarized,project="T4")
T5.seurat = CreateSeuratObject(T5.summarized,project="T5")
T6.seurat = CreateSeuratObject(T6.summarized,project="T6")
T7.seurat = CreateSeuratObject(T7.summarized,project="T7")
T8.seurat = CreateSeuratObject(T8.summarized,project="T8")

# Check raw dimensions
dim(T1.seurat)
dim(T2.seurat)
dim(T3.seurat)
dim(T4.seurat)
dim(T5.seurat)
dim(T6.seurat)
dim(T7.seurat)
dim(T8.seurat)

combined_2 = merge(x=T1.seurat, y=c(T2.seurat,T3.seurat,T4.seurat,T5.seurat,T6.seurat,T7.seurat,T8.seurat))
combined_2 = subset(combined_2, subset = nCount_RNA > 500 & nFeature_RNA > 250)
combined_2_2 = subset(combined_2, subset = nCount_RNA > 1000 & nFeature_RNA > 500)
combined_2_3 = subset(combined_2, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)
dim(combined_2) # [1]  55291 235546
dim(combined_2_2) # [1]  55291 182484
dim(combined_2_3) # [1]  55291 132692
# Look up barcodes (last 8, this orientation) to get sample name
# Searches for a hamming distance within 1, but for this dataset, it seems to be all exact matches
# If no match, cell name will be prepended with "TRASH" for subsequent filtering
# Matches will have cell name prepended with GTrep-ampBC

bclookup=read.table("/proj/zylkalab/Esther/splitseq/SSSC081621/barcode_files/CC_day2_BC1.txt",header=F,sep="\t")

ids = colnames(combined_2@assays$RNA)
ids_l8 = substr(ids,22,29)
head(ids_l8)


newNames = ids
for(i in 1:length(ids_l8)) {
  ham = stringdist(ids_l8[i],bclookup$V1,method="hamming")
  ct = length(ham[ham<2])
  if(ct==0){
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else if(ct > 1){
    print(paste0("Barcode ",ids_l8[i]," had ",ct," matches within a Hamming distance of 1"))
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else{
    bestMatch = which.min(ham)
    matchingBC = as.character(bclookup[bestMatch,1])
    newNames[i] = paste0(as.character(bclookup[bestMatch,2]),"_",newNames[i])
  }
}
combined_2 = RenameCells(combined_2,new.names = newNames)

dim(combined_2)
saveRDS(combined_2,"/proj/zylkalab/Esther/splitseq/Rmd_files/CC_splitseq/output/seurat/alevin-fry_combined2_220316.rds")
