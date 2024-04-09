# CC_snRNAseq

## Data summary
This repository contains information about the single nucleus RNA sequencing (snRNAseq) data generated from cerebral cortex (specifically, right cortical hemisphere) of the 14 collaborative cross mouse strains.
The scRNAseq data was generated using 92 animals from 14 different collaborative cross mouse strains to understand the influence of genetic background on brain cell type composition. We used Split Pool Ligation-based Transcriptome sequencing (SPLiT-seq) for generating the snRNAseq data, which utilizes a combinatorial barcoding technique to identify the cellular origin of RNA during snRNAseq analysis. This repository includes metadata, barcode sequences used for SPLiT-seq, scripts to generate the process snRNAseq seurat data (.rds file). 

## Contents
Barcode_data.txt: SPLiT-seq uses 4 sets of barcodes. In this study, the number of barcode sequences for each barcode set is 48 (BC1), 96 (BC2), 96 (BC3) and 16 (BC4; marks each sublibrary). This file contains sequences of BC1 (Round 1 v2), BC2 (Round 2 v1), and BC3 (Round 3 v1). 

Barcodesharing.txt: This file contains 48 BC1 oligodT barcodes and 48 BC1 random hexamer barcodes used in the study. Each sample was added into a well that has one oligodT barcode (BC1) and one random hexamer barcode (BC1). This barcode sharing files was used to collapse these two types of barcodes (oligodT and random hexamer).

CC_day1_BC1.txt: 92 animals were processed in two batches (day 1 and day 2). This file is for the samples processed in the first batch (day 1). It contains information about the the BC1 sequence (oligodT and random hexamer) for each animal ID and its well ID. 

CC_day2_BC1.txt: 92 animals were processed in two batches (day 1 and day 2). This file is for the samples processed in the second batch (day 2). It contains information about the the BC1 sequence (oligodT and random hexamer) for each animal ID and its well ID. 

fullList_barcode_combination_zylkalab.txt: this file contains all possible combination of barcode sequences (96 BC1 X 96 BC2 X 96 BC3 = 884726; 48 oligodT BC1 and 48 random hex BC1) for each sublibrary. This file is used to filter out any reads with unexpected barcode sequences (if more than 1 bp mismatch exists, the read is excluded). 

alevin-fry.sh: this is the script for alevin-fry to utilize *.fastq.gz file to generate gene x nucleus matrix. This script was submitted for each of the sublibrary (total 16 sublibrary). This is an example script for sublibrary 1 (S1). 

gencode.vM25.annotation.lookup.txt: this file is a reference genome annotation file (GENCODE M25) used to convert transcript id into gene names

build_seurat_object.R: this R script is to load gene x nucleus matrix into R and generate a seurat object

filter_seurat_object.R: this R script is to filter out low quality nuclei from the data

normalize_merge_dim_reduction_clustering.R: this R script is to SCTransform normalize, merge, run dimension reduction, and clustering. 


