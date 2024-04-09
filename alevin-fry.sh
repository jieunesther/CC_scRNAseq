#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem 20g
#SBATCH --job-name=alevinfry
#SBATCH --mail-user=esther_park@med.unc.edu
#SBATCH --mail-type=ALL
#SBATCH -t 99:00:00 
#SBATCH --output 2022-3-15_CCS1_alevinfry_barcodelist.%J.log  

#splitp -r /proj/zylkalab/Esther/splitseq/CC_JN_rawdata/00001_CAGATC_S1_ME_L001_R2_001.fastq.gz -b /proj/zylkalab/Esther/splitseq/SSSC081621/barcode_files/alevin-fry_barcodesharing_all48.txt -s 87 -e 94 -o | gzip > /pine/scr/j/i/jieunp/alevin-fry/CC/S1_R2_barcodecollapsed.fastq.gz

#fastp --disable_adapter_trimming -i /proj/zylkalab/Esther/splitseq/CC_JN_rawdata/00001_CAGATC_S1_ME_L001_R1_001.fastq.gz -I /pine/scr/j/i/jieunp/alevin-fry/CC/S1_R2_barcodecollapsed.fastq.gz -o /pine/scr/j/i/jieunp/alevin-fry/CC/out.S1_R1.fastq.gz -O /pine/scr/j/i/jieunp/alevin-fry/CC/out.S1_R2_barcodecollapsed.fastq.gz

salmon alevin -i /proj/jmsimon/genomeAnnotation/gencode_vM25_transcriptome_splici_fl95_idx -l A -1 /pine/scr/j/i/jieunp/alevin-fry/CC/out.S1_R2_barcodecollapsed.fastq.gz -2 /pine/scr/j/i/jieunp/alevin-fry/CC/out.S1_R1.fastq.gz \
--read-geometry 2[2-100] \
--umi-geometry 1[1-10] \
--bc-geometry 1[11-18,49-56,87-94] \
-p 32 -o /pine/scr/j/i/jieunp/alevin-fry/CC/S1_run --rad --sketch

alevin-fry generate-permit-list -d both -i /pine/scr/j/i/jieunp/alevin-fry/CC/S1_run/ --output-dir /pine/scr/j/i/jieunp/alevin-fry/CC/S1_out_permit_knee --unfiltered-pl /proj/zylkalab/Esther/splitseq/SSSC081621/barcode_files/fullList_barcode_combination_zylkalab.txt --min-reads 22

alevin-fry collate -r /pine/scr/j/i/jieunp/alevin-fry/CC/S1_run/ -t 16 -i /pine/scr/j/i/jieunp/alevin-fry/CC/S1_out_permit_knee

alevin-fry quant -m /proj/jmsimon/genomeAnnotation/gencode_vM25_transcriptome_splici_fl95_t2g_3col.tsv -i /pine/scr/j/i/jieunp/alevin-fry/CC/S1_out_permit_knee -o /pine/scr/j/i/jieunp/alevin-fry/CC/S1_counts -t 16 -r cr-like-em --use-mtx

