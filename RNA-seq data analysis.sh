#!/bin/bash
# ---
# title: "RNA-seq bash script"
# author: "Tuan Zea Tan"
# last update: "Jan 23, 2024"
# ---

# Alignment & transcript quantification
SAMPLE=$FILENAME
GNOMEVER='44'
THREAD=8

./STAR/2.7.5b/bin/Linux_x86_64/STAR --runThreadN ${THREAD} --genomeDir star.hg38.Gv${GNOMEVER} --readFilesCommand zcat --readFilesIn ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz --outFileNamePrefix ${SAMPLE} --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM

./rsem/1.3.0/bin/rsem-calculate-expression -p ${THREAD} --bam --paired-end ${SAMPLE}Aligned.toTranscriptome.out.bam hg38.Gv${GNOMEVER}.rsem ${SAMPLE}


# generate a data matrix for downstream analyses
./rsem/1.3.0/bin/rsem-generate-data-matrix *.genes.results > ${SAMPLE}_count.txt

#differential expression
Rscript gEBseq.R ${SAMPLE}_count.txt

