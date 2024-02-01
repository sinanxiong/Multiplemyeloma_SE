---
title: "1.ChIP-seq and Super Enhancer analysis"
author: "Tan Tze King"
date: "2024/1/1"
output: html_document
---


### 1 ChIP-seq alignment

```
bowtie2 --no-unal --very-sensitive -p 5 -x Bowtie_Index_hg19/hg19 -U $id.fastq.gz | samtools view -bS - | samtools sort -o $id.sort.bam -@ 5 -l 9 - 

```

### 2 Peak Calling

```
#broadPeak Calling
macs2 callpeak -t $id.sort.bam -c $id_input.sort.bam -f BAM --keep-dup 1 --outdir $id_MACS2 -n $id_broadPeak --SPMR --broad -p 1e-5 --broad-cutoff 1e-5 -B

#narrowPeak Calling
macs2 callpeak -t $id.sort.bam -c $id_input.sort.bam -f BAM --keep-dup 1 --outdir $id_MACS2 -n $id_broadPeak --SPMR -p 1e-5 -B --cutoff-analysis

#Signal track processing
macs2 bdgcmp -t $id_treat_pileup.bdg -c $id_control_lambda.bdg -o $id.bdg -m subtract
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $id.bdg | awk -F't' '$4>0.1' | sed 's/MT/M/g' | sort -k1,1 -k2,2n > $id_filter.bdg
bedGraphToBigWig $id_filter.bdg hg19/chrom.sizes $id_subtract_SPMR.bw
```

### 3 Super-enhancer calling

```
ROSE2 -i $id.broadPeak.bed -r $id.sort.bam -c $id_input.sort.bam -g HG19 -o ROSE2/$id_ROSE2 -s 12500 -t 2000 --mask hg19/blacklist.bed

#enhancer annotation
annotatePeaks.pl $id_AllEnhancer.bed hg19 > $id_AllEnhancer.annotate

```

### 4 DiffBind differential analysis

```
library(DiffBind)

sample<- read.csv("diffbind.csv", sep='\t')
dataset <- dba(sampleSheet=sample )
dataset.count <- dba.count(dataset, bUseSummarizeOverlaps=TRUE,bParallel=TRUE,summits=FALSE,filter=FALSE)
dbObj <- dba.blacklist(dataset.count, blacklist=DBA_BLACKLIST_HG19, greylist=FALSE)
dbObj <- dba.contrast(dbObj, categories=c( DBA_FACTOR) , minMembers = 2, contrast=c("Factor","MM","control"))
dba.show(dbObj, bContrasts=T)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(res_deseq)
write.table(out, file="diffBind_cellLine_MM_vs_Control_H3K27ac_NoSummitPeak.txt", sep="\t", quote=F, row.names=F)

```

