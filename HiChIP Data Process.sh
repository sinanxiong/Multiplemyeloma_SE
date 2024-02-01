---
title: "2. HiChIP-seq analysis"
author: "Tan Tze King"
date: "2024/1/1"
output: html_document
---


### 1 HiChIP alignment and Pair Process

```
for i in $id ; do  
bwa mem -5SP -T 0 -t 20 bwa_index/hg19.fa $id_1.fastq.gz  $id_2.fastq.gz 
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 20 --nproc-out 20 --chroms-path hg19/chrom.sizes 
pairtools sort --tmpdir . --nproc 20
pairtools dedup --nproc-in 20 --nproc-out 20 --mark-dups --output-stats $id_stats.txt 
pairtools split --nproc-in 20 --nproc-out 20 --output-pairs $id_mapped.pairs --output-sam -
samtools view -bS -@ 20 
samtools sort -@ 20 -o $id_mapped.PT.bam
samtools index $id_mapped.PT.bam 
done

```



### 2 Calling peaks with MACS2 on HiChIP data

```
samtools –view –h –F 0x900 $id_mapped.PT.bam  | bedtools bamtobed -i stdin > $id.primary.aln.bed
macs2 callpeak -t $id_primary.aln.bed -n $id_macs2_p05 --keep-dup 1 --outdir MACS2 --broad --SPMR -B -p 1e-5 --broad-cutoff 1e-5

```


### 3 HiChipper Interaction Calling

```
grep -v '#' $id_mapped.pairs | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > $id_allValidPairs
hichipper call -o $id_hichipper -ii $id_allValidPairs -l 150 -p $id_macs2_broadPeak.bed --make-ucsc --make-washu --skip-resfrag-pad
awk ' {{printf $1 "\t" $2 "\t" $6 "\t"$7 "\t"}  if($8*10>1000){printf "1000"} else {printf $8*10} {print  "\t" $8 "\t\.\t0\t" $1 "\t" $2 "\t" $3 "\t\.\t\.\t" $4 "\t" $5 "\t" $6 "\t\.\t\." }}' $id_hichipper_filt.bedpe > $id_tmp
sort -k1,1 -k2,2n $id_tmp > $id_tmp2
bedToBigBed -as=hg19/interaction.as -type=bed5+13 $id_tmp2 hg19/chrom.sizes $id_hichipper_filt.bb 

```

### 4 Virtual 4C Analysis

```
sambamba sort -o $id_NameSrt.bam -n -M -t 30 $id_mapped.PT.bam
pairToBed -abam $id_NameSrt.bam -b target.bed > $id_v4C.bam
sambamba sort -o $id_v4C.bam -M -t 30 $id_v4C.sort.bam
multiBamSummary bins --bamfiles $id_v4C.sort.bam --smartLabels --scalingFactors scale_factor.txt -bs 500 -r chr1:204338000:204498000 --outRawCounts raw_count.txt 
```


### 5 hic file convert

```
java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre --threads 10 $id_mapped.pairs $id_H3K27ac.hic hg19.genome 
```