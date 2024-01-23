# Multiplemyeloma_SE
RNA-seq analysis used RSEM (rsem-calculate-expression, v1.3.0) with STAR alignment (v2.5.7).
HiChIP sequencing data were processed with HiC-Pro pipeline version 2.11.1 and aligned to hg19 genome using Bowtie2 version 2.4.1. Sequencing reads were filtered and QC was performed for genome-wide signal correlation. Valid interactions were identified in each cell type-specific dataset using hichipper pipeline version 2.7.9 at 2 kb resolution by limiting the interactions between 5 kb to 2 Mb.
ChIP-seq peaks were called with MACS2 software version 2.2.7.1 with --keep-dup=1 --SPMR -p 1e-5 -f BAM –broad for H3K27ac broadPeak. The output bedGraph data were normalized by subtracting the corresponding background values using MACS2 bdgcmp -m subtract.

