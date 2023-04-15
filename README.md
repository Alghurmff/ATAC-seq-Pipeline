# ATAC-seq Pipeline

## Description 

This ATAC-seq pipeline is both straightforward and reproducible, covering everything from raw FASTQ files to transcription factor activity and generating a signal track file (bigwig). It can handle both single-end and paired-end data. Detailed error reporting is provided, and the pipeline allows for easy resumption of interrupted runs. This pipeline is run using Linux and R.

The figure below provides an overview of all the steps involved, along with the corresponding tools and programming language.

<img width="819" alt="ATACseq_pipeline" src="https://user-images.githubusercontent.com/114663482/232249040-fdc83234-737b-4942-8b00-508094ca3c0c.png">

## List of tools : 

### Bash Script

* fastqc/0.11.8
* trimmomatic/0.38

* Alignment 
bowtie2/2.3.5

* Convert sam to bam , sorting and indexing , mark duplicate QC and filtering 
samtools/1.8

* Mark duplicate 
picard/2.20.4

* Peak Calling
macs2/2.1.1.2/anaconda2-2.5.0

* Frip scores - manuel equation
bedtools/2.29.0

* Visualization (bigwig)
 deeptools/python2.7/3.3.1

* MultiQC
multiqc/1.9


### R Packages
R version 4.2.2 (2022-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.4

* Consensus peak  
rtracklayer_1.58.0 
GenomicRanges_1.50.2  
* Consensus peak similarity between cell types (Heatmap and MDS)
Jaccard

* Count table 
Multicov function from bedtools 2.29.0

* Peak annotation 
ChIPseeker v1.34.1
TxDb.Hsapiens.UCSC.hg38.knownGene
GenomicFeatures
 AnnotationDbi
GenomicRanges_1.50.2  
rtracklayer_1.58.0 
clusterProfiler v4.6.0

* Transcription factor activity 
chromVAR_1.20.2 
Motifmatchr_1.20.0
Matrix_1.5-3
SummarizedExperiment_1.28.0
BiocParallel_1.32.5
BSgenome.Hsapiens.UCSC.hg38_1.4.5 

## Quick Start

Open runATAC.sh in a text editor and make the necessary changes by editing the path according to your file, updating the tools and files or variables name. After that go to templateATAC.sh and make sure that the changes you made in runATAC.sh is updated in the template. 

### To run the script from FASTQC to peak calling for each sample do the followin:
```
bash runATAC.sh "sample_name"
```
After that all the files and reports will be generated for ech samples.

 
### Consensus peaks
First merge all narrowpeak files to use it for consensus peak and then run the  10_consensus_peaks.Rmd file.

For consensus peaks QC, run 11_Jaccard_heatmap.sh. This will perform a loop between all bed files from different cell types in pairwise_jaccard.txt file. This file will be used in 11_jaccard_samples.R  to creat the jacccard heat map and mds similarity.

### Count table
Before doing the peak annotation we need to do count table to get the peaks that overlap between clean.bam files from all samples bed file for all cell type using bash script (bedCount.sh) with the function multicov. 


### Peak annotation 

(12_peak_annotation.R) using the consensus peak for all cell types file

### Transcription factor activity

 (13_chromVAR_TF .R) which has a detailed explination of the process



