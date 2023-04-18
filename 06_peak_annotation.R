#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 12_peak_annotation.R      #
#***********************************#
### http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# Library required 
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read consensus peaks
peaks <- read.table("consensus_cell_types.bed") 
colnames(peaks) <- c("chr","start","end","","")
peaksID <- paste(peaks[,1],peaks[,2], peaks[,3], sep="_")

#compute the mean coverage
count_table <- read.table("count_table.txt", sep="\t", header=TRUE)

# Remove non-numeric columns 
#and intervals numbers that determine the start and the end of the peak
mean_cov<-rowMeans(count_table[,-c(1:6)])

#mean_cov <- rowMeans(count_table)
peaks$coverage <- mean_cov

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks[,-c(4:6)], keep.extra.columns=TRUE)
head(peak_granges)
str(peak_granges)

##Coverage plot
png("peaks_over_chr.png")
covplot(peak_granges, title = "ATAC-Seq Peaks over Chromosomes", weightCol="coverage")
dev.off()

# Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak_granges, windows=promoter)

# Heatmap of ChIP binding to TSS regions
png("binding_TSS.png")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()

# Average Profile of ChIP peaks binding to TSS region
png("average_binding_TSS.png")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

# Peak Annotation

library(org.Hs.eg.db)

peakAnno <- annotatePeak(peak_granges, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db", level="gene")
head(peakAnno)
peak_annotation <- as.data.frame(peakAnno) 

peak_annotation$annotation2 <- peak_annotation$annotation
peak_annotation$annotation2[grep("Promoter",peak_annotation$annotation)] <- "Promoter"
peak_annotation$annotation2[grep("Intron",peak_annotation$annotation)] <- "Intron"
peak_annotation$annotation2[grep("Distal Intergenic",peak_annotation$annotation)] <- "Distal Intergenic"
peak_annotation$annotation2[grep("Exon",peak_annotation$annotation)] <- "Exon"


write.csv(peak_annotation, "consensus_peaks_annotation_chipseeker.txt")

# Visualize Genomic Annotation
png("pie_annotation.png")
plotAnnoPie(peakAnno)
dev.off()

