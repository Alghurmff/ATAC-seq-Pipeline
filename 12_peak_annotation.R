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

library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)


gene <- seq2gene(peak_granges, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)


peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

## Use biomart to fill missing ensemble data
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

#test
x <- c(1, 1, 4, 5, 4, 6)
sum(is.na(x))
sum(duplicated(x))



atac_peaks<-read.csv("ATACseq_matrix_v2.txt", header=TRUE, sep=",")
sum(is.na(atac_peaks$SYMBOL)) # 549
sum(duplicated(atac_peaks$SYMBOL)) #33402

symbolnoduplicacetd<-atac_peaks$SYMBOL[!duplicated(atac_peaks$SYMBOL)]

library("biomaRt")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembleID <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                    filters= c("hgnc_symbol"), values= symbolnoduplicacetd, mart= mart)


write.table(ensembleID,file="ATACseq_ensemble.txt", col.names=TRUE, sep=";", row.names = FALSE)

ensemble<-read.csv("ATACseq_ensemble.txt", header=TRUE, sep=",") # 18225
sum(is.na(ensemble$ensembl_gene_id.hgnc_symbol)) 

atac_peaks$ensemble2 <- ensemble

##### useEnsembl isntead of useMart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembleID <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                    filters= c("hgnc_symbol"), values= symbolnoduplicacetd, mart= ensembl)

write.table(ensembleID,file="ATACseq_ensemble_v2.txt", col.names=TRUE, sep=";", row.names = FALSE)
ensemble_v2<-read.csv("ATACseq_ensemble_v2.txt", header=TRUE, sep=";")

## try to add the new ensemble to the whole matrix
## by excel use vlookup 
atac_peaks$ensemble2 <- NA
ensemble_v3 <- ensemble_v2[,-2]

atac_peaks$ensemble2 <- ensemble_v2$ensembl_gene_id

