#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 13_chromVAR_TF.R          #
#***********************************#



# https://greenleaflab.github.io/chromVAR/articles/Introduction.html
# https://bioconductor.org/packages/3.3/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment


# loading required packages
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)


# first we need to construct SummarizedExperiment object 

 ## Getting the required files
# 1- Set of OCR
peaks_file="consensus_cell_types.bed"
# Read consensus peaks
peaks <- read.table(peaks_file)
rownames(peaks) <- paste(peaks[,1],peaks[,2],peaks[,3],sep="_")
colnames(peaks) <- c("chr","start","end")

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks)
peak_granges

# 2- OCR counts and metadata
counts_file <- "count_table.txt"
my_counts_matrix <- read.delim(counts_file,
    sep = "\t",
    dec = ",",
    check.names = FALSE, 
    header = TRUE)

my_counts_matrix <- my_counts_matrix [,-c(1:6)]
my_counts_matrix <- as.matrix(my_counts_matrix)


# metadata
library(readxl)

atac_metadata <- read_excel("BMAtlas_ATACseq_metadata.xlsx")
rownames(atac_metadata) <- atac_metadata$SampleID
atac_metadata2 <- atac_metadata[colnames(my_counts_matrix),]

#order cell type labels
atac_metadata2$cell_type <- factor(as.character(atac_metadata2$cell_type),
                                   levels=c("HSC","CLP","ProB","PreB","Immature.B","neutrophils","PreProB-progenitors"))
table(atac_metadata2$cell_type)
sum(table(atac_metadata2$cell_type)) # 55 
## Creating the object 

fragment_counts <- SummarizedExperiment(assays=list(counts=my_counts_matrix),
                     rowRanges=peak_granges, colData=atac_metadata2)

# Getting GC content of peaks
fragment_counts2 <- addGCBias(fragment_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg38)

#fragment_counts2
#head(rowData(fragment_counts2))

# Filtering inputs
counts_filtered <- filterPeaks(fragment_counts2, non_overlapping = TRUE)

#Get motifs and what peaks contain motifs
library(JASPAR2022)
library(TFBSTools)

##
getJasparMotifs2022 <- function(species = "Homo sapiens",
                                collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}
motifs2022 <- getJasparMotifs2022()


motif_ix <- matchMotifs(motifs2022, fragment_counts2, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

# Compute deviations
#dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)

#Background Peaks
bg <- getBackgroundPeaks(object = fragment_counts2)

dev <- computeDeviations(object = fragment_counts2, annotations = motif_ix,
                         background_peaks = bg)

# Variability
variability <- computeVariability(dev)

plotVariability(variability, use_plotly = FALSE) 

#plotVariability(variability) # interactive plot

#Visualizing Deviations
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)
#plot(tsne_results[,1], tsne_results[,2])

#atac_metadata[atac_metadata$SampleID[rownames(tsne_results)],]
#Error in Rtsne.default(t(mat[ix2, , drop = FALSE]), perplexity = perplexity,  : 
#Remove duplicates before running TSNE.

tsne_plots <- plotDeviationsTsne(dev, tsne_results,
                                 sample_column = "cell_type", 
                                 shiny = FALSE)
tsne_plots[[1]]
#plot(tsne_results[,1], tsne_results[,2])
#tsne_plots[[2]]

write.csv(deviationScores(dev),file = "ATACseq_TF-zScores")
