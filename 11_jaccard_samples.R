#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 11_ccard_heatmap.sh       #
# script: 11_jaccard_heatmap.R      #
#***********************************#

# 11. Quality Control of Consensus Peaks - Measuring dataset similarity


###########################################
# PLOT JACCARD SCORE - SIMILARITY MEASURE #
###########################################
##https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp6--measuring-dataset-similarity

jaccard2 <- read.table("pairwise_jaccard.txt", header=TRUE, sep=" ", row.names=1)

jaccard <- jaccard2 [-1,-1] # remove  Allconsensus_cell column and row 

colnames(jaccard) <- c("CLP", "HSC", "Immature.B", "Nutrophils",
                       "ProB-Progenitors","PreB", "ProB","PC")
colnames(jaccard) <- gsub(".bed","",gsub("Allconsensus_cell_types","",colnames(jaccard)))
rownames(jaccard) <- gsub(".bed","",gsub("Allconsensus_cell_types","",colnames(jaccard)))


jaccard[,ncol(jaccard)] <- as.numeric(gsub("[\\n]","",jaccard[,ncol(jaccard)]))

jaccard <- as.matrix(jaccard)


# List of colors
# https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
library(RColorBrewer)

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:ncol(jaccard)]
names(color_types) <- colnames(jaccard)

library(ComplexHeatmap)
ha = HeatmapAnnotation(Cell_type=colnames(jaccard),
                       col = list(Cell_type = color_types),  
                       gp = gpar(col = "black"),  #black border per cell
                       annotation_legend_param = list(Cell_type = list(border="black"))
)

ra = rowAnnotation(Cell_type=colnames(jaccard),
                   col = list(Cell_type = color_types),  
                   gp = gpar(col = "black"),  #black border per cell
                   show_legend=FALSE
)


# Heatmap colors
library(circlize)
library(grid)
col_heatmap = colorRamp2(seq(0, 1, length = 3), c("white", "gray", "red"), space = "RGB")

nr <- nrow(jaccard)
png("jaccard_consensus2.png", width = 0.8*nr + 1.3150, height = 0.8*nr + 1.3150, units = "in", res = 100)
draw(Heatmap(jaccard, top_annotation = ha, show_row_names = TRUE, left_annotation = ra, col=col_heatmap, column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10), border=TRUE, heatmap_legend_param = list(title = "jaccard", border="black"), height = unit(6, "mm")*nr))
dev.off()


# MDS
library(ggplot2)

mds <- cmdscale(1-jaccard,eig=TRUE, k=2)

mplot <- data.frame(MDS1=mds$points[,1], MDS2=mds$points[,2], Cell_type=as.factor(names(color_types)))

p <- ggplot(mplot, aes(x=MDS1, y=MDS2, color = Cell_type, label=rownames(jaccard))) + 
  geom_point(size=2) +
  scale_color_manual(values=color_types) + 
  labs(title="MDS - Measure of similarity (Jaccard)", x ="Coordinate 1", y="Coordinate 2")+
  geom_text() +
  geom_label() +
  xlim(-0.5, 0.5)

nr=150
ggsave(plot = p, width = 0.06*nr, height = 0.04*nr, dpi = 100, filename = "jaccard_mds_tissue2.pdf")

