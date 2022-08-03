library(devtools)
library(Seurat)
library(BiocManager)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
library(ggraph)
library(wesanderson)
library(fgsea)
library(msigdbr)
library(ggraph)
library(ggrepel)
library(cowplot)
library(magrittr)
library(data.table)
library(stringr)
library(Matrix)
library(dplyr)
library(tidyr)
library("edgeR")


#functions
DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$Gene <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=logFC, y=-log10(FDR), color=factor(Significant), label = DEGs$Gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 19)+
    scale_color_manual(values = c("grey", color)) +
    geom_text_repel(size=4, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}


#####IMPORT DATA FROM GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131032 
# from the publication Czarnewski P, Parigi SM, Sorini C, Diaz OE et al. Conserved transcriptomic profile between mouse and human colitis 
# allows unsupervised patient stratification. Nat Commun 2019 Jun 28;10(1):2892.  PMID: 31253778 

all_samples_counts <- read.csv("~/NAS/Coco/Collaborations/Eosinophils BD/Data analysis/Final/DSS_Villablanca/GSE131032_kallisto_counts.csv")
View(all_samples_counts)
dim(all_samples_counts)
all_samples_mat <- as.matrix(all_samples_counts[2:27])
rownames(all_samples_mat)<- all_samples_counts$X
colnames(all_samples_mat)  <- c("d0_1", "d0_2", "d0_3", "d2_1", "d2_2", "d2_3", "d4_1", "d4_2", "d4_3",
                               "d6_1", "d6_2", "d6_3", "d7_1", "d7_2","d7_3", "d8_1", "d8_2", "d8_3",
                               "d10_1", "d10_2", "d12_1", "d12_2","d12_3", "d14_1", "d14_2", "d14_3")
head(all_samples_mat)

#Create DGEList object
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9))
y <- DGEList(counts=all_samples_mat,group=group)
y$samples
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes")

# Filter reads by counts: Most of the samples should have at least 10 reads, normalize the library and estimate dispersion
keep <- filterByExpr(y, min.count = 10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
y <- estimateDisp(y,design)
plotMDS(y, pch = 2, label=colnames(y)) 
logcpm <- cpm(y, log=TRUE)
head(logcpm)

#extract FPKM value for marker gene heatmap
markergenes <- c("Siglecf", "Ccr3", "Il5ra", "Cd80", "Cd274", "Il33", "Ifng",
                  "Epx", "Ear1", "Ear2", "Ear6", "Prg2","Prg3",  "Prtn3","Elane", "Tuba1b", "Cebpe", 
                 "S100a6","S100a8", "S100a9", "Gbp2", "Gbp7", "Ltf", 
                 "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7")


all_samples_df <- as.data.frame(logcpm)
all_samples_df$Gene <- rownames(all_samples_df)
markergenes_values <- filter(all_samples_df, Gene %in% markergenes)
#sort 
markergenes_values_sorted <- markergenes_values[match(markergenes, markergenes_values$Gene),]

#make numeric
markergenes_values_sorted_num            <-  data.matrix(markergenes_values_sorted[,1:26 ], rownames.force = NA)

#prepare palette for pheatmap
paletteLength   <- 50
myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
breaksList      = seq(-2, 2, by = 0.04)

##prepare annotation for pheatmap
#annotation_rows             <- data.frame(markers = rep(c("Precursors", "Immature", "Circulating", "Basal", "Active"), c(10, 6, 4, 6, 10)))
#rownames(annotation_rows)   <- rownames(markergenes_values_sorted_num)
#annotation_rows$markers     <- factor(annotation_rows$markers, levels = c("Precursors", "Immature", "Circulating", "Basal", "Active"))


library(pheatmap)
library(grid)
pdf("heatmap.pdf",  width = 5, height = 10)
p <- pheatmap(markergenes_values_sorted_num, scale="row", 
         color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         cluster_rows = F, cluster_cols = F, 
         border_color = "black", 
         legend_breaks = -2:2, 
         cellwidth = 15, cellheight = 10,
         angle_col = "45", 
         #annotation_row = annotation_rows,
         fontsize = 10)
pdf(file="heatmap.pdf")
grid.draw(p)
dev.off()

