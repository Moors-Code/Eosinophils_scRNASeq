#####PACKAGES and FUNCTIONS####
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
  DEGs$Gene <- ifelse((DEGs$FDR < p_treshold/2 & abs(DEGs$logFC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=logFC, y=-log10(FDR), color=factor(Significant), label = DEGs$Gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 20)+
    scale_color_manual(values = c("grey", color)) +
    geom_text_repel(size=2, max.overlaps = Inf, segment.size=0.3)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=10),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}

#####EGDER####
##import and format count table from GEO
result_Il33_over_ctrl <- read.csv("~/NAS/Coco/Collaborations/Eosinophils BD/bulk/result--Il33--over--ctrl.csv")
all_samples_counts <- as.data.frame(result_Il33_over_ctrl[, c(3,23:38)])
head(all_samples_counts)
all_samples_mat <- as.matrix(all_samples_counts[,c(2:17)])
head(all_samples_mat)
rownames(all_samples_mat)<- all_samples_counts$gene_name
colnames(all_samples_mat) <- c("ctrl_1", "ctrl_2", "ctrl_3", "ctrl_4", "IL33_1", "IL33_2", "IL33_3", "IL33_4",
                               "IFNy_1", "IFNy_2", "IFNy_3", "IFNy_4", "IL33+IFNy_1","IL33+IFNy_2", "IL33+IFNy_3", "IL33+IFNy_4")
head(all_samples_mat)

#Create DGEList object
group <- factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4))
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

pdf("MDSplot.pdf",  width = 5, height = 4)
pch <- c(19, 19, 19, 19)
colors <- rep(c("darkgreen", "red", "blue", "orange"), 2)
plotMDS(y, col=colors[group], pch=pch[group], ylim=c(-1, 1.5))
legend("topright", legend=c("ctrl", "IL33", "IFNy", "IL33+IFNy"), pch=pch, col=colors, ncol=2)
dev.off() 



#####DEGS and VOLCANOS#####
# IL33_vs_ctrl
IL33_vs_ctrl <- exactTest(y, pair=c(1,2)) #the first group listed in the pair is the baseline for the comparison
topTags(IL33_vs_ctrl)
summary(decideTests(IL33_vs_ctrl))

all_IL33_vs_ctrl = topTags(IL33_vs_ctrl, n = Inf)
dim(all_IL33_vs_ctrl)
head(all_IL33_vs_ctrl$table)
all_IL33_vs_ctrl$table$Gene <- rownames(all_IL33_vs_ctrl$table)
View(all_IL33_vs_ctrl$table)

write.table(all_IL33_vs_ctrl$table, file = "IL33_vs_ctrl.txt")

DEGs<-all_IL33_vs_ctrl$table
pdf("IL33volcano.pdf",  width = 5, height = 4)
DEGs_volcano(DEGs, 0.05, 2.75, "IL33 vs ctrl", "red", 100, 7)
dev.off() 

#IFNy_vs_ctrl
IFNy_vs_ctrl <- exactTest(y, pair=c(1,3)) #the first group listed in the pair is the baseline for the comparison
topTags(IFNy_vs_ctrl)
summary(decideTests(IFNy_vs_ctrl))

all_IFNy_vs_ctrl = topTags(IFNy_vs_ctrl, n = Inf)
dim(all_IFNy_vs_ctrl)
head(all_IFNy_vs_ctrl$table)
all_IFNy_vs_ctrl$table$Gene <- rownames(all_IFNy_vs_ctrl$table)
View(all_IFNy_vs_ctrl$table)
write.table(all_IFNy_vs_ctrl$table, file = "IFNy_vs_ctrl.txt")

DEGs<-all_IFNy_vs_ctrl$table
pdf("IFNyvolcano.pdf",  width = 5, height = 4)
DEGs_volcano(DEGs, 0.05, 2.75, "IFNy vs ctrl", "blue", 100, 10)
dev.off() 

# IL33IFNy_vs_IL33
IL33IFNy_vs_IL33 <- exactTest(y, pair=c(2,4)) #the first group listed in the pair is the baseline for the comparison
topTags(IL33IFNy_vs_IL33)
summary(decideTests(IL33IFNy_vs_IL33))

all_IL33IFNy_vs_IL33 = topTags(IL33IFNy_vs_IL33, n = Inf)
dim(all_IL33IFNy_vs_IL33)
head(all_IL33IFNy_vs_IL33$table)
all_IL33IFNy_vs_IL33$table$Gene <- rownames(all_IL33IFNy_vs_IL33$table)
View(all_IL33IFNy_vs_IL33$table)
write.table(all_IL33IFNy_vs_IL33$table, file = "all_IL33IFNy_vs_IL33.txt")

DEGs<-all_IL33IFNy_vs_IL33$table
pdf("IL33IFNyvolcano.pdf",  width = 5, height = 4)
DEGs_volcano(DEGs, 0.05, 3, "IL33+IFNy vs IL33", "orange", 100, 10)
dev.off() 

#IL33+IFNy vs ctrls
IL33IFNy_vs_ctrl <- exactTest(y, pair=c(1,4)) #the first group listed in the pair is the baseline for the comparison
topTags(IL33IFNy_vs_ctrl)
summary(decideTests(IL33IFNy_vs_ctrl))

all_IL33IFNy_vs_ctrl = topTags(IL33IFNy_vs_ctrl, n = Inf)
dim(all_IL33IFNy_vs_ctrl)
head(all_IL33IFNy_vs_ctrl$table)
all_IL33IFNy_vs_ctrl$table$Gene <- rownames(all_IL33IFNy_vs_ctrl$table)
View(all_IL33IFNy_vs_ctrl$table)
write.table(all_IL33IFNy_vs_ctrl$table, file = "all_IL33IFNy_vs_ctrl.txt")

DEGs<-all_IL33IFNy_vs_ctrl$table
pdf("IL33IFNy_vs_ctrl_volcano.pdf",  width = 5, height = 4)
DEGs_volcano(DEGs, 0.05, 2.75, "IL33+IFNy vs ctrl", "orange", 100, 10)
dev.off() 

#####VENN DIAGRAM####
# Load library
library(VennDiagram)
# Generate 3 sets of 200 words
positionsIL33 <- which(all_IL33_vs_ctrl$table$FDR < 0.05 & (all_IL33_vs_ctrl$table$logFC < 2 | all_IL33_vs_ctrl$table$logFC > 2))
positionsIFNy <- which(all_IFNy_vs_ctrl$table$FDR < 0.05 & (all_IFNy_vs_ctrl$table$logFC < 2 | all_IFNy_vs_ctrl$table$logFC > 2))
positionsIL33IFNy <- which(all_IL33IFNy_vs_ctrl$table$FDR < 0.05 & (all_IL33IFNy_vs_ctrl$table$logFC < 2 | all_IL33IFNy_vs_ctrl$table$logFC > 2))

IL33 <- all_IL33_vs_ctrl$table$Gene[positionsIL33]
IFNy <- all_IFNy_vs_ctrl$table$Gene[positionsIFNy]
IL33IFNy <- all_IL33IFNy_vs_ctrl$table$Gene[positionsIL33IFNy]

# Chart
pdf("venn.pdf")
temp <- venn.diagram(x = list(IL33, IFNy, IL33IFNy),
                     category.names = c("IL33" , "IFNy" , "IL33+IFNy"), 
                     lwd = 0,
                     lty = 'blank',
                     fill = c("red", "blue", "orange"),
                     filename = NULL, 
                     cex = 1,
                     area.vector = c(length(IL33)/2, length(IFNy)/2, length(IL33IFNy)/2),
                     fontface = "bold",
                     fontfamily = "sans",
                     cat.cex = 1.5,
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans")
grid.draw(temp)
pdf(file="venn.pdf")
grid.draw(temp)
dev.off()

venn.diagram(x = list(IL33, IFNy, IL33IFNy),
             category.names = c("IL33" , "IFNy" , "IL33+IFNy"), filename="venn.png")

####UPSET PLOT####
library(UpSetR)
sigDEGs   <- list("IL33"=IL33 , 
                  "IFNy" =IFNy,
                  "IL33+IFNy"= IL33IFNy)
pdf("upset.pdf")
upset(fromList(sigDEGs), 
      order.by = "freq", 
      nintersects = 10,  
      main.bar.color = "darkblue", 
      keep.order = T, 
      matrix.color = "darkred", 
      mainbar.y.label = "intersecting DEGs", 
      sets.x.label= "DEGs",
      sets.bar.color = "lightgrey", 
      point.size = 7, 
      empty.intersections=F, 
      text.scale = 2)
dev.off()

####HEATMAPS####
#extract FPKM value for marker gene heatmap
markergenes <- c("Mki67", "Tuba1b", "Epx", "Prg3", "Prg2","Ear1","Ear2", "Ear6",  "Cd63", "Cebpe",
                 "Alox15", "Aldh2", "S100a9", "S100a6", "S100a10", "Retnla", "Ccl9", "Il1rl1", 
                 "Cd24a", "Mmp9", "Icosl", "Il4", "Tgfb1", "Pirb", "Rara", "Cd80", "Cd274", "Ptgs2", "Il1rn", "Il1b", 
                 "Vegfa", "Ccl3", "Cxcl2", "Il16", "Tnf")

all_samples_df <- as.data.frame(logcpm)
all_samples_df$Gene <- rownames(all_samples_df)
markergenes_values <- filter(all_samples_df, Gene %in% markergenes)
#sort 
markergenes_values_sorted <- markergenes_values[match(markergenes, markergenes_values$Gene),]

#make numeric
markergenes_values_sorted_num            <-  data.matrix(markergenes_values_sorted[, c(1:16)], rownames.force = NA)

#prepare palette for pheatmap
paletteLength   <- 50
myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
breaksList      = seq(-2, 2, by = 0.04)

#prepare annotation for pheatmap
annotation_rows             <- data.frame(markers = rep(c("Precursors", "Immature", "Circulating", "Basal", "Active"), c(10, 5, 4, 6, 10)))
rownames(annotation_rows)   <- rownames(markergenes_values_sorted_num)
annotation_rows$markers     <- factor(annotation_rows$markers, levels = c("Precursors", "Immature", "Circulating", "Basal", "Active"))

mycolors <- c("#5BBCD6",  "#F98400", "#F2AD00","#00A08A", "#FF0000")
names(mycolors) <- unique(annotation_rows$markers)
mycolors <- list(category = mycolors)
annot_colors=list(Subset=c(Precursors="#5BBCD6", Immature= "#F98400", Circulating= "#F2AD00", Basal="#00A08A",Active= "#FF0000"))


library(pheatmap)
library(grid)
pdf("heatmap.pdf",  width = 5, height = 10)
p <- pheatmap(markergenes_values_sorted_num, scale="row", 
              color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
              breaks = breaksList,
              cluster_rows = F, cluster_cols = T, 
              border_color = "black", 
              legend_breaks = -2:2, 
              cellwidth = 15, cellheight = 10,
              angle_col = "45", 
              annotation_colors = annot_colors,
              annotation_row = annotation_rows,
              fontsize = 10)
pdf(file="heatmap.pdf")
grid.draw(p)
dev.off()

##with subset markers.. weird > reorder
markers_steadystate <- read.csv("~/NAS/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Steadystate/markers_steadystate.csv")
top50 <- markers_steadystate %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

all_samples_df <- as.data.frame(all_samples_mat)
all_samples_df$Gene <- rownames(all_samples_df)
marker50genes_values <- filter(all_samples_df, Gene %in% top50$gene)
#sort 
marker50genes_values_sorted <- marker50genes_values[match(top50$gene, marker50genes_values$Gene),]

#make numeric
marker50genes_values_sorted_num            <-  data.matrix(marker50genes_values_sorted[, c(1:16)], rownames.force = NA)

#prepare palette for pheatmap
paletteLength   <- 50
myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
breaksList      = seq(-2, 2, by = 0.04)

pdf("heatmap_markers.pdf",  width = 5, height = 10)
p <- pheatmap(marker50genes_values_sorted_num, scale="row", 
              color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
              breaks = breaksList,
              cluster_rows = F, cluster_cols = T, 
              border_color = NA, 
              legend_breaks = -2:2, 
              cellwidth = 15, 
              cellheight =3,
              angle_col = "45", 
              fontsize = 6)
pdf(file="heatmap_markers.pdf")
grid.draw(p)
dev.off()


###SIGNATURES
Antigen_processing <- c("H2.D1", "H2.Q7", "H2.K1", "H2.T23","H2.Q4","Tap1","Tapbp","B2m", "Psmb8", "Psme1",
                             "Psmb9", "Calr", "Psmb10", "Ncf1", "Fcer1g")

granules_antimicrobial <- c( "Epx", "Ear1", "Ear2", "Ear6",  "Prg2","Prg3",  "Prtn3","Elane", "Tuba1b", "Cebpe", 
                             "Ltf",  "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "S100a6", "S100a8", 
                             "S100a9", "Gbp2", "Gbp7", "Cd63")

Nfkb_list <- c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz",  "Rela", "Relb")

Regulatory_list <-  c("Cd274", "Cd80", "Cd9")

all_signatures <- c(Nfkb_list, Regulatory_list,granules_antimicrobial, Antigen_processing)

markergenes <- all_signatures
all_samples_df <- as.data.frame(logcpm)
all_samples_df$Gene <- rownames(all_samples_df)
markergenes_values <- filter(all_samples_df, Gene %in% markergenes)
#sort 
markergenes_values_sorted <- markergenes_values[match(markergenes, markergenes_values$Gene),]

#make numeric
markergenes_values_sorted_num  <-  data.matrix(markergenes_values_sorted[, c(1:16)], rownames.force = NA)

#prepare palette for pheatmap
paletteLength   <- 50
myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
breaksList      = seq(-2, 2, by = 0.04)

#prepare annotation for pheatmap
annotation_rows             <- data.frame(signatures = rep(c("Nfkb signalling", "Immune regulation", "Granules and antimicrobial peptides", "Antigen processing and presentation"), 
                                                        c(length(Nfkb_list), length(Regulatory_list), length(granules_antimicrobial), length(Antigen_processing))))
rownames(annotation_rows)   <- rownames(markergenes_values_sorted_num)
annotation_rows$signatures     <- factor(annotation_rows$signatures, levels = c("Nfkb signalling", "Immune regulation", "Granules and antimicrobial peptides", "Antigen processing and presentation"))

library(pheatmap)
library(grid)
pdf("heatmap_signatures.pdf",  width = 5, height = 10)
p <- pheatmap(markergenes_values_sorted_num, scale="row", 
              color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
              breaks = breaksList,
              cluster_rows = F, cluster_cols = T, 
              border_color = "black", 
              legend_breaks = -2:2, 
              cellwidth = 12, cellheight =6,
              angle_col = "45", 
              annotation_row = annotation_rows,
              fontsize = 6)
pdf(file="heatmap_signatures.pdf")
grid.draw(p)
dev.off()



