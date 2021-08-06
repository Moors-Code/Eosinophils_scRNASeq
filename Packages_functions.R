#PACKAGES AND FUNCTIONS
library(devtools)
library(Seurat)
library(BiocManager)
library(biomaRt)
library(sctransform)
library(fgsea)
library(msigdbr)
library(SeuratData)
library(progeny)
library(coexnet)


library(dplyr)
library(magrittr)
library(data.table)
library(stringr)
library(Matrix)
library(tidyr)
library(readr)
library(tibble)
library(ggpubr)
library(proxyC)

library(ggplot2)
library(ggraph)
library(ggraph)
library(ggrepel)
library(cowplot)
library(pheatmap)
library("Nebulosa")
library(UpSetR)
library(grid)


library(RColorBrewer) 
library(viridis)
library(wesanderson)



# import function
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Gene, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Gene")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}


DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$Gene <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=avg_log2FC, y=-log10(p_val_adj), fill=factor(Significant), label = DEGs$Gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=4, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}



DEGs_plot <- function(DEGs, p_treshold, FC_treshold, title, color) {
  DEGs$Significant <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$Gene <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(y=avg_log2FC,x=1, fill=factor(Significant), label = DEGs$Gene) ) + 
    theme_classic()+ 
    geom_point(shape = 21,color="black")+
    scale_fill_manual(values = c("black", color)) +
    geom_text_repel(size=3, max.overlaps = 100)+
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}




preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=avg_log2FC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  BP_x <- fgsea(pathways = BP, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

preranked_REACTOME <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=avg_log2FC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  REACTOME_x <- fgsea(pathways = REACTOME, 
                      stats = ranks,
                      minSize=10,
                      maxSize=500,
                      nperm=1000000)
  
  REACTOME_x$pathway<-gsub("REACTOME_","",REACTOME_x$pathway)
  REACTOME_x$pathway<-gsub("_"," ",REACTOME_x$pathway)
  return(REACTOME_x)
}
preranked_Hallmark <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=avg_log2FC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  Hallmarks_x <- fgsea(pathways = Hallmarks, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500,
                       nperm=1000000)
  
  Hallmarks_x$pathway<-gsub("HALLMARK_","",Hallmarks_x$pathway)
  Hallmarks_x$pathway<-gsub("_"," ",Hallmarks_x$pathway)
  return(Hallmarks_x)
}

GSEA_dotplot <- function(BP_cluster, positions, title) {
  selected_BP_cluster <- BP_cluster[positions,c("pathway","NES","padj", "size") ]
  selected_BP_cluster <- arrange(selected_BP_cluster, NES) 
  selected_BP_cluster$pathway <- factor(selected_BP_cluster$pathway, levels = selected_BP_cluster$pathway)
  plot <- ggplot(data = selected_BP_cluster , aes(x = NES, y = pathway, size = size, label= pathway))+ 
    geom_point(aes(color=-log10(selected_BP_cluster$padj)), size=abs(selected_BP_cluster$NES)*5)  +theme_bw()+ 
    theme(axis.text.y = element_text( size=10, colour = "black"))+ labs(title = title, y = "", x="")+
    theme(axis.text.x = element_text(face="bold", size=10, colour = "black",hjust=1))+ xlim(-3.5,3.5)+
    scale_colour_gradientn(colours = pal[3:9]) + theme(legend.position="none")
  return(plot)
}



###IMPORT DATASET: Biological process ####
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)


####COLOR PALETTE FOR PLOTS ####
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector

col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]


