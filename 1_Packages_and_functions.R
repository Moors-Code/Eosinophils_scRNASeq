######PACKAGES######
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


####FUNCTIONS#####
# import function .st data to sparse matrix
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

#volcano plot
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

####FGSEA####

#import dataset
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)

#rank genes and calculate BP enrichment
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


