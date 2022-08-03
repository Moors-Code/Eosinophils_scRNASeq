if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea", force=T)

library(fgsea)
library(msigdbr)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)

####MAGECK COUNT in terminal###
mageck count -l GWIlibrary.csv -n GWI1_active_BMSC --sample-label "active,BMSC" --fastq F5_S5_R1_001.bam F3_S3_R1_001.bam --norm-method total
mageck count -l GWIlibrary.csv -n GWI2_active_BMSC --sample-label "active,BMSC" --fastq F6_S6_R1_001.bam F4_S4_R1_001.bam --norm-method total

#import in R
GWI1_active_BMSC.count_normalized.txt <- read.delim("/media/Coco/Collaborations/Eosinophils BD/GWI screen/fastqs/bowtie_output/GWI1_active_BMSC.count_normalized.txt")
GWI2_active_BMSC.count_normalized.txt <- read.delim("/media/Coco/Collaborations/Eosinophils BD/GWI screen/fastqs/bowtie_output/GWI2_active_BMSC.count_normalized.txt")
head(GWI1_active_BMSC.count_normalized.txt)
paired_counts <- merge(GWI1_active_BMSC.count_normalized.txt, GWI2_active_BMSC.count_normalized.txt, by="sgRNA")
paired_counts <- paired_counts[,c(1,2,3,4,6,7)]
colnames(paired_counts) <- c( "sgRNA" ,"Gene", "active1","BMSC1","active2","BMSC2")
head(paired_counts)
write.table(paired_counts, "fastqs/bowtie_output/paired_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

###MAGECK TEST in terminal ###
mageck test -k paired_counts.txt -t active1,active2 -c BMSC1,BMSC2 -n paired_active_BMSC --paired

#import in R
paired_active_BMSC.gene_summary <- read.delim("/media/Coco/Collaborations/Eosinophils BD/GWI screen/fastqs/bowtie_output/paired_active_BMSC.gene_summary.txt")
View(paired_active_BMSC.gene_summary)
ggplot(data=paired_active_BMSC.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value))) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

paired_active_BMSC.gene_summary$zscore <- paired_active_BMSC.gene_summary %>% mutate(zscore = (neg.lfc - mean(neg.lfc))/sd(neg.lfc))
head(paired_active_BMSC.gene_summary)

####FGSEA####
#import dataset
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)

m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(KEGG)

#rank genes
dge2 <- read.delim("fastqs/bowtie_output/paired_active_BMSC.gene_summary.txt")
rownames(dge2) <- dge2$id
View(dge2)
#dge2_down <-subset(dge2, neg.p.value < 0.08, select = c(1:14))

#rank dge based on logFC and p-value
ranks <- dge2 %>% 
  na.omit()%>%
  mutate(ranking=neg.lfc)
#mutate(ranking=-log10(neg.p.value)/sign(neg.lfc))
ranks <- ranks$ranking
names(ranks) <- rownames(dge2)
print(head(ranks, 10))

#run gsea for pathway of interest 
KEGG_basal <- fgsea(pathways = KEGG, stats = ranks, minSize=10,maxSize=500,  nperm=1000000)
View(KEGG_basal)

Hallmarks_basal <- fgsea(pathways = Hallmarks, stats = ranks, minSize=10,maxSize=500,  nperm=1000000)
View(Hallmarks_basal)

KEGG_active$pathway<-gsub("KEGG_","",KEGG_active$pathway)
KEGG_active$pathway<-gsub("_"," ",KEGG_active$pathway)

Hallmarks_active$pathway<-gsub("HALLMARK_","",Hallmarks_active$pathway)
Hallmarks_active$pathway<-gsub("_"," ",Hallmarks_active$pathway)

pos_KEGGactive<- which(KEGG_active$pval<0.05)
pos_Hallmarksactive<- which(Hallmarks_active$pval<0.05)

merged_active <- rbind(KEGG_active[pos_KEGGactive,], Hallmarks_active[pos_Hallmarksactive,])
head(merged_active)
merged_active <- arrange(merged_active, -pval) 
merged_active$pathway <- factor(merged_active$pathway, levels = merged_active$pathway)

#plot
merged_active$pathway <- tolower(merged_active$pathway)
merged_active$pathway <- factor(merged_active$pathway, levels = merged_active$pathway)
ggplot(merged_active, aes(x=pathway, y=-NES, size = size)) + 
  geom_point(alpha = 1.0, stat = "identity",aes(colour = pval)) +   scale_size_area(max_size = 8)+
  coord_flip() +theme_classic() + scale_color_gradientn(colours = c("#FF0000", "black")) + 
  theme(legend.position="right") +
  theme(axis.text.y = element_text( size=10, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",hjust=1))+ 
  labs(title = "CD80+PDL1+ vs BMSC", y = "normalized enrichment score", x="")+
  ggsave("GSEA_active.pdf", width = 5, height = 4)


#plot genes
DEGs <- paired_active_BMSC.gene_summary
tnfa_sig <- merged_active[15,8]
tnfa_sig<- tnfa_sig$leadingEdge
tnfa_sig <- c("Bmp2"   ,  "Phlda1" ,  "Serpinb8", "Birc2"   , "Fosl1" ,   "Hes1"    , "Pde4b" ,   "Mxd1"  ,   "Fos"    ,  "Tap1"   ,  "Rhob"   ,  "Ier5"   ,  "Plk2" ,   
              "Litaf" ,   "Tnc"   ,   "Tnip1"   , "Tlr2"  ,   "Ets2"  ,   "Stat5a"   ,"Fjx1"  ,   "Plaur" ,   "Gadd45b" , "Tnfaip3",  "Nfkb1"  ,  "Nampt"   , "Trib1",   
              "Maff"  ,   "Ppp1r15a", "Gpr183" ,  "Ifit2" ,   "Birc3" ,   "Pfkfb3"  , "Tnfaip2" , "Il7r"  ,   "Inhba"  ,  "Spsb1"  ,  "Dusp4"  ,  "Klf9"    , "Panx1",   
              "Ptgs2" ,   "Atf3"  ,   "Dnajb4" ,  "Tgif1" ,   "Dusp1" ,   "Myc"   ,   "Nfkb2" ,   "Pnrc1" ,   "Nr4a3"  ,  "Klf6"  ,   "Csf1"   ,  "Dusp5"  ,  "Relb" ,   
             "Kynu"   ,  "Gadd45a" , "Il12b"   , "Ccnl1"  ,  "Il23a"  ,  "Sod2"   ,  "Fosl2"  ,  "Smad3"  ,  "Ripk2"   , "Egr2"   ,  "Ehd1"   ,  "Serpinb2", "Mcl1"  ,  
              "Irf1"  ,   "Zfp36"  ,  "Il6"    ,  "Ddx58" ,   "Map2k3" ,  "Plau"  ,   "Traf1"  ,  "Eif1"  ,   "Marcks" ,  "Cd83"  ,   "Sgk1"   ,  "Ccl4"   ,  "Fosb" ,   
              "Zbtb10" ,  "Ccl2"   ,  "Cxcl2"   , "Areg"  ,   "Ifih1"  ,  "Cxcl3" ,   "Tnfaip8" , "Bcl6"  ,   "Rnf19b" ,  "Ackr3" ,   "Tnfrsf9" , "Olr1"   ,  "Cxcl5",   
               "Jag1"  ,   "Cxcl1"   , "Irs2"  ,   "Plek"  ,   "Junb"   ,  "Pmepa1"   ,"Ninj1"  ,  "Snn"  ,    "Dram1" ,   "Efna1" ,   "Sdc4"   ,  "Il18"  ,   "Tnip2",   
              "Ier3"   ,  "Btg3" )
tnfa_sig_data <- DEGs[DEGs$id %in% tnfa_sig,]

mapk_sig <- merged_active[14,8]
mapk_sig<-mapk_sig$leadingEdge
mapk_sig <- c("Mecom"   , "Mapkapk5", "Map2k6"  , "Map3k3"   ,"Pla2g2f" , "Tgfbr1" ,  "Fos"   ,   "Elk4"   ,  "Cacnb4"  , "Mknk2"   , "Arrb2"  ,  "Dusp3"   , "Nfatc4" , 
              "Pla2g2e",  "Rps6ka4" , "Map2k4"  , "Cacna2d1", "Fgf9"  ,   "Fgf5",     "Prkacb" ,  "Pla2g4e" , "Gadd45b" , "Map4k2" ,  "Prkaca" ,  "Nfkb1"  ,  "Flnc"  ,  
              "Rac3"   ,  "Mapk11"  , "Pla2g2d" , "Cacng7"  , "Cacng2"  , "Akt2"   ,  "Mapt"   ,  "Prkca"  ,  "Egf"    ,  "Pla2g1b" , "Cacng4"  , "Rap1a"  ,  "Mapk8ip1",
              "Rasgrf2" , "Stk3"   ,  "Hspa1l"  , "Fgfr2"  ,  "Map4k1" ,  "Mapkapk3", "Rap1b" ,   "Pdgfra"  , "Casp3"  ,  "Dusp16" ,  "Map3k6" ,  "Pak1"   ,  "Fgfr4",   
             "Hspa8"   , "Rras2"   , "Max"  ,    "Ptpn5"  ,  "Map2k7"  , "Cdc42"  ,  "Dusp4"  ,  "Bdnf"    , "Ptprr"  ,  "Cacna1e",  "Traf6"   , "Dusp6"   , "Akt1"   , 
             "Taok1"   , "Map3k4" ,  "Fgf2" ,    "Mapk14"  , "Pla2g5"  , "Cacna1s" , "Dusp1"  ,  "Mapk13"  , "Myc"  ,    "Pak2"    , "Nfkb2"   , "Arrb1"   , "Map3k5"  ,
              "Fgf1"   ,  "Braf"   ,  "Ecsit"  ,  "Flnb" ,    "Dusp5"  ,  "Map3k12" , "Cacna1b",  "Relb"  ,   "Fgf18"  ,  "Fgf6"   ,  "Gadd45a" , "Rasgrf1" , "Il1r1"  , 
              "Mapk8ip2" ,"Ddit3" ,   "Dusp8"  ,  "Mapk7"  ,  "Cacna1i" , "Cacna1h" , "Rasa2"  ,  "Mapk9"  ,  "Mapk12" ,  "Daxx"   ,  "Map3k11" , "Stk4"   ,  "Cdc25b" , 
              "Sos2"  ,   "Tab1"   ,  "Rasgrp3" , "Map3k13",  "Rps6ka1" , "Map2k3"  , "Ppp3r2" ,  "Dusp9"  ,  "Ikbkg" ,   "Chuk"   ,  "Flna"   ,  "Fasl"   ,  "Rasgrp2" ,
               "Fgf8" ,    "Fgf14" ,   "Fgf23"  ,  "Cacng3" ,  "Ntf3"  ,   "Fgf17"  ,  "Taok3" ,   "Pla2g12b", "Hspa2"  ,  "Mapkapk2", "Map4k4" ,  "Srf"   ,   "Rasgrp1", 
              "Ppm1a" ,   "Cacnb2" ,  "Ntrk2"   , "Nras"    , "Fgf12"    ,"Map3k2"  , "Tgfb3"  ,  "Fgfr3"  ,  "Ppp3cb"  , "Ppp3cc" ,  "Cacng6"  , "Mknk1"  ,  "Grb2"  ,  
              "Cacng5" ,  "Nfatc2" ,  "Gadd45g"  ,"Pla2g10" , "Mapk8ip3" ,"Egfr"   ,  "Cacng1" ,  "Fgf7"   ,  "Nf1"    ,  "Tnfrsf1a", "Cacna1a" , "Ntrk1"  ,  "Fgf15",   
              "Map3k14" , "Ikbkb" ,   "Gng12"  ,  "Fgf11"   , "Pla2g6" ,  "Pla2g12a", "Tab2"  ,   "Rps6ka6" , "Pdgfa"   , "Rac1"   ,  "Pdgfb"   )
mapk_sig_data <- DEGs[DEGs$id %in% mapk_sig,]

control_data <- DEGs[DEGs$id == "negative_control",]
markers_list <- c("Cd274", "Cd80")
markers_data <- DEGs[DEGs$id %in% markers_list,]

ggplot(data=DEGs, 
       aes(x=id, y=log10(neg.score), size=(neg.lfc - mean(neg.lfc))/sd(neg.lfc))) +
  geom_point(shape = 20) +scale_y_reverse(limits = c(0,-5)) + 
  theme_bw() + scale_size_area(max_size = 6)+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  geom_point(data=tnfa_sig_data, aes(x=id, y=log10(neg.score),  size=(neg.lfc - mean(neg.lfc))/sd(neg.lfc)), colour="#EA0D03") +
  geom_point(data=mapk_sig_data, aes(x=id, y=log10(neg.score),  size=(neg.lfc - mean(neg.lfc))/sd(neg.lfc)), colour="#801B0C") +
  geom_point(data=markers_data, aes(x=id, y=log10(neg.score),  size=(neg.lfc - mean(neg.lfc))/sd(neg.lfc)), colour="#F98400") +
    ggsave("dotplot_active.pdf", width = 8, height = 5)  

 
