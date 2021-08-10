source('Packages_functions.R', echo=TRUE)
readRDS("eosinophil_pure.rds")

#####SUBSET STEADYSTATE ORGANS#####
DimPlot(eosinophil_pure)
Idents(eosinophil_pure) <- "orig.ident"
eosinophils_steadystate <- subset(eosinophil_pure,  idents = c("bonemarrow", "blood", "Spleen", "stomach", "SI", "colon"))
current.cluster.ids <- c("bonemarrow", "blood", "Spleen", "stomach", "SI", "colon")
new.cluster.ids <-  c("bonemarrow", "blood", "spleen", "stomach", "small intestine", "colon") #all lowercase
eosinophils_steadystate$orig.ident <- plyr::mapvalues(x = eosinophils_steadystate$orig.ident, from = current.cluster.ids, to = new.cluster.ids)
eosinophils_steadystate$orig.ident <- factor(x = eosinophils_steadystate$orig.ident, levels = c("bonemarrow", "blood", "spleen", "stomach", "small intestine", "colon"))
eosinophils_steadystate <- NormalizeData(eosinophils_steadystate, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophils_steadystate <- FindVariableFeatures(eosinophils_steadystate, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophils_steadystate)
eosinophils_steadystate <- ScaleData(eosinophils_steadystate, features = all.genes, vars.to.regress = "percent.mt")
eosinophils_steadystate <- RunPCA(eosinophils_steadystate, features = VariableFeatures(object = eosinophils_steadystate))
DimPlot(eosinophils_steadystate, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eosinophils_steadystate)
eosinophils_steadystate <- FindNeighbors(eosinophils_steadystate, dims = 1:20)
eosinophils_steadystate <- FindClusters(eosinophils_steadystate, resolution = 0.3)
eosinophils_steadystate <- RunUMAP(eosinophils_steadystate, dims = 1:20, return.model=TRUE)
DimPlot(eosinophils_steadystate, cols = col_vector, label = T)
FeaturePlot(eosinophils_steadystate, features = c("Mki67", "Camp", "Ltf", "Ly6a2", "Ly6g", "Epx",  "Siglece", "Retnlg", "Retnla", "Cd274"), order=T)
FeaturePlot(eosinophils_steadystate, features = c("Icosl")

#remove cluster mito high and differentiating neutrophils
eosinophils_steadystate <- subset(eosinophils_steadystate,  idents = c(0,1,2,3)) #then run the analysis from line 6

#rename clusters
current.cluster.ids <- c(0, 1, 2, 3, 4,5,6)
new.cluster.ids <-  c("basal eosinophils", "active eosinophils",  "circulating eosinophils", "immature eosinophils", "basal eosinophils",
                      "eosinophil progenitors", "basal eosinophils")
eosinophils_steadystate@meta.data$seurat_clusters <- plyr::mapvalues(x = eosinophils_steadystate@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(eosinophils_steadystate) <- "seurat_clusters"
eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = rev(c("active eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors")))
DimPlot(eosinophils_steadystate, reduction = "umap", pt.size = .5, label=F, cols = col_vector[1:5]) + 
  ggsave("Figures/UMAP.pdf", width = 8, height = 5)

#expression of key markers
length(eosinophils_steadystate@active.ident)
a<-plot_density(eosinophils_steadystate, "Siglecf", pal = "magma")
b<-plot_density(eosinophils_steadystate, "Il5ra",pal = "magma")
c<-plot_density(eosinophils_steadystate, "Ccr3", pal = "magma")
d<-plot_density(eosinophils_steadystate, "Epx", pal = "magma")
ggarrange(a, b, c, d, ncol = 4, nrow = 1) + ggsave("Figures/keymarkers.pdf", width = 18, height = 4)

#######COMPOSITIONAL ANALYSIS######
#frequencies per cluster
numberofcells         <- table(eosinophils_steadystate$orig.ident, eosinophils_steadystate$seurat_clusters)
numberofcells
totalcellsperorgan   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]), sum(numberofcells[3,]), sum(numberofcells[4,]),
                          sum(numberofcells[5,]), sum(numberofcells[6,]))
a                     <- cbind(numberofcells,totalcellsperorgan)
a
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]), sum(a[,5]), sum(a[,6]))
b                     <- rbind(a, totalcellspercluster)
b

c0 <- (b[1:6,1]/totalcellsperorgan)*100
c1 <- (b[1:6,2]/totalcellsperorgan)*100
c2 <- (b[1:6,3]/totalcellsperorgan)*100
c3 <- (b[1:6,4]/totalcellsperorgan)*100
c4 <- (b[1:6,5]/totalcellsperorgan)*100
c5 <- (b[1:6,6]/totalcellsperorgan)*100

c <- rbind(c0,c1,c2,c3,c4)
colSums(c)
rownames(c) =  rev(c("active eosinophils", "basal eosinophils", "circulating eosinophils", "immature eosinophils", "eosinophil progenitors"))
c

#plot
par(mar=c(5,8,2,14))
pdf(file="Figures/Clusterbreakdown.pdf")
barplot(c, horiz=TRUE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=180, cex=.8),
        main = "Cluster breakdown per organ", 
        las = 1, 
        col= rev(col_vector[1:5]) )
dev.off()

#LOCAL DENSITY PLOT
Idents(eosinophils_steadystate) <- "orig.ident"
col_vector[1:7]
DimPlot(eosinophils_steadystate, order=T, group.by = "orig.ident", pt.size = .2, label=F, cols = col_vector) + theme(legend.position = "right") + labs(title=" ")
a <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="black", "blood"="grey",  "spleen"="grey", "stomach"="grey", "small intestine"="grey",
                       "colon"="grey"),
             order = c("bonemarrow", "blood", "spleen", "stomach", "small intestine", "colon"))  + theme_void()+labs(title="bone marrow")+ theme(legend.position = "none")
b <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="grey", "blood"="black",  "spleen"="grey", "stomach"="grey", "small intestine"="grey",
                       "colon"="grey"),
             order = c("blood", "bonemarrow", "spleen", "stomach", "small intestine", "colon"))  +theme_void()+ labs(title="blood")+ theme(legend.position = "none")
c <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="grey", "blood"="grey",  "spleen"="black", "stomach"="grey", "small intestine"="grey",
                       "colon"="grey"),
             order = c("spleen", "blood", "bonemarrow", "stomach", "small intestine", "colon"))  + theme_void()+labs(title="spleen")+ theme(legend.position = "none")
d <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="grey", "blood"="grey",  "spleen"="grey", "stomach"="black", "small intestine"="grey",
                       "colon"="grey"),
             order = c( "stomach", "blood", "bonemarrow", "spleen", "small intestine", "colon")) + theme_void()+labs(title="stomach")+ theme(legend.position = "none")
e <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="grey", "blood"="grey",  "spleen"="grey", "stomach"="grey", "small intestine"="black",
                       "colon"="grey"),
             order = c( "small intestine", "blood", "bonemarrow", "spleen", "stomach", "colon"))+ theme_void()+labs(title="small intestine")+ theme(legend.position = "none") 
f <- DimPlot(eosinophils_steadystate, reduction = "umap", group.by = c("ident"), pt.size=.2,
             cols = c( "bonemarrow"="grey", "blood"="grey",  "spleen"="grey", "stomach"="grey", "small intestine"="grey",
                       "colon"="black"),
             order = c("colon", "blood", "bonemarrow", "spleen", "stomach", "small intestine"))+labs(title="colon") + theme_void()+ theme(legend.position = "none") 

ggarrange(a, b, c, d, e, f, ncol = 3, nrow = 2)+ ggsave("Figures/organUMAP.pdf", width = 12, height = 8)


###MARKERS####
markers_steadystate <- FindAllMarkers(object = eosinophils_steadystate, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_steadystate %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
top200 <- markers_steadystate %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(top200,"markers_steadystate.csv", row.names = TRUE)

#markers density
a<-plot_density(eosinophils_steadystate, "Mki67", pal = "magma")
b<-plot_density(eosinophils_steadystate, "S100a9", pal = "magma")
c<-plot_density(eosinophils_steadystate, "Cd24a", pal = "magma")
d<-plot_density(eosinophils_steadystate, "Siglece", pal = "magma")
e<-plot_density(eosinophils_steadystate, "Cd274", pal = "magma")
ggarrange(a, b, c,d,e, ncol = 3, nrow = 2) + ggsave("Figures/MarkersFeature.pdf", width = 17, height = 10)

#steadystate dotplot
final.markers <- c("Mki67", "Tuba1b", "Epx", "Prg3", "Prg2","Ear1","Ear2", "Ear6",  "Cd63", "Cebpe",
                   "Alox15", "Aldh2", "S100a9", "S100a6", "S100a10", "Il5", "Retnla", "Ccl9", "Il1rl1", 
                   "Cd24a", "Mmp9", "Icosl", "Il4", "Tgfb1", "Pirb", "Rara", "Cd80", "Cd274", "Ptgs2", "Il1rn", "Il1b", 
                   "Vegfa", "Ccl3", "Cxcl2", "Il16", "Tnf")

eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("active eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors"))
Idents(eosinophils_steadystate) <- "seurat_clusters"
DotPlot(eosinophils_steadystate, features = final.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")+
  ggsave("Figures/steadystate_dotplot.pdf", width = 15, height = 3.5)

###SIGNATURES######
#Cell cycle score
eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("eosinophil progenitors", "immature eosinophils", "circulating eosinophils", "basal eosinophils", "active eosinophils"))
DoHeatmap(eosinophils_steadystate, features = c("Mki67", "Cdk1", "Pcna", "Cdt1", "Fbxo5", "Spc24", "Ranbp1", "Rad21", 
                                                "Nusap1", "Cdc20", "Pmf1", "Cdc45", "Cenpf", "Smc4", "Tpx2", "Cdk2",
                                                "E2f8", "Top2a", "Stmn1","Nuf2"),    label=F, draw.lines	= T,    group.by = "seurat_clusters",
          lines.width = 100, group.colors	= rev(col_vector[1:5]))+  scale_fill_gradientn(colors = pal,  na.value = "white") + theme(axis.text.y = element_text(face = "italic", color="black") ) +
  ggsave("Figures/progenitors_heatmap.pdf", width = 7, height = 2.8)

eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("active eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors"))

tmp <- CellCycleScoring(
  object = eosinophils_steadystate,
  g2m.features = g2m.genes,
  s.features = s.genes
)

s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)

eosinophils_steadystate <- AddModuleScore(eosinophils_steadystate, features = cc.genes, name = "CC")
names(x = eosinophils_steadystate[[]])

RidgePlot(eosinophils_steadystate, features="CC1", group.by = "seurat_clusters", rev(col_vector[1:5])) +
  theme_classic() +  
  theme(text = element_text(size=25)) + labs (title = "Cell cycle score ", y = " ", x= " ") +theme(legend.position="none")+
  ggsave("Figures/CC_ridge.pdf", width = 8, height = 6)

#test
Idents(eosinophils_steadystate) <- "seurat_clusters"
eos_prog <- subset(eosinophils_steadystate, idents = c("eosinophil progenitors"))
eos_immature <- subset(eosinophils_steadystate, idents = c("immature eosinophils"))
wilcox.test(eos_prog$CC1, eos_immature$CC1, alternative = "two.sided") #p-value < 2.2e-16

#stemness score (Koeva et al, 2011)
stemness_list <- list(c("Orc1l", "Impdh2", "Cct5", "Nap1l1", "Ccnd2", "Smo", "Mcm4", "Mcm5", "Hells", "Hnrnpa2b1", "Cct8", "Col18a1", "Sfrs3", 
                        "Rrm2", "Bub1", "Ncl", "Kpna2", "Shmt1", "Ipo5", "Ruvbl1", "Shroom3", "Dnahc11", "Cdc6", "Ttk", "Cks2", "Mcm2", "Fignl1", 
                        "Dph5", "Cdt1", "Cct3", "Eya2", "Pcna", "Set", "Prps1", "Fbl", "Dtymk", "Ssbp1", "Depdc6", "Top2a", "Csrp2"))

eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= stemness_list,name = "Stemness")
names(x = eosinophils_steadystate[[]])

VlnPlot(eosinophils_steadystate, features= "Stemness1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " Stemness score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Stemness_violin.pdf", width = 8, height = 6)

#test
wilcox.test(eos_prog$Stemness1, eos_immature$Stemness1, alternative = "two.sided") #p-value < 2.2e-16

#Granules synthesis
Granules_synthesis_list <-   list(c("Prg2","Prg3",  "Epx", "Ear6", "Ear1", "Ear2"))
eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= Granules_synthesis_list,name = "GranulesSynthesis")
VlnPlot(eosinophils_steadystate, features="GranulesSynthesis1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + ylim(-1.3,5)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = " Granulopoiesis score",  title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/GranulesSynthesis_violin.pdf", width = 8, height = 6)

#test
eos_circ <- subset(eosinophils_steadystate, idents = c("circulating eosinophils"))
eos_basal <- subset(eosinophils_steadystate, idents = c("basal eosinophils"))
eos_active <- subset(eosinophils_steadystate, idents = c("active eosinophils"))
wilcox.test(eos_prog$GranulesSynthesis1, eos_immature$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_circ$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_basal$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_active$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16

           
#FGSEA#####
Idents(eosinophils_steadystate) <- "seurat_clusters"
DimPlot(eosinophils_steadystate)

prog_markers <- FindMarkers(object = eosinophils_steadystate, ident.1="eosinophil progenitors", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
prog_markers$p_val_adj[prog_markers$p_val_adj == 0] <- 2.225074e-308 #replace 0 with lowest number
BP_progenitors <- preranked_BP(prog_markers)
sig_BP_progenitors <- as.data.frame(BP_progenitors %>% filter(padj<0.05))
View(sig_BP_progenitors)

immature_markers <- FindMarkers(object = eosinophils_steadystate, ident.1="immature eosinophils", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
immature_markers$p_val_adj[immature_markers$p_val_adj == 0] <- 2.225074e-308 #replace 0 with lowest number
BP_immature <- preranked_BP(immature_markers)
sig_BP_immature <- BP_immature %>% filter(padj<0.05)
View(sig_BP_immature)

basal_markers <- FindMarkers(object = eosinophils_steadystate, ident.1="basal eosinophils", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
basal_markers$p_val_adj[basal_markers$p_val_adj == 0] <- 2.225074e-308 #replace 0 with lowest number
BP_basal <- preranked_BP(basal_markers)
View(BP_basal%>% filter(pval<0.05))
sig_BP_basal <- BP_basal %>% filter(abs(NES)>1 & padj<0.05)

active_markers <- FindMarkers(object = eosinophils_steadystate, ident.1="active eosinophils", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
active_markers$p_val_adj[active_markers$p_val_adj == 0] <- 2.225074e-308 #replace 0 with lowest number
BP_active <- preranked_BP(active_markers)
View(BP_active %>% filter(padj<0.05))
sig_BP_active <- BP_active %>% filter(abs(NES)>1 & padj<0.05)


circulating_markers <- FindMarkers(object = eosinophils_steadystate, ident.1="circulating eosinophils", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
circulating_markers$p_val_adj[circulating_markers$p_val_adj == 0] <- 2.225074e-308 #replace 0 with lowest number
BP_circulating <- preranked_BP(circulating_markers)
View(BP_circulating%>% filter( pval<0.05))
sig_BP_circulating <- BP_circulating %>% filter(abs(NES)>1 & pval<0.05) #ATTENTION PVAL

#plot GSEA as heatmap with enrichment as rows (color=NES) and cluster as columns
pos_progenitors <- which(BP_progenitors$pathway ==  "CELL DIVISION"| 
                           BP_progenitors$pathway ==  "CELL CYCLE"| 
                           BP_progenitors$pathway ==  "DNA REPLICATION"| 
                           BP_progenitors$pathway == "EXOCYTOSIS" | 
                           BP_progenitors$pathway ==  "CELL ACTIVATION"| 
                           BP_progenitors$pathway == "CYTOKINE SECRETION"| 
                           BP_progenitors$pathway == "CELL ADHESION MEDIATED BY INTEGRIN")

pos_immature <- which(BP_immature$pathway ==  "PROTEIN TARGETING TO MEMBRANE" | 
                        BP_immature$pathway == "EXOCYTOSIS"|
                        BP_immature$pathway == "PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM"|
                        BP_immature$pathway == "MYELOID LEUKOCYTE ACTIVATION"|
                        BP_immature$pathway =="SECRETION"|
                        BP_immature$pathway =="CELL MATURATION")

pos_basal<- which(BP_basal$pathway ==  "NEGATIVE REGULATION OF LYMPHOCYTE ACTIVATION"|
                    BP_basal$pathway ==  "REGULATION OF WOUND HEALING"| 
                    BP_basal$pathway ==  "PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM"|
                    BP_basal$pathway ==  "CELLULAR EXTRAVASATION"|
                    BP_basal$pathway ==  "LEUKOCYTE DIFFERENTIATION"|
                    BP_basal$pathway ==  "CELL ACTIVATION"|
                    BP_basal$pathway ==  "IMMUNE EFFECTOR PROCESS"|
                    BP_basal$pathway ==  "RESPONSE TO CYTOKINE"|
                    BP_basal$pathway ==  "SECRETION"|
                    BP_basal$pathway ==  "RESPONSE TO MOLECULE OF BACTERIAL ORIGIN"|
                    BP_basal$pathway ==  "EXOCYTOSIS"|
                    BP_basal$pathway ==  "LEUKOCYTE MIGRATION"|
                    BP_basal$pathway ==  "LEUKOCYTE CHEMOTAXIS"|
                    BP_basal$pathway ==  "GRANULOCYTE CHEMOTAXIS"|
                    BP_basal$pathway ==  "MYELOID LEUKOCYTE MIGRATION"|
                    BP_basal$pathway ==  "GRANULOCYTE MIGRATION" |
                    BP_basal$pathway == "LYMPHOCYTE MEDIATED IMMUNITY" | 
                    BP_basal$pathway == "NEGATIVE REGULATION OF MAP KINASE ACTIVITY" |
                    BP_basal$pathway == "GRANULOCYTE CHEMOTAXIS" |
                    BP_basal$pathway ==  "INFLAMMATORY RESPONSE")


pos_actuve <- which(BP_active$pathway ==   "RESPONSE TO LIPID"|
                     BP_active$pathway == "RESPONSE TO MOLECULE OF BACTERIAL ORIGIN" |
                     BP_active$pathway == "INFLAMMATORY RESPONSE"|
                     BP_active$pathway == "MYELOID LEUKOCYTE MIGRATION"|
                     BP_active$pathway == "RESPONSE TO MOLECULE OF BACTERIAL ORIGIN"|
                     BP_active$pathway == "RESPONSE TO CYTOKINE"|
                     BP_active$pathway == "RESPONSE TO TUMOR NECROSIS FACTOR"|
                     BP_active$pathway == "LEUKOCYTE CHEMOTAXIS"|
                     BP_active$pathway == "GRANULOCYTE CHEMOTAXIS"|
                     BP_active$pathway == "PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM"|
                     BP_active$pathway =="POSITIVE REGULATION OF REACTIVE OXYGEN SPECIES METABOLIC PROCESS"|
                     BP_active$pathway =="REGULATION OF WOUND HEALING"|
                     BP_active$pathway =="REGULATION OF ADAPTIVE IMMUNE RESPONSE"|
                     BP_active$pathway =="RESPONSE TO INTERLEUKIN 1"|
                     BP_active$pathway =="P38MAPK CASCADE"|
                     BP_active$pathway =="NIK NF KAPPAB SIGNALING"|
                     BP_active$pathway =="REGULATION OF CELL MATRIX ADHESION")

pos_circulating <- which(BP_circulating$pathway == "REGULATION OF CELL POPULATION PROLIFERATION" |
                           BP_circulating$pathway == "REACTIVE OXYGEN SPECIES BIOSYNTHETIC PROCESS" |
                     BP_circulating$pathway == "PROTEIN LOCALIZATION TO ENDOPLASMIC RETICULUM"|
                     BP_circulating$pathway == "RESPONSE TO MECHANICAL STIMULUS"|
                     BP_circulating$pathway == "LEUKOCYTE MIGRATION"|  
                      BP_circulating$pathway == "REGULATION OF MYELOID CELL DIFFERENTIATION")



#plot
merged <- merge(BP_progenitors[pos_progenitors,c(1,5)], BP_immature[pos_immature,c(1,5)], by = "pathway" , all = T,
                suffixes = c(".progenitors",".immature"))
merged2 <- merge(merged, BP_circulating[pos_circulating,c(1,5)], by = "pathway" , all = T,
                 suffixes = c(".progenitors",".immature"))
merged3 <- merge(merged2,  BP_basal[pos_basal,c(1,5)], by = "pathway" , all = T, 
                 suffixes = c(".circulating",".basal"))
merged4 <- merge(merged3, BP_active[pos_active,c(1,5)], by = "pathway" , all = T)
colnames(merged4) <- c("pathway","progenitors", "immature", "circulating", "basal", "active")
merged4[is.na(merged4)] <- 0
rownames(merged4) <- merged4$pathway
rownames(merged4)<- tolower(rownames(merged4))
merged4[,1] <- NULL


paletteLength <- 100
myColor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
myBreaks = c(seq(min(merged4), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(merged4)/paletteLength,
                      max(merged4),
                      length.out=floor(paletteLength/2)))
pheatmap(merged4, cluster_cols = F, cluster_rows = T, border_color="grey", 
         cellwidth =30, col=myColor, breaks =myBreaks, main = "Biological process enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10, filename = "Figures/GO.pdf")





##GATING STRATEGY####
#1. import list of surface proteins and add colnames
surface_markers <- read.delim("surface_markers.txt", header=FALSE, row.names=1, comment.char="#")
head(surface_markers)
colnames(surface_markers)[1] <- "Gene"
head(surface_markers)

#3. subset clusterN_marker dataframe by only keeping rows that are in the surface_markers df
prog_surface_markers <- prog_markers[rownames(prog_markers) %in% surface_markers$Gene,]
prog_surface_markers$Gene <- rownames(prog_surface_markers)
immature_surface_markers <- immature_markers[rownames(immature_markers) %in% surface_markers$Gene,]
immature_surface_markers$Gene <- rownames(immature_surface_markers)
circulating_surface_markers <- circulating_markers[rownames(circulating_markers) %in% surface_markers$Gene,]
circulating_surface_markers$Gene <- rownames(circulating_surface_markers)
basal_surface_markers <- basal_markers[rownames(basal_markers) %in% surface_markers$Gene,]
basal_surface_markers$Gene <- rownames(basal_surface_markers)
intestinal_surface_markers <- intestinal_markers[rownames(intestinal_markers) %in% surface_markers$Gene,]
intestinal_surface_markers$Gene <- rownames(intestinal_surface_markers)
View(intestinal_surface_markers)

#merge and make heatmap
merged <- merge(prog_surface_markers[,c(2,6)], immature_surface_markers[,c(2,6)], by = "Gene" , all = T,
      suffixes = c(".progenitors",".immature"))
View(merged4)
merged2 <- merge(merged, circulating_surface_markers[,c(2,6)], by = "Gene" , all = T,
                suffixes = c(".progenitors",".immature"))
merged3 <- merge(merged2, basal_surface_markers[,c(2,6)], by = "Gene" , all = T, 
                 suffixes = c(".circulating",".basal"))
merged4 <- merge(merged3, intestinal_surface_markers[,c(2,6)], by = "Gene" , all = T)
colnames(merged4) <- c("Gene","progenitors", "immature", "circulating", "basal", "intestinal")
merged4[is.na(merged4)] <- 0
row.names(merged4) <- merged4$Gene
merged4[1] <- NULL

pheatmap(merged4, cluster_cols = F, border_color="black",fontsize_number = 0.5, cellwidth =30, cellheight = 10, scale = "row", filename="Figures/surfacemarkers_heatmap.pdf")

##CITEK
cytek.markers <- c(  "Clec12a", "Siglece", "Itga4",    "Cd274",   "Cd80",  "Ly6a",
                     "Cd9",  "Fas", "Icam1", 
                     "Itgax",  "Lamp1", "Fcgr3",
                     "Slc3a2", "Cd74", "Cd44", "Il1rl1", "Cd24a",
                     "Itgam", "Epcam", "Pecam1",  "Sell", "Cdh1",  "Ly6c1","Spn","Cd53", "Cd63", "Itga2b") 

plot <- DotPlot(eosinophils_steadystate, features = cytek.markers)
plot + theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+ labs(title = "Cytek markers", y = "", x="")+
  ggsave("cytekmarkers.png", width = 12, height = 3)

#FACS.markers.ordered <- c( "Itgb2","Cd33", "Cd53", "Hmgb1", #gmps
                           "Pecam1", "Spn", "Cd63", "Itga2b" , "Clec12a" ,#progenitors
                           "Lamp1",  "Cdh1",  #immature
                           "Siglece", "Itgal", "Csf1r", "Ly6c1", #basal
                           "Cd74","Fcgr2b", "Epcam",  #fos jun high
                           "Icam1", "Slc3a2", "Ccr1", "Itgax", "Notch1","Notch2", "Ldlr","Fas", "Cxcr4",
                           "Cd274","Cd80","Cd9","Fcgr3", #cd274
                           "Sell", "Cd48", "Itgb1", "Itgb3", "Il1rl1", "Cd24a","Cd44") #blood

#plot <- DotPlot(eosinophils_steadystate, features = FACS.markers.ordered)
#plot + theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
#  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+ labs(title = "surface markers", y = "", x="")

###FACS MARKERS
Idents(eosinophils_steadystate) <- "seurat_clusters"
basal_intestinal <- subset(eosinophils_steadystate, idents = c("intestinal eosinophils", "basal eosinophils"))
surface.markers <- c("Cd274", "Cd80","Cd9","Pecam1", "Icam1", "Fas", "Ly6a", "Itgax", "Itga4", "Clec12a", "Siglece")
DotPlot(basal_intestinal, features = selected_markers,  dot.scale = 40)
plot + theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1,size=20), axis.text.y = element_text(face="bold", size=20)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+ labs(title = "FACS markers", y = "", x="")+
  ggsave("FACSmarkers.pdf", width = 12, height = 3)

data <- AverageExpression(basal_intestinal, features =  surface.markers)
data
pheatmap(data$RNA[,c(2,1)], cluster_columns = F, scale = "row", 
         fontsize_number = 0.5, cellwidth =30, cellheight = 30, treeheight_row=0, treeheight_col=0,	
         border_color = "white", color = colorRampPalette(rev(brewer.pal(n = 4, name =  "RdYlBu")))(100), filename = "FACSmarkers_heatmap_h.pdf")


##blended featureplot
FeaturePlot(eosinophils_steadystate, features = c("Cd80", "Cd274"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T) +
  ggsave("Figures/markers_blend.pdf", width = 16, height = 5)

FeaturePlot(eosinophils_steadystate, features = c("Cd80", "Clec12a"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T) +
  ggsave("Figures/markers_blend2.pdf", width = 16, height = 5)




