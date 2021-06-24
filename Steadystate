source('~/NAS/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Packages_functions.R', echo=TRUE)

#####SUBSET STEADYSTATE ORGANS#####
DimPlot(eosinophil_pure)
Idents(eosinophil_pure) <- "orig.ident"
eosinophils_steadystate <- subset(eosinophil_pure,  idents = c("bonemarrow", "blood", "Spleen", "stomach", "SI", "colon"))
current.cluster.ids <- c("bonemarrow", "blood", "Spleen", "stomach", "SI", "colon")
new.cluster.ids <-  c("bonemarrow", "blood", "spleen", "stomach", "small intestine", "colon")
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

#remove cluster mito high and diff neu
eosinophils_steadystate <- subset(eosinophils_steadystate,  idents = c(0,1,2,3)) #then run the analysis from line 6

#rename clusters
current.cluster.ids <- c(0, 1, 2, 3, 4,5,6)
new.cluster.ids <-  c("basal eosinophils", "intestinal eosinophils",  "circulating eosinophils", "immature eosinophils", "basal eosinophils",
                      "eosinophil progenitors", "basal eosinophils")
eosinophils_steadystate@meta.data$seurat_clusters <- plyr::mapvalues(x = eosinophils_steadystate@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(eosinophils_steadystate) <- "seurat_clusters"
eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = rev(c("intestinal eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors")))
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
rownames(c) =  rev(c("intestinal eosinophils", "basal eosinophils", "circulating eosinophils", "immature eosinophils", "eosinophil progenitors"))
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
write.csv(top200,"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Steadystate/markers_steadystate.csv", row.names = TRUE)

#steadystate dotplot
final.markers <- c("Mki67", "Tuba1b", "Epx", "Prg3", "Prg2","Ear1","Ear2", "Ear6",  "Cd63", "Cebpe",
                   "Alox15", "Aldh2", "S100a9", "S100a6", "S100a10", "Il5", "Retnla", "Ccl9", "Il1rl1", 
                   "Cd24a", "Mmp9", "Icosl", "Il4", "Tgfb1", "Cd80", "Cd274", "Ptgs2", "Il1rn", "Il1b", 
                   "Vegfa", "Ccl3", "Cxcl2", "Il16", "Tnf")

eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("intestinal eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors"))
Idents(eosinophils_steadystate) <- "seurat_clusters"
DotPlot(eosinophils_steadystate, features = final.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")+
  ggsave("Figures/steadystate_dotplot.pdf", width = 15, height = 3.5)

#top20 <- markers_steadystate %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#DoHeatmap(eosinophils_steadystate, features = top20$gene,    label=F, draw.lines	= T,    group.by = "seurat_clusters",
       #   lines.width = 100, group.colors	= col_vector)+  scale_fill_gradientn(colors = pal,  na.value = "white") + theme(axis.text.y = element_text(face = "italic") )

#selected.markers <- c("Mki67", "Spi1", "Chil3", "S100a8", "S100a9", "Prtn3", "Elane", "Tuba1b", "Cebpe", "Prg3",  "Prg2",  "Epx", 
#                      "Ear6", "Ear1", "Ear2", "Cd63", "Alox15","Gata2","Retnla","Retnlg","Ltb", "Mmp9", "Ffar2",  "Icosl","Il1rn", "Il1b", 
#                      "Cd274", "Cd80", "Ccl3", "Cxcl2", "Vegfa", "Csf2rb", "Il10rb","Ahr", "Ptafr", "Ldlr")

#selected.markers.2 <- c("Alox15", "Vegfa", "Il1b", "Tnf", "Tnfrsf1a", "Il4", "Il13", "Il6", "Osm", "Ccl3", "Cxcl2", "Ifng", "Ccl9",
#                       "Cxcr4","Il10rb","Ifngr1", "Tgfbr2", "Csf2rb","Rela", "Myd88", "Stat3", "Stat1", "Sirpa", "Cd300a", "Cd300f", "Cd300lf", 
#                        "Aldh1", "Aldh2", "Ahr")

#cytokines <- c("Il1b", "Il2", "Il3","Il4", "Il5", "Il6", "Il10", "Il11", "Il12a", "Il13","Il15", "Il16",
#               "Il18", "Il23a", "Il25", "Il33", 
#               "Csf2", "Tnf", "Ifng", "Osm", "Tgfb1")

#chemokines <- c("Ccl3", "Ccl5", "Ccl7","Ccl11", "Ccl13", "Ccl17", "Ccl22", "Ccl23", "Cxcl1", "Cxcl2", "Cxcl3",
#                "Cxcl5", "Cxcl8", "Cxcl9", "Cxcl10", "Cxcl11")

#plot1 <- DotPlot(eosinophils_steadystate, features = selected.markers , dot.scale = 10) + RotatedAxis()+ 
#  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
#  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")+ labs(title = "selected markers", y = "", x="")

#plot2 <- DotPlot(eosinophils_steadystate, features = selected.markers.2 , dot.scale = 10) + RotatedAxis()+ 
#  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")+ labs(title = "selected markers", y = "", x="")
#  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 

#plot3 <- DotPlot(eosinophils_steadystate, features = cytokines , dot.scale = 10) + RotatedAxis() + 
#  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
#  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")+ labs(title = "cytokines", y = "", x="")

#plot4 <- DotPlot(eosinophils_steadystate, features = chemokines , dot.scale = 10) + RotatedAxis() +
#  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
#  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "chemokines", y = "", x="")

#ggarrange(plot1, plot2, plot3, plot4, ncol = 1, nrow = 4)

#il10.markers <- c("Il10ra", "Il10rb", "Stat3")
#DotPlot(eosinophils_steadystate, features = il10.markers , dot.scale = 10) + RotatedAxis() +
#  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="") +
#  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
#  ggsave("Il10_dotplot.pdf", width = 8, height = 3)


###PROGENITORS######
#Cell cycle score
eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("eosinophil progenitors", "immature eosinophils", "circulating eosinophils", "basal eosinophils", "intestinal eosinophils"))
DoHeatmap(eosinophils_steadystate, features = c("Mki67", "Cdk1", "Pcna", "Cdt1", "Fbxo5", "Spc24", "Ranbp1", "Rad21", 
                                                "Nusap1", "Cdc20", "Pmf1", "Cdc45", "Cenpf", "Smc4", "Tpx2", "Cdk2",
                                                "E2f8", "Top2a", "Stmn1","Nuf2"),    label=F, draw.lines	= T,    group.by = "seurat_clusters",
          lines.width = 100, group.colors	= rev(col_vector[1:5]))+  scale_fill_gradientn(colors = pal,  na.value = "white") + theme(axis.text.y = element_text(face = "italic", color="black") ) +
  ggsave("progenitors_heatmap.pdf", width = 7, height = 2.8)

eosinophils_steadystate$seurat_clusters <- factor(x = eosinophils_steadystate$seurat_clusters, levels = c("intestinal eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors"))

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

#VlnPlot(eosinophils_steadystate, features= "CC1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
#  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() +
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
#  labs(title = "", y = "", x="") + theme(legend.position="right") +  
#  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
#  ggsave("Figures/CC_violin.pdf", width = 8, height = 6)

#test
Idents(eosinophils_steadystate) <- "seurat_clusters"
eos_prog <- subset(eosinophils_steadystate, idents = c("eosinophil progenitors"))
eos_immature <- subset(eosinophils_steadystate, idents = c("immature eosinophils"))
wilcox.test(eos_prog$CC1, eos_immature$CC1, alternative = "two.sided") #p-value < 2.2e-16

#S <- RidgePlot(tmp, features=(vars = c("S.Score")), group.by = "seurat_clusters", rev(col_vector[1:5])) +
#  theme_classic()+
#  theme(text = element_text(size=25)) + labs (title = "S phase score", y = "", x= " ")+theme(legend.position="none")

#G2M<-RidgePlot(tmp, features=(vars = c("G2M.Score")), group.by = "seurat_clusters", rev(col_vector[1:5])) +
#  theme_classic() +
#  theme(text = element_text(size=25)) + labs (title = "G2M phase score", y = "", x= " ") +theme(legend.position="none")

#ggarrange(S, G2M, ncol = 2, nrow = 1)+
#  ggsave("CellCycle.pdf", width = 12, height = 3)


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


###APOPTOSIS####
apoptosis.score <- list(c(BP[["GOBP_APOPTOTIC_PROCESS"]]))
eosinophils_steadystate <- AddModuleScore(eosinophils_steadystate, features= apoptosis.score,name = "apoptosis.score")
VlnPlot(eosinophils_steadystate, features="apoptosis.score1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Apoptosis score ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Apoptosis_violin.pdf", width = 8, height = 6)

####GRANULES AND ANTIMICROBIAL####
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
eos_intestinal <- subset(eosinophils_steadystate, idents = c("intestinal eosinophils"))
wilcox.test(eos_prog$GranulesSynthesis1, eos_immature$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_circ$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_basal$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_prog$GranulesSynthesis1, eos_intestinal$GranulesSynthesis1, alternative = "two.sided") #p-value < 2.2e-16

#antimicrobial synthesis
Antimicrobial_list <- list(c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                             "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= Antimicrobial_list,name = "Antimicrobial")
names(x = eosinophils_steadystate[[]])

VlnPlot(eosinophils_steadystate, features="Antimicrobial1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ylim(-0.5,3.4)+
  labs(y = "Antimicrobial signature", title = "", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  annotate(geom="text", x=0.5, y=3, label="(S100a8, S100a9, Gbp2, Prg2, Epx, Ear1, \nEar2, Ear6, Gbp7, Ltf, Lcn2, Lyz2, Irgm1, \nCamp, Adam17, Serpin1, H2-D1, \nH2-T23, H2-Q7)", 
           color="black", size=5,   fontface="italic",  hjust = 0)+
  ggsave("Figures/Antimicrobial_violin.pdf", width = 8, height = 6)

#heatmap
antimicrobial_data <- AverageExpression(eosinophils_steadystate, features = c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                                                                              "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))

pheatmap(antimicrobial_data$RNA, scale = 'row', cluster_rows = F, angle_col = "45", filename = "Figures/Antimicrobial_heatmap.pdf")

#ifng regulated antimicrobial
IFNg_regulated <- list(c("Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= IFNg_regulated,name = "IFNg_regulated")
names(x = eosinophils_steadystate[[]])

VlnPlot(eosinophils_steadystate, features="IFNg_regulated1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ylim(-0.5,3.4)+
  labs(title = "", y = "IFNg-regulated antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  annotate(geom="text", x=0.5, y=3, label="(Gbp2, Gbp7, Irgm1, Serpin1, H2-D1, H2-T23, H2-Q7)", 
           color="black", size=5,   fontface="italic",  hjust = 0)+
  ggsave("Figures/Ifng_Antimicrobial_violin.pdf", width = 8, height = 6)


####BASAL vs. INTESTINAL COMPARISON#####
basal_intestinal <- subset(eosinophils_steadystate, idents = c("intestinal eosinophils", "basal eosinophils"))

#receptors
selected_markers <- c("Ptafr", "Ahr", "Fcgr3", "Fgfr1", "Ccr1", "Cxcr4", "Csf2rb", "Csf2rb2", "Il10ra", "Ifngr1", "Tgfbr2")
DotPlot(basal_intestinal, features = selected_markers , dot.scale = 15) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1, size=20), axis.text.y = element_text(face="bold", size=20)) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = " ", y = "", x="") +
  ggsave("basalintestinal_dotplot.pdf", width = 7, height = 3.5)

#NFKB transcription factors
Nfkb.markers <- c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz", "Rela", "Relb")
DotPlot(eosinophils_steadystate, features = Nfkb.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="") +
  ggsave("Nfkb_dotplot.pdf", width = 7, height = 2.8)

#Nfkb_list <- list(c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz",  "Rela", "Relb"))
#eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= Nfkb_list,name = "Nfkb")
#names(x = eosinophils_steadystate[[]])
#VlnPlot(eosinophils_steadystate, features="Nfkb1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = "Nfkb activity score", x="", title=" ") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Nfkb_activity_violin.pdf", width = 8, height = 4)

#regulatory proteins
a<-plot_density(eosinophils_steadystate, "Cd80", pal = "magma")
b<-plot_density(eosinophils_steadystate, "Cd9", pal = "magma")
c<-plot_density(eosinophils_steadystate, "Cd274", pal = "magma")
ggarrange(a, b, c, ncol = 3, nrow = 1) + ggsave("Figures/CostimulatoryMarkers.pdf", width = 15, height = 4)
  
Regulatory_list <-  list(c("Icosl", "Cd274", "Cd80", "Cd9"))
eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= Regulatory_list,name = "Regulatory")
names(x = eosinophils_steadystate[[]])
VlnPlot(eosinophils_steadystate, features="Regulatory1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = " ", y = "Immune regulatory score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Regulatory_violin.pdf", width = 8, height = 4)

#test
wilcox.test(eos_intestinal$Regulatory1, eos_basal$Regulatory1, alternative = "two.sided") #p-value < 2.2e-16

