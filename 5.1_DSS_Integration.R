source('Packages_functions.R', echo=TRUE)

#####IMPORT DATASET (subsetted eosinophils from dataset from Schwarzfischer et al, 2021)###
eosinophils <- readrds("eosinophils.Rds")
DimPlot(eosinophils, split.by = "condition")
Idents(eosinophils) <- "condition"
eosinophils_DSS <- subset(eosinophils, idents = c("wtDSS"))
DSS <- CreateSeuratObject(eosinophils_DSS@assays$RNA@counts, project ="DSS")
DSS <- NormalizeData(DSS, normalization.method = "LogNormalize", scale.factor = 10000)
DSS <- FindVariableFeatures(DSS, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(DSS)
DSS[["percent.mt"]] <- PercentageFeatureSet(DSS, pattern = "^mt.")
DSS <- ScaleData(DSS, features = all.genes, vars.to.regress = "percent.mt")
DSS <- RunPCA(DSS, features = VariableFeatures(object = DSS))
DimPlot(DSS, reduction = "pca", group.by = "orig.ident")
ElbowPlot(DSS)
DSS <- FindNeighbors(DSS, dims = 1:20)
DSS <- FindClusters(DSS, resolution = 0.3)
DSS <- RunUMAP(DSS, dims = 1:20, return.model=TRUE)
DimPlot(DSS, group.by = "orig.ident")

#####INTEGRATION#####
anchors <- FindTransferAnchors(
  reference = eosinophils_steadystate,
  query = DSS,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

DSS <- MapQuery(
  anchorset = anchors,
  query = DSS,
  reference = eosinophils_steadystate,
  reference.reduction = "pca", 
  reduction.model = "umap",
  refdata = eosinophils_steadystate$orig.ident
)     


eosinophils_steadystate$id <- 'steadystate'
DSS$id <- 'DSS'
DSSquery <- merge(eosinophils_steadystate, DSS)
DSSquery[["pca"]] <- merge(eosinophils_steadystate[["pca"]], DSS[["ref.pca"]])
DSSquery <- RunUMAP(DSSquery, reduction = 'pca', dims = 1:20)
DimPlot(DSSquery, group.by = 'id', shuffle = TRUE, cols = c("darkgrey", "darkgreen"), pt.size =0.8, order= c("DSS", "steadystate") )+
  ggsave("UMAP_DSS.pdf", width = 8, height = 5)
DimPlot(DSSquery, reduction = "umap", group.by = "orig.ident")
DimPlot(eosinophils_steadystate, group.by= "seurat_clusters")
DSSquery <- FindNeighbors(DSSquery, dims = 1:50)
DSSquery <- FindClusters(DSSquery, resolution = 0.3)
DimPlot(DSSquery, reduction = "umap", group.by = "seurat_clusters")


###SUBSET INTESTINAL
Idents(DSSquery) <- "seurat_cluster"
intestinal_DSS <- subset(DSSquery, idents=1)         
DimPlot(intestinal_DSS, reduction = "umap", group.by = "orig.ident")
intestinal_DSS <- NormalizeData(intestinal_DSS, normalization.method = "LogNormalize", scale.factor = 10000)
intestinal_DSS <- FindVariableFeatures(intestinal_DSS, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(intestinal_DSS)
intestinal_DSS <- ScaleData(intestinal_DSS, features = all.genes, vars.to.regress = "percent.mt")
intestinal_DSS <- RunPCA(intestinal_DSS, features = VariableFeatures(object = intestinal_DSS))
DimPlot(intestinal_DSS, reduction = "pca", group.by = "orig.ident")
ElbowPlot(intestinal_DSS)
intestinal_DSS <- FindNeighbors(intestinal_DSS, dims = 1:20)
intestinal_DSS <- FindClusters(intestinal_DSS, resolution = 0.3)
intestinal_DSS <- RunUMAP(intestinal_DSS, dims = 1:20, return.model=TRUE)
DimPlot(intestinal_DSS, group.by = "orig.ident")
DSS_markers <- FindMarkers(intestinal_DSS, ident.1 = "DSS", ident.2="colon")
View(DSS_markers)

###MERGE WITH BACTERIAL CHALLENGE######
Idents(intestinal_DSS) <- "orig.ident"
intestinal_DSS_colon <- subset(intestinal_DSS, ident= c("colon", "DSS"))
DimPlot(refquery, group.by = "orig.ident")
intestinal_bac <- subset(refquery, ident= "intestinal eosinophils")
Idents(intestinal_bac) <- "orig.ident"
intestinal_DSS_bac_colon <- subset(intestinal_bac, ident= c("colonCR"))
intestinal_DSS_bac <- merge(intestinal_DSS_bac_colon, intestinal_DSS_colon)
intestinal_DSS_bac <- NormalizeData(intestinal_DSS_bac, normalization.method = "LogNormalize", scale.factor = 10000)
intestinal_DSS_bac <- FindVariableFeatures(intestinal_DSS_bac, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(intestinal_DSS_bac)
intestinal_DSS_bac <- ScaleData(intestinal_DSS_bac, features = all.genes, vars.to.regress = "percent.mt")
intestinal_DSS_bac <- RunPCA(intestinal_DSS_bac, features = VariableFeatures(object = intestinal_DSS_bac))
DimPlot(intestinal_DSS_bac, reduction = "pca", group.by = "orig.ident")
ElbowPlot(intestinal_DSS_bac)
intestinal_DSS_bac <- FindNeighbors(intestinal_DSS_bac, dims = 1:20)
intestinal_DSS_bac <- FindClusters(intestinal_DSS_bac, resolution = 0.3)
intestinal_DSS_bac <- RunUMAP(intestinal_DSS_bac, dims = 1:20, return.model=TRUE)
DimPlot(intestinal_DSS_bac, group.by = "orig.ident", cols=col_vector)
Idents(intestinal_DSS_bac) <- "orig.ident"
markers_intestinal_DSS_CR <- FindMarkers(intestinal_DSS_bac, ident.1 = "DSS", ident.2 = "colonCR")
View(markers_intestinal_DSS_CR)
#add module score antigen presentation
VlnPlot(intestinal_DSS_bac, features="H2-K1", pt.size = 0, cols = col_vector, log=F) + ylim(1,5)+ stat_summary(fun =mean, geom="point", color="black")+stat_summary(fun=sd)

markers_DSS_CR <- FindMarkers(intestinal_DSS_bac, ident.1 = "DSS", ident.2 = "colonCR")
markers_DSS <- FindMarkers(intestinal_DSS_bac, ident.1 = "DSS", ident.2 = "colon")

intestinal_DSS_bac$orig.ident <- factor(x = intestinal_DSS_bac$orig.ident, levels = c("colon", "colonCR", "DSS"))
Idents(intestinal_DSS_bac) <- "orig.ident"
RidgePlot(intestinal_DSS_bac, features = c("Cd274", "Cd80", "Cd9"),cols= c("grey",   "darkred", "darkgreen"))+labs(y = "")+
  ggsave("test2.pdf", width = 10, height = 4)

####SIGNATURES#####
#granules
Granules_synthesis_list <-   list(c("Prg2", "Prg3",  "Epx", "Ear6", "Ear1", "Ear2"))
intestinal_DSS_bac <-AddModuleScore(intestinal_DSS_bac, features= Granules_synthesis_list,name = "Granules")
names(x = intestinal_DSS_bac[[]])

VlnPlot(intestinal_DSS_bac, features="Granules1"), group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(Y = " Granule protein expression", title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Granules_violin.pdf", width = 8, height = 6)

#antigen processing
Antigen_processing <- list(c("H2-D1", "H2-Q7", "H2-K1", "H2-T23","H2-Q4","Tap1","Tapbp","B2m", "Psmb8", "Psme1",
                             "Psmb9", "Calr", "Psmb10", "Ncf1", "Fcer1g"))
intestinal_DSS_bac <-AddModuleScore(intestinal_DSS_bac, features= Antigen_processing, name = "Antigen_processing")
VlnPlot(intestinal_DSS_bac, features="Antigen_processing1", group.by = "orig.ident", cols= c("grey",   "darkred", "darkgreen"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(y = "Antigen processing signature", title = "", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Antigen_processing_violin.pdf", width = 6, height = 6)

Idents(intestinal_DSS_bac) <- "orig.ident"
intestinal_DSS <- subset(intestinal_DSS_bac, idents = c("DSS"))
intestinal_DSS_ss<- subset(intestinal_DSS_bac, idents = c("colon"))
intestinal_DSS_citro <- subset(intestinal_DSS_bac, idents = c("colonCR"))
wilcox.test(intestinal_DSS_ss$Antigen_processing1, intestinal_DSS$Antigen_processing1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(intestinal_DSS_ss$Antigen_processing1, intestinal_DSS_citro$Antigen_processing1, alternative = "two.sided") #p-value < 2.2e-16


#ifng
IFNg_regulated <- list(c("Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
intestinal_DSS_bac <-AddModuleScore(intestinal_DSS_bac, features= IFNg_regulated,name = "IFNg_regulated")
names(x = intestinal_DSS_bac[[]])

VlnPlot(intestinal_DSS_bac, features="IFNg_regulated1", group.by = "orig.ident", cols= c("grey",   "darkred", "darkgreen"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "IFNg-regulated antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Ifng_Antimicrobial_violin.pdf", width = 6, height = 6)

wilcox.test(intestinal_DSS_ss$IFNg_regulated1, intestinal_DSS$IFNg_regulated1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(intestinal_DSS_ss$IFNg_regulated1, intestinal_DSS_citro$IFNg_regulated1, alternative = "two.sided") #p-value < 2.2e-16


IFNg_signature <- list(c("Ccl9", "Cxcl9", "Cxcl10", "Cd274", "Gbp2", "Irf1", "H2-D1", "H2-Q4", "H2-Q7", "H2-K1", 
                         "Gbp7", "Stat1", "Irgm1", "B2m", "Irf9", "Ifitm3", "H2-T23", "Icam1", "Jak2", "Irf2", "Ifngr2"))

intestinal_DSS_bac <-AddModuleScore(intestinal_DSS_bac, features= IFNg_signature,name = "IFNg_signature")
names(x = intestinal_DSS_bac[[]])

VlnPlot(intestinal_DSS_bac, features="IFNg_signature1", group.by = "orig.ident", cols= c("grey",   "darkred", "darkgreen"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "IFNg score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("IFNg_signature_violin.pdf", width = 6, height = 6)

wilcox.test(intestinal_DSS_ss$IFNg_signature1, intestinal_DSS$IFNg_signature1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(intestinal_DSS_ss$IFNg_signature1, intestinal_DSS_citro$IFNg_signature1, alternative = "two.sided") #p-value < 2.2e-16


#antimicrobial peptides
antimicrobial <- list(c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                        "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
intestinal_DSS_bac <-AddModuleScore(intestinal_DSS_bac, features= antimicrobial,name = "antimicrobial")
names(x = intestinal_DSS_bac[[]])
VlnPlot(intestinal_DSS_bac, features="antimicrobial1", group.by = "orig.ident", cols= c("grey",   "darkred", "darkgreen"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Antimicrobial.pdf", width = 8, height = 6)

intestinal_DSS_bac$orig.ident <- factor(x = intestinal_DSS_bac$orig.ident , levels = c("colon", "colonCR","DSS"))
antimicrobial_data <- AverageExpression(intestinal_challenge, features = c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                                                                           "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
pheatmap(antimicrobial_data$RNA, scale = 'row', cluster_rows = T, cluster_cols = F, border_color="white", 
         fontsize_number = 0.5, cellwidth =30, cellheight = 20, treeheight_row=0, treeheight_col=0,	
         angle_col = "45",color = colorRampPalette(rev(brewer.pal(n = 5, name =  "RdYlBu")))(100), filename = "Antimicrobial_heatmap.pdf")





