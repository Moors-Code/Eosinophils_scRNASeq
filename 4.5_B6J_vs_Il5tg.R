#compare old and new B6
source("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Packages_functions.R")
colon_wt_1 <- data_to_sparse_matrix("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Sample_preprocessing/Expression_data/Eos_revision_ST_SampleTag06_mm_colon_wt_Expression_Data.st")
colon_wt_CR_1 <- data_to_sparse_matrix("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Sample_preprocessing/Expression_data/Eos_revision_ST_SampleTag07_mm_colon_wt_CR_Expression_Data.st")
colon_wt_1 <- CreateSeuratObject(colon_wt_1, project = "colon_wt_1")
colon_wt_CR_1 <- CreateSeuratObject(colon_wt_CR_1, project = "colon_wt_CR_1")

eosinophil_rev_B6_tworounds <- merge(colon_wt, c(colon_wt_CR, colon_wt_1, colon_wt_CR_1), add.cell.ids = c("colon_wt","colon_wt_CR", "colon_wt_1", "colon_wt_CR_1"))
eosinophil_rev_B6_tworounds[["percent.mt"]] <- PercentageFeatureSet(eosinophil_rev_B6_tworounds, pattern = "^mt.")
VlnPlot(eosinophil_rev_B6_tworounds, features = "nFeature_RNA")
VlnPlot(eosinophil_rev_B6_tworounds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
RidgePlot(eosinophil_rev_B6_tworounds, features="nFeature_RNA")
length(eosinophil_rev_B6_tworounds@active.ident)
plot1 <- FeatureScatter(eosinophil_rev_B6_tworounds, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eosinophil_rev_B6_tworounds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
eosinophil_rev_B6_tworounds <- subset(eosinophil_rev_B6_tworounds, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
length(eosinophil_rev_B6_tworounds@active.ident)
eosinophil_rev_B6_tworounds <- NormalizeData(eosinophil_rev_B6_tworounds, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophil_rev_B6_tworounds <- FindVariableFeatures(eosinophil_rev_B6_tworounds, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophil_rev_B6_tworounds)
eosinophil_rev_B6_tworounds <- ScaleData(eosinophil_rev_B6_tworounds, features = all.genes, vars.to.regress = "percent.mt")
eosinophil_rev_B6_tworounds <- RunPCA(eosinophil_rev_B6_tworounds, features = VariableFeatures(object = eosinophil_rev_B6_tworounds))
DimPlot(eosinophil_rev_B6_tworounds, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eosinophil_rev_B6_tworounds)
eosinophil_rev_B6_tworounds <- FindNeighbors(eosinophil_rev_B6_tworounds, dims = 1:20)
eosinophil_rev_B6_tworounds <- FindClusters(eosinophil_rev_B6_tworounds, resolution = 0.3)
eosinophil_rev_B6_tworounds <- RunUMAP(eosinophil_rev_B6_tworounds, dims = 1:20)
DimPlot(eosinophil_rev_B6_tworounds, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eosinophil_rev_B6_tworounds, order=T, group.by = "orig.ident", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eosinophil_rev_B6_tworounds, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
FeaturePlot(eosinophil_rev_B6_tworounds, features = c("Siglecf", "Il5ra", "Ccr3", "Epx"), cols=pal, pt.size = .2, order = T)
FeaturePlot(eosinophil_rev_B6_tworounds, features = c("Siglecf", "Il5ra", "Ccr3", "Cd80", "Cd274", "Thsb1"), cols=pal, pt.size = .2, order = T)

#####SUBSET EOSINOPHIL CLUSTER####
eosinophil_rev_B6_tworounds_markers <- FindAllMarkers(object = eosinophil_rev_B6_tworounds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eosinophil_rev_B6_tworounds_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))

eos_B6 <- subset(eosinophil_rev_B6_tworounds, idents=c("11"))
length(eos_B6@active.ident)
eos_B6 <- NormalizeData(eos_B6, normalization.method = "LogNormalize", scale.factor = 10000)
eos_B6 <- FindVariableFeatures(eos_B6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eos_B6)
eos_B6 <- ScaleData(eos_B6, features = all.genes, vars.to.regress = "percent.mt")
eos_B6 <- RunPCA(eos_B6, features = VariableFeatures(object = eos_B6), npcs = 50)
DimPlot(eos_B6, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eos_B6)
eos_B6 <- FindNeighbors(eos_B6, dims = 1:20)
eos_B6 <- FindClusters(eos_B6, resolution = 0.3)
eos_B6 <- RunUMAP(eos_B6, dims = 1:20,  return.model=TRUE)
DimPlot(eos_B6, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eos_B6, order=T, group.by = "orig.ident", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eos_B6, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
FeaturePlot(eos_B6, features = c("Siglecf", "Il5ra", "Ccr3", "Epx"), cols=pal, pt.size = .2, order = T)
FeaturePlot(eos_B6, features = c("Siglecf", "Il5ra", "Ccr3", "Cd80", "Cd274", "Thbs1"), cols=pal, pt.size = .2, order = T)

####INTEGRATION IN STEADYSTATE + CHALLENGE DATASET (refquery)####
DefaultAssay(refquery) <- "RNA" 
DefaultAssay(eos_B6) <- "RNA" 
refquery <- FindVariableFeatures(refquery, selection.method = "vst", nfeatures = 2000)

anchors <- FindTransferAnchors(
  reference = refquery,
  query = eos_B6,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

eos_B6 <- MapQuery(
  anchorset = anchors,
  query = eos_B6,
  reference = refquery,
  reference.reduction = "pca", 
  reduction.model = "umap",
  refdata = refquery$orig.ident
)     


refquery$id <- 'Il5tg'
eos_B6$id <- 'B6'
B6query <- merge(refquery, eos_B6)
B6query[["pca"]] <- merge(refquery[["pca"]], eos_B6[["ref.pca"]])
B6query <- RunUMAP(B6query, reduction = 'pca', dims = 1:50)
DimPlot(B6query, group.by = 'id', shuffle = TRUE, cols = c("darkgrey", "darkgreen"), pt.size =0.8, order= c("B6", "Il5tg") )+
  ggsave("UMAP_B6.pdf", width = 6, height = 5)
DimPlot(B6query, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")

DimPlot(B6query, reduction = "umap", group.by = c("orig.ident"), pt.size=1, 
        cols = c("bonemarrowCR"="grey" , "bonemarrow"="grey","blood"="grey", "bloodCR"="grey", "stomach"="grey", 
                 "colon"="black", "small intestine"="grey", "spleen"="grey", "colonCR"="black", 
                 "colon_wt"="darkorange", "colon_wt_1"="darkorange", "colon_wt_CR"="darkorange", "colon_wt_CR_1"="darkorange"),
        order = c("colon_wt", "colon_wt_1", "colon_wt_CR","colon_wt_CR_1", "colon", "colonCR",
                  "bonemarrow", "bonemarrowCR",  "blood", "bloodCR", "stomach", "spleen", "small intestine"))+ theme_void()+
  theme(legend.position = "none") + ggsave("UMAP_B6_new.pdf", width = 5, height = 4)

B6query <- FindNeighbors(B6query, dims = 1:20)
B6query <- FindClusters(B6query, resolution = 0.3)
B6query <- RunUMAP(B6query, dims = 1:20, return.model=TRUE)
DimPlot(B6query, reduction = "umap", group.by = "seurat_clusters", label=T)

#####CHECK EXPRESSION OF SIGNATURES####
B6_only <- subset(intestinal_Il5_and_B6, idents= c("colon_wt_CR", "colon_wt", "colon_wt_CR_1", "colon_wt_1"))
B6_only <- RenameIdents(B6_only, 'colon_wt_1' = 'colon_wt', 'colon_wt_CR_1' = 'colon_wt_CR')
Idents(B6_only) <- "orig.ident"

#antigen processing
Antigen_processing <- list(c("H2-D1", "H2-Q7", "H2-K1", "H2-T23","H2-Q4","Tap1","Tapbp","B2m", "Psmb8", "Psme1",
                             "Psmb9", "Calr", "Psmb10", "Ncf1", "Fcer1g"))
B6_only <-AddModuleScore(B6_only, features= Antigen_processing, name = "Antigen_processing")

VlnPlot(B6_only, features="Antigen_processing1", cols= c("grey", "darkorange"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(y = "Antigen processing signature", title = "", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Antigen_processing_violin_B6.pdf", width = 5, height = 6)
B6_ss <- subset(B6_only, idents=c("colon_wt"))
B6_CR <- subset(B6_only, idents=c("colon_wt_CR"))
wilcox.test(B6_ss$Antigen_processing1, B6_CR$Antigen_processing1, alternative = "two.sided") 

#IFNg signature
IFNg_signature <- list(c("Ccl9", "Cxcl9", "Cxcl10", "Cd274", "Gbp2", "Irf1", "H2-D1", "H2-Q4", "H2-Q7", "H2-K1", 
                         "Gbp7", "Stat1", "Irgm1", "B2m", "Irf9", "Ifitm3", "H2-T23", "Icam1", "Jak2", "Irf2", "Ifngr2"))

B6_only <-AddModuleScore(B6_only, features= IFNg_signature,name = "IFNg_signature")

VlnPlot(B6_only, features="IFNg_signature1",cols= c("grey", "darkorange"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "IFNg score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("IFNg_signature_violin_B6.pdf", width = 5, height = 6)
wilcox.test(B6_ss$IFNg_signature1, B6_CR$IFNg_signature1, alternative = "two.sided") 

#antimicrobial peptides
antimicrobial <- list(c( "S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                        "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
B6_only <-AddModuleScore(B6_only, features= antimicrobial,name = "antimicrobial")

VlnPlot(B6_only, features="antimicrobial1", cols= c("grey", "darkorange"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Antimicrobial_B6.pdf", width = 5, height = 6)
wilcox.test(B6_ss$antimicrobial1, B6_CR$antimicrobial1, alternative = "two.sided") 

