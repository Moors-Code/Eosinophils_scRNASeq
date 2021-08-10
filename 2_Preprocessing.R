# the dgMatrix is a valid input to create the Seurat object
stomach <- data_to_sparse_matrix("/Expression_data/Stomach.st")
colon <- data_to_sparse_matrix("/Expression_data/Colon.st")
SI <- data_to_sparse_matrix("/Expression_data/SI.st")
Spleen <- data_to_sparse_matrix("/Expression_data/Spleen.st")
blood <- data_to_sparse_matrix("/Expression_data/blood.st")
bloodCR <- data_to_sparse_matrix("/Expression_data/bloodCR.st")
bonemarrow <- data_to_sparse_matrix("/Expression_data/bonemarrow.st")
bonemarrowCR <- data_to_sparse_matrix("/Expression_data/bonemarrowCR.st")
colonCR <- data_to_sparse_matrix("/Expression_data/colonCR.st")
stomachHP <- data_to_sparse_matrix("/Expression_data/stomachHP.st")

stomach <- CreateSeuratObject(stomach, project = "stomach")
colon <- CreateSeuratObject(colon, project = "colon")
SI <- CreateSeuratObject(SI, project = "SI")
Spleen <- CreateSeuratObject(Spleen, project = "Spleen")
blood <- CreateSeuratObject(blood, project = "blood")
bloodCR <- CreateSeuratObject(bloodCR, project = "bloodCR")
bonemarrow <- CreateSeuratObject(bonemarrow, project = "bonemarrow")
bonemarrowCR <- CreateSeuratObject(bonemarrowCR, project = "bonemarrowCR")
colonCR <- CreateSeuratObject(colonCR, project = "colonCR")
stomachHP <- CreateSeuratObject(stomachHP, project = "stomachHP")


####MERGE ALL SEURAT OBJECTS INTO ONE SAMPLE AND THEN REMOVE CONTAMINANTS####
eosinophil_allsamples <- merge(blood, c(bloodCR, 
                                        bonemarrow, bonemarrowCR, 
                                        colon, colonCR,
                                        lung,
                                        neonates,
                                        SI,
                                        Spleen,
                                        stomach,
                                        stomachHP), add.cell.ids = c("blood", "bloodCR", 
                    "bonemarrow", "bonemarrowCR", 
                    "colon", "colonCR",
                    "lung",
                    "neonates",
                    "smallintestine",
                    "spleen",
                    "stomach",
                    "stomachHP"))

##### QUALITY FILTERING #####
length(eosinophil_allsamples@active.ident)
eosinophil_allsamples[["percent.mt"]] <- PercentageFeatureSet(eosinophil_allsamples, pattern = "^mt.")
VlnPlot(eosinophil_allsamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
RidgePlot(eosinophil_allsamples, features="nFeature_RNA")
plot1 <- FeatureScatter(eosinophil_allsamples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eosinophil_allsamples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
eosinophil_allsamples <- subset(eosinophil_allsamples, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
length(eosinophil_allsamples@active.ident)

##### PRE-PROCESSING #####
eosinophil_allsamples <- NormalizeData(eosinophil_allsamples, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophil_allsamples <- FindVariableFeatures(eosinophil_allsamples, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophil_allsamples)
eosinophil_allsamples <- ScaleData(eosinophil_allsamples, features = all.genes, vars.to.regress = "percent.mt")
eosinophil_allsamples <- RunPCA(eosinophil_allsamples, features = VariableFeatures(object = eosinophil_allsamples))
DimPlot(eosinophil_allsamples, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eosinophil_allsamples)
eosinophil_allsamples <- FindNeighbors(eosinophil_allsamples, dims = 1:50)
eosinophil_allsamples <- FindClusters(eosinophil_allsamples, resolution = 0.3)
eosinophil_allsamples <- RunUMAP(eosinophil_allsamples, dims = 1:20)
DimPlot(eosinophil_allsamples, order=T, group.by = "orig.ident", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eosinophil_allsamples, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=T, cols = col_vector)+
  ggsave("Figures/allsampleUMAP.pdf", width = 6, height = 5) 
FeaturePlot(eosinophil_allsamples, features = c("Siglecf", "Il5ra", "Ccr3", "Epx"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_allsamples, features = c("Mki67"), cols=pal, pt.size = .5, order = T)

#####CLUSTER MARKERS AND ANNOTATION#######
eosinophil_allsamples_markers <- FindAllMarkers(object = eosinophil_allsamples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eosinophil_allsamples_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))

#general markers
markers.to.plot <- c("Ptprc","Ccr3","Siglecf" ,"Itgax", "Il5ra","Epx", "Fcer1a","Ms4a7", "H2-Ab1", "Cst3","S100a4", 
                      "Epcam", "Cdh1",  "Ccr7", "Cd19", "Irf4", "Col1a1","Cd14", "Nkg7","Il7r",  "Cd8a")
                    
plot <- DotPlot(eosinophil_allsamples, features = markers.to.plot, dot.scale = 8, cols = col_vector)
plot + 
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+
  ggsave("Figures/allsampleDotPlot.pdf", width = 7, height = 5)


######SUBSETTING OF EOSINOPHIL CLUSTERS#######
eosinophil_pure <- subset(eosinophil_allsamples,  idents = c(0, 1, 2, 3, 4, 8))
DimPlot(eosinophil_pure, reduction = "umap", pt.size = .5, label=T, cols = col_vector)
DimPlot(eosinophil_pure, reduction = "umap",group.by = "orig.ident", pt.size = .5, label=T, cols = col_vector)
length(eosinophil_pure@active.ident)
eosinophil_pure <- NormalizeData(eosinophil_pure, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophil_pure <- FindVariableFeatures(eosinophil_pure, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophil_pure)
eosinophil_pure <- ScaleData(eosinophil_pure, features = all.genes, vars.to.regress = "percent.mt")
eosinophil_pure <- RunPCA(eosinophil_pure, features = VariableFeatures(object = eosinophil_pure))
DimPlot(eosinophil_pure, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eosinophil_pure)
eosinophil_pure <- FindNeighbors(eosinophil_pure, dims = 1:20)
eosinophil_pure <- FindClusters(eosinophil_pure, resolution = 0.2)
eosinophil_pure <- RunUMAP(eosinophil_pure, dims = 1:20)
DimPlot(eosinophil_pure, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eosinophil_pure, order=T, group.by = "orig.ident", pt.size = 0.1, label=F, cols = col_vector, split.by = "orig.ident")
DimPlot(eosinophil_pure, order=T, group.by = "seurat_clusters", pt.size = 0.2, label=T, cols = col_vector)
DimPlot(eosinophil_pure, reduction = "umap", pt.size = .5, split.by = "orig.ident", label=T)

#save the eosinophil_pure dataset for downstream analysis
saveRDS(eosinophil_pure, file = "eosinophil_pure.rds")








