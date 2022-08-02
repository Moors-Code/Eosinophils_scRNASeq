source("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Packages_functions.R")
lung <- data_to_sparse_matrix("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Sample_preprocessing/Expression_data/Eos_revision_ST_SampleTag01_mm_lung_Expression_Data.st")
adipose <- data_to_sparse_matrix("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Sample_preprocessing/Expression_data/Eos_revision_ST_SampleTag03_mm_adipose_Expression_Data.st")
uterus <- data_to_sparse_matrix("/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Sample_preprocessing/Expression_data/Eos_revision_ST_SampleTag05_mm_uterus_Expression_Data.st")

lung <- CreateSeuratObject(lung, project = "lung")
adipose <- CreateSeuratObject(adipose, project = "adipose")
uterus <- CreateSeuratObject(uterus, project = "uterus")

####MERGE SEURAT OBJECTS####
eosinophil_revision_allsamples <- merge(lung, c(adipose,uterus), add.cell.ids = c("lung","adipose","uterus"))
length(eosinophil_revision_allsamples@active.ident)

eosinophil_revision_allsamples[["percent.mt"]] <- PercentageFeatureSet(eosinophil_revision_allsamples, pattern = "^mt.")
VlnPlot(eosinophil_revision_allsamples, features = "nFeature_RNA")

VlnPlot(eosinophil_revision_allsamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
RidgePlot(eosinophil_revision_allsamples, features="nFeature_RNA")
length(eosinophil_revision_allsamples@active.ident)
plot1 <- FeatureScatter(eosinophil_revision_allsamples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eosinophil_revision_allsamples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
eosinophil_revision_allsamples <- subset(eosinophil_revision_allsamples, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
length(eosinophil_revision_allsamples@active.ident)
eosinophil_revision_allsamples <- NormalizeData(eosinophil_revision_allsamples, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophil_revision_allsamples <- FindVariableFeatures(eosinophil_revision_allsamples, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophil_revision_allsamples)
eosinophil_revision_allsamples <- ScaleData(eosinophil_revision_allsamples, features = all.genes, vars.to.regress = "percent.mt")
eosinophil_revision_allsamples <- RunPCA(eosinophil_revision_allsamples, features = VariableFeatures(object = eosinophil_revision_allsamples))
DimPlot(eosinophil_revision_allsamples, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eosinophil_revision_allsamples)
eosinophil_revision_allsamples <- FindNeighbors(eosinophil_revision_allsamples, dims = 1:50)
eosinophil_revision_allsamples <- FindClusters(eosinophil_revision_allsamples, resolution = 0.3)
eosinophil_revision_allsamples <- RunUMAP(eosinophil_revision_allsamples, dims = 1:20)
DimPlot(eosinophil_revision_allsamples, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector, split.by = "orig.ident")
DimPlot(eosinophil_revision_allsamples, order=T, group.by = "orig.ident", pt.size = 0.1, label=F, cols = col_vector)
DimPlot(eosinophil_revision_allsamples, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=F, cols = col_vector)
FeaturePlot(eosinophil_revision_allsamples, features = c("Siglecf", "Il5ra", "Ccr3", "Epx"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Mki67"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Siglecf"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Cd80"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Cd274"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Siglecf"), cols=pal, pt.size = .5, order = T)
FeaturePlot(eosinophil_revision_allsamples, features = c("Cd14"), cols=pal, pt.size = .5, order = T)
saveRDS(eosinophil_revision_allsamples, file = "/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Revision/eosinophil_revision_allsamples.rds")

#####CLUSTER MARKERS AND ANNOTATION#######
eosinophil_allsamples_revision_markers <- FindAllMarkers(object = eosinophil_revision_allsamples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eosinophil_allsamples_revision_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))

#general markers
markers.to.plot <- c("Ptprc","Ccr3","Siglecf" ,"Itgax", "Il5ra","Epx", "Fcer1a","Ms4a7", "H2-Ab1", "Cst3","S100a4", 
                     "Epcam", "Cdh1",  "Ccr7", "Cd19", "Irf4", "Col1a1","Cd14", "Nkg7","Il7r",  "Cd8a")
plot <- DotPlot(eosinophil_revision_allsamples, features = markers.to.plot, dot.scale = 8, cols = col_vector)
plot + 
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")#+
#  ggsave("allsamplerevisionDotPlot.pdf", width = 7, height = 5)

#####MERGING#####
eosinophil_allsamples$id <- 'submission'
eosinophil_revision_allsamples$id <- 'revision'
eosinophil_everything <- merge(eosinophil_allsamples, c(eosinophil_revision_allsamples))
length(eosinophil_everything@active.ident)
eosinophil_everything <- NormalizeData(eosinophil_everything, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophil_everything <- FindVariableFeatures(eosinophil_everything, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophil_everything)
eosinophil_everything <- ScaleData(eosinophil_everything, features = all.genes, vars.to.regress = "percent.mt")
eosinophil_everything <- RunPCA(eosinophil_everything, features = VariableFeatures(object = eosinophil_everything))
DimPlot(eosinophil_everything, reduction = "pca", group.by = "id")
ElbowPlot(eosinophil_everything)
eosinophil_everything <- FindNeighbors(eosinophil_everything, dims = 1:50)
eosinophil_everything <- FindClusters(eosinophil_everything, resolution = 0.3)
eosinophil_everything <- RunUMAP(eosinophil_everything, dims = 1:20)
DimPlot(eosinophil_everything, order=T, group.by = "orig.ident", pt.size = 0.1, label=T, cols = col_vector)
DimPlot(eosinophil_everything, order=T, group.by = "id", pt.size = 0.1, label=T, cols = col_vector)
saveRDS(eosinophil_everything, file = "/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Revision/eosinophil_everything.rds")

#####BATCH CORRECTION#####
# install.packages("devtools")
#devtools::install_github("immunogenomics/harmony")
library(harmony)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = eosinophil_everything, reduction = "pca", pt.size = .1, group.by = "id")
p2 <- VlnPlot(object = eosinophil_everything, features = "PC_1", group.by = "id", pt.size = .1)
plot_grid(p1,p2)

harmonized_eosinophil_everything <- RunHarmony(eosinophil_everything, "id", plot_convergence = TRUE)
harmonized_eosinophil_everything <- RunUMAP(harmonized_eosinophil_everything, reduction = "harmony")

harmony_embeddings <- Embeddings(harmonized_eosinophil_everything, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = harmonized_eosinophil_everything, reduction = "harmony", pt.size = .1, group.by = "id")
p2 <- VlnPlot(object = harmonized_eosinophil_everything, features = "harmony_1", group.by = "id", pt.size = .1)
plot_grid(p1,p2)

harmonized_eosinophil_everything <- harmonized_eosinophil_everything %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
a <- DimPlot(eosinophil_everything, reduction = "umap", group.by = "orig.ident", pt.size = .1)
b<-DimPlot(harmonized_eosinophil_everything, reduction = "umap", group.by = "orig.ident", pt.size = .1)
plot_grid(a,b, ncol=1)
DimPlot(harmonized_eosinophil_everything, reduction = "umap", group.by = "id", pt.size = .1)
DimPlot(harmonized_eosinophil_everything, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'id', label=T)
FeaturePlot(harmonized_eosinophil_everything, features="Siglecf", reduction = "umap", pt.size = .1, split.by = 'id', label=T)

harmonized_eosinophil_everything_markers <- FindAllMarkers(object = harmonized_eosinophil_everything, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(harmonized_eosinophil_everything_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
#0 eosinophils, 1 eosinophils, 2 epithelial, 3 mito, 4 eosinophils, 5 eosinophils? 6 neutrophils, 
#7 mhc II, 8 mhc ii, 9 eosinophils, 10 mesenchymal, 11 b cells, 12 goblet cells, 13 mesenchymal, 14 mhcii, 15 t cells
FeaturePlot(harmonized_eosinophil_everything, feature="Ly6g")
FeaturePlot(harmonized_eosinophil_everything, feature="Siglecf", label=T)
FeaturePlot(harmonized_eosinophil_everything, feature="Il5ra", order=T, label=T)
FeaturePlot(harmonized_eosinophil_everything, feature="Ccr3")
FeaturePlot(harmonized_eosinophil_everything, feature="Cd274")
FeaturePlot(harmonized_eosinophil_everything, feature="Cd80")
Idents(harmonized_eosinophil_everything) <- "seurat_clusters"
saveRDS(eosinophil_everything, file = "/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Revision/harmonized_eosinophil_everything.rds")

#check lung for batch
lung <- subset(all_organs_steadystate, idents = "lung")
DimPlot(lung, group.by = "id")
DimPlot(lung, group.by = "id", reduction="pca")

#####REMOVE EPITHELIAL GENES FROM ADDITIONAL ORGANS#####
all_organs_steadystate <- subset(harmonized_eosinophil_pure_no6, idents=c("bloodCR", "bonemarrowCR", "colonCR", "stomachHP"),  invert=T)

#Find out the index of the genes you want removed: check markers of epithelial cluster
DimPlot(harmonized_eosinophil_everything)
View(harmonized_eosinophil_everything_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))

#find unique markers of cluster 2 (epithelial) that are expressed in less than 5% of cells in other clusters (pct.2=0.05)
positions_un <-which(harmonized_eosinophil_everything_markers$cluster==2 & harmonized_eosinophil_everything_markers$pct.2 < 0.05)
epithelial_markers <- harmonized_eosinophil_everything_markers$gene[positions_un]
DotPlot(eosinophil_everything, features= unique(epithelial_markers))+theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1))

#extract counts of all_organ steadystate
counts <- GetAssayData(all_organs_steadystate, assay = "RNA")
#subset counts
counts <- counts[-(which(rownames(counts) %in% epithelial_markers)),]

#subset seurat object and process
all_organs_steadystate_clean <- subset(all_organs_steadystate, features = rownames(counts))
all_organs_steadystate_clean <- NormalizeData(all_organs_steadystate_clean, normalization.method = "LogNormalize", scale.factor = 10000)
all_organs_steadystate_clean <- FindVariableFeatures(all_organs_steadystate_clean, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all_organs_steadystate_clean)
all_organs_steadystate_clean <- ScaleData(all_organs_steadystate_clean, features = all.genes, vars.to.regress = "percent.mt")
all_organs_steadystate_clean <- RunPCA(all_organs_steadystate_clean, features = VariableFeatures(object = all_organs_steadystate_clean))
DimPlot(all_organs_steadystate_clean, reduction = "pca", group.by = "orig.ident")
ElbowPlot(all_organs_steadystate_clean)
all_organs_steadystate_clean <- FindNeighbors(all_organs_steadystate_clean, dims = 1:15)
all_organs_steadystate_clean <- FindClusters(all_organs_steadystate_clean, resolution = 0.55)
all_organs_steadystate_clean <- RunUMAP(all_organs_steadystate_clean, dims = 1:15, return.model=T )
DimPlot(all_organs_steadystate_clean, order=T, group.by = "seurat_clusters", label=T)
DimPlot(all_organs_steadystate_clean, order=T, group.by = "seurat_clusters", split.by="orig.ident")
DotPlot(all_organs_steadystate_clean, features=c("Cd80", "Cd274"), group.by = "seurat_clusters")
Idents(all_organs_steadystate_clean) <- "orig.ident"
table(all_organs_steadystate_clean@active.ident)

#cluster markers
Idents(all_organs_steadystate_clean) <- "seurat_clusters"
all_organs_steadystate_clean_markers <- FindAllMarkers(all_organs_steadystate_clean)
View(all_organs_steadystate_clean_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))

#remove clusters 7 (MHC2 pos) and 11 (epithelial?)
all_organs_steadystate_pure <- subset(all_organs_steadystate_clean, idents=c("7", "11"), invert=T) 
DimPlot(all_organs_steadystate_pure, cols = col_vector)
all_organs_steadystate_pure <- NormalizeData(all_organs_steadystate_pure, normalization.method = "LogNormalize", scale.factor = 10000)
all_organs_steadystate_pure <- FindVariableFeatures(all_organs_steadystate_pure, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all_organs_steadystate_pure)
all_organs_steadystate_pure <- ScaleData(all_organs_steadystate_pure, features = all.genes, vars.to.regress = "percent.mt")
all_organs_steadystate_pure <- RunPCA(all_organs_steadystate_pure, features = VariableFeatures(object = all_organs_steadystate_pure))
DimPlot(all_organs_steadystate_pure, reduction = "pca", group.by = "orig.ident")
ElbowPlot(all_organs_steadystate_pure)
all_organs_steadystate_pure <- FindNeighbors(all_organs_steadystate_pure, dims = 1:20)
all_organs_steadystate_pure <- FindClusters(all_organs_steadystate_pure, resolution = 0.6)
all_organs_steadystate_pure <- RunUMAP(all_organs_steadystate_pure, dims = 1:20, return.model=T )
DimPlot(all_organs_steadystate_pure, order=T, group.by = "seurat_clusters", label=T)
DimPlot(all_organs_steadystate_pure, order=T, group.by = "seurat_clusters", split.by="orig.ident")
DotPlot(all_organs_steadystate_pure, features=c("Cd80", "Cd274"), group.by = "seurat_clusters")
Idents(all_organs_steadystate_pure) <- "orig.ident"
table(all_organs_steadystate_pure@active.ident)

Idents(all_organs_steadystate_pure) <- "seurat_clusters"
all_organs_steadystate_pure_markers <- FindAllMarkers(all_organs_steadystate_pure)
View(all_organs_steadystate_pure_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))
cluster_markers <- all_organs_steadystate_pure_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(cluster_markers,"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Revision/additional_organs_cluster_markers.csv", row.names = TRUE)

#rename clusters
current.cluster.ids <- c(0, 1, 2, 3, 4,5,6,7,8,9)
new.cluster.ids <-  c("basal eosinophils", "active eosinophils",   "circulating eosinophils", "immature eosinophils", "basal eosinophils", "tissue eosinophils",
                      "circulating eosinophils", "eosinophil precursors", "basal eosinophils", "circulating eosinophils")
all_organs_steadystate_pure@meta.data$seurat_clusters <- plyr::mapvalues(x = all_organs_steadystate_pure@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(all_organs_steadystate_pure) <- "seurat_clusters"
all_organs_steadystate_pure$seurat_clusters <- factor(x = all_organs_steadystate_pure$seurat_clusters, levels = rev(c("tissue eosinophils", "active eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil precursors")))
DimPlot(all_organs_steadystate_pure, reduction = "umap", pt.size = .5, label=F, cols = c( "#5BBCD6", "#F98400","#F2AD00" ,"#00A08A" ,  "#FF0000" ,"#F1BB7B"   )) + 
  ggsave("UMAP_additional_organs.pdf", width = 8, height = 5)

DimPlot(all_organs_steadystate_pure, reduction = "umap", group.by = c("orig.ident"), pt.size=.5, label=F,
           cols = c( "bone marrow"="grey", "blood"="grey",  "stomach"="grey", "colon"="grey", "small intestine"="grey", "spleen"="grey",
                     "uterus" ="#BF5B17", "lung" = "#386CB0", "adipose tissue" = "#1B9E77"),
           order = rev(c( "bonemarrow", "blood", "stomach","colon", "small intestine", "spleen", "adipose tissue", "lung", "uterus")))+ 
  ggsave("UMAP_additional_organs_perorgan.pdf",width = 8, height = 5)

DimPlot(all_organs_steadystate_pure, reduction = "umap", pt.size = .5, label=F, cols = c( "#5BBCD6", "#F98400","#F2AD00" ,"#00A08A" ,  "#FF0000" ,"#F1BB7B"   ), split.by="orig.ident") 
tissue_eos <-  subset(all_organs_steadystate_pure, idents="tissue eosinophils")

##find out if small cluster of active eos was tissue eos
eosinophils_steadystate_more_clustered <- eosinophils_steadystate
eosinophils_steadystate_more_clustered <- FindClusters(eosinophils_steadystate_more_clustered, resolution = 0.5)
eosinophils_steadystate_more_clustered <- RunUMAP(eosinophils_steadystate_more_clustered, dims = 1:20, return.model=T )
DimPlot(eosinophils_steadystate_more_clustered, order=T, group.by = "seurat_clusters", label=T)
cluster7 <- subset(eosinophils_steadystate_more_clustered, idents="7")
colnames(cluster7)
intersect(colnames(cluster7), colnames(tissue_eos))

#organ markers
Idents(all_organs_steadystate_pure) <- "orig.ident"
View(all_organs_steadystate_pure_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))
markers <- all_organs_steadystate_pure_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(markers,"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Revision/additional_organs_cluster_markers.csv", row.names = TRUE)

#progenitors check
progenitor_panel <- AverageExpression(all_organs_steadystate_pure, features=c("Mki67" ,"Il5ra", "Siglecf", "Cd34", "Itgam", "Ly6a"))
pheatmap(progenitor_panel$RNA, scale = "row", cluster_cols = F)  

#active vs tissue eosinophils
active_vs_tissue <- FindMarkers(all_organs_steadystate_pure, ident.1="active eosinophils", ident.2="tissue eosinophils")
View(active_vs_tissue)
active_vs_tissue$Gene <- rownames(active_vs_tissue)
DEGs<-active_vs_tissue
pdf("active_vs_tissue.pdf",  width = 4, height = 4)
DEGs_volcano(DEGs, 0.05, 1, "active vs tissue eosinophils", "black", 65, 2)
dev.off() 

####COMPOSITIONAL ANALYSIS####
#frequencies per cluster
numberofcells         <- table(all_organs_steadystate_pure$orig.ident, all_organs_steadystate_pure$seurat_clusters)
numberofcells
numberofcells <- numberofcells[c(1,2,4,6:11),]
totalcellsperorgan   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]), sum(numberofcells[3,]), sum(numberofcells[4,]),
                          sum(numberofcells[5,]), sum(numberofcells[6,]), sum(numberofcells[7,]), sum(numberofcells[8,]),
                            sum(numberofcells[9,]))
                          
a                     <- cbind(numberofcells,totalcellsperorgan)
a
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]), sum(a[,5]), sum(a[,6]), sum(a[,7]))
b                     <- rbind(a, totalcellspercluster)
b

c0 <- (b[1:9,1]/totalcellsperorgan)*100
c1 <- (b[1:9,2]/totalcellsperorgan)*100
c2 <- (b[1:9,3]/totalcellsperorgan)*100
c3 <- (b[1:9,4]/totalcellsperorgan)*100
c4 <- (b[1:9,5]/totalcellsperorgan)*100
c5 <- (b[1:9,6]/totalcellsperorgan)*100

c <- rbind(c0,c1,c2,c3,c4,c5)
colSums(c)
rownames(c) =  rev(c("tissue eosinophils", "active eosinophils", "basal eosinophils", "circulating eosinophils", "immature eosinophils", "eosinophil precursors"))
c

#plot
#par(mar=c(5,8,2,14))
pdf(file="Clusterbreakdown.pdf")
barplot(c, horiz=TRUE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=180, cex=.8),
        main = "Cluster breakdown per organ", 
        las = 1, 
        col= c( "#5BBCD6", "#F98400","#F2AD00" ,"#00A08A" ,  "#FF0000" ,"#F1BB7B"   ))
dev.off()
