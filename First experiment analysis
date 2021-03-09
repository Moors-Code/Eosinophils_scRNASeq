cp /IMCR_shares/Moorlab/Coco/BD_data_Eosinophils/eos_prj_4organs_counts.tmp.RData /home/moorlab/mnt_rstudio/Coco/Eosinophils
load("eos_prj_4organs_counts.tmp.RData")

source("/home/mnt_rstudio/Coco/functions.R") #or Functions_Seuratv3.R ??  
library(devtools)
library(Seurat)
library(BiocManager)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
library(ggraph)
library(clustree)
library(harmony)
library(viridis)
library(sctransform)
library(wesanderson)
library(fgsea)
library(msigdbr)
library(ggraph)
library(ggrepel)
library(cowplot)
library(SeuratWrappers)
library(conos)
library(SeuratData)
library('psupertime')
library('SingleCellExperiment')
library(magrittr)
library(data.table)
library(stringr)


###CREATE SEURAT OBJECTS WITH COUNT MATRICES FROM EACH ORGAN#######
colon <- CreateSeuratObject(
  colon_counts,
  project = "colon",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

spleen <- CreateSeuratObject(
  spleen_counts,
  project = "spleen",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

smallintestine <- CreateSeuratObject(
  sint_counts,
  project = "small intestine",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

stomach <- CreateSeuratObject(
  stomach_counts,
  project = "stomach",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

#MERGE ALL ORGANS INTO ONE SEURAT OBJECT#######
eosinophil_merge <- merge(stomach, c(smallintestine, colon, spleen), add.cell.ids = c("stomach", "small intestine", "colon", "spleen"))
length(eosinophil_merge@active.ident)

#FILTERING FOR LOW QUALITY CELLS
eosinophil_merge[["percent.mt"]] <- PercentageFeatureSet(eosinophils, pattern = "^mt.")
VlnPlot(eosinophil_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
RidgePlot(eosinophil_merge, features="percent.mt", group.by = "orig.ident")
FeatureScatter(eosinophil_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
eosinophils_mito <- subset(eosinophil_merge, subset = nFeature_RNA > 100 & nFeature_RNA < 5000)
length(eosinophils_mito@active.ident)
summary(eosinophils_mito@meta.data)

#NORMALIZATION AND VARIABLE GENES
eosinophils_mito <- NormalizeData(eosinophils_mito, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophils_mito <- FindVariableFeatures(eosinophils_mito, selection.method = "vst", nfeatures = 2000)
top50 <- head(VariableFeatures(eosinophils_mito), 50)
plot1 <- VariableFeaturePlot(eosinophils_mito)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot2

#SCALING AND DIMENSIONALITY REDUCTION
all.genes <- rownames(eosinophils_mito)
eosinophils_mito <- ScaleData(eosinophils_mito, features = all.genes)
eosinophils_mito <- RunPCA(eosinophils_mito, features = VariableFeatures(object = eosinophils_mito))
DimPlot(eosinophils_mito, reduction = "pca")
ElbowPlot(eosinophils_mito)
eosinophils_mito <- FindNeighbors(eosinophils_mito, dims = 1:50)
eosinophils_mito <- FindClusters(eosinophils_mito, resolution = 0.4)
eosinophils_mito <- RunUMAP(eosinophils_mito, dims = 1:50)
DimPlot(eosinophils_mito, order=T, group.by = "orig.ident", pt.size = 0.1)
DimPlot(eosinophils_mito, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=T)
DimPlot(eosinophils_mito, reduction = "umap", pt.size = .1, split.by = "orig.ident", label=T)
FeaturePlot(eosinophils_mito, feature="percent.mt", order=T, pt.size = 1, split.by = "orig.ident")

#more clusters
eosinophils_mito_clusters <- FindClusters(eosinophils_pure, resolution = 0.8)
eosinophils_mito_clusters <- RunUMAP(eosinophils_mito_clusters, dims = 1:50)
DimPlot(eosinophils_mito_clusters, order=T, group.by = "seurat_clusters", pt.size = 0.1)

#CLUSTER MARKERS AND ANNOTATION#######
eosinophils_mito_markers <- FindAllMarkers(object = eosinophils_mito, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eosinophils_mito_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC))
top100 <-eosinophils_mito_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(top100,"/home/moorlab/mnt_rstudio/Coco/Eosinophils/top100clustermarkers.csv", row.names = FALSE)
#cp  /home/moorlab/mnt_rstudio/Coco/Eosinophils/top100clustermarkers.csv /IMCR_shares/Moorlab/Coco/BD_data_Eosinophils

#cluster 0 > eos, mostly stomach, siglecfe, siglecf
#cluster 1 > eos mostly colon > siglecF low > Vegf > resident? tissue repair
#cluster 2 > mito high eos
#cluster 3 > Md14, Ms4a7, Mrc1 > monocyte/macrophage
#cluster 4 > eosinophils, Epx, Ear2, Retnlg, Ear1, Prg2, Ear6, S100a6 > degranulating
#cluster 5 > eosinophils, express EPX, Ear1, mostly stomach, express "epithelial" genes Gkn, Tff, Car2
#cluster 6 > not eos, collagen production Ly6c > stromal cells 
#cluster 7 > B cells 
#cluster 8 > eosinophil subset > Lyz, defensins, MHCII- 
#cluster 9 > epithelial contaminant
#cluster 10 > epithelial contaminant
#cluster 11 > B cells 
#cluster 12 > eos subset
#cluster 13 > eos subset
#cluster 14 > dendritic cells

#general immune markers
markers.to.plot <- c("Nkg7","Il7r",  "Ccr7", "Cd8a",   #CD8, NK, DC
                      "Siglecf",  "Cd14",
                     "Epcam", "Cdh1", "Ptprc", #eos
                     "Lyz1",
                     "S100a4",   "Cd19", #B
                     "Col1a1",
                      "Epx",
                        "Ms4a7", #mono
                     "H2.Aa",
                      #NK
                     "Cst3", "Fcer1a", "Irf4", "Itgax", #DC > cluster 12
                     "Il5ra","Ccr3") #Platelet

Idents(eosinophils_mito) <- "seurat_clusters"


FeaturePlot(eosinophils_mito, features = "Tcrb")


#eosinophils markers
DotPlot(eosinophils_mito, features = markers.eosinophils, dot.scale = 8, split.by = ) + 
  RotatedAxis()

#cluster annotation
eosinophils_mito <- RenameIdents(eosinophils_mito, 
                                 `0` = "Eosinophils (stomach)", 
                                 `1` = "Eosinophils (colon)", 
                                 `2` = "Eosinophils (mito high)", 
                                 `3` = "Monocytes/macrophages", 
                                 `4` = "Eosinophils (granule-rich)", 
                                 `5` = "Eosinophils (epithelial-like)", 
                                 `6` = "Stromal cells", 
                                 `7` = "B cells", 
                                 `8` = "Eosinophils (defensins)", 
                                 `9` = "Epithelial contaminant", 
                                 `10` = "Epithelial contaminant 2",
                                 `11` = "B cells", 
                                 `12` = "Eosinophils", 
                                 `13` = "Eosinophils 2", 
                                 `14` = "DCs")

DimPlot(eosinophils_mito, label = TRUE, pt.size = 1)
plot <- DotPlot(eosinophils_mito, features = markers.to.plot, dot.scale = 8, split.by = )
plot + 
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust=1), axis.text.y = element_text(face="italic")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")


#SUBSETTING OF EOSINOPHIL ONLY CLUSTERS#######
eosinophils_pure <- subset(eosinophils_mito,  idents = c("Eosinophils (stomach)", "Eosinophils (colon)", "Eosinophils (mito high)", 
                                                         "Eosinophils (granule-rich)", "Eosinophils (epithelial-like)", "Eosinophils (defensins)",
                                                         "Eosinophils", "Eosinophils 2"))
DimPlot(eosinophils_pure)
eosinophils_pure <- NormalizeData(eosinophils_pure, normalization.method = "LogNormalize", scale.factor = 10000)
eosinophils_pure <- FindVariableFeatures(eosinophils_pure, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eosinophils_pure)
eosinophils_pure <- ScaleData(eosinophils_pure, features = all.genes)
eosinophils_pure <- RunPCA(eosinophils_pure, features = VariableFeatures(object = eosinophils_pure))
DimPlot(eosinophils_pure, reduction = "pca")
ElbowPlot(eosinophils_pure)
eosinophils_pure <- FindNeighbors(eosinophils_pure, dims = 1:20)
eosinophils_pure <- FindClusters(eosinophils_pure, resolution = 0.4)
eosinophils_pure <- RunUMAP(eosinophils_pure, dims = 1:20)
DimPlot(eosinophils_pure, order=T, group.by = "seurat_clusters", split.by = "orig.ident", label=T, pt.size = 0.1)
DimPlot(eosinophils_pure, reduction = "umap", pt.size = .5, split.by = "orig.ident", label=T)
DimPlot(eosinophils_pure, reduction = "umap", pt.size = .5, label=T)

#markers (unsupervised)
markers_eos_pure <- FindAllMarkers(object = eosinophils_pure, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_eos_pure %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC))
top100 <-markers_eos_pure %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
write.csv(top100,"/home/moorlab/mnt_rstudio/Coco/Eosinophils/top100eosonlyclustermarkers.csv", row.names = FALSE)
cp  /home/moorlab/mnt_rstudio/Coco/Eosinophils/top100eosonlyclustermarkers.csv /IMCR_shares/Moorlab/Coco/BD_data_Eosinophils

#candidate markers
DotPlot(eosinophils_pure, features = markers.eosinophils) + 
  RotatedAxis()

cytokines <- c("Tnf", "Il23a", "Csf2", "Csf3",  "Il33", "Il4", "Ccr3", "Ccl5", 
               "Ccl11", "Il1b", "Il12a", "Il6", "Tgfb1", "Il10")
receptors <- c("Ffar2", "Ffar3", "Hcar2", "Il4ra", "Tlr2", "Tlr4", "Ifngr1", 
               "Il10ra", "Tnfrsf1a", "Csf2ra", "Csf2rb", "Il1rl1", "Il5ra", "Cd69")
TFTs_enzymes <- c("Rela", "Irf5", "Irf4", "Aldh1a1", "Ido1", "Ptgs2", "Nos2", "Arg1", "Arg2")

surface.markers <- c("Cd24a", "Siglece","Siglecf", "Ccr3", "Clec12a", "Cd274", "Cd80", "Cd84", "Cxcr4", "Osm", "Vegfa", "Ahr", "Ccl4")

granules <- c("Epx", "Ear1","Ear2", "Ear10", "Ear14", "Elane",  "Prg2", "Gzmb", "Rnase2b", "Nos2")

DotPlot(eosinophils_pure, features = cytokines) + NoLegend()+RotatedAxis()

DotPlot(eosinophils_pure, features = receptors) + NoLegend()+  RotatedAxis()

DotPlot(eosinophils_pure, features = TFTs_enzymes) + NoLegend()+ RotatedAxis()

DotPlot(eosinophils_pure, features = surface.markers) + NoLegend()+ RotatedAxis()

DotPlot(eosinophils_pure, features = granules) + NoLegend()+ RotatedAxis()

DimPlot(eosinophils_pure)

#cluster 1 <- gastro specific, "regulatory eosinophils"
#cluster 4 <- granules high
#cluster 6 <- defensins and lysozyme, small intestine specific


#for trajectory analysis: what is more progenitor-like population/cluster?
FeaturePlot(eosinophils_pure, feature="Cd34", split.by = "orig.ident", order=T)
FeaturePlot(eosinophils_pure, feature="Lst1", split.by = "orig.ident", order=T)
FeaturePlot(eosinophils_pure, feature="Mki67", order=T)
FeaturePlot(eosinophils_pure, feature="Ly6e", split.by = "orig.ident", order=T)
FeaturePlot(eosinophils_pure, feature="Ccr3", split.by = "orig.ident", order=T)
FeaturePlot(eosinophils_pure, feature="Il5ra", order = T)
FeaturePlot(eosinophils_pure, feature="Ccnd2", order = T)


#######CLUSTER BREAKDOWN PER SAMPLE######
numberofcells         <- table(eosinophils_pure$orig.ident, eosinophils_pure$seurat_clusters)
totalcellsperorgan   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]), sum(numberofcells[3,]), sum(numberofcells[4,]))
a                     <- cbind(numberofcells,totalcellsperorgan)
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]), sum(a[,5]), sum(a[,6]), 
                           sum(a[,7]), sum(a[,8]), sum(a[,9]), sum(a[,10]))
b                     <- rbind(a, totalcellspercluster)

c0 <- (b[1:4,1]/totalcellsperorgan)*100
c1 <- (b[1:4,2]/totalcellsperorgan)*100
c2 <- (b[1:4,3]/totalcellsperorgan)*100
c3 <- (b[1:4,4]/totalcellsperorgan)*100
c4 <- (b[1:4,5]/totalcellsperorgan)*100
c5 <- (b[1:4,6]/totalcellsperorgan)*100
c6 <- (b[1:4,7]/totalcellsperorgan)*100
c7 <- (b[1:4,8]/totalcellsperorgan)*100
c8 <- (b[1:4,9]/totalcellsperorgan)*100
c9 <- (b[1:4,10]/totalcellsperorgan)*100

c <- rbind(c0,c1,c2,c3,c4,c5,c6,c7,c8)
colSums(c)

par(mar=c(2,2,2,2))

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=10)

barplot(c, horiz=TRUE,
        legend = F,
       # args.legend=list(bty = "n", x=225, cex=1.5),
        main = "Cluster breakdown per organ", 
        las = 1, 
        col=color_list,
        cex.axis=1.5, cex.names=1.5, cex.main=2)


####SAVE IMAGE####
save.image(file = "Eosinophils.RData")
cp /home/moorlab/mnt_rstudio/Coco/Eosinophils/Eosinophils.RData /IMCR_shares/Moorlab/Coco/
  


##GATING STRATEGY

#1. import list of surface proteins and add colnames
#cp /IMCR_shares/Moorlab/Coco/BD_data_Eosinophils/surface_markers.txt /home/moorlab/mnt_rstudio/Coco/Eosinophils
surface_markers <- read.delim("~/mnt_rstudio/Coco/Eosinophils/surface_markers.txt", header=FALSE, row.names=1, comment.char="#")head(surface_markers)
head(surface_markers)
colnames(surface_markers)[1] <- "Gene"
head(surface_markers)

#1. find markers of every cluster
cluster0_markers <- FindMarkers(eosinophils_pure, ident.1=0, verbose=F)
cluster1_markers <- FindMarkers(eosinophils_pure, ident.1=1, verbose=F)
cluster2_markers <- FindMarkers(eosinophils_pure, ident.1=2, verbose=F)
cluster3_markers <- FindMarkers(eosinophils_pure, ident.1=3, verbose=F)
cluster4_markers <- FindMarkers(eosinophils_pure, ident.1=4, verbose=F)
cluster5_markers <- FindMarkers(eosinophils_pure, ident.1=5, verbose=F)
cluster6_markers <- FindMarkers(eosinophils_pure, ident.1=6, verbose=F)


#3. subset clusterN_marker dataframe by only keeping rows that are in the surface_markers df
cluster0_surface_markers <- cluster0_markers[rownames(cluster0_markers) %in% surface_markers$Gene,]
cluster1_surface_markers <- cluster1_markers[rownames(cluster1_markers) %in% surface_markers$Gene,]
cluster2_surface_markers <- cluster2_markers[rownames(cluster2_markers) %in% surface_markers$Gene,]
cluster3_surface_markers <- cluster3_markers[rownames(cluster3_markers) %in% surface_markers$Gene,]
cluster4_surface_markers <- cluster4_markers[rownames(cluster4_markers) %in% surface_markers$Gene,]
cluster5_surface_markers <- cluster5_markers[rownames(cluster5_markers) %in% surface_markers$Gene,]
cluster6_surface_markers <- cluster6_markers[rownames(cluster6_markers) %in% surface_markers$Gene,]
View(cluster0_surface_markers)

FeaturePlot(eosinophils_pure, feature="Thbs1")

FACS.markers <- c("Ccr1","Cd24a","Cd274","Cd33","Cd44","Cd48","Cd53","Cd63","Cd74","Cd80",
                  "Cd9","Cdh1","Clec12a","Csf1r","Cxcr4","Epcam","Fas","Fcer1g","Fcgr2b",
                  "Fcgr3","Hmgb1","Icam1","Il1rl1","Itga2b","Itga4","Itgal","Itgax","Itgb1","Itgb2",
                  "Itgb3","Itgb7","L1cam","Lamp1","Ldlr","Notch1","Notch2","Pecam1","Sell",
                  "Siglece","Slc3a2","Spn","Treml2")

FACS.markers.ordered <- c("L1cam","Itgb7","Cd24a","Clec12a", "Il1rl1", "Cxcr4",
                          "Cd274","Cd80","Cd9","Fcgr3","Ldlr","Fas","Notch1","Notch2",
                          "Treml2",
                          "Epcam","Cdh1",
                          "Itgb3","Fcgr2b","Itga2b",
                          "Csf1r","Itgb1","Itgax","Itgal","Cd63", "Siglece","Lamp1","Spn","Pecam1", "Cd48",
                          "Cd74","Ccr1", "Fcer1g","Cd53","Cd33","Cd44","Sell","Slc3a2","Itgb2","Icam1",
                           "Hmgb1")

DoHeatmap(eosinophils_pure, features = FACS.markers, label=F, draw.lines	=T )+ 
  scale_fill_gradientn(colors = brewer.pal(11,"RdYlBu")[11:1])  +
  theme(axis.text.y = element_text(face = "bold.italic") )

plot <- DotPlot(eosinophils_pure, features = FACS.markers.ordered)
plot + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust=1), axis.text.y = element_text(face="italic")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")

FeaturePlot(eosinophils_pure, features= "Cd44")
DoHeatmap(eosinophils_pure, features = top100$gene, label=F, draw.lines	=T, lines.width =10)+ 
  scale_fill_gradientn(colors = brewer.pal(11,"RdYlBu")[11:1])  +
  theme(axis.text.y = element_text(face = "bold.italic") )



library(devtools)
library(Seurat)
library(BiocManager)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggraph)
library(monocle)

###1. CREATE CDS OBJECT
#1. extract data from Seurat object
data <- as(as.matrix(eosinophils_pure@assays$RNA@data), 'sparseMatrix') 
pd <- new('AnnotatedDataFrame', data = eosinophils_pure@meta.data) 
head(pd)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
fd <- new('AnnotatedDataFrame', data = fData) 
monocle_eos <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


#2. ESTIMATE SIZE FACTORS AND DISPERSION 
#size factors help us normalize for differences in mRNA recovered across cells, 
#"dispersion" values will help us perform differential expression analysis later.
monocle_eos <- estimateSizeFactors(monocle_eos)
monocle_eos <- estimateDispersions(monocle_eos)


#3. SET ORDERING FILTER
#use 2000 highly variable genes obtained from Seurat and used for clustering 
varianblefeatures <- eosinophils_pure@assays$RNA@var.features
monocle_eos <- setOrderingFilter(monocle_eos, ordering_genes = varianblefeatures)


#4. REDUCE DIMENSION
#reduce data dimensionality by Reversed Graph Embedding, use as many components as you used for clustering
monocle_eos <- reduceDimension(monocle_eos, 
                               max_components = 10,
                               method = 'DDRTree')


#5. ORDER CELLS ALONG TRAJECTORY
monocle_eos <- orderCells(monocle_eos)

trajectory_plot <- plot_cell_trajectory(monocle_eos,
                     show_tree=T,
                     show_branch_points = T, cell_size =0.5, color_by = "seurat_clusters")  
trajectory_plot 
trajectory_plot + facet_wrap(~orig.ident, nrow = 1)

plot_cell_trajectory(monocle_eos, color_by = "State") + facet_wrap(~orig.ident, nrow = 1)
plot_cell_trajectory(monocle_eos, color_by = "Pseudotime")  

plot_cell_trajectory(monocle_eos,show_branch_points = F, backbone_color= "red", markers = "Cd274", use_color_gradient=T, cell_size =0.5) 
plot_cell_trajectory(monocle_eos,show_branch_points = F, backbone_color= "red", markers = "Epx", use_color_gradient=T, cell_size =0.5)  
plot_cell_trajectory(monocle_eos,show_branch_points = F, backbone_color= "red", markers = "Il5ra", use_color_gradient=T, cell_size =0.5) 

plot_cell_trajectory(monocle_eos, show_tree=T,
                                        show_branch_points = T, cell_size =0.5, color_by = "seurat_clusters")  


#branch point analysis
BEAM_res1 <- BEAM(monocle_eos, branch_point = 5, cores = 4)
BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(monocle_eos[row.names(subset(BEAM_res1,
                            qval < 1e-45)),],
                            branch_point = 5,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)


regulatory_genes <- row.names(subset(fData(monocle_eos),
                               gene_short_name %in% c("Cd274", "Csf2ra", "Csf2rb")))
plot_genes_branched_pseudotime(monocle_eos[regulatory_genes,],
                               branch_point = 5,
                               color_by = "seurat_clusters",
                               ncol = 1)


#surface_proteins<- intercept(c(cluster markes), c(surface protein list))



###GSEA ON DIFFERENTIAL EXPRESSION ##
#select species and set
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
CGP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(CGP)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(KEGG)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
m_df<- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
TFT <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(TFT)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)
m_df<- msigdbr(species = "Mus musculus", category = "C6")
oncSig <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(oncSig)


Idents(eosinophils_pure) <- "orig.ident" #set original identity as active identity 
cluster4_markers <- FindMarkers(eosinophils_pure, ident.1 = 4, verbose = FALSE)
View(cluster4_markers)
cluster4_markers$p_val_adj[cluster4_markers$p_val_adj == 0] <- 8.637356e-302

ranks <- cluster4_markers %>% 
  na.omit()%>%
  mutate(ranking=-log10(p_val_adj)/sign(avg_logFC))
ranks <- ranks$ranking
names(ranks) <- rownames(cluster4_markers)
head(ranks, 10)

Hallmarks_cluster4 <- fgsea(pathways = Hallmarks, 
                            stats = ranks,
                            minSize=10,
                            maxSize=500,
                            nperm=1000000)
View(Hallmarks_cluster4)
View(Hallmarks_cluster4 %>% filter(abs(NES)>1 & padj<0.05))

BP_cluster4 <- fgsea(pathways = BP, 
                     stats = ranks,
                     minSize=10,
                     maxSize=500,
                     nperm=1000000)
View(BP_cluster4)

BP_cluster4$pathway<-gsub("GO_","",BP_cluster4$pathway)
BP_cluster4$pathway<-gsub("_"," ",BP_cluster4$pathway)

ggplot(BP_cluster4 %>% filter(abs(NES)>0.5 & padj<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=padj, width=3)) +
  scale_size_area(max_size = 10)+
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x=" ", y="Normalized Enrichment Score",
       title="GO: Biological processes enriched in cluster 4", cols="black") + 
  theme_classic()+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))


KEGG_cluster4 <- fgsea(pathways = KEGG, 
                       stats = ranks,
                       minSize=10,
                       maxSize=500,
                       nperm=1000000)
View(KEGG_cluster4 %>% filter(abs(NES)>1 & padj<0.05)) #REACTOME, for some reason



