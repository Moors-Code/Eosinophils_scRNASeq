###TRAJECTORY ALIGNMENT####
library(dtw)
library(reshape2)

#####Create BBC and BBCCR Seurat objects and export data for Monocle (run in R4)######
Idents(refquery) <- "orig.ident"
refquery$subset <- refquery$seurat_clusters
BBC <- subset(refquery,  idents = c("bonemarrow", "blood", "colon"))
BBC <- NormalizeData(BBC, normalization.method = "LogNormalize", scale.factor = 10000)
BBC <- FindVariableFeatures(BBC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BBC)
BBC <- ScaleData(BBC, features = all.genes, vars.to.regress = "percent.mt")
BBC <- RunPCA(BBC, features = VariableFeatures(object = BBC))
DimPlot(BBC, reduction = "pca", group.by = "orig.ident")
ElbowPlot(BBC)
BBC <- FindNeighbors(BBC, dims = 1:20)
BBC <- FindClusters(BBC, resolution = 0.4)
BBC <- RunUMAP(BBC, dims = 1:20)
DimPlot(BBC, cols = col_vector, label = T)   
DimPlot(BBC, cols = col_vector, label = T, group.by = "seurat.clusters")   
data_BBC <- as(as.matrix(BBC@assays$RNA@data), 'sparseMatrix') 
save(data_BBC,file="/Trajectory_alignment/BBC_data.Rdata")
pd_BBC <- new('AnnotatedDataFrame', data = BBC@meta.data) 
head(pd)
save(pd,file="/Trajectory_alignment/BBC_pd.Rdata")
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
save(fData,file="/Trajectory_alignment/BBC_fData.Rdata")
fd <- new('AnnotatedDataFrame', data = fData) 
save(fd,file="/Trajectory_alignment/BBC_fd.Rdata")
variablefeatures <- BBC@assays$RNA@var.features
save(variablefeatures, file="/Trajectory_alignment/BBC_variablefeatures.Rdata")

BBCCR <- subset(refquery,  idents = c("bonemarrowCR", "bloodCR", "colonCR"))
BBCCR <- NormalizeData(BBCCR, normalization.method = "LogNormalize", scale.factor = 10000)
BBCCR <- FindVariableFeatures(BBCCR, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BBCCR)
BBCCR <- ScaleData(BBCCR, features = all.genes, vars.to.regress = "percent.mt")
BBCCR <- RunPCA(BBCCR, features = VariableFeatures(object = BBCCR))
DimPlot(BBCCR, reduction = "pca", group.by = "orig.ident")
ElbowPlot(BBCCR)
BBCCR <- FindNeighbors(BBCCR, dims = 1:20)
BBCCR <- FindClusters(BBCCR, resolution = 0.2)
BBCCR <- RunUMAP(BBCCR, dims = 1:20)
DimPlot(BBCCR, cols = col_vector, label = T)
DimPlot(BBCCR, cols = col_vector, label = T, group.by = "seurat_clusters")
FeaturePlot(BBCCR, features = "S100a9")
FeaturePlot(BBCCR, features = "Cd274")
data <- as(as.matrix(BBCCR@assays$RNA@data), 'sparseMatrix') 
save(data,file="/Trajectory_alignment/BBCCR_data.Rdata")
pd <- new('AnnotatedDataFrame', data = BBCCR@meta.data) 
head(pd)
save(pd,file="/Trajectory_alignment/BBCCR_pd.Rdata")
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
save(fData,file="/Trajectory_alignment/BBCCR_fData.Rdata")
fd <- new('AnnotatedDataFrame', data = fData) 
save(fd,file="/Trajectory_alignment/BBCCR_fd.Rdata")
variablefeatures <- BBCCR@assays$RNA@var.features
save(variablefeatures, file="/Trajectory_alignment/BBCCR_variablefeatures.Rdata")

######CREATE STEADYSTAE AND CHALLENGE TRAJECTORIES (run in R3.6)
library(devtools)
library(Seurat)
library(BiocManager)
library(dplyr)
library(ggplot2)
library(ggraph)
library(monocle)
library(viridis)

#####steadystate bonemarrow-blood-colon#####
monocle_BBC <- newCellDataSet(data,
                             phenoData = pd,
                             featureData = fd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())

monocle_BBC <- estimateSizeFactors(monocle_BBC)
monocle_BBC <- estimateDispersions(monocle_BBC)
monocle_BBC <- setOrderingFilter(monocle_BBC, ordering_genes = variablefeatures)
monocle_BBC <- reduceDimension(monocle_BBC, 
                              max_components = 10,
                              method = 'DDRTree')
monocle_BBC <- orderCells(monocle_BBC, reverse = T)
trajectory_plot <- plot_cell_trajectory(monocle_BBC,
                                        show_tree=T,
                                        show_branch_points = T, cell_size =0.5, color_by = "orig.ident")+
  ggsave("Figures/trajectoryBBC.pdf", width = 5, height = 5)


monocle_BBCR <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_BBCR <- estimateSizeFactors(monocle_BBCR)
monocle_BBCR <- estimateDispersions(monocle_BBCR)
monocle_BBCR <- setOrderingFilter(monocle_BBCR, ordering_genes = variablefeatures)
monocle_BBCR <- reduceDimension(monocle_BBCR, 
                               max_components = 10,
                               method = 'DDRTree')
monocle_BBCR <- orderCells(monocle_BBCR, reverse = T)
plot_cell_trajectory(monocle_BBCR, show_tree=F, show_branch_points = F, cell_size =1, color_by = "seurat_clusters")  + 
  theme(legend.position = "right") + scale_color_manual(values=rev(col_vector[1:5]))+ 
  ggsave("Figures/trajectoryBBCCR.pdf", width = 7, height =4)

####ALIGN TRAJECTORY use functions from https://github.com/cole-trapnell-lab/pseudospace/blob/master/code/Pseudospace_support_functions.R####
# Identify genes that are expressed in at least 50 of cells 
expressed_genes.list <- list()
expressed_genes.list[["steadystate"]] <- row.names(fData(monocle_BBC[Matrix::rowSums(Biobase::exprs(monocle_BBC) > 0) > 50 ,]))
length(expressed_genes.list[["steadystate"]])
expressed_genes.list[["challenge"]] <- row.names(fData(monocle_BBCR[Matrix::rowSums(Biobase::exprs(monocle_BBCR) > 0) > 50 ,]))
length(expressed_genes.list[["challenge"]])

expressed_genes <- unique(union(expressed_genes.list[["steadystate"]], expressed_genes.list[["challenge"]]))
length(expressed_genes)

# Use dynamic time warping to align Mock and TGFB pseudospatial trajectories and create a cds object of aligned trajectories
challgene.to.steadystate.aligned.cds <- getDTWcds(monocle_BBCR,monocle_BBC, 
                                      ref = "steadystate", query = "challenge", 
                                      expressed_genes = expressed_genes, cores = 1)
cds.aligned.list <- list()
cds.aligned.list[["steadystate to challenge"]] <- challgene.to.steadystate.aligned.cds

for(alignment in names(cds.aligned.list)){
  
  cds.aligned.list[[alignment]] <- preprocess_cds(cds.aligned.list[[alignment]])
  
}


# Identify genes that are differentially expressed across pseudotime as a function of treatment
aligned.pseudotime.DEG.test.list <- list()

for(alignment in names(cds.aligned.list)) {
  
  aligned.pseudotime.DEG.test.list[[alignment]] <- differentialGeneTest(cds.aligned.list[[alignment]][expressed_genes], 
                                                                         fullModelFormulaStr="~sm.ns(Pseudotime, df=3)*Cell.Type", 
                                                                         reducedModelFormulaStr="~sm.ns(Pseudotime, df=3) + Cell.Type", cores = 1)
  
}


diff_test_res <-aligned.pseudotime.DEG.test.list[[alignment]]
sig_gene_names <- row.names(subset(diff_test_res, qval  <= 1e-10))
View(diff_test_res[sig_gene_names,])

compare_cell_types_in_pseudotime(cds.aligned.list[["steadystate to challenge"]][row.names(subset(fData(cds.aligned.list[["steadystate to challenge"]]), 
                                                                                                  gene_short_name %in% c("Thbs1","Ly6c2"   ,"Retnla"   ,"Retnlg"))),], 
                                  color_by="orig.ident", df=3, min_expr=0.1, cell_alpha = 0.05, line_size = 1, ncol = 2,
                                  panel_order = c("Thbs1"   ,"Ly6c2"   ,"Retnla"   ,"Retnlg")) + 
  scale_color_manual(values = c("steadystate" = "grey", "challenge" = "red", "bonemarrow"="#5BBCD6","bonemarrowCR"="#5BBCD6",
                                "blood"="#F2AD00", "bloodCR"="#F2AD00", "colon"="#FF0000", "colonCR"="#FF0000"), name = "Treatment") +
  #theme(legend.position="none", text=element_text(size=20)) +
  ggsave("Figures/aligned_trajectories.pdf", width = 9, height = 5)

Pseudospatial.aligned.sig.genes.list <- list()

for(sample in names(aligned.pseudotime.DEG.test.list)){
  
  Pseudospatial.aligned.sig.genes.list[[sample]] <- row.names(subset(aligned.pseudotime.DEG.test.list[[sample]], 
                                                                     qval <= 1e-10))
  print(sample)
  print(length(Pseudospatial.aligned.sig.genes.list[[sample]]))
}

genes_of_interest <- Pseudospatial.aligned.sig.genes.list[[sample]]



#plot selected genes
compare_cell_types_in_pseudotime(cds.aligned.list[["steadystate to challenge"]][row.names(subset(fData(cds.aligned.list[["steadystate to challenge"]]), 
                                                                                                  gene_short_name %in% c( "Ahr", "Camp", "B2m", "Cd47", "Chil3", "Ear1", "Ear2", "Ear6"))),], 
                                  color_by="orig.ident", df=3, min_expr=0.1, cell_alpha = 0.05, line_size = 1, ncol = 2,
                                  panel_order = c( "Ahr", "Camp", "B2m", "Cd47", "Chil3", "Ear1", "Ear2", "Ear6")) + 
  scale_color_manual(values = c("steadystate" = "grey86", "challenge" = "darkred", "bonemarrow"="#5BBCD6","bonemarrowCR"="#5BBCD6",
                                "blood"="#F2AD00", "bloodCR"="#F2AD00", "colon"="#FF0000", "colonCR"="#FF0000"), name = "Treatment") +
  #theme(legend.position="none", text=element_text(size=20)) +
  ggsave("_aligned_pseudotime_final_1.pdf", width = 9, height = 7)

compare_cell_types_in_pseudotime(cds.aligned.list[["steadystate to challenge"]][row.names(subset(fData(cds.aligned.list[["steadystate to challenge"]]), 
                                                                                                  gene_short_name %in% c(  "Epx", "Retnla", "Gbp7", "Ifngr1", "Irf9", "Lcn2", "Lyz2", "Ltf" ))),], 
                                  color_by="orig.ident", df=3, min_expr=0.1, cell_alpha = 0.05, line_size = 1, ncol = 2,
                                  panel_order = c(  "Epx", "Retnla", "Gbp7", "Ifngr1", "Irf9", "Lcn2", "Lyz2", "Ltf")) + 
  scale_color_manual(values = c("steadystate" = "grey86", "challenge" = "darkred", "bonemarrow"="#5BBCD6","bonemarrowCR"="#5BBCD6",
                                "blood"="#F2AD00", "bloodCR"="#F2AD00", "colon"="#FF0000", "colonCR"="#FF0000"), name = "Treatment") +
  #theme(legend.position="none", text=element_text(size=20)) +
  ggsave("_aligned_pseudotime_final_2.pdf", width = 9, height = 7)

compare_cell_types_in_pseudotime(cds.aligned.list[["steadystate to challenge"]][row.names(subset(fData(cds.aligned.list[["steadystate to challenge"]]), 
                                                                                                  gene_short_name %in% c(  "Prg2", "Prg3", "S100a8", "S100a9", "Xbp1", "Gata1"))),], 
                                  color_by="orig.ident", df=3, min_expr=0.1, cell_alpha = 0.05, line_size = 1, ncol = 2,
                                  panel_order = c(  "Prg2", "Prg3", "S100a8", "S100a9", "Xbp1", "Gata1")) + 
  scale_color_manual(values = c("steadystate" = "grey86", "challenge" = "darkred", "bonemarrow"="#5BBCD6","bonemarrowCR"="#5BBCD6",
                                "blood"="#F2AD00", "bloodCR"="#F2AD00", "colon"="#FF0000", "colonCR"="#FF0000"), name = "Treatment") +
  #theme(legend.position="none", text=element_text(size=20)) +
  ggsave("_aligned_pseudotime_final_3.pdf", width = 9, height = 5.25)

