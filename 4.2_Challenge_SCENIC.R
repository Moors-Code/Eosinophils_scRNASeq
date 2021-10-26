#####SCENIC######
source('Packages_functions.R', echo=TRUE)
library(SCENIC)

#get data
DimPlot(refquery, group.by = "seurat_clusters")
Idents(refquery) <- "seurat_clusters"
active <- subset(refquery, ident= "active eosinophils")
Idents(active) <- "orig.ident"
active <- subset(active, ident= c("colon", "stomach", "stomachHP",
                                          "colonIF", "colonCR", "colonCRIF"))
exprMat <- as.matrix(active@assays$RNA@counts)
cellInfo <- active@meta.data

dir.create("SCENIC")
setwd("SCENIC") 

### Initialize settings
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "orig.ident"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

Idents(active) <- "orig.ident"
DimPlot(active, cols = rev(col_vector[1:6]))
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(orig.ident=c("colon"="#FF0000",
                            "colonCR"="#F98400" , 
                           "colonIF"="#F1BB7B", 
                           "colonCRIF"= "#F2AD00" ,
                            "stomach"="#00A08A", 
                           "stomachHP"="#5BBCD6" ))
colVars$orig.ident <- colVars$orig.ident[intersect(names(colVars$orig.ident), cellInfo$orig.ident)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$orig.ident, legend=names(colVars$orig.ident))

org="mgi" # or hgnc, or dmel
dbDir="/SCENIC/" # RcisTarget databases location
myDatasetTitle="SCENIC Challenge" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"
scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions@settings$db_mcVersion <- "v8"
saveRDS(scenicOptions, file="scenicOptions.Rds") 


#Gene filter: Filter by the number of cells in which the gene is detected (minCountsPerGene, by default 6 UMI counts across all samples)
# and by the number of cells in which the gene is detected (minSamples, by default 1% of the cells)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
interestingGenes <- c("Siglecf", "Cd274", "Csfr2b", "Ear1")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)] #no Cd274!!! Run again with lower threshold!

#calculating correlation
runCorrelation(exprMat_filtered, scenicOptions)

#GENIE3: infer potential transcription factor targets based on the expression data
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123


exprMat_log <- log2(exprMat+1)
dim(exprMat)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

#This part I ran in science cloud.
#Clustering / dim reduction on the regulon activity.
scenicOptions <- readRDS("int/scenicOptions.Rds")
nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them

#Binarizing the network
runSCENIC_4_aucell_binarize(scenicOptions)

# # Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
 
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="orig.ident", cex=.5)
logMat <- log2(exprMat+1)

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)
colnames(aucell_regulonAUC.t) <- gsub(" ", "_", colnames(aucell_regulonAUC.t))
rownames(aucell_regulonAUC.t) <- gsub("[.]", "-", rownames(aucell_regulonAUC.t))

active.scenic <- active
active.scenic@meta.data <- cbind(active.scenic@meta.data, aucell_regulonAUC.t[rownames(active.scenic@meta.data),])
pdf("challenge_regulon_activity.pdf")
DimPlot(active.scenic, reduction = "umap")
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[35:43], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[9:12], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[13:16], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[17:20], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[21:24], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[25:28], cols = c("yellow", "red"))
FeaturePlot(active.scenic, reduction = "umap", features = colnames(active.scenic@meta.data)[29:32], cols = c("yellow", "red"))

#Regulon scores heatmap
cells.ord.cluster <- active.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)

library(gplots)
heatmap.2(regulon.scores)
library(pheatmap)
pheatmap(data_subset_norm, cluster_rows = F,cluster_cols = F, annotation_col = NA, show_colnames = F)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(regulon.scores, 1, cal_z_score))
pheatmap(data_subset_norm)

cluster.col <- data.frame(active.scenic@active.ident, row.names = names(active.scenic@active.ident))
colnames(cluster.col) <- "Seurat_cluster"
pheatmap(data_subset_norm, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)

pheatmap(regulon.scores.scaled, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F, fontsize_row=5)




#Binary regulon activity
library(AUCell)
binary.regulon.activity <- Binarize_regulon_activity(scenicOptions, skipBoxplot = FALSE, skipHeatmaps = FALSE, 
                                                   skipTsne = FALSE, exprMat = NULL)
cells.ord.cluster <- active.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
binary.regulon.activity <- binaryRegulonActivity
colnames(binary.regulon.activity) <- gsub("[.]", "-", colnames(binary.regulon.activity))
binary.regulon.activity <- binary.regulon.activity[,names(cells.ord.cluster)]
active.binary.regulon.activity <- binary.regulon.activity
save(active.binary.regulon.activity, file= "active_scenic_binaryRegulonActivity_table.RData")
binary.regulon.activity <- binary.regulon.activity[which(rowSums(binary.regulon.activity) > 80),]

cluster.col <- data.frame(active.scenic@active.ident, row.names = names(active.scenic@active.ident))
colnames(cluster.col) <- "Seurat_cluster"
pdf("Figures/active_scenic_binaryRegulonActivity.pdf")
pheatmap(binary.regulon.activity, cluster_rows = T,cluster_cols = T, annotation_col = cluster.col, show_colnames = F,  fontsize_row=7)
dev.off()

active.scenic <- active
active.scenic@meta.data <- cbind(active.scenic@meta.data, t(binary.regulon.activity)[rownames(active.scenic@meta.data),])
library(Seurat)
pdf("Figures/active_regulon_activity_on_Umap.pdf")
DimPlot(active.scenic, reduction = "umap")
FeaturePlot(active.scenic, reduction = "umap", features = "Nfkb1_(24g)", cols = c("grey", "red"), split.by="orig.ident")
FeaturePlot(active.scenic, reduction = "umap", features = c("Gata2 (20g)"), cols = c("grey", "red"))
dev.off()
FeaturePlot(active.scenic,features = "Stat1_extended_(462g)", split.by = "orig.ident")
DotPlot(active.scenic,features = "Ltf_top50")

#heatmap of TF expressions
TFs <- rownames(regulon.scores.scaled)
TFs <- as.character(lapply(TFs, function(x) unlist(strsplit(x, "_"))[1]))

DoHeatmap(active.scenic, features = TFs, size = 3)
challenge.exp.mat <- active.scenic@assays$RNA@data
TFs.exp.mat <- challenge.exp.mat[TFs,names(cells.ord.cluster)]
TFs.exp.mat <- as.numeric(TFs.exp.mat)
pheatmap(TFs.exp.mat, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)


#Average Regulon Activity 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$orig.ident),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                        cluster_rows=T, cluster_columns = T, column_names_rot = 30,
                        column_names_gp = gpar(fontsize = 11), row_names_gp = gpar(fontsize = 8))

#Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- 0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$orig.ident), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

#plot selected data
data <- regulonActivity_byCellType_Scaled[,c(1,2,5,6)]
regulons <- rownames(data)
View(regulons)
positions<- which(regulons ==  "Gata2_extended (110g)"|
                  regulons ==  "Xbp1 (40g)"|
                  regulons ==  "Gata1 (32g)"|
                  regulons ==  "Gata1 (32g)"|
                  regulons ==    "Foxo1 (17g)"  |
                  regulons ==   "Klf2_extended (89g)"|
                  regulons ==  "Klf3_extended (21g)" |
                  regulons == "Klf7 (65g)" |
                  regulons == "Klf10_extended (11g)"|
                  regulons == "Klf13 (14g)" |
                  regulons == "Rel (184g)"|
                  regulons == "Rela (208g)"| 
                  regulons == "Relb (68g)"|
                  regulons == "Nfkb1 (24g)"| 
                  regulons == "Nfkb2 (111g)"|
                  regulons == "Ahr (15g)" |
                  regulons == "Egr1 (89g)"|
                  regulons ==  "Hmgb1 (11g)"|
                  regulons == "Irf1 (46g)" |
                  regulons ==  "Irf2 (76g)" | 
                  regulons ==  "Irf5 (17g)"|
                  regulons ==    "Irf7 (23g)"|
                  regulons ==   "Irf9 (50g)" |
                  regulons ==  "Rara_extended (21g)" |
                  regulons ==  "Rxra_extended (128g)"|
                  regulons ==   "Rxrb_extended (17g)" |
                  regulons == "Nfya (120g)" |
                  regulons ==   "Nfyb (95g)" | 
                  regulons ==  "Nfyc_extended (17g)"|
                  regulons ==  "Nfe2_extended (419g)"|
                  regulons ==  "Mafg (10g)"|
                  regulons ==  "Stat1 (130g)"|
                  regulons ==  "Stat3_extended (254g)"|
                  regulons ==  "Stat4_extended (21g)" |
                  regulons ==  "Stat5b_extended (93g)"|
                  regulons == "Stat6 (27g)")

pdf("regulonActivity_byCellType_Scaled.pdf", width = 5, height = 5)
ComplexHeatmap::Heatmap(data[positions,], name="Regulon activity",
                        cluster_columns = F, column_names_rot = 30, border_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize = 11), row_names_gp = gpar(fontsize = 9), 
                        row_title_gp = gpar(fill = col_vector[c(1,5,4,2,3)]))
dev.off()




#Cell-type specific regulators (based on the Regulon Specificity Score (RSS) proposed by Suo et al. for the Mouse Cell Atlas in 2018). 
#Useful for big analysis with many cell types, to identify the cell-type specific regulons.
# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
## Showing regulons and cell types with any RSS > 0.01 (dim: 6x5)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "basal eosinophils")
plotRSS_oneSet(rss, setName = "active eosinophils")
plotRSS_oneSet(rss, setName = "circulating eosinophils")
plotRSS_oneSet(rss, setName = "immature eosinophils")
plotRSS_oneSet(rss, setName = "eosinophil progenitors")


# load everything from SCENIC run
setwd("SCENIC")

GENIE3_linkList <- readRDS("int/1.4_GENIE3_linkList.Rds")
tfModules <- readRDS("int/1.6_tfModules_asDF.Rds")
tfModules.MotifEnrichmet <- readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds") # THIS contains top50 targets of each TF

motif_AUC <- readRDS("int/2.2_motifs_AUC.Rds")
regulonTargetInfo <- readRDS("int/2.5_regulonTargetsInfo.Rds")
regulons_asGeneSet <- readRDS("int/2.6_regulons_asGeneSet.Rds")
regulons_asIncidMat <- readRDS("int/2.6_regulons_asIncidMat.Rds")
regulons_forAUCell <- readRDS("int/3.1_regulons_forAUCell.Rds")

aucellRankings <- readRDS("int/3.3_aucellRankings.Rds")
regulonAUC <- readRDS("int/3.4_regulonAUC.Rds")


