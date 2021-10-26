##### SCENIC ###
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
library(AUCell)
## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")

options(timeout = max(300, getOption("timeout")))

for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

library(SCENIC)
exprMat <- as.matrix(eosinophils_steadystate@assays$RNA@counts)
cellInfo <- eosinophils_steadystate@meta.data

dir.create("SCENIC")
setwd("SCENIC") 

### Initialize settings
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "CellType"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

Idents(eosinophils_steadystate) <- "seurat_clusters"
DimPlot(eosinophils_steadystate, cols = rev(col_vector[1:5]))

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("eosinophil progenitors"="#5BBCD6", 
                           "immature eosinophils"="#F98400", 
                           "circulating eosinophils"= "#F2AD00" ,
                            "basal eosinophils"="#00A08A", 
                           "active eosinophils"= "#FF0000"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

org="mgi" # or hgnc, or dmel
dbDir="/SCENIC/" # RcisTarget databases location
myDatasetTitle="SCENIC Steady State Organs" # choose a name for your analysis
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
interestingGenes <- c("Siglecf", "Cd274", "Csfr2b")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)] 

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
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
logMat <- log2(exprMat+1)

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)
colnames(aucell_regulonAUC.t) <- gsub(" ", "_", colnames(aucell_regulonAUC.t))
rownames(aucell_regulonAUC.t) <- gsub("[.]", "-", rownames(aucell_regulonAUC.t))

steadystate.scenic <- eosinophils_steadystate
steadystate.scenic@meta.data <- cbind(steadystate.scenic@meta.data, aucell_regulonAUC.t[rownames(steadystate.scenic@meta.data),])
pdf("Figures/steadystate_regulon_activity.pdf")
DimPlot(steadystate.scenic, reduction = "umap")

#Regulon scores heatmap
cells.ord.cluster <- steadystate.scenic@active.ident
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

colnames(cluster.col) <- "Seurat_cluster"
pheatmap(data_subset_norm, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)
pheatmap(regulon.scores.scaled, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F, fontsize_row=5)

#Binary regulon activity
library(AUCell)
binary.regulon.activity <- Binarize_regulon_activity(scenicOptions, skipBoxplot = FALSE, skipHeatmaps = FALSE, 
                                                   skipTsne = FALSE, exprMat = NULL)
cells.ord.cluster <- steadystate.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
binary.regulon.activity <- binaryRegulonActivity
colnames(binary.regulon.activity) <- gsub("[.]", "-", colnames(binary.regulon.activity))
binary.regulon.activity <- binary.regulon.activity[,names(cells.ord.cluster)]
steadystate.binary.regulon.activity <- binary.regulon.activity
save(steadystate.binary.regulon.activity, file= "/scenic_binaryRegulonActivity_table.RData")
binary.regulon.activity <- binary.regulon.activity[which(rowSums(binary.regulon.activity) > 80),]

cluster.col <- data.frame(steadystate.scenic@active.ident, row.names = names(steadystate.scenic@active.ident))
colnames(cluster.col) <- "Seurat_cluster"
pdf("scenic_binaryRegulonActivity.pdf")
pheatmap(binary.regulon.activity, cluster_rows = T,cluster_cols = T, annotation_col = cluster.col, show_colnames = F,  fontsize_row=7)
dev.off()

steadystate.scenic@meta.data <- cbind(steadystate.scenic@meta.data, t(binary.regulon.activity)[rownames(steadystate.scenic@meta.data),])
library(Seurat)
DimPlot(steadystate.scenic, reduction = "umap")
FeaturePlot(steadystate.scenic, reduction = "umap", features = c("E2f1 (1388g)","Nfe2 (109g)", "Mafg_extended (112g)", "Xbp1_extended (15g)", "Batf3_extended (12g)",
                                                                 "Meis1_extended (19g)",
                                                                 "Nfkb1 (107g)","Irf5 (57g)",  "Gata2 (20g)"), cols = c("grey", "brown3"), ncol=3, order=T)+ 
  ggsave("RegulonActivityUMAP.pdf", width = 17, height = 15)

FeaturePlot(steadystate.scenic, reduction = "umap", features = c("Gata1 (13g)", "Rara (208g)", "Egr1 (434g)", "Stat6 (43g)",
                                                                 "Stat3_extended (226g)", "Hif1a (150g)", "Ahr_extended (55g)"), cols = c("grey", "brown3"), ncol=3, order=T)+ 
  ggsave("RegulonActivityUMAP2.pdf", width = 17, height = 15)


FeaturePlot(steadystate.scenic, reduction = "umap", features = c("Nfkb1 (107g)",  "Rela (185g)", "Stat3_extended (226g)", "Irf5 (57g)", "Ahr_extended (55g)", "Hif1a (150g)"), cols = c("grey", "brown3"), ncol=3, order=T)+ 
  ggsave("RegulonActivityUMAP_main.pdf", width = 19, height = 12)

FeaturePlot(steadystate.scenic, reduction = "umap", features = c("E2f1 (1388g)",  "Mafg_extended (112g)","Nfe2 (109g)", "Xbp1_extended (15g)","Rara (208g)", "Meis1_extended (19g)"), cols = c("grey", "brown3"), ncol=3, order=T)+ 
  ggsave("RegulonActivityUMAP_sup.pdf", width = 19, height = 12)


FeaturePlot(steadystate.scenic, reduction = "umap", features = c("Gata2 (20g)"), cols = c("grey", "darkred"))
dev.off()
DotPlot(steadystate.scenic,features = c("Stat1_extended_(481g)", "Stat1_(145g)"," Nfkb2_extended_(237g)", "Nfil3_(31g)"))

#heatmap of TF expressions
TFs <- rownames(regulon.scores.scaled)
TFs <- as.character(lapply(TFs, function(x) unlist(strsplit(x, "_"))[1]))

DoHeatmap(steadystate.scenic, features = TFs, size = 3)
steadystate.exp.mat <- steadystate.scenic@assays$RNA@data
TFs.exp.mat <- steadystate.exp.mat[TFs,names(cells.ord.cluster)]
TFs.exp.mat <- as.numeric(TFs.exp.mat)
pheatmap(TFs.exp.mat, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)

#Average Regulon Activity 
library(grid)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pdf("regulonActivity_byCellType_Scaled.pdf", width = 10, height = 20)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                         cluster_columns = F, column_names_rot = 30, border_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize = 11), row_names_gp = gpar(fontsize = 9), row_km = 5, row_gap = unit(3, "mm"),
                        row_title_gp = gpar(fill = col_vector[c(1,5,4,2,3)]))
dev.off()
#Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- .3
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"),
                        cluster_rows=T, cluster_columns = F, column_names_rot = 30,
                        column_names_gp = gpar(fontsize = 11), row_names_gp = gpar(fontsize = 9))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Rela" & highConfAnnot==TRUE]
write.csv(tableSubset, "Rela.csv")
viewMotifs(tableSubset, options=list(pageLength=5)) 

#Cell-type specific regulators (based on the Regulon Specificity Score (RSS) proposed by Suo et al. for the Mouse Cell Atlas in 2018). 
#Useful for big analysis with many cell types, to identify the cell-type specific regulons.
# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
## Showing regulons and cell types with any RSS > 0.01 (dim: 6x5)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "basal eosinophils")
plotRSS_oneSet(rss, setName = "intestinal eosinophils")
plotRSS_oneSet(rss, setName = "circulating eosinophils")
plotRSS_oneSet(rss, setName = "immature eosinophils")
plotRSS_oneSet(rss, setName = "eosinophil progenitors")


# load everything from SCENIC run
GENIE3_linkList <- readRDS("int/1.4_GENIE3_linkList.Rds")
tfModules <- readRDS("int/1.6_tfModules_asDF.Rds")
tfModules.MotifEnrichmet <- readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds") # THIS contains top50 targets of each TF
View(tfModules.MotifEnrichmet)

motif_AUC <- readRDS("int/2.2_motifs_AUC.Rds")
regulonTargetInfo <- readRDS("int/2.5_regulonTargetsInfo.Rds")
View(regulonTargetInfo)
regulons_asGeneSet <- readRDS("int/2.6_regulons_asGeneSet.Rds")
regulons_asIncidMat <- readRDS("int/2.6_regulons_asIncidMat.Rds")
regulons_forAUCell <- readRDS("int/3.1_regulons_forAUCell.Rds")

aucellRankings <- readRDS("int/3.3_aucellRankings.Rds")
regulonAUC <- readRDS("int/3.4_regulonAUC.Rds")



###### PROGENY ####
CellsClusters <- data.frame(Cell = names(Idents(eosinophils_steadystate)),
                            CellType = as.character(Idents(eosinophils_steadystate)),
                            stringsAsFactors = FALSE)
Idents(eosinophils_steadystate) <- "seurat_clusters"
eosinophils_steadystate <- progeny(eosinophils_steadystate, scale=FALSE, organism="Mouse", top=500, perm=1, return_assay = TRUE)
eosinophils_steadystate <- Seurat::ScaleData(eosinophils_steadystate, assay = "progeny")

progeny_scores_df <-
  as.data.frame(t(GetAssayData(eosinophils_steadystate, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(summarized_progeny_scores_df[,-1],fontsize=14,
                        fontsize_row = 14, cluster_rows=F, cluster_cols = T,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy: pathway activity analysis", angle_col = 45, cellwidth = 30, cellheight = 20,
                        treeheight_col = 0,  border_color = "black", filename = "Figures/progeny.pdf")

