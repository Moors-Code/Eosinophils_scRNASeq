library(devtools)
library(Seurat)
library(BiocManager)
library(dplyr)
library(ggplot2)
library(ggraph)
library(monocle)
library(viridis)

###EXPORT MONOCLE INPUT FROM SEURAT (run in R4####
data <- as(as.matrix(eosinophils_steadystate@assays$RNA@data), 'sparseMatrix') 
save(data,file="/Monocle/data.Rdata")

pd <- new('AnnotatedDataFrame', data = eosinophils_steadystate@meta.data) 
head(pd)
save(pd,file="/Monocle/pd.Rdata")

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
save(fData,file="/Monocle/fData.Rdata")

fd <- new('AnnotatedDataFrame', data = fData) 
save(fd,file="/Monocle/fd.Rdata")

#use 2000 highly variable genes obtained from Seurat and used for clustering as ordering filter
variablefeatures <- eosinophils_steadystate@assays$RNA@var.features
save(variablefeatures, file="/Monocle/variablefeatures.Rdata")

#####MONOCLE (run in R studio 3.6.1)#####
monocle_ss <- newCellDataSet(data,
                             phenoData = pd,
                             featureData = fd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())

monocle_ss <- estimateSizeFactors(monocle_ss)
monocle_ss <- estimateDispersions(monocle_ss)

#use 2000 highly variable genes obtained from Seurat and used for clustering as ordering filter
monocle_ss <- setOrderingFilter(monocle_ss, ordering_genes = variablefeatures)
monocle_ss <- reduceDimension(monocle_ss, 
                              max_components = 10,
                              method = 'DDRTree')
monocle_ss <- orderCells(monocle_ss)
trajectory_plot <- plot_cell_trajectory(monocle_ss,
                                        show_tree=T,
                                        show_branch_points = T, cell_size =0.5, color_by = "seurat_clusters")  
trajectory_plot 
trajectory_plot + facet_wrap(~orig.ident, nrow = 1)

plot_cell_trajectory(monocle_ss, color_by = "State") 
plot_cell_trajectory(monocle_ss, color_by = "Pseudotime")  

plot_cell_trajectory(monocle_ss,show_branch_points = F, backbone_color= "red", markers = "Mki67", cell_size =0.5) 
plot_cell_trajectory(monocle_ss,show_branch_points = F, backbone_color= "red", markers = "Cd274", use_color_gradient=T, cell_size =0.5) 
plot_cell_trajectory(monocle_ss,show_branch_points = F, backbone_color= "red", markers = "Epx", use_color_gradient=T, cell_size =0.5)  
plot_cell_trajectory(monocle_ss,show_branch_points = F, backbone_color= "red", markers = "Il5ra", use_color_gradient=T, cell_size =0.5) 
plot_cell_trajectory(monocle_ss,show_branch_points = F, backbone_color= "red", markers = "Cd34", use_color_gradient=T, cell_size =0.5) 

plot_cell_trajectory(monocle_ss, show_tree=F, show_branch_points = F, cell_size =1, color_by = "seurat_clusters")  + 
  theme(legend.position = "right") + scale_color_manual(values=rev(col_vector[1:5]))+ 
  ggsave("Figures/Trajectory.png", width = 10, height =5)

timevarying_genes1 <- row.names(subset(fData(monocle_ss), gene_short_name %in% c( "Mki67", "Epx", "Ear1", "S100a6" )))
cds_subset <- monocle_ss[timevarying_genes1, ]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters", cell_size = 1, ncol=1, min_expr=0.5)  + xlim(0,18) + ylim(0,5) +
  scale_color_manual(values=rev(col_vector[1:5])) + theme(legend.position = "none") + 
  theme(axis.line.y = element_line(size=1, color="black")) + 
  theme(axis.line.x = element_line(size=1, color="black")) +
  theme(strip.text = element_text(size=10, face="italic")) + 
  theme(axis.title = element_text(size=11,face="bold")) + 
  theme(axis.text.x = element_text(size=9, color="black")) + 
  theme(axis.text.y = element_text(size=9, color="black")) + 
  theme(axis.line.y = element_line(size=1, color="black")) + 
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())+ 
  ggsave("Figures/Jitter1.pdf", width = 5, height =8)

timevarying_genes2 <- row.names(subset(fData(monocle_ss), gene_short_name %in% c( "Cd274", "Cd80", "Cd9")))
cds_subset <- monocle_ss[timevarying_genes2, ]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters", cell_size = 1, ncol=1, min_expr=0.5)  + xlim(0,18) + ylim(0,5) +
  scale_color_manual(values=rev(col_vector[1:5])) + theme(legend.position = "none") + 
  theme(axis.line.y = element_line(size=1, color="black")) + 
  theme(axis.line.x = element_line(size=1, color="black")) +
  theme(strip.text = element_text(size=10, face="italic")) + 
  theme(axis.title = element_text(size=11,face="bold")) + 
  theme(axis.text.x = element_text(size=9, color="black")) + 
  theme(axis.text.y = element_text(size=9, color="black")) + 
  theme(axis.line.y = element_line(size=1, color="black")) + 
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())+ 
  ggsave("Figures/Jitter2.pdf", width = 5, height =7.5)


diff_test_res_ss <- differentialGeneTest(monocle_ss[variablefeatures,],
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res_ss, qval < 0.01))
View(diff_test_res_ss[sig_gene_names,])
top100_gene_names <- row.names(subset(diff_test_res_ss, qval<4.655098e-103))

plot_pseudotime_heatmap(monocle_ss[sig_gene_names,], num_clusters = 3, hmcols = magma(300), 
                        show_rownames = T, return_heatmap=T, cores=5) + 
                        ggsave("Figures/Heatmap1.png", width = 5, height =10) 

plot_pseudotime_heatmap(monocle_ss[top100_gene_names,], num_clusters = 3, hmcols = magma(300), 
                        show_rownames = T, return_heatmap=T, cores=5) + 
  ggsave("Figures/Heatmap1.png", width = 5, height =10) 


