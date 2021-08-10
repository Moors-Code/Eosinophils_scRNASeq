##GATING STRATEGY####
#1. import list of surface proteins and add colnames
surface_markers <- read.delim("surface_markers.txt", header=FALSE, row.names=1, comment.char="#")
head(surface_markers)
colnames(surface_markers)[1] <- "Gene"
head(surface_markers)

#3. subset clusterN_marker dataframe by only keeping rows that are in the surface_markers df
prog_surface_markers <- prog_markers[rownames(prog_markers) %in% surface_markers$Gene,]
prog_surface_markers$Gene <- rownames(prog_surface_markers)
immature_surface_markers <- immature_markers[rownames(immature_markers) %in% surface_markers$Gene,]
immature_surface_markers$Gene <- rownames(immature_surface_markers)
circulating_surface_markers <- circulating_markers[rownames(circulating_markers) %in% surface_markers$Gene,]
circulating_surface_markers$Gene <- rownames(circulating_surface_markers)
basal_surface_markers <- basal_markers[rownames(basal_markers) %in% surface_markers$Gene,]
basal_surface_markers$Gene <- rownames(basal_surface_markers)
active_surface_markers <- active_markers[rownames(active_markers) %in% surface_markers$Gene,]
active_surface_markers$Gene <- rownames(active_surface_markers)
View(active_surface_markers)

#merge and make heatmap
merged <- merge(prog_surface_markers[,c(2,6)], immature_surface_markers[,c(2,6)], by = "Gene" , all = T,
      suffixes = c(".progenitors",".immature"))
View(merged4)
merged2 <- merge(merged, circulating_surface_markers[,c(2,6)], by = "Gene" , all = T,
                suffixes = c(".progenitors",".immature"))
merged3 <- merge(merged2, basal_surface_markers[,c(2,6)], by = "Gene" , all = T, 
                 suffixes = c(".circulating",".basal"))
merged4 <- merge(merged3, active_surface_markers[,c(2,6)], by = "Gene" , all = T)
colnames(merged4) <- c("Gene","progenitors", "immature", "circulating", "basal", "active")
merged4[is.na(merged4)] <- 0
row.names(merged4) <- merged4$Gene
merged4[1] <- NULL

pheatmap(merged4, cluster_cols = F, border_color="black",fontsize_number = 0.5, cellwidth =30, cellheight = 10, scale = "row", filename="Figures/surfacemarkers_heatmap.pdf")

##CITEK
cytek.markers <- c(  "Clec12a", "Siglece", "Itga4",    "Cd274",   "Cd80",  "Ly6a",
                     "Cd9",  "Fas", "Icam1", 
                     "Itgax",  "Lamp1", "Fcgr3",
                     "Slc3a2", "Cd74", "Cd44", "Il1rl1", "Cd24a",
                     "Itgam", "Epcam", "Pecam1",  "Sell", "Cdh1",  "Ly6c1","Spn","Cd53", "Cd63", "Itga2b") 

plot <- DotPlot(eosinophils_steadystate, features = cytek.markers)
plot + theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+ labs(title = "Cytek markers", y = "", x="")+
  ggsave("Figures/cytekmarkers.png", width = 12, height = 3)


### SURFACER MARKERS ###
Idents(eosinophils_steadystate) <- "seurat_clusters"
basal_active <- subset(eosinophils_steadystate, idents = c("active eosinophils", "basal eosinophils"))
surface.markers <- c("Cd274", "Cd80","Cd9","Pecam1", "Icam1", "Fas", "Ly6a", "Itgax", "Itga4", "Clec12a", "Siglece")
DotPlot(basal_active, features = selected_markers,  dot.scale = 40)
plot + theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1,size=20), axis.text.y = element_text(face="bold", size=20)) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")+ labs(title = "FACS markers", y = "", x="")+
  ggsave("Figures/FACSmarkers.pdf", width = 12, height = 3)

data <- AverageExpression(basal_active, features =  surface.markers)
data
pheatmap(data$RNA[,c(2,1)], cluster_columns = F, scale = "row", 
         fontsize_number = 0.5, cellwidth =30, cellheight = 30, treeheight_row=0, treeheight_col=0,	
         border_color = "white", color = colorRampPalette(rev(brewer.pal(n = 4, name =  "RdYlBu")))(100), filename = "Figures/FACSmarkers_heatmap_h.pdf")


##blended featureplot
FeaturePlot(eosinophils_steadystate, features = c("Cd80", "Cd274"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T) +
  ggsave("Figures/markers_blend.pdf", width = 16, height = 5)

FeaturePlot(eosinophils_steadystate, features = c("Cd80", "Clec12a"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T) +
  ggsave("Figures/markers_blend2.pdf", width = 16, height = 5)



