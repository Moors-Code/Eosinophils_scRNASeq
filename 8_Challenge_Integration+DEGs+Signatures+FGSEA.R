source('Packages_functions.R', echo=TRUE)

###INTEGRATION####
Idents(eosinophil_pure) <- "seurat_clusters"
DimPlot(eosinophil_pure)
eosinophil_pure <- subset(eosinophil_pure, idents=c(0,1,2,3,5))
Idents(eosinophil_pure) <- "orig.ident"
DimPlot(eosinophil_pure)
eos_challenge <- subset(eosinophil_pure, idents=c("colonCR", "bloodCR", "bonemarrowCR", "stomachHP"))
DimPlot(eos_challenge, group.by = "orig.ident")
eos_challenge <- NormalizeData(eos_challenge, normalization.method = "LogNormalize", scale.factor = 10000)
eos_challenge <- FindVariableFeatures(eos_challenge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(eos_challenge)
eos_challenge <- ScaleData(eos_challenge, features = all.genes, vars.to.regress = "percent.mt")
eos_challenge <- RunPCA(eos_challenge, features = VariableFeatures(object = eos_challenge))
DimPlot(eos_challenge, reduction = "pca", group.by = "orig.ident")
ElbowPlot(eos_challenge)
eos_challenge <- FindNeighbors(eos_challenge, dims = 1:20)
eos_challenge <- FindClusters(eos_challenge, resolution = 0.3)
eos_challenge <- RunUMAP(eos_challenge, dims = 1:20, return.model=TRUE)
DimPlot(eos_challenge, group.by = "orig.ident")


anchors <- FindTransferAnchors(
  reference = eosinophils_steadystate,
  query = eos_challenge,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

eos_challenge <- MapQuery(
  anchorset = anchors,
  query = eos_challenge,
  reference = eosinophils_steadystate,
  reference.reduction = "pca", 
  reduction.model = "umap",
  refdata = eosinophils_steadystate$orig.ident
  )

DimPlot(eos_challenge, reduction = "ref.umap", group.by = "orig.ident")
DimPlot(eos_challenge, reduction = "umap", group.by = "orig.ident")
DimPlot(eosinophils_steadystate, reduction = "umap", group.by = "orig.ident")
DimPlot(eosinophil_pure, reduction = "umap", group.by = "orig.ident")


eos_challenge <- FindNeighbors(eos_challenge, dims = 1:50)
eosinophils_steadystate$id <- 'steadystate'
eos_challenge$id <- 'bacterial challenge'
refquery <- merge(eosinophils_steadystate, eos_challenge)
refquery[["pca"]] <- merge(eosinophils_steadystate[["pca"]], eos_challenge[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:20)
DimPlot(refquery, group.by = 'id',  pt.size =.8, cols= c("darkred", "darkgray"))+
  ggsave("Figures/UMAP_challenge.pdf", width = 8, height = 5)
DimPlot(refquery, reduction = "umap", group.by = "orig.ident")
refquery <- FindNeighbors(refquery, dims = 1:50)
refquery <- FindClusters(refquery, resolution = 2)
refquery <- RunUMAP(refquery, dims = 1:20, return.model=TRUE)
DimPlot(refquery, reduction = "umap", group.by = "seurat_clusters", label=T)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16)
new.cluster.ids <-  c("basal eosinophils", 
                      "basal eosinophils", 
                      "active eosinophils", 
                      "active eosinophils", 
                      "circulating eosinophils",
                      "immature eosinophils", 
                      "immature eosinophils", 
                      "basal eosinophils", 
                      "basal eosinophils", 
                      "circulating eosinophils",
                      "eosinophil progenitors",
                      "basal eosinophils",
                      "immature eosinophils",
                      "basal eosinophils",
                      "eosinophil progenitors",
                      "immature eosinophils",
                      "eosinophil progenitors")
refquery@meta.data$seurat_clusters <- plyr::mapvalues(x = refquery@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(refquery) <- "seurat_clusters"
refquery$seurat_clusters <- factor(x = refquery$seurat_clusters, levels = rev(c("active eosinophils", "basal eosinophils","circulating eosinophils","immature eosinophils",  "eosinophil progenitors")))
DimPlot(refquery, reduction = "umap", pt.size = .5, label=F, cols = rev(col_vector[1:5]))
markers <- FindAllMarkers(refquery, min.pct = 0.25, logfc.threshold = 0.25)
VlnPlot(eos_challenge, features = "predicted.id.score", group.by = "orig.ident")

View(markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))

#######COMPOSITIONAL ANALYSIS######
#frequencies per cluster
numberofcells         <- table(refquery$orig.ident, refquery$seurat_clusters)
numberofcells
totalcellsperorgan   <- c(sum(numberofcells[1,]), sum(numberofcells[2,]), sum(numberofcells[3,]), sum(numberofcells[4,]),
                          sum(numberofcells[5,]), sum(numberofcells[6,]), sum(numberofcells[7,]), sum(numberofcells[8,]),
                          sum(numberofcells[9,]), sum(numberofcells[10,]))
a                     <- cbind(numberofcells,totalcellsperorgan)
a
totalcellspercluster  <- c(sum(a[,1]), sum(a[,2]), sum(a[,3]), sum(a[,4]), sum(a[,5]), sum(a[,6]))
b                     <- rbind(a, totalcellspercluster)
b

c0 <- (b[c(1:6,9:10),1]/totalcellsperorgan[c(1:6,9:10)])*100
c1 <- (b[c(1:6,9:10),2]/totalcellsperorgan[c(1:6,9:10)])*100
c2 <- (b[c(1:6,9:10),3]/totalcellsperorgan[c(1:6,9:10)])*100
c3 <- (b[c(1:6,9:10),4]/totalcellsperorgan[c(1:6,9:10)])*100
c4 <- (b[c(1:6,9:10),5]/totalcellsperorgan[c(1:6,9:10)])*100
c5 <- (b[c(1:6,9:10),6]/totalcellsperorgan[c(1:6,9:10)])*100
c6 <- (b[c(1:6,9:10),7]/totalcellsperorgan[c(1:6,9:10)])*100
c7 <- (b[c(1:6,9:10),8]/totalcellsperorgan[c(1:6,9:10)])*100


c <- rbind(c0,c1,c2,c3,c4)
colSums(c)
rownames(c) =  rev(c("active eosinophils","basal eosinophils", "circulating eosinophils","immature eosinophils",  "eosinophil progenitors"))
c

par(mar=c(16,8,2,14))
pdf(file="Figures/Cluster_breakdown.pdf")
barplot(c, horiz=TRUE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=180, cex=.8),
        main = "Cluster breakdown per organ", 
        las = 1, 
        col= rev(col_vector[1:5]))
dev.off()

#hypergeometric test to test significance of enrichment in colon CR
set1        <- 1030+179 #total active eosinophils colon + colon CR
set2        <- 207 #cells in colonCR
overlap     <- 179 #active eos in colonCR
allterms    <- 1030+179+207 #union
phyper(overlap, set1, allterms-set1, set2, lower.tail=F) #7.365844e-25

#hypergeometric test to test significance of enrichment in stomach HP
set1        <- 73+284 #total active eosinophils stomach + stomach HP
set2        <- 359 #cells in stomach HP
overlap     <- 284 #active eos in colonCR
allterms    <- 73+284+359 #union
phyper(overlap, set1, allterms-set1, set2, lower.tail=F) #7.365844e-25


#LOCAL DENSITY
a <- DimPlot(refquery, order=T, group.by = "orig.ident", pt.size = .2, label=F, cols = col_vector) + theme(legend.position = "none") + labs(title=" ")

b <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
             cols = c( "bonemarrow"="black","bloodCR"="grey", "blood"="grey", "bonemarrowCR"= "grey",  "stomach"="grey", 
                       "stomachHP "="grey", "colonCR"="grey", "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
             order = c("bonemarrow", "bonemarrowCR", "bloodCR", "blood", "stomach","colon", "colonCR", "stomachHP"))+ theme_void() + 
  theme(legend.position = "none") + labs(title="bone marrow")

b1 <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
              cols = c( "bloodCR"="grey", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "black", "stomach"="grey", 
                        "stomachHP "="grey", "colonCR"="grey",  "colon"="grey","small intestine"="grey",  "bone marrow"="grey"),
              order = c("bonemarrowCR","bonemarrow",  "bloodCR", "blood", "stomach", "colon", "colonCR", "stomachHP"))+ theme_void() + 
  theme(legend.position = "none") + labs(title="bone marrow CR")

c <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
             cols = c( "blood"="black",  "bloodCR"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey", "stomach"="grey", 
                       "stomachHP "="grey", "colonCR"="grey",  "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
             order = c("blood","bonemarrowCR", "bloodCR", "bonemarrow",   "stomach", "colon", "colonCR", "stomachHP")) + theme_void()+ 
  theme(legend.position = "none") + labs(title="blood")

c1 <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
              cols = c( "bloodCR"="black", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey", "stomach"="grey", 
                        "stomachHP "="grey","colonCR"="grey", "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
              order = c("bloodCR", "bonemarrow", "blood", "stomach", "colon", "colonCR", "bonemarrowCR", "stomachHP"))+ theme_void() +
  theme(legend.position = "none") + labs(title="blood CR")

d <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
             cols = c( "colon"="black",  "bloodCR"="grey", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey","stomach"="grey", 
                       "stomachHP "="grey", "colonCR"="grey", "small intestine"="grey",  "bone marrow"="grey"),
             order = c("colon", "bloodCR", "bonemarrow", "blood", "stomach", "colonCR", "bonemarrowCR", "stomachHP"))+ theme_void() +
  theme(legend.position = "none") + labs(title="colon")

d1 <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
              cols = c( "bloodCR"="grey", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey", "stomach"="grey", 
                        "stomachHP "="grey", "colonCR"="black",  "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
              order = c("colonCR","bloodCR", "bonemarrow", "blood", "stomach", "colon",  "bonemarrowCR", "stomachHP")) + theme_void()+
  theme(legend.position = "none") + labs(title="colon CR")


e <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
             cols = c("stomach"="black", "bloodCR"="grey", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey", 
                      "stomachHP"="grey","colonCR"="grey",   "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
             order = c( "stomach", "colonCR","bloodCR", "bonemarrow", "blood",  "colon", "bonemarrowCR", "stomachHP")) + theme_void()+
  theme(legend.position = "none") + labs(title="stomach")


e1 <- DimPlot(refquery, reduction = "umap", group.by = c("orig.ident"), pt.size=.2,
              cols = c( "bloodCR"="grey", "blood"="grey",  "bonemarrow"= "grey",  "bonemarrowCR"= "grey", "stomach"="grey", 
                        "stomachHP"="black", "colonCR"="grey",  "colon"="grey", "small intestine"="grey",  "bone marrow"="grey"),
              order = c( "stomachHP", "bloodCR", "bonemarrow", "blood","stomach",  "colon", "colonCR", "bonemarrowCR"))+ theme_void() +
  theme(legend.position = "none") + labs(title="stomach HP")

ggarrange( b, b1, c, c1, d, d1, e, e1, ncol = 4, nrow = 2)+ ggsave("Figures/conditionUMAP.pdf", width = 16, height = 8)



###DIFFERENTIAL GENE EXPRESSION####
#active cluster only
DimPlot(refquery, group.by = "seurat_clusters")
Idents(refquery) <- "seurat_clusters"
active <- subset(refquery, ident= "active eosinophils")
Idents(active) <- "orig.ident"
DimPlot(active, group.by = "orig.ident")
active <- subset(active, ident= c("colon", "stomach", "stomachHP",  "colonCR"))
Idents(active) <- "orig.ident"
stomachHP_markers <- FindMarkers(active, ident.1 = "stomachHP", ident.2 = "stomach", verbose = FALSE, only.pos = F)
write.csv(stomachHP_markers,"stomachHP_markers.csv", row.names = TRUE)
sig_stomachHP_markers <- stomachHP_markers %>% filter(p_val_adj <0.05)

colonCR_markers <- FindMarkers(active, ident.1 = "colonCR", ident.2 = "colon", verbose = FALSE, only.pos = F)
write.csv(colonCR_markers,"colonCR_markers.csv", row.names = TRUE)
sig_colonCR_markers <-colonCR_markers %>% filter(p_val_adj <0.05)

a <- DEGs_volcano(stomachHP_markers, 0.05, 1, "StomachHP vs stomach", "grey89", upperylim = 45, xlim = 4)
b <- DEGs_volcano(colonCR_markers, 0.05, 1, "ColonCR vs colon", "grey89", upperylim = 50, xlim = 4)
c<- DEGs_volcano(blood_CR_markers, 0.05, 1, "BloodCR vs blood", "grey89", upperylim = 65, xlim = 3)
d<- DEGs_volcano(bonemarrow_CR_markers, 0.05, 1, "bonemarrowCR vs bonemarrow", "grey89", upperylim = 60, xlim = 2)
ggarrange(a, b, c,d,ncol = 4, nrow = 1) +
  ggsave("Figures/Volcanos.pdf", width = 20, height = 5)


#UPSET PLOT
commonbacterial <- sig_stomachHP_markers[rownames(sig_stomachHP_markers) %in% sig_colonCR_markers$Gene,]
nrow(commonbacterial)
View(commonbacterial)
write.csv(commonbacterial,"commonbacterial.csv", row.names = TRUE)
seta <- sig_stomachHP_markers$Gene
setb <- sig_colonCR_markers$Gene
a<- setdiff(seta,setb)
only_stomachHP <- sig_stomachHP_markers[a,]
write.csv(only_stomachHP,"only_stomachHP.csv", row.names = TRUE)
b<- setdiff(setb,seta)
only_colonCR <- sig_colonCR_markers[b,]
write.csv(only_colonCR,"only_colonCR.csv", row.names = TRUE)

#plot
sig_colonCR_markers$Gene <- rownames(sig_colonCR_markers) 
sig_stomachHP_markers$Gene <- rownames(sig_stomachHP_markers)
bacterialDEGS <- list("colonCR vs colon" = sig_colonCR_markers$Gene,
                      "stomachHP vs stomach" = sig_stomachHP_markers$Gene)

pdf(file="Figures/Upset.pdf", width = 7, height =6) 
upset(fromList(bacterialDEGS), group.by = "freq", keep.order = F,
      main.bar.color = "#3B9AB2",  matrix.color =  "#E8C520",
      sets.bar.color = "#F21A00", point.size = 8, empty.intersections=T, text.scale = 2)
dev.off()


####SIGNATURES####
#####proliferation and stemness
refquery <- AddModuleScore(refquery, features = cc.genes, name = "CC")
names(x = refquery[[]])
RidgePlot(refquery_sub, features="CC2", group.by = "orig.ident", col= c("gray45", "khaki2", "gray86", "orange2", "grey", "red", "grey16", "darkred")) +
  theme_classic() +  
  theme(text = element_text(size=25)) + labs (title = "", y = " ", x= "Cell cycle score") +theme(legend.position="none")+
  ggsave("Figures/CC_ridge_organs.pdf", width = 8, height = 6)

stemness_list <- list(c("Orc1l", "Impdh2", "Cct5", "Nap1l1", "Ccnd2", "Smo", "Mcm4", "Mcm5", "Hells", "Hnrnpa2b1", "Cct8", "Col18a1", "Sfrs3", 
                        "Rrm2", "Bub1", "Ncl", "Kpna2", "Shmt1", "Ipo5", "Ruvbl1", "Shroom3", "Dnahc11", "Cdc6", "Ttk", "Cks2", "Mcm2", "Fignl1", 
                        "Dph5", "Cdt1", "Cct3", "Eya2", "Pcna", "Set", "Prps1", "Fbl", "Dtymk", "Ssbp1", "Depdc6", "Top2a", "Csrp2"))

Idents(refquery) <- "orig.ident"
refquery_sub <- subset(refquery, idents= c("bonemarrow", "bonemarrowCR", "blood", "bloodCR", "colon", "colonCR", "stomach", "stomachHP"))
refquery_sub$orig.ident <- factor(x = refquery_sub$orig.ident , levels = c("bonemarrow", "bonemarrowCR", "blood", "bloodCR", "colon", "colonCR", "stomach", "stomachHP"))
refquery_sub <-AddModuleScore(refquery_sub, features= stemness_list,name = "Stemness")
names(x = refquery[[]])
VlnPlot(refquery_sub, features= "Stemness1", group.by = "orig.ident",col= c("gray45", "khaki2", "gray86", "orange2", "grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " Stemness score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Stemness_violin.pdf", width = 10, height = 6)

#test
wilcox.test(eos_prog$Stemness1, eos_immature$Stemness1, alternative = "two.sided") #p-value < 2.2e-16

colononly <- subset(refquery, idents="colon")
colononly <- AddModuleScore(colononly, features = cc.genes, name = "CC")
names(x = colononly[[]])
a <- RidgePlot(colononly, features="CC2", group.by = "seurat_clusters", rev(col_vector[1:4])) +
  theme_classic() +  xlim(-0.2,1.2)+
  theme(text = element_text(size=25)) + labs (title = "Colon", y = " ", x= "cell cycle score") +theme(legend.position="none")

colonCRonly <- subset(refquery, idents="colonCR")
Idents(colonCRonly) <- "seurat_clusters"
colonCRonly <- subset(colonCRonly, idents=c("immature eosinophils", "circulating eosinophils", "basal eosinophils", "active eosinophils"))
colonCRonly <- AddModuleScore(colonCRonly, features = cc.genes, name = "CC")
names(x = colonCRonly[[]])
b <- RidgePlot(colonCRonly, features="CC2", group.by = "seurat_clusters", rev(col_vector[1:4])) +
  theme_classic() +  xlim(-0.2,1.2)+
  theme(text = element_text(size=25)) + labs (title = "ColonCR", y = " ", x= "cell cycle score") +theme(legend.position="none")

stomachonly <- subset(refquery, idents="stomach")
stomachonly <- AddModuleScore(stomachonly, features = cc.genes, name = "CC")
names(x = stomachonly[[]])
c <- RidgePlot(stomachonly, features="CC2", group.by = "seurat_clusters", rev(col_vector[1:5])) +
  theme_classic() +  xlim(-0.2,1.2)+
  theme(text = element_text(size=25)) + labs (title = "Stomach", y = " ", x= "cell cycle score") +theme(legend.position="none")

stomachHPonly <- subset(refquery, idents="stomachHP")
stomachHPonly <- AddModuleScore(stomachHPonly, features = cc.genes, name = "CC")
names(x = stomachHPonly[[]])
d <- RidgePlot(stomachHPonly, features="CC2", group.by = "seurat_clusters", rev(col_vector[1:4])) +
  theme_classic() +  xlim(-0.2,1.2)+
  theme(text = element_text(size=25)) + labs (title = "StomachHP", y = " ", x= "cell cycle score") +theme(legend.position="none")

ggarrange(a,b,c,d, ncol = 2, nrow = 2)+ ggsave("Figures/CC_ridge.pdf", width = 16, height =16)

#compute signatures in stomach and colon only
active_challenge <- subset(active, idents=c("stomach", "stomachHP", "colon", "colonCR"))
active_challenge$orig.ident <- factor(x = active_challenge$orig.ident , levels = c("colon", "colonCR","stomach", "stomachHP"))
Idents(active_challenge) <-"orig.ident"

#granules
Granules_synthesis_list <-   list(c("Prg2", "Prg3",  "Epx", "Ear6", "Ear1", "Ear2"))
active_challenge <-AddModuleScore(active_challenge, features= Granules_synthesis_list,name = "Granules")
names(x = active_challenge[[]])

VlnPlot(active_challenge, features="Granules1", group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(Y = " Granule protein expression", title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Granules_violin.pdf", width = 8, height = 6)

wilcox.test(eos_colonCR$Granules1, eos_colon$Granules1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_stomachHP$Granules1, eos_stomach$Granules1, alternative = "two.sided") #p-value < 2.2e-16

#antigen processing
Antigen_processing <- list(c("H2-D1", "H2-Q7", "H2-K1", "H2-T23","H2-Q4","Tap1","Tapbp","B2m", "Psmb8", "Psme1",
                             "Psmb9", "Calr", "Psmb10", "Ncf1", "Fcer1g"))
active_challenge <-AddModuleScore(active_challenge, features= Antigen_processing, name = "Antigen_processing")
VlnPlot(active_challenge, features="Antigen_processing1", group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(y = "Antigen processing signature", title = "", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Antigen_processing_violin.pdf", width = 8, height = 6)

wilcox.test(eos_colonCR$Antigen_processing1, eos_colon$Antigen_processing1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_stomachHP$Antigen_processing1, eos_stomach$Antigen_processing1, alternative = "two.sided") #p-value < 2.2e-16

#ifng-regulated antimicrobial signature
IFNg_regulated <- list(c("Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
active_challenge <-AddModuleScore(active_challenge, features= IFNg_regulated,name = "IFNg_regulated")
names(x = active_challenge[[]])

VlnPlot(active_challenge, features="IFNg_regulated1", group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "IFNg-regulated antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Ifng_Antimicrobial_violin.pdf", width = 8, height = 6)

wilcox.test(eos_colonCR$IFNg_regulated1, eos_colon$IFNg_regulated1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_stomachHP$IFNg_regulated1, eos_stomach$IFNg_regulated1, alternative = "two.sided") #p-value < 2.2e-16

#ifng signature
IFNg_signature <- list(c("Ccl9", "Cxcl9", "Cxcl10", "Cd274", "Gbp2", "Irf1", "H2-D1", "H2-Q4", "H2-Q7", "H2-K1", 
                         "Gbp7", "Stat1", "Irgm1", "B2m", "Irf9", "Ifitm3", "H2-T23", "Icam1", "Jak2", "Irf2", "Ifngr2"))

active_challenge <-AddModuleScore(active_challenge, features= IFNg_signature,name = "IFNg_signature")
names(x = active_challenge[[]])

VlnPlot(active_challenge, features="IFNg_signature1", group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "IFNg score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/IFNg_signature_violin.pdf", width = 8, height = 6)

wilcox.test(eos_colonCR$IFNg_signature1, eos_colon$IFNg_signature1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_stomachHP$IFNg_signature1, eos_stomach$IFNg_signature1, alternative = "two.sided") #p-value < 2.2e-16

#antimicrobial peptides
antimicrobial_data <- AverageExpression(active_challenge, features = c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                                                                              "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))

antimicrobial <- list(c("S100a8", "S100a9", "Gbp2", "Prg2", "Epx", "Ear1", "Ear2", "Ear6", "Gbp7", "Ltf", 
                        "Lcn2", "Lyz2", "Irgm1", "Camp", "Adam17", "Serpine1", "H2-D1", "H2-T23", "H2-Q7"))
active_challenge <-AddModuleScore(active_challenge, features= antimicrobial,name = "antimicrobial")
names(x = active_challenge[[]])
VlnPlot(active_challenge, features="antimicrobial1", group.by = "orig.ident", cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Antimicrobial signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Antimicrobial.pdf", width = 8, height = 6)

Idents(active_challenge) <- "orig.ident"
eos_colon <- subset(active_challenge, idents = c("colon"))
eos_colonCR <- subset(active_challenge, idents = c("colonCR"))
wilcox.test(eos_colonCR$antimicrobial1, eos_colon$antimicrobial1, alternative = "two.sided") #p-value < 2.2e-16

eos_stomach <- subset(active_challenge, idents = c("stomach"))
eos_stomachHP <- subset(active_challenge, idents = c("stomachHP"))
wilcox.test(eos_stomachHP$antimicrobial1, eos_stomach$antimicrobial1, alternative = "two.sided") #p-value < 2.2e-16

pheatmap(antimicrobial_data$RNA, scale = 'row', cluster_rows = T, cluster_cols = F, border_color="white", 
         fontsize_number = 0.5, cellwidth =30, cellheight = 20, treeheight_row=0, treeheight_col=0,	
         angle_col = "45",color = colorRampPalette(rev(brewer.pal(n = 5, name =  "RdYlBu")))(100), filename = "Figures/Antimicrobial_heatmap.pdf")

#apoptosis
apoptosis.score <- list(c(BP[["GOBP_APOPTOTIC_PROCESS"]]))
active_challenge <- AddModuleScore(active_challenge, features= apoptosis.score,name = "apoptosis.score")
VlnPlot(active_challenge, features="apoptosis.score1", group.by = "orig.ident",cols= c("grey", "red", "grey16", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Apoptosis score ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Apoptosis_violin_challenge.pdf", width = 8, height = 6)
wilcox.test(eos_colonCR$apoptosis.score1, eos_colon$apoptosis.score1, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_stomachHP$apoptosis.score1, eos_stomach$apoptosis.score1, alternative = "two.sided") #p-value < 2.2e-16


######BONE MARROW AND BLOOD C.Rodentium#####
Idents(refquery) <- "orig.ident"
bonemarrow_CR_markers <- FindMarkers(refquery, ident.1 = "bonemarrowCR", ident.2 = "bonemarrow", verbose = FALSE)
View(bonemarrow_CR_nomarkers)
write.csv(bonemarrow_CR_markers %>% top_n(n = 200, wt = avg_log2FC),"bonemarrow_CR_markers.csv", row.names = TRUE)


blood_CR_markers <- FindMarkers(refquery, ident.1 = "bloodCR", ident.2 = "blood", verbose = FALSE)
View(blood_CR_markers)
write.csv(blood_CR_markers %>% top_n(n = 200, wt = avg_log2FC),"blood_CR_markers.csv", row.names = TRUE)


BMBlood_challenge <- subset(refquery, idents=c("bonemarrow", "bonemarrowCR", "blood", "bloodCR"))
BMBlood_challenge$orig.ident <- factor(x = BMBlood_challenge$orig.ident, levels = c("bonemarrow", "bonemarrowCR", "blood", "bloodCR"))
BMBlood_challenge <-AddModuleScore(BMBlood_challenge, features=Granules_synthesis_list ,name = "Granules")
VlnPlot(BMBlood_challenge, features="Granules1", group.by = "orig.ident", cols= c("gray45", "khaki2", "gray86", "orange2"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = "Granule signature", title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Granules_violin_BMblood.pdf", width = 8, height = 6)g

wilcox.test(eos_BMCR$Granules1, eos_BM$Granules1, alternative = "two.sided")
wilcox.test(eos_BLCR$Granules1, eos_BL$Granules1, alternative = "two.sided") 

antimicrobial2 <- list(c("S100a6", "S100a8", "S100a9", "Lcn2", "Lyz2", "Adam17", "Camp", "Ltf"))
BMBlood_challenge <-AddModuleScore(BMBlood_challenge, features= antimicrobial2,name = "Antimicrobial2")
VlnPlot(BMBlood_challenge, features="Antimicrobial21", group.by = "orig.ident", cols= c("gray45", "khaki2", "gray86", "orange2"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = "Non-gamma antimicrobial signature", title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Antimicrobial_nongamma_nongranule_violin_BMblood.pdf", width = 8, height = 6)

Idents(BMBlood_challenge) <- "orig.ident"
eos_BM <- subset(BMBlood_challenge, idents = c("bonemarrow"))
eos_BMCR <- subset(BMBlood_challenge, idents = c("bonemarrowCR"))
eos_BL <- subset(BMBlood_challenge, idents = c("blood"))
eos_BLCR <- subset(BMBlood_challenge, idents = c("bloodCR"))

wilcox.test(eos_BMCR$Antimicrobial21, eos_BM$Antimicrobial21, alternative = "two.sided") #p-value < 2.2e-16
wilcox.test(eos_BLCR$Antimicrobial21, eos_BL$Antimicrobial21, alternative = "two.sided") #p-value < 2.2e-16


BMBlood_challenge <-AddModuleScore(BMBlood_challenge, features= IFNg_signature,name = "IFNg_signature")
VlnPlot(BMBlood_challenge, features="IFNg_signature1", group.by = "orig.ident", cols= c("gray45", "khaki2", "gray86", "orange2"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = "IFNg signature", title = " ", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/IFNg_signature_violin_BMblood.pdf", width = 8, height = 6)

#proliferation score
BMBlood_challenge <- AddModuleScore(BMBlood_challenge, features = cc.genes, name = "CC")
names(x = BMBlood_challenge[[]])
RidgePlot(BMBlood_challenge, features="CC2", group.by = "orig.ident", cols= c("gray45", "khaki2", "gray86", "orange2")) +
  theme_classic() +  
  theme(text = element_text(size=25)) + labs (title = "Cell cycle score ", y = " ", x= " ") +theme(legend.position="none")+
  ggsave("Figures/CC_ridgeBMBLOOD.pdf", width = 8, height = 6)


#####CIRCULATING EOSINOPHILS DURING CHALLENGE#####
DimPlot(refquery, group.by = "orig.ident")
eos_circulating <- subset(refquery, idents = "circulating eosinophils")
DimPlot(eos_circulating,  group.by = "orig.ident")
FeaturePlot(eos_circulating, feature=c("Thbs1"), split.by = "orig.ident")
Idents(eos_circulating) <- "orig.ident"
circulating_markers_HP <- FindMarkers(eos_circulating, ident.1 = "stomachHP", ident.2 = "blood", verbose = FALSE)
View(circulating_markers_HP)
circulating_markers_CR <- FindMarkers(eos_circulating, ident.1 = "bloodCR", ident.2 = "blood", verbose = FALSE) 
View(circulating_markers_CR)

circulating_markers_CR_colonblood <- FindMarkers(eos_circulating, ident.1 = "colonCR", ident.2 = "bloodCR", verbose = FALSE) 
View(circulating_markers_CR_colonblood)

circulating_markers_CR$Gene <- rownames(circulating_markers_CR)
commoncirculating_HP_CR <- circulating_markers_HP[rownames(circulating_markers_HP) %in% circulating_markers_CR$Gene,]
commoncirculating_HP_CR
a <- DEGs_volcano(circulating_markers_CR, 0.05, 1, "circulating eos: bloodCR vs blood", "grey89", upperylim = 45, xlim = 3)
b <- DEGs_volcano(circulating_markers_HP, 0.05, 1, "circulating eos: stomachHP vs blood SS", "grey89", upperylim = 40, xlim = 5)
ggarrange(a, b, c, ncol = 3, nrow=1) 
DEGs_volcano(circulating_markers_CR_colonblood, 0.05, 1, "circulating eos: colonCR vs bloodCR", "grey89", upperylim = 10, xlim = 5)+ggsave("Figures/Circulating_volcanos.pdf", width = 5, height = 5)

write.csv(only_colonCR,"only_colonCR.csv", row.names = TRUE)


#####FGSEA####
BP_bonemarrowCR <- preranked_BP(bonemarrow_CR_markers)
View(BP_bonemarrowCR %>% filter(padj <0.05))

BP_bloodCR <- preranked_BP(blood_CR_markers)
View(BP_bloodCR %>% filter(padj <0.05))

BP_colonCR <- preranked_BP(colonCR_markers)
View(BP_colonCR %>% filter(padj <0.05))


BP_stomachHP <- preranked_BP(stomachHP_markers)
View(BP_stomachHP %>% filter(padj <0.05))

pos_bonemarrowCR <-  which(BP_bonemarrowCR$pathway == "ANTIMICROBIAL HUMORAL RESPONSE"|
                           BP_bonemarrowCR$pathway =="APOPTOTIC SIGNALING PATHWAY"|
                           BP_bonemarrowCR$pathway =="BIOLOGICAL PROCESS INVOLVED IN INTERACTION WITH HOST"|
                           BP_bonemarrowCR$pathway =="BLOOD VESSEL MORPHOGENESIS"|
                           BP_bonemarrowCR$pathway =="CELL ACTIVATION"|
                             BP_bonemarrowCR$pathway =="CYTOKINE PRODUCTION"|
                          BP_bonemarrowCR$pathway =="DEFENSE RESPONSE TO BACTERIUM"|
                        BP_bonemarrowCR$pathway =="GRANULOCYTE CHEMOTAXIS"|
                        BP_bonemarrowCR$pathway =="HUMORAL IMMUNE RESPONSE"|
                        BP_bonemarrowCR$pathway =="INNATE IMMUNE RESPONSE"|
                          BP_bonemarrowCR$pathway =="EXOCYTOSIS"|
                        BP_bonemarrowCR$pathway =="INTERLEUKIN 6 PRODUCTION"|
                        BP_bonemarrowCR$pathway =="INTRACELLULAR RECEPTOR SIGNALING PATHWAY"|
                        BP_bonemarrowCR$pathway =="MYELOID CELL DIFFERENTIATION"|
                          BP_bonemarrowCR$pathway ==  "RESPONSE TO BACTERIUM"|
                        BP_bonemarrowCR$pathway == "SECRETION"|
                        BP_bonemarrowCR$pathway =="PATTERN RECOGNITION RECEPTOR SIGNALING PATHWAY"|
                        BP_bonemarrowCR$pathway =="POSITIVE REGULATION OF NF KAPPAB TRANSCRIPTION FACTOR ACTIVITY"|
                        BP_bonemarrowCR$pathway =="RESPONSE TO TUMOR NECROSIS FACTOR"|
                        BP_bonemarrowCR$pathway =="WOUND HEALING"|
                        BP_bonemarrowCR$pathway =="SPROUTING ANGIOGENESIS")


pos_bloodCR <-  which(BP_bloodCR$pathway == "ANTIMICROBIAL HUMORAL RESPONSE"|
                        BP_bloodCR$pathway == "IMMUNE EFFECTOR PROCESS"|
                      BP_bloodCR$pathway == "LEUKOCYTE MEDIATED IMMUNITY"|
                      BP_bloodCR$pathway == "MYELOID LEUKOCYTE ACTIVATION"|
                        BP_bloodCR$pathway =="ESTABLISHMENT OF PROTEIN LOCALIZATION TO MEMBRANE"|
                        BP_bloodCR$pathway == "CELL ACTIVATION"|
                        BP_bloodCR$pathway == "EXOCYTOSIS"|
                        BP_bloodCR$pathway == "SECRETION")


pos_colonCR <- which(BP_colonCR$pathway ==   "ANTIGEN PROCESSING AND PRESENTATION"|
                       BP_colonCR$pathway == "ANTIGEN PROCESSING AND PRESENTATION OF EXOGENOUS PEPTIDE ANTIGEN VIA MHC CLASS I" |
                       BP_colonCR$pathway == "ANTIMICROBIAL HUMORAL RESPONSE"|
                       BP_colonCR$pathway == "CELL ACTIVATION"|
                       BP_colonCR$pathway == "TOLL LIKE RECEPTOR SIGNALING PATHWAY"|
                       BP_colonCR$pathway == "CYTOKINE PRODUCTION"|
                       BP_colonCR$pathway == "DEFENSE RESPONSE TO BACTERIUM"|
                       BP_colonCR$pathway == "DEFENSE RESPONSE TO GRAM POSITIVE BACTERIUM"|
                       BP_colonCR$pathway == "EXOCYTOSIS"|
                       BP_colonCR$pathway == "SECRETION"|
                       BP_colonCR$pathway =="INTERFERON GAMMA MEDIATED SIGNALING PATHWAY"|
                       BP_colonCR$pathway =="LEUKOCYTE MEDIATED CYTOTOXICITY"|
                       BP_colonCR$pathway =="PATTERN RECOGNITION RECEPTOR SIGNALING PATHWAY"|
                       BP_colonCR$pathway =="POSITIVE REGULATION OF ADAPTIVE IMMUNE RESPONSE"|
                       BP_colonCR$pathway =="RESPONSE TO BACTERIUM"|
                       BP_colonCR$pathway =="RESPONSE TO INTERFERON GAMMA"|
                       BP_colonCR$pathway =="T CELL MEDIATED IMMUNITY")

pos_stomachHP <- which(BP_stomachHP$pathway == "ANTIMICROBIAL HUMORAL RESPONSE" |
                         BP_stomachHP$pathway == "EPITHELIAL CELL DIFFERENTIATION" |
                         BP_stomachHP$pathway == "APOPTOTIC SIGNALING PATHWAY"|
                         BP_stomachHP$pathway == "GRANULOCYTE CHEMOTAXIS"|
                         BP_stomachHP$pathway == "MYELOID LEUKOCYTE ACTIVATION"|  
                         BP_stomachHP$pathway == "POSITIVE REGULATION OF INFLAMMATORY RESPONSE"|
                         BP_stomachHP$pathway == "CELL ACTIVATION" |
                         BP_stomachHP$pathway == "EXOCYTOSIS"|
                         BP_stomachHP$pathway == "SECRETION")


#plot
merged <- merge(BP_bonemarrowCR[pos_bonemarrowCR,c(1,5)], BP_bloodCR[pos_bloodCR,c(1,5)], by = "pathway" , all = T,
                suffixes = c(".BMCR",".bloodCR"))
merged2 <- merge(merged, BP_colonCR[pos_colonCR,c(1,5)], by = "pathway" , all = T)
merged3 <- merge(merged2,  BP_stomachHP[pos_stomachHP,c(1,5)], by = "pathway" , all = T, 
                 suffixes = c(".colonCR",".stomachHP"))
colnames(merged3) <- c("pathway","bonemarrowCR","bloodCR", "colonCR", "stomachHP")

merged3[is.na(merged3)] <- 0
rownames(merged3) <- merged3$pathway
rownames(merged3)<- tolower(rownames(merged3))
merged3[,1] <- NULL


paletteLength <- 20
myColor = colorRampPalette(c( "white", "red"))(paletteLength)
pheatmap(merged3, cluster_cols = F, cluster_rows = T, border_color="grey", 
         cellwidth =30, col=myColor,main = "Biological process enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10, filename = "GO.pdf")




