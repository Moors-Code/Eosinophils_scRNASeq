source('~/NAS/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Packages_functions.R', echo=TRUE)

ifg <- subset(eosinophil_pure, idents=c("colonIF", "colonCR", "colon"))
DimPlot(ifg)
ifg <- NormalizeData(ifg, normalization.method = "LogNormalize", scale.factor = 10000)
ifg <- FindVariableFeatures(ifg, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ifg)
ifg <- ScaleData(ifg, features = all.genes, vars.to.regress = "percent.mt")
ifg <- RunPCA(ifg, features = VariableFeatures(object = ifg))
DimPlot(ifg, reduction = "pca", group.by = "orig.ident")
ElbowPlot(ifg)
ifg <- FindNeighbors(ifg, dims = 1:20)
ifg <- FindClusters(ifg, resolution = 0.3)
ifg <- RunUMAP(ifg, dims = 1:20, return.model=TRUE)
DimPlot(ifg, group.by = "orig.ident")

DimPlot(ifg, group.by = "seurat_clusters")
ifg_active <- subset(ifg, idents=c(0,1))
DimPlot(ifg_active, reduction = "umap", pt.size = .5, label=F, group.by = "orig.ident", cols=c("grey", "darkblue", "darkred")) + 
  ggsave("UMAP_IF.pdf", width = 4.5, height = 2.5)


###DIFFERENTIAL GENE EXPRESSION###
colonIF_active_markers <- FindMarkers(ifg_active, ident.1 = "colonIF", ident.2 = "colon", only.pos = F, verbose = FALSE, min.pct = 0, logfc.threshold = 0.25)
View(colonIF_active_markers)
colonIF_active_markers$gene <- rownames(colonIF_active_markers)
write.csv(colonIF_active_markers %>% top_n(n = 200, wt = avg_log2FC),"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Challenge/colonIF_vs_colon.csv", row.names = TRUE)
BP_colonIF <- preranked_BP(colonIF_active_markers)
View(BP_colonIF %>% filter(abs(NES)>.5 & padj<0.05))

colonCR_active_markers <- FindMarkers(ifg_active, ident.1 = "colonCR", ident.2 = "colon", only.pos = F, verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
View(colonCR_active_markers)
colonCR_active_markers$gene <- rownames(colonCR_active_markers)

colonCR_vs_IF_active_markers <- FindMarkers(ifg_active, ident.1 = "colonCR", ident.2 = "colonIF", only.pos = F, verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
View(colonCR_vs_IF_active_markers)

intersect <- intersect(colonIF_active_markers$gene, colonCR_active_markers$gene)
intersect

onlyCR <- setdiff(colonCR_active_markers$gene, colonIF_active_markers$gene)
onlyCR

onlyIF <- setdiff(colonIF_active_markers$gene,colonCR_active_markers$gene)
onlyIF

DotPlot(ifg_active, features=c("Cd274", "Stat1", "Irf1"), group.by = "orig.ident", cols=c("grey", "darkblue", "darkred"), dot.scale = 10) + 
  RotatedAxis()+  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right") +
   ggsave("IFNy_IF.pdf", width = 5, height = 2)


###SIGNATURES###
Idents(ifg_active) <- "orig.ident"
Granules_list <-   list(c("S100a6","S100a8", "S100a9", "Prtn3", "Elane", "Tuba1b", "Cebpe", "Prg3",  "Prg2",  "Epx", 
                          "Ear6", "Ear1", "Ear2", "Cd63", "Ltf"  ,    "Lcn2"   ,  "Lyz2" ))
ifg_active <-AddModuleScore(ifg_active_2, features= Granules_list,name = "Granules")
VlnPlot(ifg_active, features="Granules1", group.by = "orig.ident", cols=c("grey", "darkblue", "darkred"), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ylim(-0.5,2.4)+
  labs(title = "Granule and antimicrobial protein expression", y = "Expression", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Granules_violin_IF.pdf", width = 8, height = 6)




