####ACTIVE PERCENTAGE####
d <- rbind(c[5,], colSums(c[c(1,2,3,4),]))
d
pdf(file="Figures/ActivePercentage.pdf")
barplot(d, horiz=FALSE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=180, cex=.8),
        las = 1, 
        col= c(col_vector[1], "gray"))
dev.off()

#hypergeometric test to test significance of enrichment in colon
set1        <- 1430 #total active eosinophils 
set2        <- 1291 #cells in colon
overlap     <- 977 #active eos in colon 
allterms    <- 1430+1291 #union
phyper(overlap, set1, allterms-set1, set2, lower.tail=F) #1.255248e-121



####BASAL vs. ACTIVE COMPARISON#####
basal_active <- subset(eosinophils_steadystate, idents = c("active eosinophils", "basal eosinophils"))

#receptors
selected_markers <- c("Ptafr", "Ahr", "Fcgr3", "Fgfr1", "Ccr1", "Cxcr4", "Csf2rb", "Csf2rb2", "Il10ra", "Ifngr1", "Tgfbr2")
DotPlot(basal_active, features = selected_markers , dot.scale = 15) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1, size=20), axis.text.y = element_text(face="bold", size=20)) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = " ", y = "", x="") +
  ggsave("Figures/basalactive_dotplot.pdf", width = 7, height = 3.5)

#NFKB transcription factors
Nfkb.markers <- c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz", "Rela", "Relb")
DotPlot(eosinophils_steadystate, features = Nfkb.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="") +
  ggsave("Figures/Nfkb_dotplot.pdf", width = 7, height = 2.8)

#VlnPlot(eosinophils_steadystate, features="Nfkb1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(y = "Nfkb activity score", x="", title=" ") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Nfkb_activity_violin.pdf", width = 8, height = 4)

#regulatory proteins
a<-plot_density(eosinophils_steadystate, "Cd80", pal = "magma")
b<-plot_density(eosinophils_steadystate, "Cd9", pal = "magma")
c<-plot_density(eosinophils_steadystate, "Cd274", pal = "magma")
ggarrange(a, b, c, ncol = 3, nrow = 1) + ggsave("Figures/CostimulatoryMarkers.pdf", width = 15, height = 4)
  
Regulatory_list <-  list(c("Icosl", "Cd274", "Cd80", "Cd9"))
eosinophils_steadystate <-AddModuleScore(eosinophils_steadystate, features= Regulatory_list,name = "Regulatory")
names(x = eosinophils_steadystate[[]])
VlnPlot(eosinophils_steadystate, features="Regulatory1", group.by = "seurat_clusters", rev(col_vector[1:5]), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = " ", y = "Immune regulatory score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("Figures/Regulatory_violin.pdf", width = 8, height = 4)

#test
wilcox.test(eos_active$Regulatory1, eos_basal$Regulatory1, alternative = "two.sided") #p-value < 2.2e-16


#####DIFFERENTIALLY EXPRESSED GENES IN BASAL AND ACTIVE CLUSTER COLON VS SMALL INTESTINE#####
basal <- subset(eosinophils_steadystate, idents= "basal eosinophils")
Idents(basal) <- "orig.ident"
basal_CO_vs_SI <- FindMarkers(basal, ident.1 = "colon", ident.2 = "small intestine", only.pos = F)
View(basal_SPL_vs_SI)
basal_STO_vs_SI <- FindMarkers(basal, ident.1 = "stomach", ident.2 = "small intestine", only.pos = F)
basal_SPL_vs_SI <- FindMarkers(basal, ident.1 = "spleen", ident.2 = "small intestine", only.pos = F)


active <- subset(eosinophils_steadystate, idents= "active eosinophils")
Idents(active) <- "orig.ident"
active_CO_vs_SI <- FindMarkers(active, ident.1 = "colon", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
active_STO_vs_SI <- FindMarkers(active, ident.1 = "stomach", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
sig_active_CO_vs_SI <- active_CO_vs_SI %>% filter(p_val_adj<0.05)
write.csv(sig_active_CO_vs_SI,"sig_active_CO_vs_SI.csv", row.names = TRUE)

View(basal_CO_vs_SI)
BP_active_CO_vs_SI <- preranked_BP(active_CO_vs_SI)
sig_BP_active_CO_vs_SI <- BP_active_CO_vs_SI%>% filter(padj<0.05)
write.csv(sig_BP_active_CO_vs_SI,"BP_active_CO_vs_SI.csv", row.names = TRUE)

View(basal_STO_vs_SI)
a <- DEGs_volcano(active_CO_vs_SI, 0.05, 0, "Active eos: colon vs. SI","grey89", upperylim = 20, xlim = 1.5)
b <- DEGs_volcano(active_STO_vs_SI, 0.05, 0, "Active eos: stomach vs. SI", "grey89", upperylim = 20, xlim = 2)
c <- DEGs_volcano(basal_CO_vs_SI, 0.05, 0, "Basal eos: colon vs. SI","grey89", upperylim =20, xlim = 1.5)
d <- DEGs_volcano(basal_STO_vs_SI, 0.05, 0, "Basal eos: stomach vs. SI","grey89", upperylim = 20, xlim = 1.5)
e <- DEGs_volcano(basal_SPL_vs_SI, 0.05, 0, "Basal eos: spleen vs. SI","grey89", upperylim = 20, xlim = 1.5)

ggarrange(a, b, c,d, e,  ncol = 5, nrow = 1) +
  ggsave("Figures/Volcanos_orig.pdf", width = 25, height = 5)

