####BASAL vs. ACTIVE COMPARISON#####
basal_active <- subset(eosinophils_steadystate, idents = c("active eosinophils", "basal eosinophils"))

#receptors
selected_markers <- c("Ptafr", "Ahr", "Fcgr3", "Fgfr1", "Ccr1", "Cxcr4", "Csf2rb", "Csf2rb2", "Il10ra", "Ifngr1", "Tgfbr2")
DotPlot(basal_active, features = selected_markers , dot.scale = 15) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1, size=20), axis.text.y = element_text(face="bold", size=20)) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = " ", y = "", x="") +
  ggsave("basalactive_dotplot.pdf", width = 7, height = 3.5)

#NFKB transcription factors
Nfkb.markers <- c("Nfkb1", "Nfkbib", "Nfkbia", "Nfkbie", "Nfkb2", "Nfkbiz", "Rela", "Relb")
DotPlot(eosinophils_steadystate, features = Nfkb.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="") +
  ggsave("Nfkb_dotplot.pdf", width = 7, height = 2.8)

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

 
