####INTESTINAL PERCENTAGE####
d <- rbind(c[5,], colSums(c[c(1,2,3,4),]))
d
pdf(file="Figures/IntestinalPercentage.pdf")
barplot(d, horiz=FALSE,
        legend = T, border=NA,
        args.legend=list(bty = "n",x=180, cex=.8),
        las = 1, 
        col= c(col_vector[1], "gray"))
dev.off()

#hypergeometric test to test significance of enrichment in colon
set1        <- 1430 #total intestinal eosinophils 
set2        <- 1291 #cells in colon
overlap     <- 977 #intestinal eos in colon 
allterms    <- 1430+1291 #union
phyper(overlap, set1, allterms-set1, set2, lower.tail=F) #1.255248e-121


#####UPSET PLOT#####
Idents(eosinophils_steadystate) <- "orig.ident"
colonmarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "colon")
SImarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "small intestine")
spleenmarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "spleen")
bonemarrowmarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "bonemarrow")
bloodmarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "blood")
stomachmarkers <- FindMarkers(object = eosinophils_steadystate, ident.1 = "stomach")

colonmarkers$Gene <- rownames(colonmarkers) 
SImarkers$Gene <- rownames(SImarkers)
spleenmarkers$Gene <- rownames(spleenmarkers)
bonemarrowmarkers$Gene <- rownames(bonemarrowmarkers)
bloodmarkers$Gene <- rownames(bloodmarkers)
stomachmarkers$Gene <- rownames(stomachmarkers)


steadystateDEGS <- list("bone marrow" = bonemarrowmarkers$Gene,
                        "blood" = bloodmarkers$Gene,
                        "stomach" = stomachmarkers$Gene,
                        "spleen" = spleenmarkers$Gene,
                        "small intestine" = SImarkers$Gene,
                        "colon" = colonmarkers$Gene)


upset(fromList(steadystateDEGS), nsets = 7, group.by = "sets", keep.order = F, cutoff = 1,
      main.bar.color = "#3B9AB2",  matrix.color =  "#E8C520", mainbar.y.label = "organ-specific DEGs", sets.x.label= "total number of DEGs",
      sets.bar.color = "#F21A00", point.size = 8, empty.intersections=T, text.scale = 2)

pdf(file="Figures/Upset.pdf", width = 10, height =6) 
upset(fromList(steadystateDEGS),  nsets = 6, text.scale = 1.5)
dev.off()

#####DIFFERENTIALLY EXPRESSED GENES IN BASAL AND INTESTINAL CLUSTER COLON VS SMALL INTESTINE#####
basal <- subset(eosinophils_steadystate, idents= "basal eosinophils")
Idents(basal) <- "orig.ident"
basal_CO_vs_SI <- FindMarkers(basal, ident.1 = "colon", ident.2 = "small intestine", only.pos = F)
View(basal_SPL_vs_SI)
basal_STO_vs_SI <- FindMarkers(basal, ident.1 = "stomach", ident.2 = "small intestine", only.pos = F)
basal_SPL_vs_SI <- FindMarkers(basal, ident.1 = "spleen", ident.2 = "small intestine", only.pos = F)


intestinal <- subset(eosinophils_steadystate, idents= "intestinal eosinophils")
Idents(intestinal) <- "orig.ident"
intestinal_CO_vs_SI <- FindMarkers(intestinal, ident.1 = "colon", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
intestinal_STO_vs_SI <- FindMarkers(intestinal, ident.1 = "stomach", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
sig_intestinal_CO_vs_SI <- intestinal_CO_vs_SI %>% filter(p_val_adj<0.05)
write.csv(sig_intestinal_CO_vs_SI,"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Steadystate/sig_intestinal_CO_vs_SI.csv", row.names = TRUE)

View(basal_CO_vs_SI)
BP_intestinal_CO_vs_SI <- preranked_BP(intestinal_CO_vs_SI)
sig_BP_intestinal_CO_vs_SI <- BP_intestinal_CO_vs_SI%>% filter(padj<0.05)
write.csv(sig_BP_intestinal_CO_vs_SI,"/media/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Steadystate/BP_intestinal_CO_vs_SI.csv", row.names = TRUE)

View(basal_STO_vs_SI)
a <- DEGs_volcano(intestinal_CO_vs_SI, 0.05, 0, "Intestinal eos: colon vs. SI","grey89", upperylim = 20, xlim = 1.5)
b <- DEGs_volcano(intestinal_STO_vs_SI, 0.05, 0, "Intestinal eos: stomach vs. SI", "grey89", upperylim = 20, xlim = 2)
c <- DEGs_volcano(basal_CO_vs_SI, 0.05, 0, "Basal eos: colon vs. SI","grey89", upperylim =20, xlim = 1.5)
d <- DEGs_volcano(basal_STO_vs_SI, 0.05, 0, "Basal eos: stomach vs. SI","grey89", upperylim = 20, xlim = 1.5)
e <- DEGs_volcano(basal_SPL_vs_SI, 0.05, 0, "Basal eos: spleen vs. SI","grey89", upperylim = 20, xlim = 1.5)

ggarrange(a, b, c,d, e,  ncol = 5, nrow = 1) +
  ggsave("Figures/Volcanos_orig.pdf", width = 25, height = 5)

###COEFFICIENT OF VARIATION COLON SMALL INTESTINE#####
Idents(basal_intestinal) <- "orig.ident"
DEGScolon <- FindMarkers(basal_intestinal, ident.1 = "colon", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
sig_DEGScolon <- DEGScolon %>% filter(p_val_adj<0.05)
sig_DEGScolon$Gene <- rownames(sig_DEGScolon)
DEGSstomach <- FindMarkers(basal_intestinal, ident.1 = "stomach", ident.2 = "small intestine", only.pos = F, logfc.threshold = 0.25)
sig_DEGSstomach <- DEGSstomach %>% filter(p_val_adj<0.05)
sig_DEGSstomach$Gene <- rownames(sig_DEGSstomach)

#create normalized expression matrix of basal and intestinal cluster for colon and SI, and calculate CV
basal_data <- as(as.matrix(basal@assays$RNA@data), 'sparseMatrix') 
basal_data_CV <- cofVar(basal_data, complete = TRUE, treatment = NULL, type = NULL)
basal_data_CV$Gene <- rownames(basal_data_CV)
basal_CV_colon <- basal_data_CV[sig_DEGScolon$Gene, c("Gene", "cv")]
basal_CV_stomach <- basal_data_CV[sig_DEGSstomach$Gene, c("Gene", "cv")]

intestinal_data <- as(as.matrix(intestinal@assays$RNA@data), 'sparseMatrix') 
intestinal_data_CV <- cofVar(intestinal_data, complete = TRUE, treatment = NULL, type = NULL)
intestinal_data_CV$Gene <- rownames(intestinal_data_CV)
intestinal_CV_colon <- intestinal_data_CV[sig_DEGScolon$Gene, c("Gene", "cv")]
intestinal_CV_stomach <- intestinal_data_CV[sig_DEGSstomach$Gene, c("Gene", "cv")]

basal_CV_colon$Cluster <- "basal colon"
basal_CV_stomach$Cluster <- "basal stomach"
intestinal_CV_colon$Cluster <- "intestinal colon"
intestinal_CV_stomach$Cluster <- "intestinal stomach"
DEGS_CV <- rbind(basal_CV_colon,intestinal_CV_colon, basal_CV_stomach, intestinal_CV_stomach)

DEGS_CV$Cluster <- factor(DEGS_CV$Cluster, levels = c('basal colon','intestinal colon' ,'basal stomach','intestinal stomach'),ordered = TRUE)

ggplot(DEGS_CV, aes(x=Cluster, y=cv, fill=factor(Cluster))) + 
  geom_violin(color="black", width =.8) + theme_classic() + labs(fill = "cluster") + ylim(0,7)+
  ylab("Coefficent of variation") + ggtitle(" ")+
  scale_fill_manual(values = c(col_vector[2], col_vector[1],col_vector[2],col_vector[1])) + 
  stat_summary(fun = "mean", geom = "point", 
               size = 2) +
  stat_summary(fun.data = "mean_cl_normal",  #standard error of the mean
               geom = "errorbar",
               width = .1)+
  ggsave("Figures/CV.pdf", width = 9, height = 8)

meanb <- mean(basal_CV_colon$cv)
meani <- mean(intestinal_CV_colon$cv)
wilcox.test(intestinal_CV_colon$`cv`, basal_CV_colon$`cv`, alternative = "two.sided") #p-value  9.691e-07
wilcox.test(intestinal_CV_stomach$`cv`, basal_CV_stomach$`cv`, alternative = "two.sided") #p-value  5.121e-09


###COEFFICIENT OF VARIATION COLON STOMACH#####
Idents(eosinophils_steadystate) <- "orig.ident"
DEGS <- FindMarkers(eosinophils_steadystate, ident.1 = "colon", ident.2 = "stomach", only.pos = F)
DEGS$Gene <- rownames(DEGS)

#create normalized expression matrix of basal and intestinal cluster for colon and SI, and calculate CV
basal_data <- as(as.matrix(basal@assays$RNA@data), 'sparseMatrix') 
basal_data_CV <- cofVar(basal_data, complete = TRUE, treatment = NULL, type = NULL)
basal_data_CV$Gene <- rownames(basal_data_CV)
basal_CV <- basal_data_CV[DEGS$Gene, c("Gene", "cv")]
nrow(basal_CV)

intestinal_data <- as(as.matrix(intestinal@assays$RNA@data), 'sparseMatrix') 
intestinal_data_CV <- cofVar(intestinal_data, complete = TRUE, treatment = NULL, type = NULL)
intestinal_data_CV$Gene <- rownames(intestinal_data_CV)
intestinal_CV <- intestinal_data_CV[DEGS$Gene, c("Gene", "cv")]
nrow(intestinal_CV)

basal_CV$Cluster <- "basal"
intestinal_CV$Cluster <- "intestinal"
DEGS_CV <- rbind(basal_CV,intestinal_CV)

ggplot(DEGS_CV, aes(x=Cluster, y=cv, fill=factor(Cluster))) + 
  geom_violin(color="black", width =.5) + theme_classic() + labs(fill = "cluster") + ylim(0,7)+
  ylab("CV of colon vs small intestine DEGs") + ggtitle(" ")+
  scale_fill_manual(values = rev(col_vector[1:2])) + 
  stat_summary(fun = "mean", geom = "point", 
               size = 3) +
  stat_summary(fun.data = "mean_cl_normal",  #standard error of the mean
               geom = "errorbar",
               width = .1)+
  ggsave("Figures/CV.pdf", width = 9, height = 8)

meanb <- mean(basal_CV$cv)
meani <- mean(intestinal_CV$cv)
wilcox.test(intestinal_CV$`cv`, basal_CV$`cv`, alternative = "two.sided") #p-value  0.0009678



#####UPSET PLOT#####
#basal
Idents(basal) <- "orig.ident"
colonmarkers <- FindMarkers(object = basal, ident.1 = "colon")
sig_colonmarkers <- colonmarkers %>% filter(p_val_adj<0.05)
SImarkers <- FindMarkers(object = basal, ident.1 = "small intestine")
sig_SImarkers <- SImarkers %>% filter(p_val_adj<0.05)
spleenmarkers <- FindMarkers(object = basal, ident.1 = "spleen")
sig_spleenmarkers <- spleenmarkers %>% filter(p_val_adj<0.05)
stomachmarkers <- FindMarkers(object = basal, ident.1 = "stomach")
sig_stomachmarkers <- stomachmarkers %>% filter(p_val_adj<0.05)


sig_colonmarkers$Gene <- rownames(sig_colonmarkers) 
sig_SImarkers$Gene <- rownames(sig_SImarkers)
sig_spleenmarkers$Gene <- rownames(sig_spleenmarkers)
sig_stomachmarkers$Gene <- rownames(sig_stomachmarkers)

basalDEGS <- list("stomach" = sig_stomachmarkers$Gene,
                        "spleen" = sig_spleenmarkers$Gene,
                        "small intestine" = sig_SImarkers$Gene,
                        "colon" = sig_colonmarkers$Gene)

pdf(file="Figures/Upset_basal.pdf", width = 8, height =6) 
upset(fromList(basalDEGS), group.by = "freq", keep.order = F, mainbar.y.max = 25,
      mainbar.y.label = "intersection size", sets.x.label= "DEGs",
     point.size = 5, text.scale = 2)
dev.off()


#intestinal
Idents(intestinal) <- "orig.ident"
colonmarkers <- FindMarkers(object = intestinal, ident.1 = "colon")
sig_colonmarkers <- colonmarkers %>% filter(p_val_adj<0.05)
SImarkers <- FindMarkers(object = intestinal, ident.1 = "small intestine")
sig_SImarkers <- SImarkers %>% filter(p_val_adj<0.05)
spleenmarkers <- FindMarkers(object = intestinal, ident.1 = "spleen")
sig_spleenmarkers <- spleenmarkers %>% filter(p_val_adj<0.05)
stomachmarkers <- FindMarkers(object = intestinal, ident.1 = "stomach")
sig_stomachmarkers <- stomachmarkers %>% filter(p_val_adj<0.05)


sig_colonmarkers$Gene <- rownames(sig_colonmarkers) 
sig_SImarkers$Gene <- rownames(sig_SImarkers)
sig_spleenmarkers$Gene <- rownames(sig_spleenmarkers)
sig_stomachmarkers$Gene <- rownames(sig_stomachmarkers)

intestinalDEGS <- list("stomach" = sig_stomachmarkers$Gene,
                  "spleen" = sig_spleenmarkers$Gene,
                  "small intestine" = sig_SImarkers$Gene,
                  "colon" = sig_colonmarkers$Gene)

pdf(file="Figures/Upset_intestinal.pdf", width = 8, height =6) 
upset(fromList(intestinalDEGS), group.by = "freq", keep.order = F, mainbar.y.max = 25,
      mainbar.y.label = "intersection size", sets.x.label= "DEGs",
      point.size = 5, text.scale = 2)
dev.off()
