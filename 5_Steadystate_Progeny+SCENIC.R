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

