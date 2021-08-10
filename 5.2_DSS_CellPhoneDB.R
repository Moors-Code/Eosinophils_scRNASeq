#########################CellPhoneDB (K. Handler)################################

###load packages 
library(Seurat)
library(biomaRt)
library(ktplots)

###load annotated Seurat object from Schwarzfischer et al. 2021 (in preparation) and extract cells from WT DSS lamina propria cells 
merged_ann <- readRDS(file = ".../merged_ann.Rds")
Idents(merged_ann) <- "condition"
merged_ann_wtDSS <- subset(merged_ann, idents = "wtDSS")

###Generate input files for CellPhoneDB

###Count.txt file containing converted gene symbols from mouse to human and noramlized counts 
#Conversion of mouse gene symbols to human gene symbols  
#Code adapted partially from: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md 
#and https://www.cellphonedb.org/faq-and-troubleshooting
#take raw count data from "counts" slot of Seurat object and store in a matrix 
count_raw_meta <- GetAssayData(object = merged_ann_wtDSS, slot = "counts")[,colnames(x = merged_ann_wtDSS)]
#normalize raw counts 
count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
#load ensembl datasets from mouse and human and match to convert mouse to human symbols 
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(count_norm_meta) , 
                 mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=F),]
matrixB <- matrixA
matrixB$gene <- genesV2$Gene.stable.ID
rownames(matrixA) <- matrixB$gene
#save matrix as text file
write.table(matrixA, '.../cellphonedb_count_wtDSS.txt', sep='\t', quote=F, row.names = T)

##Meta.txt file containing information of cell annotation  
meta_data_meta <- cbind(rownames(merged_ann_wtDSS@meta.data), merged_ann_wtDSS@meta.data[,'annotation', drop=F])   
write.table(meta_data_meta, '.../cellphonedb_meta_wtDSS.txt', sep='\t', quote=F, row.names=F)


###Run CellPhoneDB in terminal 
#activate conda environment 
source /mnt/khandler/cell_phone_db/bin/activate
#generate new directory for CellPhoneDB output 
mkdir /.../ouputs_wtDSS 
#run CellPhoneDB 
cellphonedb method statistical_analysis /.../cellphonedb_meta_wtDSS.txt /.../cellphonedb_count_wtDSS.txt --output-path=/.../ouputs_wtDSS --project-name="IBD_wtDSS" --counts-data=ensembl 


###Plotting of significant interactions using ktplots 
#load output files from CellPhoneDB 
means <- read.delim(file = "/.../outputs_wtDSS/means.txt", check.names = FALSE)
pvals <- read.delim(file = ".../outputs_wtDSS/pvalues.txt", check.names = FALSE)

a <- plot_cpdb(cell_type1 = 'Eosinophils', cell_type2 = 'CD4+ T', scdata = merged_ann_WT,
               idents = 'annotation',
               means = meansAll, pvals = pvalsAll,keep_significant_only = TRUE,
               genes = c("THBS1", "CCL3", "CSF1", "CCL4","CCL4L2","PLAUR", "ICAM1","SELL",
                         "SIRPA", "CD80","THY1", "CTL4A", "ADORA2A", "IDE", "CD28", "PDCD1", "CD47", "TNFRSF1B",
                         "FAS", "CCL3L1", "SEMA7A"
               ), standard_scale = FALSE, highlight_size = 0.5, max_size = 20) + 
  small_axis(fontsize = 10) + small_grid() + small_guide() + small_legend(keysize = 2, fontsize = 10) 
a + ggtitle("L-R interactions between Eosinophils and CD4+ T cells significant interactions") + theme(title = element_text(size = 10)) + 
  ggsave("/.../Eos_CD4_sig_genes.pdf", width = 8, height = 10)


a <- plot_cpdb(cell_type1 = 'Eosinophils', cell_type2 = 'CD8+ T', scdata = merged_ann_WT,
               idents = 'annotation', 
               means = meansAll, pvals = pvalsAll,keep_significant_only = TRUE,
               genes = c("THBS1", "CCL3", "CSF1", "CCL4","CCL4L2","PLAUR", "ICAM1","SELL",
                         "SIRPA", "CD80","THY1", "CTL4A", "ADORA2A", "IDE", "CD28", "PDCD1", "CD47", "TNFRSF1B",
                         "FAS", "CCL3L1", "SEMA7A"
               ), standard_scale = FALSE, highlight_size = 0.5, max_size = 20) + 
  small_axis(fontsize = 10) + small_grid() + small_guide() + small_legend(keysize = 2, fontsize = 10) 
a + ggtitle("L-R interactions between Eosinophils and CD8+ T cells significant interactions") + theme(title = element_text(size = 10)) + 
  ggsave(".../Eos_CD8_sig_genes_CD4_and_CD8.pdf", width = 8, height = 10)


