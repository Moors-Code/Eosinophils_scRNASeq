####GENERATE LOOM FILES####

#get cellIDs from eosinophil_pure Seurat object
library(stringr)
library(Seurat)
load("~/NAS/Coco/Collaborations/Eosinophils BD/Data analysis/Final/Steadystate/eosinophil_pure.RData")
Cells(eosinophil_pure)
cellID <- str_split(Cells(eosinophil_pure), pattern = "_", simplify = T)[,2]
cellID
length(cellID)
write.table(cellID, file= "/media/Coco/Collaborations/Eosinophils BD/Data analysis/CellRank/cellIDs_eospure.tsv", quote=F, col.names = F, row.names = F, sep = "\t")


#merge multiple bams
samtools merge merged.bam bam1.bam bam2.bam #merge multiple runs if needed

#extract cellIDs from bam file
samtools sort merged.bam > sorted_merged.bam #sort bam file
samtools index sorted_merged.bam  #index bam file
wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux #download tool from 10x page
./subset-bam_linux -b sorted_merged.bam -c cellIDs_eospure.tsv -o eospure.bam --cores 15 #cellID.tsv single column with 1 cellID per line

#from bam to loom
#install velocyto into samtools conda env
conda install numpy scipy cython numba matplotlib scikit-learn h5py click 
pip install velocyto

#change the tags of filtered bam file
cat    <(samtools view -HS eospure_fast.bam) <(samtools view eospure_fast.bam | grep "MA:Z:*" | sed "s/MA:Z:/UB:Z:/" ) | samtools view -Sb -@6 > eospure_fast.bam.for_velocyto.bam 
samtools index steadystate.for_velocyto.bam --threads 50

#run velocyto
conda activate samtools
velocyto run -b cellIDs_eospure.tsv -o . -m eospure_fast.bam.for_velocyto.bam gencode.vM25.basic.annotation.gtf --samtools-memory 2000 --samtools-threads 20

