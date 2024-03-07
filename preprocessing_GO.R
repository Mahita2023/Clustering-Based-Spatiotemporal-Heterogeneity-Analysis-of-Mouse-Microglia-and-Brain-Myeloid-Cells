
# Author:Songtao Wei
# Name:preprocessing_GO

library(dplyr)
library(Seurat)
library(patchwork)

dataFileName='Gene_data_afterraw.txt';
GE_QC.data<- read.delim(dataFileName,header = TRUE,sep = "\t",row.names ="geneNames", check.names=FALSE)
GE_QC <- CreateSeuratObject(counts = GE_QC.data, min.cells = 5, min.features = 50);
GE_QC
VlnPlot(GE_QC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot1 <- FeatureScatter(GE_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GE_QC <- NormalizeData(GE_QC, normalization.method = "LogNormalize", scale.factor = 10000)
Out_data <- as.data.frame(as.matrix(GE_QC@assays$RNA@data))
write.csv(Out_data, file = "GO_before_trans.csv")

dataFileName_2='GO_data_beforeselect.txt';
GE_QC_2.data<- read.delim(dataFileName_2,header = TRUE,sep = "\t",row.names ="GOIDs", check.names=FALSE)
GE_QC_2 <- CreateSeuratObject(counts = GE_QC_2.data, min.cells = 5, min.features = 50);
GE_QC_2
GE_QC_2 <- FindVariableFeatures(GE_QC_2, selection.method = "vst", nfeatures = 2000)
all.genes_2 <- rownames(GE_QC_2)
GE_QC_2 <- ScaleData(GE_QC_2, features = all.genes_2)
Out_data_2 <- as.data.frame(as.matrix(GE_QC_2@assays$RNA@scale.data))
write.csv(Out_data_2, file = "GO_after_scaling.csv")
Out_Gene_2 <- as.matrix(as.matrix(GE_QC_2@assays$RNA@var.features))
write.csv(Out_Gene_2, file = "GO_List_selected.csv")